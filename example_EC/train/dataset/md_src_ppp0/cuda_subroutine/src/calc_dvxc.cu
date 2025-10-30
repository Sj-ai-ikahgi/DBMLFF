#include "../include/util.h"
#include "../include/calc_dvxc.h"

#define BLOCKSIZE 1024

__global__ void cudakernel_calc_dvxc_loop(double *rho, double *rho_p, double *dvxc, int n1, int n2, int n3, double pi, double f23, double f53, double fact)
{
    int k = blockIdx.x + 1;
    int iteration_counts = n1 * n2;
    int tidx = threadIdx.x;
    double uxc = 0.0;
    double uxc_tmp = 0.0;

    // int iter = tidx;
    // while (iter < iteration_counts)
    // {
    //     int j = iter / n1 + 1;
    //     int i = iter % n1 + 1;
    //     double vxc_tmp1 = UxcCA(rho(i, j, k), &uxc);
    //     vxc_tmp1 += fact * f53 * pow(abs(rho(i, j, k)), f23);
    //     double vxc_tmp2 = UxcCA(rho(i, j, k) + 1.0e-4, &uxc_tmp);
    //     vxc_tmp2 += fact * f53 * pow(abs(rho(i, j, k) + 1.0e-4), f23);
    //     dvxc[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = (vxc_tmp2 - vxc_tmp1) / 1.0e-4; 
    //     dvxc[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] *= rho_p[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2];
    //     iter += blockDim.x;
    // }


    int iter = tidx;
    while (iter < iteration_counts)
    {
        int j = iter / n1 + 1;
        int i = iter % n1 + 1;
        double vxc_tmp1;
        double vxc_tmp2;
        // double delt = 1.0e-4;
        double delt = 1.0e-3;
        double tderi;
        double sign;

        if (rho(i,j,k) < 0.0) {
            vxc_tmp1 = UxcCA(rho(i, j, k) - delt, &uxc);
            vxc_tmp2 = UxcCA(rho(i, j, k), &uxc_tmp);
            vxc_tmp1 += -1.0 * fact * f53 * pow(abs(rho(i, j, k) - delt), f23);
            vxc_tmp2 += -1.0 * fact * f53 * pow(abs(rho(i, j, k)), f23);
        } else {
            vxc_tmp1 = UxcCA(rho(i, j, k), &uxc);
            vxc_tmp2 = UxcCA(rho(i, j, k) + delt, &uxc_tmp);
            vxc_tmp1 += fact * f53 * pow(abs(rho(i, j, k)), f23);
            vxc_tmp2 += fact * f53 * pow(abs(rho(i, j, k) + delt), f23);
        }

        tderi = (vxc_tmp2 - vxc_tmp1) / delt;
        if (tderi < 0.0) {
            sign = -1.0;
        } else {
            sign = 1.0;
        }
        // tderi = sign * min(abs(tderi),100.0);
        tderi = sign * min(abs(tderi),1000.0);

        dvxc[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = tderi; 
        dvxc[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] *= rho_p[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2];
        iter += blockDim.x;
    }

}

void loop(double *rho, double *rho_p, double *dvxc, int n1, int n2, int n3, double pi, double f23, double f53, double fact)
{
    dim3 gridDim(n3, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_calc_dvxc_loop<<<gridDim, blockDim, 0, 0>>>(rho, rho_p, dvxc, n1, n2, n3, pi, f23, f53, fact);
}

void fcc_calc_dvxc(double *rho, double *rho_p, double *dvxc, int n1, int n2, int n3, double pi)
{
    double f53 = 5.0 / 3.0;
    double f23 = 2.0 / 3.0;
    double fact = 3 / 10.0 * pow((3 * pow(pi, 2)), (2.0 / 3));
    // cudaMemset(dvxc, 0, n1 * n2 * n3 * sizeof(double));
    loop(rho, rho_p, dvxc, n1, n2, n3, pi, f23, f53, fact);
}
