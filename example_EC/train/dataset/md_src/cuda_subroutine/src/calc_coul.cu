#include "../include/util.h"
#include "../include/calc_coul.h"

#define BLOCKSIZE 1024

extern void *buffer_h, *buffer_d;
extern int default_stream;
extern cudaStream_t stream1, stream2, stream3, stream4;
extern int bytes_of_bufferh, bytes_of_bufferd;

__global__ void cudakernel_calc_coul_loop(int n1, int n2, int n3, cuDoubleComplex *rho_q, double *ALI, double pi)
{
    int k = blockIdx.x + 1;
    int k1 = k - 1;
    if (k > n3 / 2)
    {
        k1 = k - n3 - 1;
    }
    int tidx = threadIdx.x;
    int iteration_counts = n1 * n2;

    int iter = tidx;
    while (iter < iteration_counts)
    {
        int j = iter / n1 + 1;
        int i = iter % n1 + 1;
        int j1 = j - 1;
        int i1 = i - 1;
        if (j > n2 / 2)
        {
            j1 = j - n2 - 1;
        }
        if (i > n1 / 2)
        {
            i1 = i - n1 - 1;
        }
        double gkx = 2 * pi * (ALI[0] * i1 + ALI[3] * j1 + ALI[6] * k1);
        double gky = 2 * pi * (ALI[1] * i1 + ALI[4] * j1 + ALI[7] * k1);
        double gkz = 2 * pi * (ALI[2] * i1 + ALI[5] * j1 + ALI[8] * k1);
        double q2 = pow(gkx, 2) + pow(gky, 2) + pow(gkz, 2);
        if (q2 > 1.0e-20)
        {
            rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] =  cuCmul(make_cuDoubleComplex(4 * pi / q2, 0), rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2]);
        }
        else
        {
            rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = make_cuDoubleComplex(0.0, 0.0);
        }
        iter += blockDim.x;
    }
}

void loop(int n1, int n2, int n3, cuDoubleComplex *rho_q, double *ALI, double pi)
{
    dim3 gridDim(n3, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_calc_coul_loop<<<gridDim, blockDim, 0, 0>>>(n1, n2, n3, rho_q, ALI, pi);
}

void fcc_calc_coul(double *rho, double *vcoul, int n1, int n2, int n3, double *AL, double pi)
{
    // cudaMemset(buffer_d, 0, bytes_of_bufferd);
    // memset(buffer_h, 0, bytes_of_bufferh);
    cuDoubleComplex *rho_q_d = reinterpret_cast<cuDoubleComplex *>(buffer_d);
    cuDoubleComplex *fftbuffer_d = reinterpret_cast<cuDoubleComplex *>(rho_q_d + n1 * n2 * n3);
    double *zero_buffer_d = reinterpret_cast<double *>(fftbuffer_d + n1 * n2 * n3);
    double *ALI_d = reinterpret_cast<double *>(zero_buffer_d + n1 * n2 * n3);

    // double *AL_h = (double *)malloc(3 * 3 * sizeof(double));
    double *AL_h = reinterpret_cast<double *>(buffer_h);
    cudaMemcpy(AL_h, AL, 3 * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    // double *ALI_h = (double *)malloc(3 * 3 * sizeof(double));
    double *ALI_h = AL_h + 3 * 3;

    // double *ALI_d = reinterpret_cast<double *>(buffer_d);
    // cudaMalloc((void **)&ALI_d, 3 * 3 * sizeof(double));

    get_ALI(AL_h, ALI_h);  
    cudaMemcpy(ALI_d, ALI_h, 3 * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cufftHandle plan = get_plan(n1, n2, n3);
    // cuDoubleComplex *fftbuffer_d = reinterpret_cast<cuDoubleComplex *>(ALI_d + 3 * 3);
    // double *zero_buffer_d = reinterpret_cast<double *>(fftbuffer_d + n1 * n2 * n3);
    // cudaMalloc((void **)&fftbuffer_d, n1 * n2 * n3 * sizeof(cuDoubleComplex));
    // cudaMalloc((void **)&zero_buffer_d, n1 * n2 * n3 * sizeof(double));
    // cudaMemset(zero_buffer_d, 0.0, n1 * n2 * n3 * sizeof(double));
    make_dcmplx(fftbuffer_d, rho, 0.0, n1 * n2 * n3);   
    // debug
    // cuDoubleComplex *fftbuffer_h = (cuDoubleComplex *)malloc(n1 * n2 * n3 * sizeof(cuDoubleComplex));
    // cudaMemcpy(fftbuffer_h, fftbuffer_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < n1 * n2 * n3; i++){
    //     printf("%20.15E %20.15E \r\n", fftbuffer_h[i].x, fftbuffer_h[i].y);
    // }
    // debug    
    //cufftPlan3d(&plan, n3, n2, n1, CUFFT_Z2Z);  
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_INVERSE); 
    double factor = 1.0 / (double)(n1 * n2 * n3);
    factorByConstant(fftbuffer_d, n1 * n2 * n3, factor);
    // cuDoubleComplex *rho_q_d = reinterpret_cast<cuDoubleComplex *>(zero_buffer_d + n1 * n2 * n3);
    // cudaMalloc((void **)&rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex));
    cudaMemcpy(rho_q_d, fftbuffer_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    // debug
    // cuDoubleComplex *rho_q_h = (cuDoubleComplex *)malloc(n1 * n2 * n3 * sizeof(cuDoubleComplex));
    // cudaMemcpy(rho_q_h, rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < 4; i++)
    //     printf("%lf %lf \r\n", rho_q_h[i].x, rho_q_h[i].y);
    // debug
    loop(n1, n2, n3, rho_q_d, ALI_d, pi);
    cudaMemcpy(fftbuffer_d, rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_FORWARD);
    getdreal(fftbuffer_d, vcoul, n1 * n2 * n3);
    // cudaFree(ALI_d);
    // cudaFree(fftbuffer_d);
    // cudaFree(zero_buffer_d);
    // cudaFree(rho_q_d);
    // free(AL_h);
    // free(ALI_h);
    // cufftDestroy(plan);
}