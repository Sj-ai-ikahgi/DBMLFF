#include "../include/util.h"
#include "../include/nonbond.h"

#define BLOCKSIZE 512
#define maxThreadsPerBlock BLOCKSIZE
#define minBlocksPerMultiprocessor 1

extern void *buffer_h, *buffer_d;
extern cudaStream_t stream1, stream2, stream3, stream4;

__global__ void cudakernel_loop_82to96(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int k1 = bidx - id3;
    int iteration_counts = (2 * id2 + 1) * (2 * id1 + 1);

    double dx3 = k1 * 1.0 / n3;
    int iter = tidx;
    while (iter < iteration_counts)
    {
        int j1 = iter / (2 * id1 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        double dx1 = i1 * 1.0 / n1;
        double dx2 = j1 * 1.0 / n2;
        double dx = AL(1, 1) * dx1 + AL(1, 2) * dx2 + AL(1, 3) * dx3;
        double dy = AL(2, 1) * dx1 + AL(2, 2) * dx2 + AL(2, 3) * dx3;
        double dz = AL(3, 1) * dx1 + AL(3, 2) * dx2 + AL(3, 3) * dx3;
        dxyz_box(1, i1, j1, k1) = dx;
        dxyz_box(2, i1, j1, k1) = dy;
        dxyz_box(3, i1, j1, k1) = dz;
        iter += blockDim.x;
    }
}

__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_134to149(int *icent, double *w_cent, int *nat_cent, double *xatom_m, int itype_mol, int imol, double *xatom_cent, int ncent_itype_mol, int natom_mm, int nmolm)
{
    int tidx = threadIdx.x;
    int iteration_counts = ncent_itype_mol * 3;
    int iter = tidx;

    while (iter < iteration_counts)
    {
        int ii = iter / 3 + 1;
        int ixyz = iter % 3 + 1;

        double x1 = xatom_m(ixyz, icent(1, ii, itype_mol), imol, itype_mol);
        double xc = x1 * w_cent(1, ii, itype_mol);
        double w_sum = w_cent(1, ii, itype_mol);

        for (int jj = 2; jj < nat_cent(ii, itype_mol) + 1; jj++)
        {
            double x2 = xatom_m(ixyz, icent(jj, ii, itype_mol), imol, itype_mol);
            if (abs(x2 + 1 - x1) < abs(x2 - x1))
                x2 = x2 + 1;
            if (abs(x2 - 1 - x1) < abs(x2 - x1))
                x2 = x2 - 1;
            xc = xc + x2 * w_cent(jj, ii, itype_mol);
            w_sum += w_cent(jj, ii, itype_mol);
        }
        xc /= w_sum;
        xatom_cent(ixyz, ii) = xc;
        iter += blockDim.x;
    }
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_149to238(double *xatom_cent, double *z_cent, int *itype_cent, double *funcr2, double *AL,int nr, double Rm2, 
                                         double *dxyz_box, double *rho, double *rho_z, int id3, int id2, int id1, int n3, int n2, int n1, 
                                         double Rbox2, int itype_mol, double pi, int* ion_type_cent, int* imax_ion, double *r_ion, double *rho_ion, 
                                         double vol, double* box, double *box2, double *box3, double *Q_type, double *z_ion, int imax_nr, int imax_ntype_cent)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id3 + 1) * (2 * id2 + 1) * (2 * id1 + 1);
    int ia = bidx + 1;
    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;

    double dx00 = AL(1, 1) * dx10 + AL(1, 2) * dx20 + AL(1, 3) * dx30;
    double dy00 = AL(2, 1) * dx10 + AL(2, 2) * dx20 + AL(2, 3) * dx30;
    double dz00 = AL(3, 1) * dx10 + AL(3, 2) * dx20 + AL(3, 3) * dx30;

    double fact11 = nr / Rm2;
    double fact22 = 0;
    double fact = 0;
    double fact2 = 0;

    __shared__ double shared_mem[3 * blocksize];
    double partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;
    if (itype_ion != 0) 
    {
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }
    int iter = tidx;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00;
        double dd = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
        if (dd < Rbox2)
        {
            double d = sqrt(dd);
            double yy = d * fact11;
            int ir = yy;
            double x = yy - ir;
            double f1 = 1 - x - 0.5 * x * (1 - x);
            double f2 = x + x * (1 - x);
            double f3 = -0.5 * x * (1 - x);
            ir += 1;
            box(i1, j1, k1, ia) = funcr2(ir, itype, itype_mol) * f1 + funcr2(ir + 1, itype, itype_mol) * f2 + funcr2(ir + 2, itype, itype_mol) * f3;
            box2(i1, j1, k1, ia) = exp(-pow((d / 1.5), 2));
            if (itype_ion != 0) 
            {
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x);
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1 - x);
                ir = ir + 1;
                box3(i1, j1, k1, ia) = rho_ion(ir, itype_ion) * f1 + rho_ion(ir + 1, itype_ion) * f2 + rho_ion(ir + 2, itype_ion) * f3;
            }
        }
        partial_sum1 += box(i1, j1, k1, ia);
        partial_sum2 += box2(i1, j1, k1, ia);
        partial_sum3 += box3(i1, j1, k1, ia);
        iter += blockDim.x;
    }
    __syncthreads();
    reduce_sharedmem<double, blocksize>(shared_mem + 0 * blocksize, partial_sum1, tidx, shared_mem + 0 * blocksize);
    reduce_sharedmem<double, blocksize>(shared_mem + 1 * blocksize, partial_sum2, tidx, shared_mem + 1 * blocksize);
    reduce_sharedmem<double, blocksize>(shared_mem + 2 * blocksize, partial_sum3, tidx, shared_mem + 2 * blocksize);
    if (tidx == 0)
    {
        shared_mem[0 * blocksize] = shared_mem[0 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[1 * blocksize] = shared_mem[1 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[2 * blocksize] = shared_mem[2 * blocksize] * vol / (n1 * n2 * n3);
    }
    __syncthreads();
    fact = (Q_type(itype, itype_mol) - shared_mem[0 + 0 * blocksize]) / shared_mem[0 + 1 * blocksize];
    if (itype_ion != 0)
    {
        fact2 = z_ion(itype_ion) / shared_mem[2 * blocksize];
    }
    iter = tidx;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (k0 + k1 - 1 + n3) % n3 + 1;
        int j = (j0 + j1 - 1 + n2) % n2 + 1;
        int i = (i0 + i1 - 1 + n1) % n1 + 1;   
        atomicAdd(&rho(i, j, k), (box(i1, j1, k1, ia) + fact * box2(i1, j1, k1, ia)));
        atomicAdd(&rho_z(i, j, k), (-fact2 * box3(i1, j1, k1, ia))); 
        iter += blockDim.x;
    }
}
template<unsigned int blocksize>
__global__ void cudakernel_loop_149to238_step_1(double *xatom_cent, double *z_cent, int *itype_cent, double *funcr2, double *AL,int nr, double Rm2, 
                                         double *dxyz_box, double *rho, double *rho_z, int id3, int id2, int id1, int n3, int n2, int n1, 
                                         double Rbox2, int itype_mol, double pi, int* ion_type_cent, int* imax_ion, double *r_ion, double *rho_ion, 
                                         double vol, double* box, double *box2, double *box3, double *Q_type, double *z_ion, int imax_nr, int imax_ntype_cent,
                                         double *block_reduce_sum)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int tidx = threadIdx.x;
    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;

    float dx10 = x1 - (i0 - 1.0) / n1;
    float dx20 = x2 - (j0 - 1.0) / n2;
    float dx30 = x3 - (k0 - 1.0) / n3;

    float dx00 = (float)AL(1, 1) * dx10 + (float)AL(1, 2) * dx20 + (float)AL(1, 3) * dx30;
    float dy00 = (float)AL(2, 1) * dx10 + (float)AL(2, 2) * dx20 + (float)AL(2, 3) * dx30;
    float dz00 = (float)AL(3, 1) * dx10 + (float)AL(3, 2) * dx20 + (float)AL(3, 3) * dx30;

    float fact11 = nr / Rm2;
    float fact22 = 0;

    __shared__ float shared_mem[blocksize];

    float partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;

    if (itype_ion != 0) 
    {
        fact22 = (imax_ion(itype_ion) - 1) / (float)r_ion(imax_ion(itype_ion), itype_ion);
    }

    while (loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        box(i1, j1, k1, ia) = 0;
        // box2(i1, j1, k1, ia) = 0;
        box3(i1, j1, k1, ia) = 0;

        float dx = (float)dxyz_box(1, i1, j1, k1) - dx00;
        float dy = (float)dxyz_box(2, i1, j1, k1) - dy00;
        float dz = (float)dxyz_box(3, i1, j1, k1) - dz00;
        float dd = powf(dx, 2) + powf(dy, 2) + powf(dz, 2);
        if (dd < Rbox2)
        {
            float d = sqrtf(dd);
            float yy = d * fact11;
            int ir = yy;
            float x = yy - ir;
            float f1 = 1 - x - 0.5 * x * (1 - x);
            float f2 = x + x * (1 - x);
            float f3 = -0.5 * x * (1 - x);
            ir += 1;
            box(i1, j1, k1, ia) = (float)funcr2(ir, itype, itype_mol) * f1 + (float)funcr2(ir + 1, itype, itype_mol) * f2 + (float)funcr2(ir + 2, itype, itype_mol) * f3;
            // box2(i1, j1, k1, ia) = exp(-pow((d / 1.5), 2));
            if (itype_ion != 0) 
            {
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x);
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1 - x);
                ir = ir + 1;
                box3(i1, j1, k1, ia) = (float)rho_ion(ir, itype_ion) * f1 + (float)rho_ion(ir + 1, itype_ion) * f2 + (float)rho_ion(ir + 2, itype_ion) * f3;
            }
        }
        // partial_sum1 += box(i1, j1, k1, ia);
        // partial_sum2 += box2(i1, j1, k1, ia);
        // partial_sum3 += box3(i1, j1, k1, ia);
        loop_id += loop_stride;
    }
    // int stride = gridDim.x * gridDim.y * gridDim.z;
    // int wrt_idx = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum1, tidx, &block_reduce_sum[wrt_idx]);
    // __syncthreads();
    // wrt_idx += stride;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum2, tidx, &block_reduce_sum[wrt_idx]); 
    // __syncthreads(); 
    // wrt_idx += stride;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum3, tidx, &block_reduce_sum[wrt_idx]);  
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_149to238_step_2(double *input, double *output, size_t reduce_size)
{
    int bidx = blockIdx.x;
    int bidy = blockIdx.y;
    int tidx = threadIdx.x;
    int stride = reduce_size * (bidx + bidy * gridDim.x);

    double *src_val = input + stride;

    __shared__ double reduce_shared_mem[blocksize];
    double tmp = 0;
    int iter = tidx;
    while (iter < reduce_size) 
    {
        tmp += src_val[iter];
        iter += blockDim.x;
    }
    reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, src_val);
}

__global__ void cudakernel_loop_149to238_step_3(double *xatom_cent, double *z_cent, int *itype_cent, double *funcr2, double *AL,int nr, double Rm2, 
                                         double *dxyz_box, double *rho, double *rho_z, int id3, int id2, int id1, int n3, int n2, int n1, 
                                         double Rbox2, int itype_mol, double pi, int* ion_type_cent, int* imax_ion, double *r_ion, double *rho_ion, 
                                         double vol, double* box, double *box2, double *box3, double *Q_type, double *z_ion, int imax_nr, int imax_ntype_cent,
                                         double *sum_res, int reduce_block_size)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;

    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    // *sum1 *= vol / (n1 * n2 * n3);
    // *sum2 *= vol / (n1 * n2 * n3);
    // *sum3 *= vol / (n1 * n2 * n3);
    double sum1 = sum_res[reduce_block_size * blockIdx.z + 0 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);
    double sum2 = sum_res[reduce_block_size * blockIdx.z + 1 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);
    double sum3 = sum_res[reduce_block_size * blockIdx.z + 2 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);

    // double fact = (Q_type(itype, itype_mol) - sum1) / sum2;
    double fact = 0;
    double fact2 = 0;
    if (itype_ion != 0)
    {
        // fact2 = z_ion(itype_ion) / sum3;
        fact2 = 1.0;
    }
    
    while (loop_id <= iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;
        int k = (k0 + k1 - 1 + n3) % n3 + 1;
        int j = (j0 + j1 - 1 + n2) % n2 + 1;
        int i = (i0 + i1 - 1 + n1) % n1 + 1;   
        atomicAdd(&rho(i, j, k), box(i1, j1, k1, ia));
        atomicAdd(&rho_z(i, j, k), (fact2 * box3(i1, j1, k1, ia)));
        loop_id += loop_stride;
    }
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_400(int *itype_cent, double *xatom_cent, double *z_cent, int *icorner, double *dbox_c, double *dbox3_c, double *dxyz_box, double *AL, 
                double *rho_m, double *rho_mz, double *funcr2, int nr, double Rm2, int id3, int id2, int id1, int nm3, int nm2, int nm1,  int n3, int n2, int n1, double Rbox2, 
                int itype_mol, double pi, int *ion_type_cent, int *imax_ion, double *r_ion, double *rho_ion, double *box,  double *box2, double *box3, double vol, 
                double *Q_type, double *z_ion, double *dbox, double *dbox2, double *dbox3, int imax_nr, int imax_ntype_cent)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id3 + 1) * (2 * id2 + 1) * (2 * id1 + 1);
    int ia = bidx + 1;
    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;
    double dx00 = AL(1, 1) * dx10 + AL(1, 2) * dx20 + AL(1, 3) * dx30;
    double dy00 = AL(2, 1) * dx10 + AL(2, 2) * dx20 + AL(2, 3) * dx30;
    double dz00 = AL(3, 1) * dx10 + AL(3, 2) * dx20 + AL(3, 3) * dx30;

    double fact11 = nr / Rm2;
    double fact22 = 0;
    double fact = 0;
    double fact2 = 0;  
    if (itype_ion != 0)
    {
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }

    __shared__ double shared_mem[3 * blocksize];
    int iter = tidx;
    double partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00;
        double dd = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
        // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
        //     printf("debugg dd : %lf \r\n", dd);
        //     printf("Rbox2 : %lf \r\n", Rbox2);
        // }
        // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
        //     printf("debugg box_value : %10.5E \r\n", box(i1, j1, k1, ia));
        // }
        if (dd < Rbox2)
        {
            double d = sqrt(dd);
            double yy = d * fact11;
            int ir = yy;
            double x = yy - ir;
            double f1 = 1 - x - 0.5 * x * (1 - x);
            double f2 = x + x * (1 - x);
            double f3 = -0.5 * x * (1 - x);
            double f11 = (x - 1.5) * fact11;
            double f22 = (2 - 2 * x) * fact11;
            double f33 = (x - 0.5) * fact11;
            ir += 1;   
            box(i1, j1, k1, ia) = funcr2(ir, itype, itype_mol) * f1 + funcr2(ir + 1, itype, itype_mol) * f2 + funcr2(ir + 2, itype, itype_mol) * f3;
            // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
            //     printf("debugg box_value : %10.5E \r\n", box(i1, j1, k1, ia));
            // }
            box2(i1, j1, k1, ia) = exp(-pow((d / 1.5), 2));
            double df = 0;
            double df2 = 0;
            double df3 = 0;
            if (d > 1.0e-10)
            {
                df = (funcr2(ir, itype, itype_mol) * f11 + funcr2(ir + 1, itype, itype_mol) * f22 + funcr2(ir + 2, itype, itype_mol) * f33) / d;
            }
            df2 = -2.0 / pow(1.5, 2) * exp(-pow((d / 1.5), 2));
            if (itype_ion != 0)
            {   
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x); 
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1 - x);
                ir += 1;
                f11 = (x - 1.5) * fact22;
                f22 = (2 - 2 * x) * fact22;
                f33 = (x - 0.5) * fact22;
                box3(i1, j1, k1, ia) = rho_ion(ir, itype_ion) * f1 + rho_ion(ir + 1, itype_ion) * f2 + rho_ion(ir + 2, itype_ion) * f3;  
                if (d > 1.0e-10)
                {
                    df3=(rho_ion(ir, itype_ion) * f11 + rho_ion(ir + 1, itype_ion) * f22 + rho_ion(ir + 2, itype_ion) * f33) / d;
                }            
            }
            dbox(i1, j1, k1, 1, ia) = -df * dx;
            dbox(i1, j1, k1, 2, ia) = -df * dy;
            dbox(i1, j1, k1, 3, ia) = -df * dz;
            dbox2(i1, j1, k1, 1, ia) = -df2 * dx;
            dbox2(i1, j1, k1, 2, ia) = -df2 * dy;
            dbox2(i1, j1, k1, 3, ia) = -df2 * dz;
            dbox3(i1, j1, k1, 1, ia) = -df3 * dx;
            dbox3(i1, j1, k1, 2, ia) = -df3 * dy;
            dbox3(i1, j1, k1, 3, ia) = -df3 * dz;          
        }     
        partial_sum1 += box(i1, j1, k1, ia);
        partial_sum2 += box2(i1, j1, k1, ia);
        partial_sum3 += box3(i1, j1, k1, ia);           
        iter += blockDim.x;
    }
    __syncthreads();
    reduce_sharedmem<double, blocksize>(shared_mem + 0 * blocksize, partial_sum1, tidx, shared_mem + 0 * blocksize);
    reduce_sharedmem<double, blocksize>(shared_mem + 1 * blocksize, partial_sum2, tidx, shared_mem + 1 * blocksize);
    reduce_sharedmem<double, blocksize>(shared_mem + 2 * blocksize, partial_sum3, tidx, shared_mem + 2 * blocksize);
    if (tidx == 0)
    {
        shared_mem[0 * blocksize] = shared_mem[0 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[1 * blocksize] = shared_mem[1 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[2 * blocksize] = shared_mem[2 * blocksize] * vol / (n1 * n2 * n3);
        // if (bidx == 0) {
        //     printf("sum1: %lf \r\n", shared_mem[0 * blocksize]);
        //     printf("sum2: %lf \r\n", shared_mem[1 * blocksize]);
        //     printf("sum3: %lf \r\n", shared_mem[2 * blocksize]);
        // }
    }
    __syncthreads();
    fact = (Q_type(itype, itype_mol) - shared_mem[0 + 0 * blocksize]) / shared_mem[0 + 1 * blocksize];
    if (itype_ion != 0)
    {
        fact2 = z_ion(itype_ion) / shared_mem[2 * blocksize];
    }
    iter = tidx;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1; 
        dbox_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact * dbox2(i1, j1, k1, 1, ia);
        dbox_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact * dbox2(i1, j1, k1, 2, ia);
        dbox_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact * dbox2(i1, j1, k1, 3, ia);
        dbox3_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact * dbox2(i1, j1, k1, 1, ia) - fact2 * dbox3(i1, j1, k1, 1, ia);
        dbox3_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact * dbox2(i1, j1, k1, 2, ia) - fact2 * dbox3(i1, j1, k1, 2, ia);
        dbox3_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact * dbox2(i1, j1, k1, 3, ia) - fact2 * dbox3(i1, j1, k1, 3, ia);  
        box(i1, j1, k1, ia) = box(i1, j1, k1, ia) + fact * box2(i1, j1, k1, ia);
        box3(i1, j1, k1, ia) = box(i1, j1, k1, ia) - fact2 * box3(i1, j1, k1, ia); 
        atomicAdd(&rho_m(i, j, k), box(i1, j1, k1, ia));
        atomicAdd(&rho_mz(i, j, k), box3(i1, j1, k1, ia)); 
        iter += blockDim.x;
    }
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_400_step_1(int *itype_cent, double *xatom_cent, double *z_cent, int *icorner, double *dbox_c, double *dbox3_c, double *dxyz_box, double *AL, 
                double *rho_m, double *rho_mz, double *funcr2, int nr, double Rm2, int id3, int id2, int id1, int nm3, int nm2, int nm1,  int n3, int n2, int n1, double Rbox2, 
                int itype_mol, double pi, int *ion_type_cent, int *imax_ion, double *r_ion, double *rho_ion, double *box,  double *box2, double *box3, double vol, 
                double *Q_type, double *z_ion, double *dbox, double *dbox2, double *dbox3, int imax_nr, int imax_ntype_cent, double *block_reduce_sum)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int loop_stride = blockDim.x * gridDim.x;
    int tidx = threadIdx.x;
    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);
    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    float dx10 = x1 - (i0 - 1.0) / n1;
    float dx20 = x2 - (j0 - 1.0) / n2;
    float dx30 = x3 - (k0 - 1.0) / n3;
    float dx00 = (float)AL(1, 1) * dx10 + (float)AL(1, 2) * dx20 + (float)AL(1, 3) * dx30;
    float dy00 = (float)AL(2, 1) * dx10 + (float)AL(2, 2) * dx20 + (float)AL(2, 3) * dx30;
    float dz00 = (float)AL(3, 1) * dx10 + (float)AL(3, 2) * dx20 + (float)AL(3, 3) * dx30;

    float fact11 = nr / Rm2;
    float fact22 = 0;
    // double fact = 0;
    // double fact2 = 0;  
    float partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;

    __shared__ float shared_mem[blocksize];
    
    if (itype_ion != 0)
    {
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }
    while (loop_id < iteration_counts) 
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        box(i1, j1, k1, ia) = 0;
        box2(i1, j1, k1, ia) = 0;
        box3(i1, j1, k1, ia) = 0;
        dbox(i1, j1, k1, 1, ia) = 0;
        dbox(i1, j1, k1, 2, ia) = 0;
        dbox(i1, j1, k1, 3, ia) = 0;
        dbox2(i1, j1, k1, 1, ia) = 0;
        dbox2(i1, j1, k1, 2, ia) = 0;
        dbox2(i1, j1, k1, 3, ia) = 0;
        dbox3(i1, j1, k1, 1, ia) = 0;
        dbox3(i1, j1, k1, 2, ia) = 0;
        dbox3(i1, j1, k1, 3, ia) = 0;

        float dx = (float)dxyz_box(1, i1, j1, k1) - dx00;
        float dy = (float)dxyz_box(2, i1, j1, k1) - dy00;
        float dz = (float)dxyz_box(3, i1, j1, k1) - dz00;
        float dd = powf(dx, 2) + powf(dy, 2) + powf(dz, 2);  
        // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
        //     printf("debugg dd : %lf \r\n", dd);
        // }
        // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
        //     printf("debugg box_value : %10.5E\r\n", box(i1, j1, k1, ia));
        // }
        if (dd < Rbox2)
        {
            float d = sqrtf(dd);
            float yy = d * fact11;
            int ir = yy;
            float x = yy - ir;
            float f1 = 1 - x - 0.5 * x * (1 - x);
            float f2 = x + x * (1 - x);
            float f3 = -0.5 * x * (1 - x);
            float f11 = (x - 1.5) * fact11;
            float f22 = (2 - 2 * x) * fact11;
            float f33 = (x - 0.5) * fact11;
            ir += 1;   
            box(i1, j1, k1, ia) = (float)funcr2(ir, itype, itype_mol) * f1 + (float)funcr2(ir + 1, itype, itype_mol) * f2 + (float)funcr2(ir + 2, itype, itype_mol) * f3;
            // if (i1 == -3 && j1 == 31 && k1 == -23 && ia == 1) {
            //     printf("debugg box_value : %10.5E \r\n", box(i1, j1, k1, ia));
            // }
            // box2(i1, j1, k1, ia) = exp(-pow((d / 1.5), 2));
            float df = 0;
            float df2 = 0;
            float df3 = 0;
            if (d > 1.0e-10)
            {
                df = ((float)funcr2(ir, itype, itype_mol) * f11 + (float)funcr2(ir + 1, itype, itype_mol) * f22 + (float)funcr2(ir + 2, itype, itype_mol) * f33) / d;
            }
            // df2 = -2.0 / pow(1.5, 2) * exp(-pow((d / 1.5), 2));
            if (itype_ion != 0)
            {   
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x); 
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1 - x);
                ir += 1;
                f11 = (x - 1.5) * fact22;
                f22 = (2 - 2 * x) * fact22;
                f33 = (x - 0.5) * fact22;
                box3(i1, j1, k1, ia) = (float)rho_ion(ir, itype_ion) * f1 + (float)rho_ion(ir + 1, itype_ion) * f2 + (float)rho_ion(ir + 2, itype_ion) * f3;  
                if (d > 1.0e-10)
                {
                    df3=((float)rho_ion(ir, itype_ion) * f11 + (float)rho_ion(ir + 1, itype_ion) * f22 + (float)rho_ion(ir + 2, itype_ion) * f33) / d;
                }            
            }
            dbox(i1, j1, k1, 1, ia) = -df * dx;
            dbox(i1, j1, k1, 2, ia) = -df * dy;
            dbox(i1, j1, k1, 3, ia) = -df * dz;
            // dbox2(i1, j1, k1, 1, ia) = -df2 * dx;
            // dbox2(i1, j1, k1, 2, ia) = -df2 * dy;
            // dbox2(i1, j1, k1, 3, ia) = -df2 * dz;
            dbox3(i1, j1, k1, 1, ia) = -df3 * dx;
            dbox3(i1, j1, k1, 2, ia) = -df3 * dy;
            dbox3(i1, j1, k1, 3, ia) = -df3 * dz;          
        }  
        // partial_sum1 += box(i1, j1, k1, ia);
        // partial_sum2 += box2(i1, j1, k1, ia);
        // partial_sum3 += box3(i1, j1, k1, ia);
        loop_id += loop_stride;
    }
    // int stride = gridDim.x * gridDim.y * gridDim.z;
    // int wrt_idx = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum1, tidx, &block_reduce_sum[wrt_idx]);
    // __syncthreads();
    // wrt_idx += stride;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum2, tidx, &block_reduce_sum[wrt_idx]); 
    // __syncthreads(); 
    // wrt_idx += stride;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum3, tidx, &block_reduce_sum[wrt_idx]);  
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_400_step_2(double *input, double *output, size_t reduce_size)
{
    int bidx = blockIdx.x;
    int bidy = blockIdx.y;
    int tidx = threadIdx.x;
    int stride = reduce_size * (bidx + bidy * gridDim.x);

    double *src_val = input + stride;

    __shared__ double reduce_shared_mem[blocksize];
    double tmp = 0;
    int iter = tidx;
    while (iter < reduce_size) 
    {
        tmp += src_val[iter];
        iter += blockDim.x;
    }
    reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, src_val);
}

__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_400_step_3(int *itype_cent, int itype_mol, double *Q_type, double *sum_res, double *dbox_c, double *dbox3_c, double *dbox, double *dbox2, double *dbox3, 
                        double *box, double *box2, double *box3, int *ion_type_cent, int nm1, int nm2, int nm3, double *xatom_cent, int n1, int n2, int n3, int id1, 
                        int id2, int id3, int *icorner, int imax_ntype_cent, double *z_ion, double vol, int reduce_block_size, double *rho_m, double *rho_mz)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int loop_stride = blockDim.x * gridDim.x;

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);

    // *sum1 *= vol / (n1 * n2 * n3);
    // *sum2 *= vol / (n1 * n2 * n3);
    // *sum3 *= vol / (n1 * n2 * n3);
    // double sum1 = sum_res[reduce_block_size * blockIdx.z + 0 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);
    // double sum2 = sum_res[reduce_block_size * blockIdx.z + 1 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);
    // double sum3 = sum_res[reduce_block_size * blockIdx.z + 2 * (reduce_block_size * gridDim.z)] * vol / (n1 * n2 * n3);

    double fact = 0;
    // double fact = (Q_type(itype, itype_mol) - sum1) / sum2;
    double fact2 = 0;
    if (itype_ion != 0)
    {
        fact2 = 1.0;
        // fact2 = z_ion(itype_ion) / sum3;
    }
    
    while (loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1; 
        dbox_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) ;
        dbox_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) ;
        dbox_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) ;
        dbox3_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact2 * dbox3(i1, j1, k1, 1, ia);
        dbox3_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact2 * dbox3(i1, j1, k1, 2, ia);
        dbox3_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact2 * dbox3(i1, j1, k1, 3, ia);  
        box(i1, j1, k1, ia) = box(i1, j1, k1, ia);
        box3(i1, j1, k1, ia) = box(i1, j1, k1, ia) + fact2 * box3(i1, j1, k1, ia); 
        // int i = (im0 + i1 - 1 + nm1) % nm1 + 1; 
        // dbox_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact * dbox2(i1, j1, k1, 1, ia);
        // dbox_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact * dbox2(i1, j1, k1, 2, ia);
        // dbox_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact * dbox2(i1, j1, k1, 3, ia);
        // dbox3_c(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact * dbox2(i1, j1, k1, 1, ia) - fact2 * dbox3(i1, j1, k1, 1, ia);
        // dbox3_c(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact * dbox2(i1, j1, k1, 2, ia) - fact2 * dbox3(i1, j1, k1, 2, ia);
        // dbox3_c(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact * dbox2(i1, j1, k1, 3, ia) - fact2 * dbox3(i1, j1, k1, 3, ia);  
        // box(i1, j1, k1, ia) = box(i1, j1, k1, ia) + fact * box2(i1, j1, k1, ia);
        // box3(i1, j1, k1, ia) = box(i1, j1, k1, ia) - fact2 * box3(i1, j1, k1, ia); 
        atomicAdd(&rho_m(i, j, k), box(i1, j1, k1, ia));
        atomicAdd(&rho_mz(i, j, k), box3(i1, j1, k1, ia));
        loop_id += loop_stride;
    }
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_500(int *itype_cent, double *xatom_cent, double *force_cent, int *icorner, double *vxc, double *vxcm, double *vcoul, double *vcoulm, double *dbox_c,
                                    double *dbox3_c, double vol_n, int id3, int id2, int id1, int n3, int n2, int n1, int nm3, int nm2, int nm1, int itype_mol)
{
    __shared__ double partial_F1[blocksize], partial_F2[blocksize], partial_F3[blocksize];

    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    double partial_F1_reg = 0.0, partial_F2_reg = 0.0, partial_F3_reg = 0.0;

    int ia = bidx + 1;
    int iteration_counts = (2 * id3 + 1) * (2 * id2 + 1) * (2 * id1 + 1);

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);    

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    int iter = tidx;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;
        double dvxc = vxc(i2, j2, k2) - vxcm(i, j, k);
        double dvcoul = vcoul(i2, j2, k2) - vcoulm(i, j, k);   
        partial_F1_reg += dvxc * dbox_c(i1, j1, k1, 1, ia) + dvcoul * dbox3_c(i1, j1, k1, 1, ia);
        partial_F2_reg += dvxc * dbox_c(i1, j1, k1, 2, ia) + dvcoul * dbox3_c(i1, j1, k1, 2, ia);
        partial_F3_reg += dvxc * dbox_c(i1, j1, k1, 3, ia) + dvcoul * dbox3_c(i1, j1, k1, 3, ia); 
        // ------debug
        // if (ia == 1) {
        //     if(tidx==32)
        //         printf("tidx : %d \r\n", tidx);
        //     __syncthreads();
        //     if(tidx==32)
        //         printf("tidx iter iteration_counts partial_F1 :%d %d, %d, %10.5E \r\n", tidx, iter, iteration_counts, dvxc * dbox_c(i1, j1, k1, 1, ia) + dvcoul * dbox3_c(i1, j1, k1, 1, ia));
        // }
        // ------debug 
        iter += blockDim.x;
    }
    reduce_sharedmem<double, blocksize>(partial_F1, partial_F1_reg, tidx, partial_F1);
    reduce_sharedmem<double, blocksize>(partial_F2, partial_F2_reg, tidx, partial_F2);
    reduce_sharedmem<double, blocksize>(partial_F3, partial_F3_reg, tidx, partial_F3);
    if (tidx == 0)
    {
        force_cent(1, ia) = partial_F1[0] * vol_n;
        force_cent(2, ia) = partial_F2[0] * vol_n;
        force_cent(3, ia) = partial_F3[0] * vol_n;
    }
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_500_step_1(int *itype_cent, double *xatom_cent, double *force_cent, int *icorner, double *vxc, double *vxcm, double *vcoul, double *vcoulm, double *dbox_c,
                                    double *dbox3_c, double vol_n, int id3, int id2, int id1, int n3, int n2, int n1, int nm3, int nm2, int nm1, int itype_mol, double *block_reduce_sum)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int loop_stride = blockDim.x * gridDim.x;
    int tidx = threadIdx.x;

    __shared__ double partial_F[blocksize];

    float partial_F1_reg = 0.0, partial_F2_reg = 0.0, partial_F3_reg = 0.0;

    double x1 = fmod(xatom_cent(1, ia) + 1.0, 1.0);
    double x2 = fmod(xatom_cent(2, ia) + 1.0, 1.0);
    double x3 = fmod(xatom_cent(3, ia) + 1.0, 1.0);    

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    while (loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;
        // int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        // int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        // int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;
        // double dvxc = vxc(i2, j2, k2) - vxcm(i, j, k);
        // double dvcoul = vcoul(i2, j2, k2) - vcoulm(i, j, k);   

        float dvxc = (float)vxc(i2, j2, k2);
        float dvcoul = (float)vcoul(i2, j2, k2);   

        partial_F1_reg += dvxc * (float)dbox_c(i1, j1, k1, 1, ia) + dvcoul * (float)dbox3_c(i1, j1, k1, 1, ia);
        partial_F2_reg += dvxc * (float)dbox_c(i1, j1, k1, 2, ia) + dvcoul * (float)dbox3_c(i1, j1, k1, 2, ia);
        partial_F3_reg += dvxc * (float)dbox_c(i1, j1, k1, 3, ia) + dvcoul * (float)dbox3_c(i1, j1, k1, 3, ia); 
        loop_id += loop_stride; 
        // ------debug
        // if (blockIdx.z == 0) {
        //     printf("partial_F1 : %10.5E \r\n", dvxc * dbox_c(i1, j1, k1, 1, ia) + dvcoul * dbox3_c(i1, j1, k1, 1, ia));
        // }
        // ------debug
    }
    int stride = gridDim.x * gridDim.y * gridDim.z;
    int wrt_idx = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    reduce_sharedmem<double, blocksize>(partial_F, (double)partial_F1_reg, tidx, &block_reduce_sum[wrt_idx]);
    __syncthreads();
    wrt_idx += stride;
    reduce_sharedmem<double, blocksize>(partial_F, (double)partial_F2_reg, tidx, &block_reduce_sum[wrt_idx]);
    __syncthreads();
    wrt_idx += stride;
    reduce_sharedmem<double, blocksize>(partial_F, (double)partial_F3_reg, tidx, &block_reduce_sum[wrt_idx]);
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_500_step_2(double *input, double *force_cent, size_t reduce_size, double vol_n)
{
    int bidx = blockIdx.x;
    int bidy = blockIdx.y;
    int tidx = threadIdx.x;
    int stride = reduce_size * (bidx + bidy * gridDim.x);

    double *src_val = input + stride;

    __shared__ double reduce_shared_mem[blocksize];
    double tmp = 0;
    int iter = tidx;
    while (iter < reduce_size) 
    {
        tmp += src_val[iter];
        iter += blockDim.x;
    }
    reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &force_cent(bidy + 1, bidx + 1), vol_n);
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_510(int itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int n3, int n2, int n1, int nm3, int nm2, int nm1, int *icorner, double *AL,
                    double vol_n, int id3, int id2, int id1, double *dxyz_box, double Rbox2, double *pxyz, double *dpxyz, double *vxc2, double *vxc2_m, double *vcoul, double *vcoul_m, int imol,
                    int natom_mm, int nmolm)
{
    int ia = blockIdx.x + 1;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id3 + 1) * (2 * id2 + 1) * (2 * id1 + 1);

    int itype_ion_t = ion_type_atomp(ia, itype_mol);
    double a1 = atom_charge_param(1, itype_ion_t, itype_mol);
    double dw1 = atom_charge_param(2, itype_ion_t, itype_mol);
    double a2 = atom_charge_param(3, itype_ion_t, itype_mol);
    double dw2 = atom_charge_param(4, itype_ion_t, itype_mol);

    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);
    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;

    double dx00 = AL(1, 1) * dx10 + AL(1, 2) * dx20 + AL(1, 3) * dx30;
    double dy00 = AL(2, 1) * dx10 + AL(2, 2) * dx20 + AL(2, 3) * dx30;
    double dz00 = AL(3, 1) * dx10 + AL(3, 2) * dx20 + AL(3, 3) * dx30;

    double px1[3] = {0};
    double px2[3] = {0};
    double dpx1[3][3] = {0};
    double dpx2[3][3] = {0};
    double dxt[3] = {0};
    double dwx1[3] = {0};
    double dwx2[3] = {0};
    int iter = tidx;
    // int count = 0;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;

        dxt[0] = dxyz_box(1, i1, j1, k1) - dx00;
        dxt[1] = dxyz_box(2, i1, j1, k1) - dy00;
        dxt[2] = dxyz_box(3, i1, j1, k1) - dz00;
        double d = sqrt(pow(dxt[0], 2) + pow(dxt[1], 2) + pow(dxt[2], 2));
        if (pow(d, 2) < Rbox2)
        {
            double w1 = exp(-pow(d, 2) / pow(0.5, 2));
            double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            dwx1[0] = -2 * dxt[0] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[0] / pow(d_tmp, 3);
            dwx1[1] = -2 * dxt[1] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[1] / pow(d_tmp, 3);
            dwx1[2] = -2 * dxt[2] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[2] / pow(d_tmp, 3);
            dwx2[0] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[0] / d_tmp - w2 * dxt[0] / pow(d_tmp, 3);
            dwx2[1] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[1] / d_tmp - w2 * dxt[1] / pow(d_tmp, 3);
            dwx2[2] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[2] / d_tmp - w2 * dxt[2] / pow(d_tmp, 3);
            double dvxc2 = vxc2(i2, j2, k2) - vxc2_m(i, j, k);
            double dvcoul = vcoul(i2, j2, k2) - vcoul_m(i, j, k);
            double dv = dvcoul + dvxc2;
            for (int it = 0; it < 3; it++) {
                px1[it] += w1 * dxt[it] / d_tmp * dv;
                px2[it] += w2 * dxt[it] / d_tmp * dv;
            }
            for (int it1 = 0; it1 < 3; it1++) {
                for (int it2 = 0; it2 < 3; it2++) {
                    dpx1[it1][it2] -= dwx1[it2] * dxt[it1] * dv;
                    dpx2[it1][it2] -= dwx2[it2] * dxt[it1] * dv;
                }
            }
            for (int it = 0; it < 3; it++) {
                dpx1[it][it] -= w1 / d_tmp * dv;
                dpx2[it][it] -= w2 / d_tmp * dv;
            }
        }
        iter += blockDim.x;
    }
    __shared__ double shared_mem[blocksize];
    for (int it = 0; it < 3; it++) {
        reduce_sharedmem<double, blocksize>(shared_mem, px1[it], tidx, shared_mem);
        if (tidx == 0) {
            shared_mem[0] = shared_mem[0] * vol_n;
            pxyz[it + (2 * ia - 1 - 1) * 3] = shared_mem[0];
        }
    }
    for (int it1 = 0; it1 < 3; it1++) {
        for (int it2 = 0; it2 < 3; it2++) {
            reduce_sharedmem<double, blocksize>(shared_mem, dpx1[it1][it2], tidx, shared_mem);
            if (tidx == 0) {
                shared_mem[0] = shared_mem[0] * vol_n;
                dpxyz[it1 + it2 * 3 + (2 * ia - 1 - 1) * 3 * 3] = shared_mem[0];
            }
        }
    }
    for (int it = 0; it < 3; it++) {
        reduce_sharedmem<double, blocksize>(shared_mem, px2[it], tidx, shared_mem);
        if (tidx == 0) {
            shared_mem[0] = shared_mem[0] * vol_n;
            pxyz[it + (2 * ia - 1) * 3] = shared_mem[0];
        }
    }
    for (int it1 = 0; it1 < 3; it1++) {
        for (int it2 = 0; it2 < 3; it2++) {
            reduce_sharedmem<double, blocksize>(shared_mem, dpx2[it1][it2], tidx, shared_mem);
            if (tidx == 0) {
                shared_mem[0] = shared_mem[0] * vol_n;
                dpxyz[it1 + it2 * 3 + (2 * ia - 1) * 3 * 3] = shared_mem[0];
            }
        }
    }
}

template<unsigned int blocksize> 
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_510_step_1(int itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int n3, int n2, int n1, int nm3, int nm2, int nm1, int *icorner, double *AL,
                    double vol_n, int id3, int id2, int id1, double *dxyz_box, double Rbox2, double *pxyz, double *dpxyz, double *vxc2, double *vxc2_m, double *vcoul, double *vcoul_m, int imol,
                    int natom_mm, int nmolm, double *block_reduce_sum)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int tidx = threadIdx.x;

    int itype_ion_t = ion_type_atomp(ia, itype_mol);
    float a1 = (float)atom_charge_param(1, itype_ion_t, itype_mol);
    float dw1 = (float)atom_charge_param(2, itype_ion_t, itype_mol);
    float a2 = (float)atom_charge_param(3, itype_ion_t, itype_mol);
    float dw2 = (float)atom_charge_param(4, itype_ion_t, itype_mol);

    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);
    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0 - icorner(1) + 1;
    int jm0 = j0 - icorner(2) + 1;
    int km0 = k0 - icorner(3) + 1;
    im0 = (im0 + n1) % n1;   
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    float dx10 = x1 - (i0 - 1.0) / n1;
    float dx20 = x2 - (j0 - 1.0) / n2;
    float dx30 = x3 - (k0 - 1.0) / n3;

    float dx00 = (float)AL(1, 1) * dx10 + (float)AL(1, 2) * dx20 + (float)AL(1, 3) * dx30;
    float dy00 = (float)AL(2, 1) * dx10 + (float)AL(2, 2) * dx20 + (float)AL(2, 3) * dx30;
    float dz00 = (float)AL(3, 1) * dx10 + (float)AL(3, 2) * dx20 + (float)AL(3, 3) * dx30;

    float px1[3] = {0};
    float px2[3] = {0};
    float dpx1[3][3] = {0};
    float dpx2[3][3] = {0};
    float dxt[3] = {0};
    float dwx1[3] = {0};
    float dwx2[3] = {0};

    float wfactor = 1.5;
    float dw;

    dw = 0.5*wfactor;
    dw1 = dw1*wfactor;
    dw2 = dw2*wfactor;

    float dwpow2 = dw*dw, dw2pow2 = dw2*dw2;
    // int count = 0;
    while (loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        // int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        // int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        // int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;

        dxt[0] = (float)dxyz_box(1, i1, j1, k1) - dx00;
        dxt[1] = (float)dxyz_box(2, i1, j1, k1) - dy00;
        dxt[2] = (float)dxyz_box(3, i1, j1, k1) - dz00;
        float d = sqrtf(dxt[0]*dxt[0] + dxt[1]*dxt[1] + dxt[2]*dxt[2]);
        float dpow2 = d*d;
        if (dpow2 < (float)Rbox2)
        {
            float w1 = expf(-dpow2 / dwpow2);
            float d_tmp = sqrtf(dpow2 + 0.09);
            float w2 = a1 * expf(-d_tmp / dw1) * (1.0 - a2 * expf(-dpow2/dw2pow2));
            float d_tmp2 = d_tmp*d_tmp;
            float tmppara1 = w1 / d_tmp * (-2 / dwpow2 - 1 / d_tmp2);
            dwx1[0] = dxt[0] * tmppara1;
            dwx1[1] = dxt[1] * tmppara1;
            dwx1[2] = dxt[2] * tmppara1;
            float tmppara2 = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / dw2pow2 * expf(-d_tmp / dw1 - dpow2 / dw2pow2) - w2 / d_tmp2) / d_tmp;
            dwx2[0] = dxt[0] * tmppara2;
            dwx2[1] = dxt[1] * tmppara2;
            dwx2[2] = dxt[2] * tmppara2;

            float dvxc2 = (float)vxc2(i2, j2, k2);
            float dvcoul = vcoul(i2, j2, k2);

            // double w1 = exp(-pow(d, 2) / pow(0.5, 2));
            // double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            // double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            // dwx1[0] = -2 * dxt[0] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[0] / pow(d_tmp, 3);
            // dwx1[1] = -2 * dxt[1] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[1] / pow(d_tmp, 3);
            // dwx1[2] = -2 * dxt[2] / pow(0.5, 2) * w1 / d_tmp - w1 * dxt[2] / pow(d_tmp, 3);
            // dwx2[0] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[0] / d_tmp - w2 * dxt[0] / pow(d_tmp, 3);
            // dwx2[1] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[1] / d_tmp - w2 * dxt[1] / pow(d_tmp, 3);
            // dwx2[2] = (-w2 / dw1 / d_tmp + 2 * a1 * a2 / pow(dw2, 2) * exp(-d_tmp / dw1 - pow((d / dw2), 2))) * dxt[2] / d_tmp - w2 * dxt[2] / pow(d_tmp, 3);
            // double dvxc2 = vxc2(i2, j2, k2) - vxc2_m(i, j, k);
            // double dvcoul = vcoul(i2, j2, k2) - vcoul_m(i, j, k);
            // double dvxc2 = vxc2(i2, j2, k2);
            // double dvcoul = vcoul(i2, j2, k2);

            float dv = dvcoul + dvxc2;
#pragma unroll
            for (int it = 0; it < 3; it++) {
                px1[it] += w1 * dxt[it] / d_tmp * dv;
                px2[it] += w2 * dxt[it] / d_tmp * dv;
            }
#pragma unroll
            for (int it1 = 0; it1 < 3; it1++) {
#pragma unroll
                for (int it2 = 0; it2 < 3; it2++) {
                    dpx1[it1][it2] -= dwx1[it2] * dxt[it1] * dv;
                    dpx2[it1][it2] -= dwx2[it2] * dxt[it1] * dv;
                }
            }
#pragma unroll
            for (int it = 0; it < 3; it++) {
                dpx1[it][it] -= w1 / d_tmp * dv;
                dpx2[it][it] -= w2 / d_tmp * dv;
            }
        }
        loop_id += stride;
    }

    __shared__ double shared_mem[blocksize];

    int wrt_idx = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int wrt_stride = gridDim.x * gridDim.y * gridDim.z;

#pragma unroll
    for (int it = 0; it < 3; it++) {
        reduce_sharedmem<double, blocksize>(shared_mem, (double)px1[it], tidx, &block_reduce_sum[wrt_idx]);
        __syncthreads();
        wrt_idx += wrt_stride;
    }

#pragma unroll
    for (int it1 = 0; it1 < 3; it1++) {
#pragma unroll
        for (int it2 = 0; it2 < 3; it2++) {
            reduce_sharedmem<double, blocksize>(shared_mem, (double)dpx1[it1][it2], tidx, &block_reduce_sum[wrt_idx]);
            __syncthreads();
            wrt_idx += wrt_stride;            
        }
    }

#pragma unroll
    for (int it = 0; it < 3; it++) {
        reduce_sharedmem<double, blocksize>(shared_mem, (double)px2[it], tidx, &block_reduce_sum[wrt_idx]);
        __syncthreads();
        wrt_idx += wrt_stride; 
    }

#pragma unroll
    for (int it1 = 0; it1 < 3; it1++) {
#pragma unroll
        for (int it2 = 0; it2 < 3; it2++) {
            reduce_sharedmem<double, blocksize>(shared_mem, (double)dpx2[it1][it2], tidx, &block_reduce_sum[wrt_idx]);
            __syncthreads();
            wrt_idx += wrt_stride; 
        }
    }
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_510_step_2(double *input, double *pxyz, double *dpxyz, size_t reduce_size, double vol_n)
{
    int bidx = blockIdx.x;
    int bidy = blockIdx.y;
    int tidx = threadIdx.x;
    int stride = reduce_size * (bidx + bidy * gridDim.x);

    double *src_val = input + stride;

    __shared__ double reduce_shared_mem[blocksize];
    double tmp = 0;
    int iter = tidx;
    while (iter < reduce_size) 
    {
        tmp += src_val[iter];
        iter += blockDim.x;
    }
    if (bidy < 3) 
    {
        int it = bidy % 3;
        reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &pxyz[it + (2 * (bidx + 1) - 1 - 1) * 3], vol_n);
    }
    else if (bidy >= 3 && bidy < 12)
    {
        int it1 = (bidy - 3) / 3;
        int it2 = (bidy - 3) % 3;
        reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &dpxyz[it1 + it2 * 3 + (2 * (bidx + 1) - 1 - 1) * 3 * 3], vol_n);
    }
    else if (bidy >= 12 && bidy < 15)
    {
        int it = (bidy - 12) % 3;
        reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &pxyz[it + (2 * (bidx + 1) - 1) * 3], vol_n);
    }
    else if (bidy >= 15 && bidy < 24)
    {
        int it1 = (bidy - 15) / 3;
        int it2 = (bidy - 15) % 3;
        reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &dpxyz[it1 + it2 * 3 + (2 * (bidx + 1) - 1) * 3 * 3], vol_n);
    }
}

void fcc_nonbond_loop_82to96(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1)
{
    dim3 gridDim(2 * id3 + 1, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_82to96<<<gridDim, blockDim, 0, 0>>>(dxyz_box, AL, id3, id2, id1, n3, n2, n1);
    // gpuErrchk(cudaPeekAtLastError());
    // gpuErrchk(cudaDeviceSynchronize());
}

void fcc_nonbond_loop_134to149(int *icent, double *w_cent, int *nat_cent, double *xatom_m, int itype_mol, int imol, double *xatom_cent, int ncent_itype_mol, int natom_mm, int nmolm)
{
    dim3 gridDim(1, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_134to149<<<gridDim, blockDim, 0, 0>>>(icent, w_cent, nat_cent, xatom_m, itype_mol, imol, xatom_cent, ncent_itype_mol, natom_mm, nmolm);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_nonbond_loop_149to238(double *xatom_cent, double *z_cent, int *itype_cent, double *funcr2, double *AL, int ncent_itype_mol,
                        int nr, double Rm2, double *dxyz_box, double *rho, double *rho_z, int id3, int id2, int id1,
                        int n3, int n2, int n1, double Rbox2, int itype_mol, double pi, int* ion_type_cent, int* imax_ion,
                        double *r_ion, double *rho_ion, double *box, double *box2, double *box3, double vol, double *Q_type, double *z_ion,
                        int imax_nr, int imax_ntype_cent)
{
    // dim3 gridDim(ncent_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_149to238<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(xatom_cent, z_cent, itype_cent, funcr2, AL, nr, Rm2, dxyz_box, rho, rho_z, id3, id2, id1, 
    //                                                                 n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, r_ion, rho_ion, vol, box, box2, 
    //                                                                 box3, Q_type, z_ion, imax_nr, imax_ntype_cent);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), ncent_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_149to238_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(xatom_cent, z_cent, itype_cent, funcr2, AL, nr, Rm2, dxyz_box, rho, rho_z, id3, id2, id1, 
                                                                    n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, r_ion, rho_ion, vol, box, box2, 
                                                                    box3, Q_type, z_ion, imax_nr, imax_ntype_cent, block_reduce_sum);
    cudakernel_loop_149to238_step_2<1024><<<dim3(ncent_itype_mol, 3, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, block_reduce_sum, reduce_block_size);
    cudakernel_loop_149to238_step_3<<<gridDim, blockDim, 0, 0>>>(xatom_cent, z_cent, itype_cent, funcr2, AL, nr, Rm2, dxyz_box, rho, rho_z, id3, id2, id1, 
                                                                    n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, r_ion, rho_ion, vol, box, box2, 
                                                                    box3, Q_type, z_ion, imax_nr, imax_ntype_cent, block_reduce_sum, reduce_block_size);

                                                                    
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_nonbond_loop_400(int *itype_cent, double *xatom_cent, double *z_cent, int *icorner, double *dbox_c, double *dbox3_c, double *dxyz_box, double *AL, 
                double *rho_m, double *rho_mz, double *funcr2, int nr, double Rm2, int id3, int id2, int id1, int nm3, int nm2, int nm1,  int n3, int n2, int n1, 
                double Rbox2, int itype_mol, int ncent_itype_mol, double pi, int *ion_type_cent, int *imax_ion, double *r_ion, double *rho_ion, double *box, 
                double *box2, double *box3, double vol, double *Q_type, double *z_ion, double *dbox, double *dbox2, double *dbox3, int imax_nr, int imax_ntype_cent)
{
    // dim3 gridDim(ncent_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_400<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, z_cent, icorner, dbox_c, dbox3_c, dxyz_box, AL, rho_m, rho_mz, funcr2, nr, Rm2, id3, id2, id1,
    //                                                         nm3, nm2, nm1, n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, r_ion, rho_ion, box, box2, box3, vol, Q_type, z_ion,
    //                                                         dbox, dbox2, dbox3, imax_nr, imax_ntype_cent);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), ncent_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_400_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, z_cent, icorner, dbox_c, dbox3_c, dxyz_box, AL, rho_m, rho_mz, funcr2, nr, Rm2, id3, id2, id1,
                                                            nm3, nm2, nm1, n3, n2, n1, Rbox2, itype_mol, pi, ion_type_cent, imax_ion, r_ion, rho_ion, box, box2, box3, vol, Q_type, z_ion,
                                                            dbox, dbox2, dbox3, imax_nr, imax_ntype_cent, block_reduce_sum);
    // ------debug
    // int icount = 0;
    // double *box_h = (double *)malloc(sizeof(double) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1));
    // cudaMemcpy(box_h, box, sizeof(double) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1), cudaMemcpyDeviceToHost);
    // for (int i = 0; i < (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1); i++) {
    //     int i1 = i % (2 * id1 + 1) - id1;
    //     int j1 = (i / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
    //     int k1 = i / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
    //     if (box_h[i] != 0) {
    //         printf("i1,j1,k1,box_value : %d %d %d %10.5E \r\n", i1, j1, k1, box_h[i]);
    //         icount++;
    //     }
    // } 
    // cudaDeviceSynchronize();
    // printf("---------- debug \r\n");
    // ------debug
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
    cudakernel_loop_400_step_2<1024><<<dim3(ncent_itype_mol, 3, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, block_reduce_sum, reduce_block_size);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
    // --------debug
    // double *sum = (double *)malloc(sizeof(double) * reduce_block_size * ncent_itype_mol * 3);
    // cudaMemcpy(sum, block_reduce_sum, sizeof(double) * reduce_block_size * ncent_itype_mol * 3, cudaMemcpyDeviceToHost);
    // for (int i = 0; i < reduce_block_size * ncent_itype_mol * 3; i += reduce_block_size * ncent_itype_mol) {
    //     printf("sum: %lf \r\n", sum[i] * vol / (n1 * n2 * n3));
    // }
    // --------debug
    cudakernel_loop_400_step_3<<<gridDim, blockDim, 0, 0>>>(itype_cent, itype_mol, Q_type, block_reduce_sum, dbox_c, dbox3_c, dbox, dbox2, dbox3, box, box2, box3, ion_type_cent, nm1, nm2,
                                                        nm3, xatom_cent, n1, n2, n3, id1, id2, id3, icorner, imax_ntype_cent, z_ion, vol, reduce_block_size, rho_m, rho_mz);

    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
    // debug
    // double *rho_mz_h = (double *)malloc(nm1 * nm2 * nm3 * sizeof(double));
    // cudaMemcpy(rho_mz_h, rho_mz, nm1 * nm2 * nm3 * sizeof(double), cudaMemcpyDeviceToHost);
    // printf("rho_mz_h(1, 41, 53): %20.15E \r\n", rho_mz_h[(1 - 1) + (41 - 1) * nm1 + (53 - 1) * nm1 * nm2]);
    // debug
}

void fcc_nonbond_loop_500(int *itype_cent, double *xatom_cent, double *force_cent, int *icorner, double *vxc, double *vxcm, double *vcoul, double *vcoulm, double *dbox_c,
                double *dbox3_c, double vol_n, int id3, int id2, int id1, int n3, int n2, int n1, int nm3, int nm2, int nm1, int itype_mol, int ncent_itype_mol)
{
    // printf("id1, id2, id3 : %d %d %d \r\n ", id1, id2, id3);
    // printf("loop_counts : %d \r\n ", (2 * id1 + 1) * (2 * id1 + 1) * (2 * id3 + 1));
    // dim3 gridDim(ncent_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_500<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, force_cent, icorner, vxc, vxcm, vcoul, vcoulm, dbox_c, dbox3_c, vol_n, id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), ncent_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_500_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, force_cent, icorner, vxc, vxcm, vcoul, vcoulm, dbox_c, dbox3_c, vol_n, id3, id2, id1, n3, n2, n1, nm3, nm2, nm1, itype_mol, block_reduce_sum);
    // // printf("vol_n : %lf \r\n", vol_n);
    cudakernel_loop_500_step_2<1024><<<dim3(ncent_itype_mol, 3, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, force_cent, reduce_block_size, vol_n);
    // cudaDeviceSynchronize();
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_nonbond_loop_510(int natom_m_itype_mol, int itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int n3, int n2, int n1, int nm3, int nm2, int nm1, int *icorner, double *AL,
                    double vol_n, int id3, int id2, int id1, double *dxyz_box, double Rbox2, double *pxyz, double *dpxyz, double *vxc2, double *vxc2_m, double *vcoul, double *vcoul_m, int imol, int natom_mm,
                    int nmolm)
{
//     const int blocksize = 256;
//     dim3 gridDim(natom_m_itype_mol, 1, 1);
//     dim3 blockDim(blocksize, 1, 1);
//     cudakernel_loop_510<blocksize><<<gridDim, blockDim, 0, 0>>>(itype_mol, ion_type_atomp, atom_charge_param, xatom_m, n3, n2, n1, nm3, nm2, nm1, icorner, AL, vol_n, id3, id2, id1, dxyz_box, Rbox2, pxyz, dpxyz, vxc2, vxc2_m, vcoul, vcoul_m, imol, natom_mm, nmolm);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), natom_m_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_510_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_mol, ion_type_atomp, atom_charge_param, xatom_m, n3, n2, n1, nm3, nm2, nm1, icorner, AL, vol_n, id3, id2, id1, dxyz_box, Rbox2, pxyz, dpxyz, vxc2, vxc2_m, vcoul, vcoul_m, imol, natom_mm, nmolm, block_reduce_sum);
    cudakernel_loop_510_step_2<1024><<<dim3(natom_m_itype_mol, 3 + 9 + 3 + 9, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, pxyz, dpxyz, reduce_block_size, vol_n);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}