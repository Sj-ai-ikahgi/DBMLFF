#include "../include/util.h"
#include "../include/calc_polar_force.h"

#define BLOCKSIZE 512
#define maxThreadsPerBlock BLOCKSIZE
#define minBlocksPerMultiprocessor 2
__global__ void cudakernel_loop_113to129(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    
    int k1 = bidx - id3;
    double dx3 = k1 * 1.0 / n3;
    int iter = tidx;
    while(iter < iteration_counts)
    {
        int j1 = iter / (2 * id1 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        double dx2 = j1 * 1.0 / n2;
        double dx1 = i1 * 1.0 / n1;
        double dx = AL(1,1) * dx1 + AL(1,2) * dx2 + AL(1,3) * dx3;
        double dy = AL(2,1) * dx1 + AL(2,2) * dx2 + AL(2,3) * dx3;
        double dz = AL(3,1) * dx1 + AL(3,2) * dx2 + AL(3,3) * dx3;
        dxyz_box(1, i1, j1, k1) = dx;
        dxyz_box(2, i1, j1, k1) = dy;
        dxyz_box(3, i1, j1, k1) = dz;
        iter += blockDim.x;
    }
}

__global__ void cudakernel_loop_500(int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int imol, int itype_mol, int n3, int n2, int n1, int id3, int id2, int id1, 
                                    int natom_mm, int nmolm, int ntype_mm, double *AL, double *dxyz_box, double *CC_pol, double *rhop_tot, double Rbox2)
{
    int ia = blockIdx.x + 1;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1);

    int ion_type_t = ion_type_atomp(ia, itype_mol);
    double a1 = atom_charge_param(1, ion_type_t, itype_mol);
    double dw1 = atom_charge_param(2, ion_type_t, itype_mol);
    double a2 = atom_charge_param(3, ion_type_t, itype_mol);
    double dw2 = atom_charge_param(4, ion_type_t, itype_mol);
    
    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;

    double dx00 = AL(1,1) * dx10 + AL(1,2) * dx20 + AL(1,3) * dx30;
    double dy00 = AL(2,1) * dx10 + AL(2,2) * dx20 + AL(2,3) * dx30;
    double dz00 = AL(3,1) * dx10 + AL(3,2) * dx20 + AL(3,3) * dx30;

    int iter = tidx;
    while(iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1; 
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;
        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00;
        double d = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        if (pow(d, 2) < Rbox2)
        {
            double w = exp((-pow(d, 2)) / pow(0.5, 2));
            double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            double rhop_tmp = w * dx / d_tmp * CC_pol(1, ia, imol, itype_mol) + 
                                w * dy / d_tmp * CC_pol(2, ia, imol, itype_mol) + 
                                w * dz / d_tmp * CC_pol(3, ia, imol, itype_mol) + 
                                w2 * dx / d_tmp * CC_pol(4, ia, imol, itype_mol) + 
                                w2 * dy / d_tmp * CC_pol(5, ia, imol, itype_mol) + 
                                w2 * dz / d_tmp * CC_pol(6, ia, imol, itype_mol);
            atomicAdd(&rhop_tot(i2, j2, k2), rhop_tmp);
        }
        iter += blockDim.x;
    }
}

__global__ void cudakernel_loop_500_step_1(int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int imol, int itype_mol, int n3, int n2, int n1, int id3, int id2, int id1, 
                                    int natom_mm, int nmolm, int ntype_mm, double *AL, double *dxyz_box, double *CC_pol, double *rhop_tot, double Rbox2)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int loop_stride = blockDim.x * gridDim.x;

    int ion_type_t = ion_type_atomp(ia, itype_mol);
    float a1 = (float)atom_charge_param(1, ion_type_t, itype_mol);
    float dw1 = (float)atom_charge_param(2, ion_type_t, itype_mol);
    float a2 = (float)atom_charge_param(3, ion_type_t, itype_mol);
    float dw2 = (float)atom_charge_param(4, ion_type_t, itype_mol);
    
    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;

    float dx10 = x1 - (i0 - 1.0) / n1;
    float dx20 = x2 - (j0 - 1.0) / n2;
    float dx30 = x3 - (k0 - 1.0) / n3;

    float dx00 = (float)AL(1,1) * dx10 + (float)AL(1,2) * dx20 + (float)AL(1,3) * dx30;
    float dy00 = (float)AL(2,1) * dx10 + (float)AL(2,2) * dx20 + (float)AL(2,3) * dx30;
    float dz00 = (float)AL(3,1) * dx10 + (float)AL(3,2) * dx20 + (float)AL(3,3) * dx30;

    float wfactor = 1.5;
    float dw;

    dw = 0.5*wfactor;
    dw1 = dw1*wfactor;
    dw2 = dw2*wfactor;



    while(loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1; 
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;
        float dx = (float)dxyz_box(1, i1, j1, k1) - dx00;
        float dy = (float)dxyz_box(2, i1, j1, k1) - dy00;
        float dz = (float)dxyz_box(3, i1, j1, k1) - dz00;
        float d = sqrtf(dx*dx + dy*dy + dz*dz);
        float dpow2 = d*d, dwpow2 = dw*dw;
        if (dpow2 < (float)Rbox2)
        {
            float dtdw22 = dpow2/dwpow2;
            float w = expf(-dtdw22);
            float d_tmp = sqrtf(dpow2 + 0.09);
            float w2 = a1 * expf(-d_tmp / dw1) * (1.0 - a2 * expf(-dpow2/(dw2*dw2)));
            float wtod_tmp = w / d_tmp, w2tod_tmp = w2 / d_tmp;
            float rhop_tmp = wtod_tmp * dx * CC_pol(1, ia, imol, itype_mol) + 
                                wtod_tmp * dy * CC_pol(2, ia, imol, itype_mol) + 
                                wtod_tmp * dz * CC_pol(3, ia, imol, itype_mol) + 
                                w2tod_tmp * dx * CC_pol(4, ia, imol, itype_mol) + 
                                w2tod_tmp * dy * CC_pol(5, ia, imol, itype_mol) + 
                                w2tod_tmp * dz * CC_pol(6, ia, imol, itype_mol);
            atomicAdd(&rhop_tot(i2, j2, k2), rhop_tmp);

            // double w = exp((-pow(d, 2)) / pow(0.5, 2));
            // double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            // double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            // double rhop_tmp = w * dx / d_tmp * CC_pol(1, ia, imol, itype_mol) + 
            //                     w * dy / d_tmp * CC_pol(2, ia, imol, itype_mol) + 
            //                     w * dz / d_tmp * CC_pol(3, ia, imol, itype_mol) + 
            //                     w2 * dx / d_tmp * CC_pol(4, ia, imol, itype_mol) + 
            //                     w2 * dy / d_tmp * CC_pol(5, ia, imol, itype_mol) + 
            //                     w2 * dz / d_tmp * CC_pol(6, ia, imol, itype_mol);
            // atomicAdd(&rhop_tot(i2, j2, k2), rhop_tmp);
        }
        loop_id += loop_stride;
    }
}

__global__ void cudakernel_loop312to339(int ncent_itype_mol, int *icent, int imol, int itype_mol, double *w_cent, double *xc_cent, double *xatom_cent, int *icorner, int n3, int n2, int n1, 
                                        int nm3, int nm2, int nm1, double *xatom_m, int natom_mm, int nmolm, int *nat_cent)
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
    __syncthreads();
    iter = tidx;
    while (iter < iteration_counts)
    {   
        int ii = iter / 3 + 1;
        int ixyz = iter % 3 + 1;        
        double dx = xatom_cent(ixyz, ii) - xatom_cent(ixyz, 1);
        if(abs(dx + 1) < abs(dx)) 
            dx = dx + 1;  
        if(abs(dx - 1) < abs(dx)) 
            dx = dx - 1;
        atomicAdd(&xc_cent(ixyz), dx + xatom_cent(ixyz,1));
        iter += blockDim.x;
    }
    __syncthreads();
    if (tidx == 0)
    {
        xc_cent(1) = xc_cent(1) / ncent_itype_mol;
        xc_cent(2) = xc_cent(2) / ncent_itype_mol;
        xc_cent(3) = xc_cent(3) / ncent_itype_mol;
        icorner(1) = fmod(xc_cent(1) + 1.0, 1.0) * n1 - (int)(0.6 + nm1 / 2.0) + 1;
        icorner(2) = fmod(xc_cent(2) + 1.0, 1.0) * n2 - (int)(0.6 + nm2 / 2.0) + 1;
        icorner(3) = fmod(xc_cent(3) + 1.0, 1.0) * n3 - (int)(0.6 + nm3 / 2.0) + 1;
    }
}   

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_600(int *itype_cent, int itype_mol, double *xatom_cent, int *icorner, double *AL, int *ion_type_cent, int *imax_ion, double *r_ion, 
                    int id3, int id2, int id1, double *box, double *box2, double *box3, double *dxyz_box, double *funcr2, double *fact_store, double *fact2_store,
                    int n3, int n2, int n1, int nr, double Rm2, double Rbox2, double *rho_ion, double vol, int nm3, int nm2, int nm1, double *rho_m, double *Q_type,
                    double *z_ion, int imax_nr, int imax_ntype_cent)
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

    if (itype_ion != 0)
    {
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }
    
    __shared__ double shared_mem[3 * blocksize];
    double partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;
    int iter = tidx;
    while(iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;

        // box(i1, j1, k1, ia) = 0;
        // box2(i1, j1, k1, ia) = 0;
        // box3(i1, j1, k1, ia) = 0;

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
    reduce_sharedmem<double ,blocksize>(shared_mem + 1 * blocksize, partial_sum2, tidx, shared_mem + 1 * blocksize);
    reduce_sharedmem<double, blocksize>(shared_mem + 2 * blocksize, partial_sum3, tidx, shared_mem + 2 * blocksize);
    if (tidx == 0)
    {
        shared_mem[0 * blocksize] = shared_mem[0 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[1 * blocksize] = shared_mem[1 * blocksize] * vol / (n1 * n2 * n3);
        shared_mem[2 * blocksize] = shared_mem[2 * blocksize] * vol / (n1 * n2 * n3);
        fact_store[ia - 1] =(Q_type(itype, itype_mol) - shared_mem[0 * blocksize]) / shared_mem[1 * blocksize];
        // if (bidx == 0)
        // {
        //     printf("%20.15E %d\r\n", Q_type[itype - 1], itype);
        // }
        fact2_store[ia - 1] = 0.0;
        if (itype_ion != 0)
        {    
            fact2_store[ia - 1] = z_ion(itype_ion) / shared_mem[2 * blocksize];
        }
    }
    __syncthreads();
    iter = tidx;
    while (iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1; 
        box(i1, j1, k1, ia) += fact_store[ia - 1] * box2(i1, j1, k1, ia);
        atomicAdd(&rho_m(i, j, k), box(i1, j1, k1, ia));      
        // if (bidx == 0)
        //     printf("%20.15E %20.15E \r\n", box(i1, j1, k1, ia), fact_store[ia - 1]);
        iter += blockDim.x;
    }   
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_600_step_1(int *itype_cent, int itype_mol, double *xatom_cent, int *icorner, double *AL, int *ion_type_cent, int *imax_ion, double *r_ion, 
                                    int id3, int id2, int id1, double *box, double *box2, double *box3, double *dxyz_box, double *funcr2, double *fact_store, double *fact2_store,
                                    int n3, int n2, int n1, int nr, double Rm2, double Rbox2, double *rho_ion, double vol, int nm3, int nm2, int nm1, double *rho_m, double *Q_type,
                                    double *z_ion, int imax_nr, int imax_ntype_cent, double *block_reduce_sum)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int itype = itype_cent(ia, itype_mol);
    int itype_ion = ion_type_cent(ia, itype_mol);
    int tidx = threadIdx.x;

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

    if (itype_ion != 0)
    {
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }
    
    __shared__ double shared_mem[blocksize];
    double partial_sum1 = 0, partial_sum2 = 0, partial_sum3 = 0;

    while(loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        box(i1, j1, k1, ia) = 0;

        // box2(i1, j1, k1, ia) = 0;
        box3(i1, j1, k1, ia) = 0;

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
            // line 413
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
                box3(i1, j1, k1, ia) = rho_ion(ir, itype_ion) * f1 + rho_ion(ir + 1, itype_ion) * f2 + rho_ion(ir + 2, itype_ion) * f3;
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
    // reduce_sharedmem<double ,blocksize>(shared_mem, partial_sum2, tidx, &block_reduce_sum[wrt_idx]);
    // __syncthreads();
    // wrt_idx += stride;
    // reduce_sharedmem<double, blocksize>(shared_mem, partial_sum3, tidx, &block_reduce_sum[wrt_idx]); 
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_600_step_2(double *input, double *output, size_t reduce_size, double scale)
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
    reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, src_val, scale);
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_600_step_3(int *itype_cent, int itype_mol, double *xatom_cent, int *icorner, double *AL, int *ion_type_cent, int *imax_ion, double *r_ion, 
                                    int id3, int id2, int id1, double *box, double *box2, double *box3, double *dxyz_box, double *funcr2, double *fact_store, double *fact2_store,
                                    int n3, int n2, int n1, int nr, double Rm2, double Rbox2, double *rho_ion, double vol, int nm3, int nm2, int nm1, double *rho_m, double *Q_type,
                                    double *z_ion, int imax_nr, int imax_ntype_cent, double *sum_res, int reduce_block_size)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int bidx = blockIdx.x;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
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

    // double sum1 = sum_res[reduce_block_size * blockIdx.z + 0 * (reduce_block_size * gridDim.z)];
    // double sum2 = sum_res[reduce_block_size * blockIdx.z + 1 * (reduce_block_size * gridDim.z)];
    // double sum3 = sum_res[reduce_block_size * blockIdx.z + 2 * (reduce_block_size * gridDim.z)];

    // double fact1 =(Q_type(itype, itype_mol) - sum1) / sum2;
    double fact2 = 0;
    if (itype_ion != 0)
    {    fact2 = 1.0;
        // fact2 = z_ion(itype_ion) / sum3;
    }
    // fact_store[ia - 1] = fact1;
    fact2_store[ia - 1] = fact2;

    while (loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1; 
        // box(i1, j1, k1, ia) += fact1 * box2(i1, j1, k1, ia);
        atomicAdd(&rho_m(i, j, k), box(i1, j1, k1, ia));      
        // if (bidx == 0)
        //     // printf("%20.15E %20.15E \r\n", box(i1, j1, k1, ia), fact_store[ia - 1]);
        //     printf("%20.15E %20.15E \r\n", box(i1, j1, k1, ia), fact_store[ia - 1]);
        loop_id += loop_stride;
    }  

}

__global__ void cudakernel_loop_700(int *ion_type_atomp, int imol, int itype_mol, double *atom_charge_param, double *xatom_m, int *icorner, int n3, int n2, int n1, 
                        int id3, int id2, int id1, double *CC_pol, double *rhop_m, double*dxyz_box, double Rbox2, int natom_mm, int nmolm, double *AL, int nm3,
                        int nm2, int nm1)
{   
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1);
    int ia =  bidx + 1;

    int ion_type_t = ion_type_atomp(ia, itype_mol);
    double a1 = atom_charge_param(1, ion_type_t, itype_mol);
    double dw1 = atom_charge_param(2, ion_type_t, itype_mol);
    double a2 = atom_charge_param(3, ion_type_t, itype_mol);
    double dw2 = atom_charge_param(4, ion_type_t, itype_mol);

    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0-icorner(1) + 1; 
    int jm0 = j0-icorner(2) + 1;
    int km0 = k0-icorner(3) + 1;
    im0 = (im0 + n1) % n1;
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;

    double dx00 = AL(1,1) * dx10 + AL(1,2) * dx20 + AL(1,3) * dx30;
    double dy00 = AL(2,1) * dx10 + AL(2,2) * dx20 + AL(2,3) * dx30;
    double dz00 = AL(3,1) * dx10 + AL(3,2) * dx20 + AL(3,3) * dx30;

    int iter = tidx;
    while(iter < iteration_counts)
    {  
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00;
        double d = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));        
        if (pow(d, 2) < Rbox2)
        {
            double w = exp(-pow(d, 2) / pow(0.5, 2));
            double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            double rhop_tmp = w * dx / d_tmp * CC_pol(1, ia, imol, itype_mol) + 
                            w * dy / d_tmp * CC_pol(2, ia, imol, itype_mol) + 
                            w * dz / d_tmp * CC_pol(3, ia, imol, itype_mol) + 
                            w2 * dx / d_tmp * CC_pol(4, ia, imol, itype_mol) + 
                            w2 * dy / d_tmp * CC_pol(5, ia, imol, itype_mol) + 
                            w2 * dz / d_tmp * CC_pol(6, ia, imol, itype_mol);
            atomicAdd(&rhop_m(i, j, k), rhop_tmp); 
        }
        iter += blockDim.x;
    }
}

__global__ void cudakernel_loop_700_step_1(int *ion_type_atomp, int imol, int itype_mol, double *atom_charge_param, double *xatom_m, int *icorner, int n3, int n2, int n1, 
                        int id3, int id2, int id1, double *CC_pol, double *rhop_m, double*dxyz_box, double Rbox2, int natom_mm, int nmolm, double *AL, int nm3,
                        int nm2, int nm1)
{   
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);

    int ion_type_t = ion_type_atomp(ia, itype_mol);
    double a1 = atom_charge_param(1, ion_type_t, itype_mol);
    double dw1 = atom_charge_param(2, ion_type_t, itype_mol);
    double a2 = atom_charge_param(3, ion_type_t, itype_mol);
    double dw2 = atom_charge_param(4, ion_type_t, itype_mol);

    double x1 = fmod(xatom_m(1, ia, imol, itype_mol) + 1.0, 1.0);
    double x2 = fmod(xatom_m(2, ia, imol, itype_mol) + 1.0, 1.0);
    double x3 = fmod(xatom_m(3, ia, imol, itype_mol) + 1.0, 1.0);

    int i0 = x1 * n1 + 1;
    int j0 = x2 * n2 + 1;
    int k0 = x3 * n3 + 1;
    int im0 = i0-icorner(1) + 1; 
    int jm0 = j0-icorner(2) + 1;
    int km0 = k0-icorner(3) + 1;
    im0 = (im0 + n1) % n1;
    jm0 = (jm0 + n2) % n2;
    km0 = (km0 + n3) % n3;

    double dx10 = x1 - (i0 - 1.0) / n1;
    double dx20 = x2 - (j0 - 1.0) / n2;
    double dx30 = x3 - (k0 - 1.0) / n3;

    double dx00 = AL(1,1) * dx10 + AL(1,2) * dx20 + AL(1,3) * dx30;
    double dy00 = AL(2,1) * dx10 + AL(2,2) * dx20 + AL(2,3) * dx30;
    double dz00 = AL(3,1) * dx10 + AL(3,2) * dx20 + AL(3,3) * dx30;

    while(loop_id < iteration_counts)
    {  
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;
        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00;
        double d = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));        
        if (pow(d, 2) < Rbox2)
        {
            double w = exp(-pow(d, 2) / pow(0.5, 2));
            double d_tmp = sqrt(pow(d, 2) + pow(0.3, 2));
            double w2 = a1 * exp(-d_tmp / dw1) * (1.0 - a2 * exp(-pow((d / dw2), 2)));
            double rhop_tmp = w * dx / d_tmp * CC_pol(1, ia, imol, itype_mol) + 
                            w * dy / d_tmp * CC_pol(2, ia, imol, itype_mol) + 
                            w * dz / d_tmp * CC_pol(3, ia, imol, itype_mol) + 
                            w2 * dx / d_tmp * CC_pol(4, ia, imol, itype_mol) + 
                            w2 * dy / d_tmp * CC_pol(5, ia, imol, itype_mol) + 
                            w2 * dz / d_tmp * CC_pol(6, ia, imol, itype_mol);
            atomicAdd(&rhop_m(i, j, k), rhop_tmp); 
        }
        loop_id += loop_stride;
    }
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_800(int * itype_cent, double *xatom_cent, int *icorner, double *AL, int n3, int n2, int n1, int *ion_type_cent, int *imax_ion,  double *r_ion, double *dxyz_box, 
                                    double *funcr2, double *box, double *box2, double *box3, double *dbox, double *dbox2, double *dbox3, double *fact_store, double *fact2_store, double Rbox2, 
                                    int nr, int nm3, int nm2, int nm1, double *vrhop_tot_coul, double *vrhop_m_coul, double *dvxc_tot, double *dvxc_m, double vol_n, int id3, int id2, int id1,
                                    double *force_cent, double Rm2, int itype_mol, double *rho_ion, int imax_nr, int imax_ntype_cent)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1);
    int ia =  bidx + 1;
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

    // double fact2 = 0;
    if (itype_ion != 0)
    {
        // fact2 = 1.0;
        fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    }
    // fact2_store[ia-1] = fact2;

    // if (itype_ion != 0)
    // {
    //     fact22 = (imax_ion(itype_ion) - 1) / r_ion(imax_ion(itype_ion), itype_ion);
    // }

    int iter = tidx;
    double sum2 = 0.0;
    double dsum3_1 = 0.0;
    double dsum3_2 = 0.0;
    double dsum3_3 = 0.0;
    double dsum_1 = 0.0;
    double dsum_2 = 0.0;
    double dsum_3 = 0.0;
    while(iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
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

        double dx = dxyz_box(1, i1, j1, k1) - dx00;
        double dy = dxyz_box(2, i1, j1, k1) - dy00;
        double dz = dxyz_box(3, i1, j1, k1) - dz00; 
        double dd = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);  
        // if (i1 == 0 && j1 == 0 && k1 == 0 && ia == 1)
        // {
        //     printf("%20.15E %20.15E %20.15E \r\n", dx00, AL(1, 1), dx10);
        //     printf("%20.15E %20.15E \r\n",  AL(1, 2), dx20);
        //     printf("%20.15E %20.15E \r\n", AL(1, 3), dx30);
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
            box2(i1, j1, k1, ia) = exp(-pow((d/1.5), 2));
            double df = 0;
            if (d > 1e-10)
            {
                df = (funcr2(ir, itype, itype_mol) * f11 + funcr2(ir + 1, itype, itype_mol) * f22 + funcr2(ir + 2, itype, itype_mol) * f33) / d;
            }
            double df2 = -2 / pow(1.5, 2) * exp(-pow((d/1.5), 2));
            double df3 = 0;
            dbox(i1, j1, k1, 1, ia) = -df * dx;
            // if (ia == 1)
            // {
            //     printf("%20.15E %20.15E %20.15E\r\n", dbox(i1, j1, k1, 1, ia), df, dx);
            // }
            dbox(i1, j1, k1, 2, ia) = -df * dy;
            dbox(i1, j1, k1, 3, ia) = -df * dz;
            dbox2(i1, j1, k1, 1, ia) = -df2 * dx;
            dbox2(i1, j1, k1, 2, ia) = -df2 * dy;
            dbox2(i1, j1, k1, 3, ia) = -df2 * dz;   
            if (itype_ion != 0)
            {
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x);
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1-x);
                f11 = (x - 1.5) * fact22;
                f22 = (2 - 2 * x) * fact22;
                f33 = (x - 0.5) * fact22;
                ir += 1;                
                box3(i1, j1, k1, ia) = rho_ion(ir, itype_ion) * f1 + rho_ion(ir + 1, itype_ion) * f2 + rho_ion(ir + 2, itype_ion) * f3;
                if (d > 1e-10)
                {
                    df3 = (rho_ion(ir, itype_ion) * f11 + rho_ion(ir + 1, itype_ion) * f22 + rho_ion(ir + 2, itype_ion) * f33) / d;
                }
                dbox3(i1, j1, k1, 1, ia) = -df3 * dx;
                // if (ia == 1)
                // {
                //     printf("%20.15E \r\n", yy);
                // }
                dbox3(i1, j1, k1, 2, ia) = -df3 * dy;
                dbox3(i1, j1, k1, 3, ia) = -df3 * dz;
            }        
        }
        // if (ia == 1)
        // {
        //     printf("%20.15E %20.15E %20.15E\r\n", dbox(i1, j1, k1, 1, ia), df, dx);
        // }
        dbox(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) + fact_store[ia - 1] * dbox2(i1, j1, k1, 1, ia);
        dbox(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) + fact_store[ia - 1] * dbox2(i1, j1, k1, 2, ia);
        dbox(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) + fact_store[ia - 1] * dbox2(i1, j1, k1, 3, ia);
        dbox3(i1, j1, k1, 1, ia) = dbox(i1, j1, k1, 1, ia) - fact2_store[ia - 1] * dbox3(i1, j1, k1, 1, ia);
        dbox3(i1, j1, k1, 2, ia) = dbox(i1, j1, k1, 2, ia) - fact2_store[ia - 1] * dbox3(i1, j1, k1, 2, ia);
        dbox3(i1, j1, k1, 3, ia) = dbox(i1, j1, k1, 3, ia) - fact2_store[ia - 1] * dbox3(i1, j1, k1, 3, ia);
        sum2 += box2(i1, j1, k1, ia);
        dsum3_1 += dbox3(i1, j1, k1, 1, ia);
        dsum3_2 += dbox3(i1, j1, k1, 2, ia);
        dsum3_3 += dbox3(i1, j1, k1, 3, ia);
        dsum_1 += dbox(i1, j1, k1, 1, ia);
        dsum_2 += dbox(i1, j1, k1, 2, ia);
        dsum_3 += dbox(i1, j1, k1, 3, ia);
        iter += blockDim.x;
    }

    __shared__ double shared_mem[blocksize];

    __syncthreads();
    reduce_sharedmem<double, blocksize>(shared_mem, sum2, tidx, shared_mem);
    __syncthreads();
    sum2 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum3_1, tidx, shared_mem);
    __syncthreads();
    dsum3_1 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum3_2, tidx, shared_mem);
    __syncthreads();
    dsum3_2 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum3_3, tidx, shared_mem);
    __syncthreads();
    dsum3_3 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum_1, tidx, shared_mem);
    __syncthreads();
    dsum_1 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum_2, tidx, shared_mem);
    __syncthreads();
    dsum_2 = shared_mem[0];
    reduce_sharedmem<double, blocksize>(shared_mem, dsum_3, tidx, shared_mem);
    __syncthreads();

    // --debug
    // if (ia == 1)
    // {
    //     printf("sum2 : %lf \r\n", sum2);
    //     printf("dsum3_1 : %lf \r\n", dsum3_1);
    //     printf("dsum3_2 : %lf \r\n", dsum3_2);
    //     printf("dsum3_3 : %lf \r\n", dsum3_3);
    //     printf("dsum_1 : %lf \r\n", dsum_1);
    //     printf("dsum_2 : %lf \r\n", dsum_2);
    //     printf("dsum_3 : %lf \r\n", dsum_3);
    // }
    // --debug

    dsum_3 = shared_mem[0];
    dsum3_1 = dsum3_1 / sum2;
    dsum3_2 = dsum3_2 / sum2;
    dsum3_3 = dsum3_3 / sum2;
    dsum_1 = dsum_1 / sum2;
    dsum_2 = dsum_2 / sum2;
    dsum_3 = dsum_3 / sum2;

    // if (tidx == 0 && bidx == 0)
    //     printf("%20.15E, %20.15E, %20.15E, %20.15E, %20.15E, %20.15E, %20.15E \r\n", dsum_1, dsum_2, dsum_3, dsum3_1, dsum3_2, dsum3_3, sum2);
    double partial_force_cent1 = 0, partial_force_cent2 = 0, partial_force_cent3 = 0;
    iter = tidx;
    while(iter < iteration_counts)
    {
        int k1 = iter / ((2 * id1 + 1) * (2 * id2 + 1)) - id3;
        int j1 = (iter / (2 * id1 + 1)) % (2 * id2 + 1) - id2;
        int i1 = iter % (2 * id1 + 1) - id1;
        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;  
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;  
        dbox3(i1, j1, k1, 1, ia) += -dsum3_1 * box2(i1, j1, k1, ia);
        dbox3(i1, j1, k1, 2, ia) += -dsum3_2 * box2(i1, j1, k1, ia);
        dbox3(i1, j1, k1, 3, ia) += -dsum3_3 * box2(i1, j1, k1, ia);
        dbox(i1, j1, k1, 1, ia) += -dsum_1 * box2(i1, j1, k1, ia);
        dbox(i1, j1, k1, 2, ia) += -dsum_2 * box2(i1, j1, k1, ia);
        dbox(i1, j1, k1, 3, ia) += -dsum_3 * box2(i1, j1, k1, ia);
        double vtmp1 = vrhop_tot_coul(i2, j2, k2) - vrhop_m_coul(i, j, k);
        double vtmp2 = dvxc_tot(i2, j2, k2) - dvxc_m(i, j, k);
        partial_force_cent1 += vtmp1 * dbox3(i1, j1, k1, 1, ia) + vtmp2 * dbox(i1, j1, k1, 1, ia);
        partial_force_cent2 += vtmp1 * dbox3(i1, j1, k1, 2, ia) + vtmp2 * dbox(i1, j1, k1, 2, ia);
        partial_force_cent3 += vtmp1 * dbox3(i1, j1, k1, 3, ia) + vtmp2 * dbox(i1, j1, k1, 3, ia);
        iter += blockDim.x;
    }
    __syncthreads();
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent1, tidx, shared_mem);
    if (tidx == 0) {
        force_cent(1, ia) = shared_mem[0] * vol_n;
    }
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent2, tidx, shared_mem);
    if (tidx == 0) {
        force_cent(2, ia) = shared_mem[0] * vol_n;
    }   
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent3, tidx, shared_mem);
    if (tidx == 0) {
        force_cent(3, ia) = shared_mem[0] * vol_n;
    }      
}

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_800_step_1(int * itype_cent, double *xatom_cent, int *icorner, double *AL, int n3, int n2, int n1, int *ion_type_cent, int *imax_ion,  double *r_ion, double *dxyz_box, 
                        double *funcr2, double *box, double *box2, double *box3, double *dbox, double *dbox2, double *dbox3, double *fact_store, double *fact2_store, double Rbox2, 
                        int nr, int nm3, int nm2, int nm1, double *vrhop_tot_coul, double *vrhop_m_coul, double *dvxc_tot, double *dvxc_m, double vol_n, int id3, int id2, int id1,
                        double *force_cent, double Rm2, int itype_mol, double *rho_ion, int imax_nr, int imax_ntype_cent, double *block_reduce_sum)
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

    float fact2 = 0;
    if (itype_ion != 0)
    {
        fact2 = 1.0;
        fact22 = (imax_ion(itype_ion) - 1) / (float)r_ion(imax_ion(itype_ion), itype_ion);
    }
    fact2_store[ia-1] = fact2;

    while(loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        box3(i1, j1, k1, ia) = 0;
        dbox(i1, j1, k1, 1, ia) = 0;
        dbox(i1, j1, k1, 2, ia) = 0;
        dbox(i1, j1, k1, 3, ia) = 0;

        dbox3(i1, j1, k1, 1, ia) = 0;
        dbox3(i1, j1, k1, 2, ia) = 0;
        dbox3(i1, j1, k1, 3, ia) = 0;

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
            float f11 = (x - 1.5) * fact11;
            float f22 = (2 - 2 * x) * fact11;
            float f33 = (x - 0.5) * fact11;
            ir += 1;
            float df = 0;
            if (d > 1e-10)
            {
                df = ((float)funcr2(ir, itype, itype_mol) * f11 + (float)funcr2(ir + 1, itype, itype_mol) * f22 + (float)funcr2(ir + 2, itype, itype_mol) * f33) / d;
            }
            float df2 = -2 / powf(1.5, 2) * expf(-powf((d/1.5), 2));

            float df3 = 0;
            dbox(i1, j1, k1, 1, ia) = -df * dx;
            dbox(i1, j1, k1, 2, ia) = -df * dy;
            dbox(i1, j1, k1, 3, ia) = -df * dz;

            if (itype_ion != 0)
            {
                yy = d * fact22;
                ir = yy;
                x = yy - ir;
                f1 = 1 - x - 0.5 * x * (1 - x);
                f2 = x + x * (1 - x);
                f3 = -0.5 * x * (1-x);
                f11 = (x - 1.5) * fact22;
                f22 = (2 - 2 * x) * fact22;
                f33 = (x - 0.5) * fact22;
                ir += 1;                
                box3(i1, j1, k1, ia) = (float)rho_ion(ir, itype_ion) * f1 + (float)rho_ion(ir + 1, itype_ion) * f2 + (float)rho_ion(ir + 2, itype_ion) * f3;
                if (d > 1e-10)
                {
                    df3 = ((float)rho_ion(ir, itype_ion) * f11 + (float)rho_ion(ir + 1, itype_ion) * f22 + (float)rho_ion(ir + 2, itype_ion) * f33) / d;
                }

                dbox3(i1, j1, k1, 1, ia) = -df3 * dx;
                dbox3(i1, j1, k1, 2, ia) = -df3 * dy;
                dbox3(i1, j1, k1, 3, ia) = -df3 * dz;
            }        
        }

        // // 这里修改 减号变加号
        dbox3(i1, j1, k1, 1, ia) = (float)dbox(i1, j1, k1, 1, ia) + (float)fact2_store[ia - 1] * (float)dbox3(i1, j1, k1, 1, ia);
        dbox3(i1, j1, k1, 2, ia) = (float)dbox(i1, j1, k1, 2, ia) + (float)fact2_store[ia - 1] * (float)dbox3(i1, j1, k1, 2, ia);
        dbox3(i1, j1, k1, 3, ia) = (float)dbox(i1, j1, k1, 3, ia) + (float)fact2_store[ia - 1] * (float)dbox3(i1, j1, k1, 3, ia);

        loop_id += loop_stride;
    }


}

template<unsigned int blocksize>
__global__ void cudakernel_loop_800_step_2(double *input, double *output, size_t reduce_size)
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

template<unsigned int blocksize>
__global__ void __launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
cudakernel_loop_800_step_3(int * itype_cent, double *xatom_cent, int *icorner, double *AL, int n3, int n2, int n1, int *ion_type_cent, int *imax_ion,  double *r_ion, double *dxyz_box, 
                        double *funcr2, double *box, double *box2, double *box3, double *dbox, double *dbox2, double *dbox3, double *fact_store, double *fact2_store, double Rbox2, 
                        int nr, int nm3, int nm2, int nm1, double *vrhop_tot_coul, double *vrhop_m_coul, double *dvxc_tot, double *dvxc_m, double vol_n, int id3, int id2, int id1,
                        double *force_cent, double Rm2, int itype_mol, double *rho_ion, int imax_nr, int imax_ntype_cent, double *sum_res, double *block_reduce_sum, int reduce_block_size)
{
    int ia = blockIdx.z + 1;
    int k1 = blockIdx.y - id3;
    int loop_id = blockIdx.x * blockDim.x + threadIdx.x;
    int loop_stride = blockDim.x * gridDim.x;
    int iteration_counts = (2 * id1 + 1) * (2 * id2 + 1);
    int tidx = threadIdx.x;

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

    // double sum2 = sum_res[reduce_block_size * blockIdx.z + 0 * (reduce_block_size * gridDim.z)];
    // double dsum3_1 = sum_res[reduce_block_size * blockIdx.z + 1 * (reduce_block_size * gridDim.z)];
    // double dsum3_2 = sum_res[reduce_block_size * blockIdx.z + 2 * (reduce_block_size * gridDim.z)];
    // double dsum3_3 = sum_res[reduce_block_size * blockIdx.z + 3 * (reduce_block_size * gridDim.z)];
    // double dsum_1 = sum_res[reduce_block_size * blockIdx.z + 4 * (reduce_block_size * gridDim.z)];
    // double dsum_2 = sum_res[reduce_block_size * blockIdx.z + 5 * (reduce_block_size * gridDim.z)];
    // double dsum_3 = sum_res[reduce_block_size * blockIdx.z + 6 * (reduce_block_size * gridDim.z)];

    // --debug
    // if (ia == 1)
    // {
    //     printf("sum2 : %lf \r\n", sum2);
    //     printf("dsum3_1 : %lf \r\n", dsum3_1);
    //     printf("dsum3_2 : %lf \r\n", dsum3_2);
    //     printf("dsum3_3 : %lf \r\n", dsum3_3);
    //     printf("dsum_1 : %lf \r\n", dsum_1);
    //     printf("dsum_2 : %lf \r\n", dsum_2);
    //     printf("dsum_3 : %lf \r\n", dsum_3);
    // }
    // --debug
    // dsum3_1 = dsum3_1 / sum2;
    // dsum3_2 = dsum3_2 / sum2;
    // dsum3_3 = dsum3_3 / sum2;
    // dsum_1 = dsum_1 / sum2;
    // dsum_2 = dsum_2 / sum2;
    // dsum_3 = dsum_3 / sum2;
    
    double partial_force_cent1 = 0, partial_force_cent2 = 0, partial_force_cent3 = 0;
    __shared__ double shared_mem[blocksize];
    while(loop_id < iteration_counts)
    {
        int j1 = loop_id / (2 * id1 + 1) - id2;
        int i1 = loop_id % (2 * id1 + 1) - id1;

        int k = (km0 + k1 - 1 + nm3) % nm3 + 1;
        int j = (jm0 + j1 - 1 + nm2) % nm2 + 1;
        int i = (im0 + i1 - 1 + nm1) % nm1 + 1;  
        int k2 = (k0 + k1 - 1 + n3) % n3 + 1;
        int j2 = (j0 + j1 - 1 + n2) % n2 + 1;
        int i2 = (i0 + i1 - 1 + n1) % n1 + 1;  
        // dbox3(i1, j1, k1, 1, ia) += -dsum3_1 * box2(i1, j1, k1, ia);
        // dbox3(i1, j1, k1, 2, ia) += -dsum3_2 * box2(i1, j1, k1, ia);
        // dbox3(i1, j1, k1, 3, ia) += -dsum3_3 * box2(i1, j1, k1, ia);
        // dbox(i1, j1, k1, 1, ia) += -dsum_1 * box2(i1, j1, k1, ia);
        // dbox(i1, j1, k1, 2, ia) += -dsum_2 * box2(i1, j1, k1, ia);
        // dbox(i1, j1, k1, 3, ia) += -dsum_3 * box2(i1, j1, k1, ia);
        // double vtmp1 = vrhop_tot_coul(i2, j2, k2) - vrhop_m_coul(i, j, k);
        // double vtmp2 = dvxc_tot(i2, j2, k2) - dvxc_m(i, j, k);

        double vtmp1 = vrhop_tot_coul(i2, j2, k2);
        double vtmp2 = dvxc_tot(i2, j2, k2);

        partial_force_cent1 += vtmp1 * dbox3(i1, j1, k1, 1, ia) + vtmp2 * dbox(i1, j1, k1, 1, ia);
        partial_force_cent2 += vtmp1 * dbox3(i1, j1, k1, 2, ia) + vtmp2 * dbox(i1, j1, k1, 2, ia);
        partial_force_cent3 += vtmp1 * dbox3(i1, j1, k1, 3, ia) + vtmp2 * dbox(i1, j1, k1, 3, ia);
        loop_id += loop_stride;
    }
    int wrt_idx = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;
    int wrt_stride = gridDim.x * gridDim.y * gridDim.z;
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent1, tidx, &block_reduce_sum[wrt_idx]);
    __syncthreads();
    wrt_idx += wrt_stride;
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent2, tidx, &block_reduce_sum[wrt_idx]);
    __syncthreads();
    wrt_idx += wrt_stride;
    reduce_sharedmem<double, blocksize> (shared_mem, partial_force_cent3, tidx, &block_reduce_sum[wrt_idx]);
}

template<unsigned int blocksize>
__global__ void cudakernel_loop_800_step_4(double *input, double *force_cent, size_t reduce_size, double scale)
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
    reduce_sharedmem<double, blocksize>(reduce_shared_mem, tmp, tidx, &force_cent(bidy + 1, bidx + 1), scale);
}

void fcc_polar_loop_113to129(double *dxyz_box, double *AL, int id3, int id2, int id1, int n3, int n2, int n1)
{
    dim3 gridDim(2 * id3 + 1, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_113to129<<<gridDim, blockDim, 0, 0>>>(dxyz_box, AL, id3, id2, id1, n3, n2, n1);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_polar_loop_500(int natom_m_itype_mol, int *ion_type_atomp, double *atom_charge_param, double *xatom_m, int imol, int itype_mol, int n3, int n2, int n1, int id3, 
                        int id2, int id1, int natom_mm, int nmolm, int ntype_mm, double *AL, double *dxyz_box, double *CC_pol, double *rhop_tot, double Rbox2)
{
    // dim3 gridDim(natom_m_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_500<<<gridDim, blockDim, 0, 0>>>(ion_type_atomp, atom_charge_param, xatom_m, imol, itype_mol, n3, n2, n1, id3, id2, id1, natom_mm, nmolm, ntype_mm, AL, dxyz_box, CC_pol, rhop_tot, Rbox2);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), natom_m_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_500_step_1<<<gridDim, blockDim, 0, 0>>>(ion_type_atomp, atom_charge_param, xatom_m, imol, itype_mol, n3, n2, n1, id3, id2, id1, natom_mm, nmolm, ntype_mm, AL, dxyz_box, CC_pol, rhop_tot, Rbox2);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_polar_loop_312to339(int ncent_itype_mol, int *icent, int imol, int itype_mol, double *w_cent, double *xc_cent, double *xatom_cent, int *icorner, int n3, int n2, int n1, 
                            int nm3, int nm2, int nm1, double *xatom_m, int natom_mm, int nmolm, int *nat_cent)
{
    dim3 gridDim(1, 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop312to339<<<gridDim, blockDim, 0, 0>>>(ncent_itype_mol, icent, imol, itype_mol, w_cent, xc_cent, xatom_cent, icorner, n3, n2, n1, nm3, nm2, nm1, xatom_m, natom_mm, nmolm, nat_cent);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_polar_loop_600(int ncent_itype_mol, int *itype_cent, int itype_mol, double *xatom_cent, int *icorner, double *AL, int *ion_type_cent, int *imax_ion, double *r_ion, 
                        int id3, int id2, int id1, double *box, double *box2, double *box3, double *dxyz_box, double *funcr2, double *fact_store, double *fact2_store, int n3, int n2, int n1,
                        int nr, double Rm2, double Rbox2, double *rho_ion, double vol, int nm3, int nm2, int nm1, double *rho_m, double *Q_type, double *z_ion, int imax_nr, int imax_ntype_cent)
{
    // dim3 gridDim(ncent_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_600<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, itype_mol, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, id3, id2, id1, box, box2, box3, dxyz_box, funcr2, fact_store, fact2_store, 
    //                                                         n3, n2, n1, nr, Rm2, Rbox2, rho_ion, vol, nm3, nm2, nm1, rho_m, Q_type, z_ion, imax_nr, imax_ntype_cent);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), ncent_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_600_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, itype_mol, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, id3, id2, id1, box, box2, box3, dxyz_box, funcr2, fact_store, fact2_store, 
                                                            n3, n2, n1, nr, Rm2, Rbox2, rho_ion, vol, nm3, nm2, nm1, rho_m, Q_type, z_ion, imax_nr, imax_ntype_cent, block_reduce_sum);
    cudakernel_loop_600_step_2<1024><<<dim3(ncent_itype_mol, 3, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, block_reduce_sum, reduce_block_size, vol / (n1 * n2 * n3));
    cudakernel_loop_600_step_3<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, itype_mol, xatom_cent, icorner, AL, ion_type_cent, imax_ion, r_ion, id3, id2, id1, box, box2, box3, dxyz_box, funcr2, fact_store, fact2_store, 
                                                            n3, n2, n1, nr, Rm2, Rbox2, rho_ion, vol, nm3, nm2, nm1, rho_m, Q_type, z_ion, imax_nr, imax_ntype_cent, block_reduce_sum, reduce_block_size);
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );
}

void fcc_polar_loop700(int natom_m_itype_mol, int *ion_type_atomp, int imol, int itype_mol, double *atom_charge_param, double *xatom_m, int *icorner, int n3, int n2, int n1,
                        int id3, int id2, int id1, double *CC_pol, double *rhop_m, double *dxyz_box, double Rbox2, int natom_mm, int nmolm, double *AL, int nm3, int nm2, int nm1)
{
    // dim3 gridDim(natom_m_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_700<<<gridDim, blockDim, 0, 0>>>(ion_type_atomp, imol, itype_mol, atom_charge_param, xatom_m, icorner, n3, n2, n1, id3, id2, id1, CC_pol, rhop_m, dxyz_box, Rbox2, natom_mm, nmolm, AL, nm3, nm2, nm1);  
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), natom_m_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_700_step_1<<<gridDim, blockDim, 0, 0>>>(ion_type_atomp, imol, itype_mol, atom_charge_param, xatom_m, icorner, n3, n2, n1, id3, id2, id1, CC_pol, rhop_m, dxyz_box, Rbox2, natom_mm, nmolm, AL, nm3, nm2, nm1);  
    // gpuErrchk( cudaPeekAtLastError() );
    // gpuErrchk( cudaDeviceSynchronize() );    
}

void fcc_polar_loop800(int ncent_itype_mol, int * itype_cent, double *xatom_cent, int *icorner, double *AL, int n3, int n2, int n1, int *ion_type_cent, int *imax_ion, double *r_ion, 
                        double *dxyz_box, double *funcr2, double *box, double *box2, double *box3, double *dbox, double *dbox2, double *dbox3, double *fact_store, double *fact2_store, 
                        double Rbox2, int nr, int nm3, int nm2, int nm1, double *vrhop_tot_coul, double *vrhop_m_coul, double *dvxc_tot, double *dvxc_m, double vol_n, int id3, int id2,
                        int id1, double *force_cent, double Rm2, int itype_mol, double *rho_ion, int imax_nr, int imax_ntype_cent)
{
    // dim3 gridDim(ncent_itype_mol, 1, 1);
    // dim3 blockDim(BLOCKSIZE, 1, 1);
    // cudakernel_loop_800<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, icorner, AL, n3, n2, n1, ion_type_cent, imax_ion, r_ion, dxyz_box, funcr2, box, box2, box3, dbox, dbox2, dbox3, fact_store, fact2_store, 
    //                                                         Rbox2, nr, nm3, nm2, nm1, vrhop_tot_coul, vrhop_m_coul, dvxc_tot, dvxc_m, vol_n, id3, id2, id1, force_cent, Rm2, itype_mol, rho_ion, imax_nr, imax_ntype_cent);
    unsigned int reduce_block_size = ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)) * (2 * id3 + 1);
    double *block_reduce_sum = reinterpret_cast<double *>(buffer_d);
    dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1) / (BLOCKSIZE * 8)), (2 * id3 + 1), ncent_itype_mol);
    // dim3 gridDim(ceil((double)(2 * id2 + 1) * (2 * id1 + 1)), (2 * id3 + 1), ncent_itype_mol);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_loop_800_step_1<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, icorner, AL, n3, n2, n1, ion_type_cent, imax_ion, r_ion, dxyz_box, funcr2, box, box2, box3, dbox, dbox2, dbox3, fact_store, fact2_store, 
                                                            Rbox2, nr, nm3, nm2, nm1, vrhop_tot_coul, vrhop_m_coul, dvxc_tot, dvxc_m, vol_n, id3, id2, id1, force_cent, Rm2, itype_mol, rho_ion, imax_nr, imax_ntype_cent, block_reduce_sum);
    cudakernel_loop_800_step_2<1024><<<dim3(ncent_itype_mol, 7, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum, block_reduce_sum, reduce_block_size);
    cudakernel_loop_800_step_3<BLOCKSIZE><<<gridDim, blockDim, 0, 0>>>(itype_cent, xatom_cent, icorner, AL, n3, n2, n1, ion_type_cent, imax_ion, r_ion, dxyz_box, funcr2, box, box2, box3, dbox, dbox2, dbox3, fact_store, fact2_store, 
                                                            Rbox2, nr, nm3, nm2, nm1, vrhop_tot_coul, vrhop_m_coul, dvxc_tot, dvxc_m, vol_n, id3, id2, id1, force_cent, Rm2, itype_mol, rho_ion, imax_nr, imax_ntype_cent, block_reduce_sum, 
                                                            block_reduce_sum + 7 * ncent_itype_mol * reduce_block_size, reduce_block_size);
    cudakernel_loop_800_step_4<1024><<<dim3(ncent_itype_mol, 3, 1), dim3(1024, 1, 1), 0, 0>>>(block_reduce_sum + 7 * ncent_itype_mol * reduce_block_size, force_cent, reduce_block_size, vol_n);
    // gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}



