#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <cufft.h>
#include <iomanip>
#include <string.h>
#include <string>
#include<fstream>
#include<iostream>
#include <map>
#include <vector>
#include "init.h"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg=cooperative_groups;

#define xatom_cent(x, y) xatom_cent[(x - 1) + (y - 1) * 3]
#define z_cent(x, y) z_cent[(x - 1) + (y - 1) * 1000]
#define itype_cent(x, y) itype_cent[(x - 1) + (y - 1) * 1000]
#define funcr2(x, y, z) funcr2[(x - 1) + (y - 1) * imax_nr + (z - 1) * imax_nr * imax_ntype_cent]
#define AL(x, y) AL[(x - 1) + (y - 1) * 3]
#define AL_h(x, y) AL_h[(x - 1) + (y - 1) * 3]
#define dxyz_box(x, y, z, w) dxyz_box[(x - 1) + (y + id1) * 3 + (z + id2) * 3 * (2 * id1 + 1) + (w + id3) * 3 * (2 * id1 + 1) * (2 * id2 + 1)]
#define rho(x, y, z) rho[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define rho_m(x, y, z) rho_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define rho_z(x, y, z) rho_z[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define rho_mz(x, y, z) rho_mz[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define xatom_m(x, y, z, w) xatom_m[(x - 1) + (y - 1) * 3 + (z - 1) * 3 * natom_mm + (w - 1) * 3 * natom_mm * nmolm]
#define icent(x, y, z) icent[(x - 1) + (y - 1) * 2 + (z - 1) * 2 * 1000]
#define w_cent(x, y, z) w_cent[(x - 1) + (y - 1) * 2 + (z - 1) * 2 * 1000]
#define nat_cent(x, y) nat_cent[(x - 1) + (y - 1) * 1000]
#define icorner(x) icorner[x - 1]
#define dbox(x, y, z, w, c) dbox[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) + \
                                    (c - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) * 3]
#define dbox2(x, y, z, w, c) dbox2[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) + \
                                    (c - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) * 3]
#define dbox3(x, y, z, w, c) dbox3[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) + \
                                    (c - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) * 3]
#define dbox_c(x, y, z, w, c) dbox_c[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) + \
                                    (c - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) * 3]
#define dbox3_c(x, y, z, w, c) dbox3_c[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) + \
                                    (c - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1) * 3]
#define force_cent(x, y) force_cent[(x - 1) + (y - 1) * 3]
#define vxc(x, y, z) vxc[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define vxc2(x, y, z) vxc2[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define vxc2_m(x, y, z) vxc2_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define vcoul(x, y, z) vcoul[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define vxcm(x, y, z) vxcm[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define vxc_m(x, y, z) vxc_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define vcoulm(x, y, z) vcoulm[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define vcoul_m(x, y, z) vcoul_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define psi(x, y, z) psi[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define imax_ion(x) imax_ion[x - 1]
#define r_ion(x, y) r_ion[(x - 1) + (y - 1) * 5000]
#define rho_ion(x, y) rho_ion[(x - 1) + (y - 1) * 5000]
#define box(x, y, z, w) box[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1)]
#define box2(x, y, z, w) box2[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1)]
#define box3(x, y, z, w) box3[(x + id1) + (y + id2) * (2 * id1 + 1) + (z + id3) * (2 * id1 + 1) * (2 * id2 + 1) + (w - 1) * (2 * id1 + 1) * (2 * id2 + 1) * (2 * id3 + 1)]
#define Q_type(x, y) Q_type[(x - 1) + (y - 1) * imax_ntype_cent]
#define z_ion(x) z_ion[x - 1]
#define ion_type_cent(x, y) ion_type_cent[(x - 1) + (y - 1) * 1000]
#define ion_type_atomp(x, y) ion_type_atomp[(x - 1) + (y - 1) * 1000]
#define atom_charge_param(x, y, z) atom_charge_param[(x - 1) + (y - 1) * 4 + (z - 1) * 4 * 20]
#define CC_pol(x, y, z, w) CC_pol[(x - 1) + (y - 1) * 6 + (z - 1) * 6 * natom_mm + (w - 1) * 6 * natom_mm * nmolm]
#define rhop_tot(x, y, z) rhop_tot[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define xc_cent(x) xc_cent[x - 1]
#define rhop_m(x, y, z) rhop_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define vrhop_tot_coul(x, y, z) vrhop_tot_coul[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define vrhop_m_coul(x, y, z) vrhop_m_coul[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]
#define dvxc_tot(x, y, z) dvxc_tot[(x - 1) + (y - 1) * n1 + (z - 1) * n1 * n2]
#define dvxc_m(x, y, z) dvxc_m[(x - 1) + (y - 1) * nm1 + (z - 1) * nm1 * nm2]

#ifdef __cplusplus
extern "C" {
#endif

void zero_gpu_mem(void *ptr_device, int len);
void print_gpu_array(void *ptr_device, int len);

#ifdef __cplusplus
}
#endif


void gaussj(double *A, int n, int np, double *B, int m, int mp);
void get_ALI(double *AL, double *ALI);
void make_dcmplx(cuDoubleComplex *complx, double *real, double *imag, size_t n);
void make_dcmplx(cuDoubleComplex *complx, double real, double imag, size_t n);
void make_dcmplx(cuDoubleComplex *complx, double *real, double imag, size_t n);
void make_dcmplx(cuDoubleComplex *complx, double real, double *imag, size_t n);
void dsqrt(double *compute_array, double *result_array, size_t n);
void setdval(double *set_array, double val, size_t n);
void getdreal(cuDoubleComplex *complx, double *real, size_t n);
void getdimag(cuDoubleComplex *complx, double *imag, size_t n);
void factorByConstant(cuDoubleComplex *complx, size_t n, const double factor);
cufftHandle get_plan(int n1, int n2, int n3);
__device__ double UxcCA(double rho_tmp, double *uxc2);

extern "C" {
    void fcc_vec_add(const double *src1, const double *src2, double *dst, size_t count);
    void nccl_mpi_allreduce(const void* sendbuff, void* recvbuff, size_t count, int datatype, int op);
}

template <unsigned int blocksize, typename T>
__device__ __forceinline__ void warpReduce(volatile T *sdata, int tid)
{
    if (blocksize >= 64)
        sdata[tid] += sdata[tid + 32];
    if (blocksize >= 32)
        sdata[tid] += sdata[tid + 16];
    if (blocksize >= 16)
        sdata[tid] += sdata[tid + 8];
    if (blocksize >= 8)
        sdata[tid] += sdata[tid + 4];
    if (blocksize >= 4)
        sdata[tid] += sdata[tid + 2];
    if (blocksize >= 2)
        sdata[tid] += sdata[tid + 1];
}

// template<typename FP_TYPE, unsigned int blocksize>
// __device__ void reduce_sharedmem(FP_TYPE *shared_mem, FP_TYPE val, int tidx)
// {
//     shared_mem[tidx] = val;
//     __syncthreads();
//     if (blocksize >= 1024)
//     {
//         if (tidx < 512)
//         {
//             shared_mem[tidx] += shared_mem[tidx + 512];
//         }
//         __syncthreads();
//     }
//     if (blocksize >= 512)
//     {
//         if (tidx < 256)
//         {
//             shared_mem[tidx] += shared_mem[tidx + 256];
//         }
//         __syncthreads();
//     }
//     if (blocksize >= 256)
//     {
//         if (tidx < 128)
//         {
//             shared_mem[tidx] += shared_mem[tidx + 128];
//         }
//         __syncthreads();
//     }
//     if (blocksize >= 128)
//     {
//         if (tidx < 64)
//         {
//             shared_mem[tidx] += shared_mem[tidx + 64];
//         }
//         __syncthreads();
//     }
//     if (tidx < 32)
//     {
//         warpReduce<blocksize, double>(shared_mem, tidx);
//     }
// }

template<typename FP_TYPE, unsigned int blocksize>
__inline__ __device__ void reduce_sharedmem(FP_TYPE *shared_mem, FP_TYPE val, int tidx, FP_TYPE *res)
{
    cg::thread_block cta = cg::this_thread_block();
    cg::thread_block_tile<32> tile = cg::tiled_partition<32>(cta);
    shared_mem[tidx] = cg::reduce(tile, val, cg::plus<FP_TYPE>());
    cg::sync(cta);
    if (cta.thread_rank() == 0) {
        val = 0;
        for (int i = 0; i < blocksize; i += tile.size()) {
            val += shared_mem[i];
        }
    *res = val;
    }
}

template<typename FP_TYPE, unsigned int blocksize>
__inline__ __device__ void reduce_sharedmem(FP_TYPE *shared_mem, FP_TYPE val, int tidx, FP_TYPE *res, FP_TYPE scale)
{
    cg::thread_block cta = cg::this_thread_block();
    cg::thread_block_tile<32> tile = cg::tiled_partition<32>(cta);
    shared_mem[tidx] = cg::reduce(tile, val, cg::plus<FP_TYPE>());
    cg::sync(cta);
    if (cta.thread_rank() == 0) {
        val = 0;
        for (int i = 0; i < blocksize; i += tile.size()) {
            val += shared_mem[i];
        }
    *res = val * scale;
    }
}

template<typename T>
__global__ void cudakernel_blockreduce(T *input, T *result, size_t n)
{
    int tidx = threadIdx.x;
    int stride = blockDim.x;
    __shared__ double share_mem[1024];
    int it = tidx;
    share_mem[tidx] = 0;
    while(it < n)
    {
        share_mem[tidx] += input[it];
        it += stride;
    }
    __syncthreads();
    if (tidx < 512)
    {
        share_mem[tidx] += share_mem[tidx + 512];
    }
    __syncthreads();
    if (tidx < 256)
    {
        share_mem[tidx] += share_mem[tidx + 256];
    }
    __syncthreads();
    if (tidx < 128)
    {
        share_mem[tidx] += share_mem[tidx + 128];
    }
    __syncthreads();
    if (tidx < 64)
    {
        share_mem[tidx] += share_mem[tidx + 64];
    }
    __syncthreads();
    if (tidx < 32)
    {
        warpReduce<1024, double>(share_mem, tidx);
    }
    if (tidx == 0)
    {
        *result = share_mem[0];
    }
}

#endif
