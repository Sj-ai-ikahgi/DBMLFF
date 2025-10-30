#ifndef INIT_H
#define INIT_H

#include <cufft.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <string.h>
#include "mpi.h"
#include <assert.h>
#include <string.h>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include "nccl.h"

extern int first_time;
extern void *buffer_h;
extern void *buffer_d;
extern int default_stream;
extern cudaStream_t stream1, stream2, stream3, stream4;
extern int bytes_of_bufferh, bytes_of_bufferd;
extern std::map<std::string, cufftHandle> plan_map;
extern ncclComm_t nccl_comm;

#ifdef __cplusplus
extern "C" {
#endif

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define MPICHECK(cmd) do {                          \
  int e = cmd;                                      \
  if( e != MPI_SUCCESS ) {                          \
    printf("Failed: MPI error %s:%d '%d'\n",        \
        __FILE__,__LINE__, e);   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)


#define CUDACHECK(cmd) do {                         \
  cudaError_t e = cmd;                              \
  if( e != cudaSuccess ) {                          \
    printf("Failed: Cuda error %s:%d '%s'\n",             \
        __FILE__,__LINE__,cudaGetErrorString(e));   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)


#define NCCLCHECK(cmd) do {                         \
  ncclResult_t r = cmd;                             \
  if (r!= ncclSuccess) {                            \
    printf("Failed, NCCL error %s:%d '%s'\n",             \
        __FILE__,__LINE__,ncclGetErrorString(r));   \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)

void init_resources(long long h_bytes, long long d_bytes, int n1, int n2, int n3, int nm1_all[], int nm2_all[], int nm3_all[], int ntype_m);
void destroy_resources();
void get_device();

#ifdef __cplusplus
}
#endif

#endif