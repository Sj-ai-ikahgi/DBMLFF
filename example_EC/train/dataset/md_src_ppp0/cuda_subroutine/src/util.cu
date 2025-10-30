#include "../include/util.h"

#define BLOCKSIZE 1024

extern std::map<std::string, cufftHandle> plan_map;

ncclDataType_t get_nccl_datatype(int type_idx) {
    switch (type_idx)
    {
        case 0:
            return ncclInt8;
        case 1:
            return ncclChar;
        case 2:
            return ncclUint8;
        case 3:
            return ncclInt32;
        case 4:
            return ncclInt;
        case 5:
            return ncclUint32;
        case 6:
            return ncclInt64;
        case 7:
            return ncclUint64;
        case 8:
            return ncclFloat16;
        case 9:
            return ncclFloat;
        case 10:
            return ncclFloat64;
        case 11:
            return ncclDouble;
        case 12:
            return ncclBfloat16;
        default:
            perror("invalid get_nccl_datatype type_idx");
            exit(0);
    }
}

ncclRedOp_t get_nccl_RedOp(int type_idx) {
    switch (type_idx)
    {
        case 0:
            return ncclSum;
        case 1:
            return ncclProd;
        case 2:
            return ncclMin;
        case 3:
            return ncclMax;
        case 4:
            return ncclAvg;
        default:
            perror("invalid get_nccl_RedOp type_idx");
            exit(0);
    }
}

cufftHandle get_plan(int n1, int n2, int n3)
{
    char str_dst[30] = {0};
    char str1[10] = {0};
    char str2[10] = {0};
    char str3[10] = {0};
    snprintf(str1, sizeof(str1), "%d", n1);
    snprintf(str2, sizeof(str2), "%d", n2);
    snprintf(str3, sizeof(str3), "%d", n3);
    strcat(str_dst, str1);
    strcat(str_dst, str2);
    strcat(str_dst, str3); 
    std::string key_str(str_dst); 
    std::map<std::string, int>::iterator iter = plan_map.find(key_str);
    if (iter != plan_map.end())
    {
        // printf("%s \r\n", iter->first.c_str());
        return iter->second;
    }
    else
    {
        printf("error : not find plan %d %d %d \r\n", n1, n2, n3);
        exit(1);
    }
}

void zero_gpu_mem(void *ptr_device, int len)
{
    cudaMemset(ptr_device, 0, len);
}

void print_gpu_array(void *ptr_device, int len)
{
    // cudaMemset(ptr_device, 0, len);
    // Todo print Array on GPU
    double *ptr_host;
    int bytes = len*sizeof(double);
    cudaMallocHost((void **)&ptr_host,bytes);
    cudaMemcpy(ptr_host, ptr_device, bytes, cudaMemcpyDeviceToHost);
    // for(int i=0;i<len;i++){
    //     printf("%f ",ptr_host[i]);
    // };

    // std::cout<<"hello "<<len<<std::endl;
    // std::fstream f;
    // f.open("data.txt",std::ios::out);
	// //输入你想写入的内容 
    // for(int i=0;i<len;i++){
    //     f<<i<<" "<<ptr_host[i]<<std::endl;
    // }
	// f.close();
    // printf("\n");
    cudaFreeHost(ptr_host);
}

__device__ double UxcCA(double rho_tmp, double *uxc2)
{
    const double thrd = 0.333333333333333;
    const double pi = 3.141592654;
    const double tft = 0.75;
    const double gamma = -0.2846;
    const double beta1 = 1.0529;
    const double beta2 = 0.3334;
    const double beta11 = 1.228383333;
    const double beta22 = 0.44453333333;
    const double cex = -1.969490099;
    const double A = 0.0622, B = -0.096, C = 0.0040, D = -0.0232;
    const double B1 = -0.11673333, C1 = 0.0026666667, D1 = -0.0168;

    // double rh3 = max(rho_tmp, 1.0e-16);

    // double rho = rh3;
    // rh3 = pow(rh3, thrd);
    // double rs = pow((pi * rho / tft), -thrd);
    // double vc = 0, ec = 0;
    // if (rs >= 1)
    // {
    //     double rootrs = sqrt(rs);
    //     vc = gamma * (1 + beta11 * rootrs + beta22 * rs) / pow((1 + beta1 * rootrs + beta2 * rs), 2);
    //     ec = gamma / (1 + beta1 * rootrs + beta2 * rs);
    // }
    // else
    // {
    //     double xlnrs = log(rs);
    //     vc = A * xlnrs + B1 + C1 * rs * xlnrs + D1 * rs;
    //     ec = A * xlnrs + B + C * rs * xlnrs + D * rs;
    // }
    // double vexcor = cex * rh3 + vc;

    // *uxc2 = (tft * cex * rh3 + ec) / 2;

    // return vexcor / 2;

    // double rho = max(rho_tmp, 1.0e-16);
    double rho = max(abs(rho_tmp), 1.0e-16);
    // double rho = sqrt(pow(rho_tmp,2.0) + 1.0e-16);

    double vexcor;
    double rh3;
    double rs = pow((pi * rho / tft), -thrd);
    double vc = 0, ec = 0;

    if (rs >= 1)
    {
        double rootrs = sqrt(rs);
        vc = gamma * (1 + beta11 * rootrs + beta22 * rs) / pow((1 + beta1 * rootrs + beta2 * rs), 2);
        ec = gamma / (1 + beta1 * rootrs + beta2 * rs);
    }
    else
    {
        double xlnrs = log(rs);
        vc = A * xlnrs + B1 + C1 * rs * xlnrs + D1 * rs;
        ec = A * xlnrs + B + C * rs * xlnrs + D * rs;
    }

    rh3 = pow(rho, thrd);

    // vexcor = (cex * rh3 + vc)*0.5;
    if (rho_tmp < 0.0) {
        vexcor = -1.0*(cex * rh3 + vc)*0.5;
    } else {
        vexcor = (cex * rh3 + vc)*0.5;
    }

    *uxc2 = (tft * cex * rh3 + ec)*0.5;

    return vexcor;
}

__global__ void cudakernel_make_dcmplx(cuDoubleComplex *complx, double *real, double *imag, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        complx[i] = make_cuDoubleComplex(real[i], imag[i]);
        i += stride;
    }
}

__global__ void cudakernel_make_dcmplx(cuDoubleComplex *complx, double real, double imag, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        complx[i] = make_cuDoubleComplex(real, imag);
        i += stride;
    }   
}

__global__ void cudakernel_make_dcmplx(cuDoubleComplex *complx, double *real, double imag, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        complx[i] = make_cuDoubleComplex(real[i], imag);
        i += stride;
    }   
}

__global__ void cudakernel_make_dcmplx(cuDoubleComplex *complx, double real, double *imag, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        complx[i] = make_cuDoubleComplex(real, imag[i]);
        i += stride;
    }   
}

__global__ void cudakernel_dsqrt(double *compute_array, double *result_array, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;

    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        result_array[i] = sqrt(abs(compute_array[i]));
        i += stride;
    }
}

__global__ void cudakernel_setval(double *set_array, double val, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        set_array[i] = val;
        i += stride;
    }
}

__global__ void cudakernel_getdreal(cuDoubleComplex *complx, double *real, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        real[i] = cuCreal(complx[i]);
        i += stride;
    }
}

__global__ void cudakernel_getdimag(cuDoubleComplex *complx, double *imag, size_t n)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        imag[i] = cuCimag(complx[i]);
        i += stride;
    }
}

__global__ void cudakernel_factorByConstant(cuDoubleComplex *complx, size_t n, const double factor)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < n)
    {
        complx[i].x = complx[i].x * factor;
        complx[i].y = complx[i].y * factor;
        i += stride;
    }
}

void make_dcmplx(cuDoubleComplex *complx, double *real, double *imag, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_make_dcmplx<<<gridDim, blockDim, 0, 0>>>(complx, real, imag, n);
}

void make_dcmplx(cuDoubleComplex *complx, double real, double imag, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_make_dcmplx<<<gridDim, blockDim, 0, 0>>>(complx, real, imag, n);
}

void make_dcmplx(cuDoubleComplex *complx, double *real, double imag, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_make_dcmplx<<<gridDim, blockDim, 0, 0>>>(complx, real, imag, n);
}

void make_dcmplx(cuDoubleComplex *complx, double real, double *imag, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_make_dcmplx<<<gridDim, blockDim, 0, 0>>>(complx, real, imag, n);
}

void dsqrt(double *compute_array, double *result_array, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_dsqrt<<<gridDim, blockDim, 0, 0>>>(compute_array, result_array, n);
}

void setdval(double *set_array, double val, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_setval<<<gridDim, blockDim, 0, 0>>>(set_array, val, n);
}

void getdreal(cuDoubleComplex *complx, double *real, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_getdreal<<<gridDim, blockDim, 0, 0>>>(complx, real, n);
}

void getdimag(cuDoubleComplex *complx, double *imag, size_t n)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_getdimag<<<gridDim, blockDim, 0, 0>>>(complx, imag, n);
}

void factorByConstant(cuDoubleComplex *complx, size_t n, const double factor)
{
    dim3 gridDim(ceil((double)n / (BLOCKSIZE)), 1, 1);
    dim3 blockDim(BLOCKSIZE, 1, 1);
    cudakernel_factorByConstant<<<gridDim, blockDim, 0, 0>>>(complx, n, factor);
}

__global__ void cudakernel_vec_add(const double *src1, const double *src2, double *dst, size_t count)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx * blockDim.x + tidx;
    int stride = gridDim.x * blockDim.x;
    while (i < count)
    {
        dst[i] = src1[i] + src2[i];
        i += stride;
    }
}

void gaussj(double *A, int n, int np, double *B, int m, int mp)
{
    const int NMAX = 200;
    int ipiv[NMAX], indxr[NMAX], indxc[NMAX];
    double dum = 0;
    for (int i = 0; i < n; i++)
    {
        ipiv[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
        double big = 0;
        int irow, icol;
        for (int j = 0; j < n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (int k = 0; k < n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (abs(A[j + k * np]) >= big)
                        {
                            big = abs(A[j + k * np]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                    {
                        printf("Singular matrix 1\r\n");
                        system("pause");
                    }
                }
            }
        }
        ipiv[icol] = ipiv[icol] + 1;
        if (irow != icol)
        {
            for (int l = 0; l < n; l++)
            {
                dum = A[irow + l * np];
                A[irow + l * np] = A[icol + l * np];
                A[icol + l * np] = dum;
            }
            for (int l = 0; l < m; l++)
            {
                dum = B[irow + l * np];
                B[irow + l * np] = B[icol + l * np];
                B[icol + l * np] = dum;
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (A[icol + icol * np] == 0)
        {
            printf("Singular matrix 2\r\n");
            system("pause");
        }
        double pivinv = 1.0 / A[icol + icol * np];
        A[icol + icol * np] = 1.0;
        for (int l = 0; l < n; l++)
        {
            A[icol + l * np] = A[icol + l * np] * pivinv;
        }
        for (int l = 0; l < m; l++)
        {
            B[icol + l * np] = B[icol + l * np] * pivinv;
        }
        for (int ll = 0; ll < n; ll++)
        {
            if (ll != icol)
            {
                dum = A[ll + icol * np];
                A[ll + icol * np] = 0;
                for (int l = 0; l < n; l++)
                {
                    A[ll + l * np] = A[ll + l * np] - A[icol + l * np] * dum;
                }
                for (int l = 0; l < m; l++)
                {
                    B[ll + l * np] = B[ll + l * np] - B[icol + l * np] * dum;
                }
            }
        }
    }
    for (int l = n - 1; l >= 0; l--)
    {
        if (indxr[l] != indxc[l])
        {
            for (int k = 0; k < n; k++)
            {
                dum = A[k + indxr[l] * np];
                A[k + indxr[l] * np] = A[k + indxc[l] * np];
                A[k + indxc[l] * np] = dum;
            }
        }
    }
}

void get_ALI(double *AL, double *ALI)
{
    double *tmp = (double *)malloc(3 * sizeof(double));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            ALI[j + i * 3] = AL[i + j * 3];
        }
        tmp[i] = 1;
    }
    gaussj(ALI, 3, 3, tmp, 1, 1);
    free(tmp);
}

void fcc_vec_add(const double *src1, const double *src2, double *dst, size_t count)
{
    int blocksize = 1024;
    int gridsize = ceil((double)count / blocksize);
    cudakernel_vec_add<<<gridsize, blocksize, 0, 0>>>(src1, src2, dst, count);
}

void nccl_mpi_allreduce(const void* sendbuff, void* recvbuff, size_t count, int datatype, int op)
{
    NCCLCHECK(ncclAllReduce(sendbuff, recvbuff, count, get_nccl_datatype(datatype), get_nccl_RedOp(op), nccl_comm, 0));
}