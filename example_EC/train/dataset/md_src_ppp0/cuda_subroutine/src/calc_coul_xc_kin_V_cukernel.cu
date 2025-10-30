#include "../include/util.h"
#include "../include/calc_coul_xc_kin_V.h"

extern void *buffer_h, *buffer_d;
extern int default_stream;
extern cudaStream_t stream1, stream2, stream3, stream4;
extern int bytes_of_bufferh, bytes_of_bufferd;

__global__ void cudakernel_calc_coul_xc_kin_V_loop_first(cuDoubleComplex *rho_q, cuDoubleComplex *psi_q, double *ALI, int n1, int n2, int n3, double pi)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int i = bidx + 1;
    int i1 = i - 1;
    if (i > n1 / 2)
        i1 = i - n1 - 1;
    int iteration_counts = n2 * n3;
    int it = tidx;

    // double q2_c2 = pow(2.0*pi/3.0, 2);
    double q2_c2 = pow(2.0*pi/2.0, 2);
    // double q2_c2 = 2.0*pi/0.05;
    double q2_c1 = 0.8*q2_c2;

    while (it < iteration_counts)
    {
        int j = it / n3 + 1;
        int k = it % n3 + 1;
        int j1 = j - 1;
        if (j > n2 / 2)
            j1 = j - n2 - 1;
        int k1 = k - 1;
        if (k > n3 / 2)
            k1 = k - n3 - 1;
        double gkx = 2 * pi * (ALI[0] * i1 + ALI[3] * j1 + ALI[6] * k1);
        double gky = 2 * pi * (ALI[1] * i1 + ALI[4] * j1 + ALI[7] * k1);
        double gkz = 2 * pi * (ALI[2] * i1 + ALI[5] * j1 + ALI[8] * k1);
        double q2 = pow(gkx, 2) + pow(gky, 2) + pow(gkz, 2);
        if (q2 > 1.0e-20)
        {
            rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = cuCmul(make_cuDoubleComplex(4 * pi / q2, 0.0), rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2]);
        }
        else
        {
            rho_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = make_cuDoubleComplex(0.0, 0.0);
        }

        // psi_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = cuCmul(make_cuDoubleComplex(0.5 * q2, 0.0), psi_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2]);

        double fact=1.0;
        if ( (q2_c1 < q2) && (q2 < q2_c2) ) {
            fact = pow(cos((q2 - q2_c1)*0.5*pi/(q2_c2-q2_c1)), 2.0);
        }
        if (q2_c2 <= q2) fact = 0.0;

        psi_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] = cuCmul(make_cuDoubleComplex(0.5 * q2 * fact, 0.0), psi_q[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2]);

        it += blockDim.x;
    }
}

__global__ void cudakernel_calc_coul_xc_kin_V_loop_second(double *vcoul, double *rho, double *rho_z, double *vxc, double *vxc2, double *psi, int n1, int n2, int n3, double fact, double fact_kin2, double f53, double f23, double *E_coul_Breduce_res, double *E_xc_Breduce_res, double *E_kin1_Breduce_res, double *E_kin2_Breduce_res)
{
    int bidx = blockIdx.x;
    int tidx = threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int it = tidx + bidx * blockDim.x;
    int iteration_counts = n1 * n2 * n3;
    __shared__ double reduce_E_coul[1024], reduce_E_xc[1024], reduce_E_kin[1024], reduce_E_kin2[1024];
    double E_coul = 0.0;
    double E_xc = 0.0;
    double E_kin = 0.0;
    double E_kin2 = 0.0;
    double uxc = 0.0;
    while (it < iteration_counts)
    {
        int k = it / (n1 * n2) + 1;
        int j = (it / n1) % n2 + 1;
        int i = it % n1 + 1;

        vxc(i, j, k) = 0;
        vxc2(i, j, k) = 0;

        E_coul += vcoul(i, j, k) * rho_z(i, j, k);
        vxc(i, j, k) = UxcCA(rho(i, j, k), &uxc);
        E_xc += uxc * rho(i, j, k);
        E_kin += fact * pow(abs(rho(i, j, k)), f53);
        E_kin2 += sqrt(abs(rho(i, j, k))) * psi(i, j, k);

        // vxc(i, j, k) += fact * f53 * pow(abs(rho(i, j, k)), f23);
        if (rho(i, j, k) < 0.0) {
            vxc(i, j, k) += -1.0 * fact * f53 * pow(abs(rho(i, j, k)), f23);
        } else {
            vxc(i, j, k) += fact * f53 * pow(abs(rho(i, j, k)), f23);
        }

        double rhom = abs(rho(i, j, k));
        if (rhom < 0.0001)
            rhom = 0.0001;
        vxc2(i, j, k) = vxc(i, j, k);
        vxc(i, j, k) += fact_kin2 * psi(i, j, k) / sqrt(rhom);
        it += stride;
    }
    reduce_E_coul[tidx] = E_coul;
    reduce_E_xc[tidx] = E_xc;
    reduce_E_kin[tidx] = E_kin;
    reduce_E_kin2[tidx] = E_kin2;
    __syncthreads();
    if (tidx < 512)
    {
        reduce_E_coul[tidx] += reduce_E_coul[tidx + 512];
        reduce_E_xc[tidx] += reduce_E_xc[tidx + 512];
        reduce_E_kin[tidx] += reduce_E_kin[tidx + 512];
        reduce_E_kin2[tidx] += reduce_E_kin2[tidx + 512];
    }
    __syncthreads();
    if (tidx < 256)
    {
        reduce_E_coul[tidx] += reduce_E_coul[tidx + 256];
        reduce_E_xc[tidx] += reduce_E_xc[tidx + 256];
        reduce_E_kin[tidx] += reduce_E_kin[tidx + 256];
        reduce_E_kin2[tidx] += reduce_E_kin2[tidx + 256];
    }
    __syncthreads();
    if (tidx < 128)
    {
        reduce_E_coul[tidx] += reduce_E_coul[tidx + 128];
        reduce_E_xc[tidx] += reduce_E_xc[tidx + 128];
        reduce_E_kin[tidx] += reduce_E_kin[tidx + 128];
        reduce_E_kin2[tidx] += reduce_E_kin2[tidx + 128];
    }
    __syncthreads();
    if (tidx < 64)
    {
        reduce_E_coul[tidx] += reduce_E_coul[tidx + 64];
        reduce_E_xc[tidx] += reduce_E_xc[tidx + 64];
        reduce_E_kin[tidx] += reduce_E_kin[tidx + 64];
        reduce_E_kin2[tidx] += reduce_E_kin2[tidx + 64];
    }
    __syncthreads();
    if (tidx < 32)
    {
        warpReduce<1024, double>(reduce_E_coul, tidx);
        warpReduce<1024, double>(reduce_E_xc, tidx);
        warpReduce<1024, double>(reduce_E_kin, tidx);
        warpReduce<1024, double>(reduce_E_kin2, tidx);
    }
    if (tidx == 0)
    {
        E_coul_Breduce_res[bidx] = reduce_E_coul[0];
        E_xc_Breduce_res[bidx] = reduce_E_xc[0];
        E_kin1_Breduce_res[bidx] = reduce_E_kin[0];
        E_kin2_Breduce_res[bidx] = reduce_E_kin2[0];
    }
}

void loop_first(cuDoubleComplex *rho_q, cuDoubleComplex *psi_q, double *ALI, int n1, int n2, int n3, double pi)
{
    dim3 gridDim(n1, 1, 1);
    dim3 blockDim(1024, 1, 1);
    cudakernel_calc_coul_xc_kin_V_loop_first<<<gridDim, blockDim, 0, 0>>>(rho_q, psi_q, ALI, n1, n2, n3, pi);
}

void loop_second(double *vcoul, double *rho, double *rho_z, double *vxc, double *vxc2, double *psi, int n1, int n2, int n3, double fact, double fact_kin2, double f53, double f23, double *results, void *buffer_d)
{
    double *E_coul_reduce_d = reinterpret_cast<double *>(buffer_d);
    double *E_xc_reduce_d = E_coul_reduce_d + static_cast<int>(ceil((double)(n1 * n2 * n3) / (1024 * 8))); 
    double *E_kin1_reduce_d = E_xc_reduce_d + static_cast<int>(ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    double *E_kin2_reduce_d = E_kin1_reduce_d + static_cast<int>(ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    // cudaMalloc((void **)&E_coul_reduce_d, ceil((n1 * n2 * n3) / (1024 * 8)) * sizeof(double));
    // cudaMalloc((void **)&E_xc_reduce_d, ceil((n1 * n2 * n3) / (1024 * 8)) * sizeof(double));
    // cudaMalloc((void **)&E_kin1_reduce_d, ceil((n1 * n2 * n3) / (1024 * 8)) * sizeof(double));
    // cudaMalloc((void **)&E_kin2_reduce_d, ceil((n1 * n2 * n3) / (1024 * 8)) * sizeof(double));
    dim3 gridDim(ceil((double)(n1 * n2 * n3) / (1024 * 8)), 1, 1);
    dim3 blockDim(1024, 1, 1);

    cudakernel_calc_coul_xc_kin_V_loop_second<<<gridDim, blockDim, 0, 0>>>(vcoul, rho, rho_z, vxc, vxc2, psi, n1, n2, n3, fact, fact_kin2, f53, f23, E_coul_reduce_d, E_xc_reduce_d, E_kin1_reduce_d, E_kin2_reduce_d);
    // cudaDeviceSynchronize();
    // gpuErrchk( cudaPeekAtLastError() );
    // cudaStream_t stream1, stream2, stream3, stream4;
    // cudaStreamCreate(&stream1); 
    // cudaStreamCreate(&stream2);
    // cudaStreamCreate(&stream3);
    // cudaStreamCreate(&stream4);
    cudakernel_blockreduce<<<1, 1024, 0, stream1>>>(E_coul_reduce_d, &results[0], (size_t)ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    cudakernel_blockreduce<<<1, 1024, 0, stream2>>>(E_xc_reduce_d, &results[1], (size_t)ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    cudakernel_blockreduce<<<1, 1024, 0, stream3>>>(E_kin1_reduce_d, &results[2], (size_t)ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    cudakernel_blockreduce<<<1, 1024, 0, stream4>>>(E_kin2_reduce_d, &results[3], (size_t)ceil((double)(n1 * n2 * n3) / (1024 * 8)));
    // cudaStreamDestroy(stream1);
    // cudaStreamDestroy(stream2);
    // cudaStreamDestroy(stream3);
    // cudaStreamDestroy(stream4);
    // cudaFree(E_coul_reduce_d);
    // cudaFree(E_xc_reduce_d);
    // cudaFree(E_kin1_reduce_d);
    // cudaFree(E_kin2_reduce_d);
}

void fcc_calc_coul_xc_kin_V(double *rho, double *rho_z, double *vxc, double *vxc2, double *vcoul, int n1, int n2, int n3, double *AL,
                            double *E_coul, double *E_xc, double *E_kin1, double *E_kin2, double fact_kin2, double pi)
{   
    // double *AL_h = (double *)malloc(3 * 3 * sizeof(double));
    // cudaMemset(buffer_d, 0, bytes_of_bufferd);
    // memset(buffer_h, 0, bytes_of_bufferh);
    cuDoubleComplex *fftbuffer_d = reinterpret_cast<cuDoubleComplex *>(buffer_d);
    cuDoubleComplex *rho_q_d = reinterpret_cast<cuDoubleComplex *>(fftbuffer_d + n1 * n2 * n3);
    cuDoubleComplex *psi_q_d = reinterpret_cast<cuDoubleComplex *>(rho_q_d + n1 * n2 * n3);
    double *zero_buffer_d = reinterpret_cast<double *>(psi_q_d + n1 * n2 * n3);
    double *psi_d = reinterpret_cast<double *>(zero_buffer_d + n1 * n2 * n3);
    double *ALI_d = reinterpret_cast<double *>(psi_d + n1 * n2 * n3);
    double *parameter_lists_d = reinterpret_cast<double *>(ALI_d + 3 * 3);

    double *AL_h = reinterpret_cast<double *>(buffer_h);
    cudaMemcpy(AL_h, AL, 3 * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    double vol = AL_h(3, 1) * (AL_h(1, 2) * AL_h(2, 3) - AL_h(1, 3) * AL_h(2, 2)) + 
                AL_h(3, 2) * (AL_h(1, 3) * AL_h(2, 1) - AL_h(1, 1) * AL_h(2, 3)) + 
                AL_h(3, 3) * (AL_h(1, 1) * AL_h(2, 2) - AL_h(1, 2) * AL_h(2, 1));
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
    make_dcmplx(fftbuffer_d, rho_z, 0.0, n1 * n2 * n3);
    //debug
    //double *rho_z_h = (double *)malloc(n1 * n2 * n3 * sizeof(double));
    //cudaMemcpy(rho_z_h, rho_z, n1 * n2 * n3 * sizeof(double), cudaMemcpyDeviceToHost);
    //int icount = 0;
    // std::cout<<"Check rho_z\n";
    // std::fstream f;
	// f.open("rho_z.dat",std::ios::out);
    // for (int i = 1; i <= n1; i++){
    //     for (int j = 1; j <= n2; j++){
    //         for (int k = 1; k <= n3; k++){
    //             if (rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] != 0 && icount < 100){
    //                 // printf("%20.15E %d %d %d \r\n", rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2], i, j, k);
    //                 f<<rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2]<<" "<<i<<" "<<j<<" "<<k<<"\n";
    //                 // icount++;
    //             }
    //         }
    //     }
    // }
    // f.close();
    // debug
    // cufftPlan3d(&plan, n3, n2, n1, CUFFT_Z2Z);
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_INVERSE);
    double factor = 1.0 / (double)(n1 * n2 * n3);
    factorByConstant(fftbuffer_d, n1 * n2 * n3, factor);

    // cuDoubleComplex *rho_q_d = reinterpret_cast<cuDoubleComplex *>(zero_buffer_d + n1 * n2 * n3);
    // cudaMalloc((void **)&rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex));

    cudaMemcpy(rho_q_d, fftbuffer_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    //debug
    // cuDoubleComplex *rho_q_h = (cuDoubleComplex *)malloc(n1 * n2 * n3 * sizeof(cuDoubleComplex));
    // cudaMemcpy(rho_q_h, rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    // int icount = 0;
    // for (int i = 1; i <= n1; i++){
    //     for (int j = 1; j <= n2; j++){
    //         for (int k = 1; k <= n3; k++){
    //             if (rho_q_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2].x != 0 && icount < 5){
    //                 printf("%20.15E %d %d %d \r\n", rho_q_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2].x, i, j, k);
    //                 icount++;
    //             }
    //         }
    //     }
    // }
    //debug

    // double *psi_d = reinterpret_cast<double *>(rho_q_d + n1 * n2 * n3);
    // cudaMalloc((void **)&psi_d, n1 * n2 * n3 * sizeof(double));

    dsqrt(rho, psi_d, n1 * n2 * n3);
    make_dcmplx(fftbuffer_d, psi_d, 0.0, n1 * n2 * n3);
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_INVERSE);
    factorByConstant(fftbuffer_d, n1 * n2 * n3, factor);

    // cuDoubleComplex *psi_q_d = reinterpret_cast<cuDoubleComplex *>(psi_d + n1 * n2 * n3);
    // cudaMalloc((void **)&psi_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex));

    cudaMemcpy(psi_q_d, fftbuffer_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cuDoubleComplex cci = make_cuDoubleComplex(0.0, 1.0);
    cuDoubleComplex gkx = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex gky = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex gkz = make_cuDoubleComplex(0.0, 0.0);
    loop_first(rho_q_d, psi_q_d, ALI_d, n1, n2, n3, pi);
    //debug
    // cuDoubleComplex *rho_q_h = (cuDoubleComplex *)malloc(n1 * n2 * n3 * sizeof(cuDoubleComplex));
    // cudaMemcpy(rho_q_h, rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    // int icount = 0;
    // std::cout<<"check rho_q\n";
    // for (int i = 1; i <= n1; i++){
    //     for (int j = 1; j <= n2; j++){
    //         for (int k = 1; k <= n3; k++){
    //             if (rho_q_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2].x != 0 && icount < 5){
    //                 printf("%20.15E %d %d %d \r\n", rho_q_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2].x, i, j, k);
    //                 icount++;
    //             }
    //         }
    //     }
    // }
    //debug
    cudaMemcpy(fftbuffer_d, rho_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_FORWARD);
    getdreal(fftbuffer_d, vcoul, n1 * n2 * n3);
    cudaMemcpy(fftbuffer_d, psi_q_d, n1 * n2 * n3 * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);
    cufftExecZ2Z(plan, fftbuffer_d, fftbuffer_d, CUFFT_FORWARD);
    getdreal(fftbuffer_d, psi_d, n1 * n2 * n3);
    double f53 = 5.0 / 3.0;
    double f23 = 2.0 / 3.0;
    double fact = 3.0 / 10.0 * pow((3.0 * pow(pi, 2)), (2.0 / 3.0));
    // parameter_lists_d[0 - 4] E_coul, E_xc, E_kin1, E_kin2

    // double *parameter_lists_d = reinterpret_cast<double *>(psi_q_d + n1 * n2 * n3);
    // double *parameter_lists = (double *)malloc(4 * sizeof(double));
    double *parameter_lists = reinterpret_cast<double *>(ALI_h + 3 * 3);
    // cudaMalloc((void **)&parameter_lists_d, 4 * sizeof(double));

    memset(parameter_lists, 0, 4 * sizeof(double));
    cudaMemset(parameter_lists_d, 0, 4 * sizeof(double));

    // cudaMemset(vxc, 0, n1 * n2 * n3 * sizeof(double));
    // cudaMemset(vxc2, 0, n1 * n2 * n3 * sizeof(double));
    //debug
    // double *vcoul_h = (double *)malloc(n1 * n2 * n3 * sizeof(double));
    // double *rho_z_h = (double *)malloc(n1 * n2 * n3 * sizeof(double));
    // cudaMemcpy(rho_z_h, rho_z, n1 * n2 * n3 * sizeof(double), cudaMemcpyDeviceToHost);
    // int icount = 0;
    // for (int i = 1; i <= n1; i++){
    //     for (int j = 1; j <= n2; j++){
    //         for (int k = 1; k <= n3; k++){
    //             if (rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] != 0 && icount < 5){
    //                 printf("%20.15E %d %d %d \r\n", rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2], i, j, k);
    //                 icount++;
    //             }
    //         }
    //     }
    // }
    //debug
    // debug liuweijian
    // double *vcoul_h = (double *)malloc(n1 * n2 * n3 * sizeof(double));
    // int icount = 0;
    // printf("check vcoul\n");
    // cudaMemcpy(vcoul_h, vcoul, n1 * n2 * n3 * sizeof(double), cudaMemcpyDeviceToHost);
    // for (int i = 1; i <= n1; i++){
    //     for (int j = 1; j <= n2; j++){
    //         // printf("%20.15E %d %d %d \r\n", rho_z_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2], i, j, k);
    //         for (int k = 1; k <= n3; k++){
    //             if (vcoul_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2] != 0 && icount < 5){
    //                 printf("%20.15E %d %d %d \r\n", vcoul_h[(i - 1) + (j - 1) * n1 + (k - 1) * n1 * n2], i, j, k);
    //                 icount=0;
    //             }
    //         }
    //     }
    // }
    // debug liuweijian
    cudaMemcpy(parameter_lists, parameter_lists_d, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    // cudaMemcpy(parameter_lists, parameter_lists_d, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    // std::cout << "Check E_coul: "<<parameter_lists[0]<<std::endl;
    loop_second(vcoul, rho, rho_z, vxc, vxc2, psi_d, n1, n2, n3, fact, fact_kin2, f53, f23, parameter_lists_d, parameter_lists_d + 4);
    cudaMemcpy(parameter_lists, parameter_lists_d, 4 * sizeof(double), cudaMemcpyDeviceToHost);
    *E_coul = parameter_lists[0];
    *E_xc = parameter_lists[1];
    *E_kin1 = parameter_lists[2];
    *E_kin2 = parameter_lists[3];
    // 这里有问题
    *E_coul = 0.5 * (*E_coul) * vol / (n1 * n2 * n3);
    // std::cout<<std::setprecision(9) << "E_coul: "<<*E_coul<<" vol:"<<vol<<std::endl;
    *E_xc = (*E_xc) * vol / (n1 * n2 * n3);
    *E_kin1 = (*E_kin1) * vol / (n1 * n2 * n3);
    *E_kin2 = (*E_kin2) * vol / (n1 * n2 * n3);

    // cudaFree(ALI_d);
    // cudaFree(fftbuffer_d);
    // cudaFree(zero_buffer_d);
    // cudaFree(rho_q_d);
    // cudaFree(psi_d);
    // cudaFree(psi_q_d);
    // cudaFree(parameter_lists_d);
    // free(parameter_lists);
    // free(AL_h);
    // free(ALI_h);
    // cufftDestroy(plan);
}
