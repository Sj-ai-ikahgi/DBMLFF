#include "../include/init.h"

std::map<std::string, cufftHandle> plan_map;
int first_time = 1;
void *buffer_h;
void *buffer_d;
int default_stream = 0;
cudaStream_t stream1, stream2, stream3, stream4;
int bytes_of_bufferh, bytes_of_bufferd;
ncclComm_t nccl_comm;

static int stringCmp(const void *a, const void* b)
{
    char* m = (char*)a;
    char* n = (char*)b;
    int i, sum = 0;
    
    for( i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
        if (m[i] == n[i])
            continue;
        else
        {
            sum = m[i] - m[i];
            break;
        }       
    
    return sum;
}

void init_resources(long long h_bytes, long long d_bytes, int n1, int n2, int n3, int nm1_all[], int nm2_all[], int nm3_all[], int ntype_m)
{
    get_device();
    buffer_h = malloc(h_bytes);
    cudaMalloc((void **)&buffer_d, d_bytes);
    cudaStreamCreate(&stream1); 
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);
    cudaStreamCreate(&stream4);
    bytes_of_bufferd = d_bytes;
    bytes_of_bufferh = h_bytes;
    // generate cufft_plan n1, n2, n3
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
    cufftHandle plan;
    cufftPlan3d(&plan, n3, n2, n1, CUFFT_Z2Z);
    plan_map.insert(std::make_pair(key_str, plan));
    //generate cufft_plan nm1, nm2, nm3 
    for (int i = 0; i < ntype_m; i++) {
        int nm1 = nm1_all[i];
        int nm2 = nm2_all[i];
        int nm3 = nm3_all[i];
        memset(str_dst, 0, sizeof(str_dst));
        memset(str1, 0, sizeof(str1));
        memset(str2, 0, sizeof(str2));
        memset(str3, 0, sizeof(str3));
        snprintf(str1, sizeof(str1), "%d", nm1);
        snprintf(str2, sizeof(str2), "%d", nm2);
        snprintf(str3, sizeof(str3), "%d", nm3);
        strcat(str_dst, str1);
        strcat(str_dst, str2);
        strcat(str_dst, str3);
        std::string key_str(str_dst);
        if (plan_map.find(key_str) == plan_map.end()) {
            cufftHandle plan;
            cufftPlan3d(&plan, nm3, nm2, nm1, CUFFT_Z2Z);
            plan_map.insert(std::make_pair(key_str, plan));
        }
    }
}

void destroy_resources()
{
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
    cudaStreamDestroy(stream4);
    cudaFree(buffer_d);
    free(buffer_h);
    std::map<std::string, cufftHandle>::iterator it;
    for (it = plan_map.begin(); it != plan_map.end(); it++){
        cufftDestroy(it->second);
    }
    ncclCommDestroy(nccl_comm);
}

void get_device()
{
    char host_name[MPI_MAX_PROCESSOR_NAME];
    char (*host_names)[MPI_MAX_PROCESSOR_NAME];
    int n, namelen, color, rank, nprocs, myrank;
    size_t bytes;
    MPI_Comm nodeComm;
    
    int dev, err1;
    struct cudaDeviceProp deviceProp;
    
    /* Check if the device has been alreasy assigned */
    if(first_time)
    {
        first_time=0;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Get_processor_name(host_name,&namelen);
        
        bytes = nprocs * sizeof(char[MPI_MAX_PROCESSOR_NAME]);
        host_names = (char (*)[MPI_MAX_PROCESSOR_NAME]) malloc(bytes);
       
        strcpy(host_names[rank], host_name);
        
        /* get all the hostnames on all the nodes */
        for (n=0; n<nprocs; n++)
            MPI_Bcast(&(host_names[n]),MPI_MAX_PROCESSOR_NAME, MPI_CHAR, n, MPI_COMM_WORLD); 

        qsort(host_names, nprocs,  sizeof(char[MPI_MAX_PROCESSOR_NAME]), stringCmp);
        
        color = 0;

        for (n=0; n<nprocs-1; n++)
        {
            if(strcmp(host_name, host_names[n]) == 0)
            {
                break;
            }
            if(strcmp(host_names[n], host_names[n+1])) 
            {
                color++;
            }
        }

        MPI_Comm_split(MPI_COMM_WORLD, color, 0, &nodeComm);
        MPI_Comm_rank(nodeComm, &myrank);

        /* Find out how many GT200 GPUs are in the system and their device number */
        int deviceCount,slot = 0;
        int *devloc;
        cudaGetDeviceCount(&deviceCount);
        devloc=(int *)malloc(deviceCount*sizeof(int));
        
        for (dev = 0; dev < deviceCount; ++dev)
        {
            cudaGetDeviceProperties(&deviceProp, dev);
            // if(deviceProp.minor <= 10 )
            // {
                devloc[slot]=dev;
                slot++;
            // };
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
//      if( (devloc[myrank] >= deviceCount) )
//      {
//          printf ("Error:::Assigning device %d  to process on node %s rank %d \n",devloc[myrank],  host_name, rank );
//          exit(-1);
//          MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TYPE);
//          MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TYPE);
//      }
        int looprank=myrank%deviceCount;
        //printf (" Assigning device %d  to process on node %s rank %d, OK\n",devloc[myrank],  host_name, rank );
        printf (" Assigning device %d  to process on node %s rank %d, OK\n",looprank,  host_name, rank );
        cudaSetDevice(looprank);
        //cudaSetDevice(devloc[myrank]);
        //cudaSetDevice(0);
        //cudaGetDevice(&dev);
        //cudaGetDeviceProperties(&deviceProp, dev);
        fflush(stdout);
        free(devloc);
        free(host_names);

        //CUDA could have "lazy context" when you call cudaMalloc for the first time.
        //Some suggest that we have a cudaFree here to oliminate the initialization time.
        //It is quite useless, because we still have it anyway. Weile, 2014-11-26
        //assert(cudaFree(0) == cudaSuccess);
        //init_cuda_device();
        ncclUniqueId id;
        if (rank == 0) ncclGetUniqueId(&id);
        MPICHECK(MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD));
        NCCLCHECK(ncclCommInitRank(&nccl_comm, nprocs, id, rank));
    }
}