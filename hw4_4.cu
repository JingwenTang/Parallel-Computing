/*
========================================================================================
Name: cuda_demo.cu
Author: Mingran Peng
Class: EECS 587, Fall 2020
Description : Demo program for HW4
P.S. Fell free to use or modify this code for future terms of EECS 587 or other courses
Add you name if you modify it and preserve all author names
========================================================================================
*/
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <numeric>
#include <iostream>
#include <cstring>
#include <sstream>
#include <string>

using namespace std;

cudaError_t addWithCuda(double *c, unsigned int size, int n, double *Along, int t);

__global__ void addKernel(double *c, int size, int n, double *dev_Along, int t)
{
    double arr[4];
    double secondSmallest=0;
    double temp;
    double update;
    // calculate index here
    int localro = threadIdx.x;
    int localco = threadIdx.y;
    int globalro = blockIdx.x * blockDim.x + threadIdx.x;
    int globalco = blockIdx.y * blockDim.y + threadIdx.y;
    int globalind = globalro * n + globalco;
    // Synchronize all threads in the block to make sure copy is done
    __syncthreads();
    if ((globalro>0) && (globalco>0) && (globalro<(n-1)) && (globalco<(n-1))){
        int il = globalro+1;
        int is = globalro-1;
        int jl = globalco+1;
        int js = globalco-1;
        arr[0] = dev_Along[il*n+jl];
        arr[1] = dev_Along[il*n+js];
        arr[2] = dev_Along[is*n+jl];
        arr[3] = dev_Along[is*n+js];
        for(int i=0;i<4;i++)
            {
                for(int j=i+1;j<4;j++)
                {
                    if(arr[i]>arr[j])
                    {
                        temp  =arr[i];
                        arr[i]=arr[j];
                        arr[j]=temp;
                    }
                }
            }
        secondSmallest = arr[1];
        update = dev_Along[globalind] + secondSmallest;
    }
    else if ((globalro==0)||(globalco==0)||(globalro==(n-1))||(globalco==(n-1))){
        update = dev_Along[globalind];
    }
    else{}
    if((globalro<n)&&(globalco<n))
        c[globalind] = update;
}

int main(int argc, char* argv[])
{
    int n;
    int t;
    n = atoi(argv[1]);
    t = atoi(argv[2]);

    int size = n*n;//number of elements
    double *c; // returned array
    //initiallize
    double *Along;
    Along = new double [size];
    c = new double [size];
    for (int i=0;i<size;i++){
        int ro = floor(i/n);
        int co = floor(i%n);
        Along[i] = pow((1+cos(2*ro)+sin(co)),2);
    }


    cudaError_t cudaStatus = addWithCuda(c, size, n, Along, t);
    if (cudaStatus != cudaSuccess) {
        cout<<"addWithCuda failed!"<<endl;
        return -1;
    }
    
    //here we get the c array then we can do the sum and check the certain element
    // examine
    double initial_sum = 0;
    double sumc = accumulate(c, c+size, initial_sum);
    cout<<"Sum: "<<sumc<<endl;
    cout<<"A(37,47): "<<c[37*n+47]<<endl;
    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(double *c, unsigned int size, int n, double *Along, int t)
{
    double *dev_Along = 0;
    //dev_Along = new double[size];
    double *dev_c = 0;

    dev_c = new double[size];
    for (int i =0;i<size;i++){
        dev_c[i]=0;
    }
    cudaError_t cudaStatus;
    cudaEvent_t start, stop;
    float gpu_time = 0.0f;
    dim3 gridSize(ceil(n/32)+1,ceil(n/32)+1,1);
    dim3 blockSize(32,32,1);


    // Choose which GPU to run on, 0 if you have only one GPU
    // on-chip GPU does not count
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaSetDevice failed!"<<endl;
        goto Error;
    }
    // Malloc memory on GPU
    
    
    cudaStatus = cudaMalloc((void**)&dev_Along, size * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaMalloc failed!"<<endl;
        goto Error;
    }
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaMalloc failed!"<<endl;
        goto Error;
    }
    
    // Copy memory from Host to Device
    cudaStatus = cudaMemcpy(dev_Along, Along, size * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaMemcpy failed!"<<endl;
        goto Error;
    }
    // Set up timing
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    // Launch a kernel on the GPU with one thread for each element
    //cout<<"allocating "<<gridSize<<" blocks, "<<blockSize<<" threads per block"<<endl;
    for (int tt=0;tt<t;tt++){
        addKernel<<<gridSize, blockSize>>>(dev_c, size, n,dev_Along,t);

        if (cudaStatus != cudaSuccess) {
            cout<<"cudaMemcpy failed!"<<endl;
            goto Error;
        }
        cudaStatus = cudaMemcpy(dev_Along, dev_c, size * sizeof(double), cudaMemcpyDeviceToDevice);
    }
    
    
    
      // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        cout<<"addKernel failed: "<<cudaGetErrorString(cudaStatus)<<endl;
        goto Error;
    }
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaDeviceSynchronize failed: "<<cudaGetErrorString(cudaStatus)<<endl;
        goto Error;
    }
    // Close timing
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpu_time, start, stop);
    cout<<"Time spent: "<<gpu_time<<"ms"<<endl;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Copy memory from devide to host
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        cout<<"cudaMemcpy failed!"<<endl;
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_Along);


    return cudaStatus;
}
