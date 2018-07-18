#include <stddef.h>
//#include <stdio.h>
//#include <cuda.h>

#define NUM_THREADS 512
// cuda (gpu) reduce_pixels! for MeasuredSymmetric source and sink
__global__ void rpg2ker(double *sr1, double *sr2, long long *srCln, long long *Nsr, long long *Ncol,
                        long long *idx, long long *ttl,double *sn1, double *sn2, long long *snCln,
                        long long *snTar, long long *Ntar)
{
  long long i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < *Ntar)
  {
    long long j,k,l,m,mytarget=snTar[i]; // C-indexed target sink pixel
    for (j=0; j < *Nsr; j++) // loop through all source pixel indicies
    {
      if (idx[j]==mytarget) // if the source pixel index matches our target
      {
        for (k=0; k < *Ncol; k++) //loop over the columns for that pixel and
        {
          l=snCln[k]+mytarget; // the index of the target's kth column
          m=srCln[k]+j;        // the index of the source's kth column
          sn1[l]+=sr1[m]; // add to the target sink1 pixel columns
          sn2[l]+=sr2[m]; // add to the target sink2 pixel columns
          ttl[l]+=1LL; // add one to the number of contributing source pixels
        }
      }
    }
  }
}

extern "C" void rpg2(double *sr1, double *sr2, long long *srCln, long long Nsr, long long Ncol,
          long long *idx, long long *ttl, double *sn1, double *sn2, long long *snCln, long long Nsn,
          long long *snTar, long long Ntar)
{
  double *d_sr1, *d_sr2, *d_sn1, *d_sn2;
  long long *d_srCln, *d_idx, *d_ttl, *d_snCln, *d_snTar;
  long long *d_Nsr, *d_Ncol, *d_Ntar;

  int sized = sizeof(double);
  int sizel = sizeof(long long);
  int Nthreads = NUM_THREADS;
  // determine how many blocks are necessary to visit all Ntar target
  // indexes using Nthreads threads.
  int Nblocks = (Ntar%Nthreads > 0) ? 1+Ntar/Nthreads : Ntar/Nthreads;

  //int gpuCount, gpu;
  //CUresult res;
  //CUdevice dev;
  //CUcontext ctx;
  //size_t dev_free, dev_total;
  
  //cuInit(0);
  //cuDeviceGetCount(&gpuCount);
  //printf( (gpuCount>1)?"There are %d GPUs present.\n":"There is %d GPU present.\n" ,gpuCount);
  //for (gpu=0; gpu<gpuCount; gpu++){
  //  cuDeviceGet(&dev,gpu);
  //  cuCtxCreate(&ctx,0,dev);
  //  res=cuMemGetInfo(&dev_free,&dev_total);
  //  if (res != CUDA_SUCCESS)
  //    printf("  cuMemGetInfo failed for GPU %d! (status=%x)\n",gpu+1,res);
  //  printf(" GPU %d has %lu bytes free of %lu total bytes.\n",gpu+1,dev_free,dev_total);
  //  cuCtxDetach(ctx);
  //}

 // printf("starting allocation of %d bytes on device", 2*Nsr*Ncol*sized+2*Nsn*Ncol*sized+2*Ncol*sizel+Nsn*Ncol*sizel+Nsr*sizel+Ntar*sizel+3*sizel);

  cudaMalloc( (void**)&d_sr1, Nsr*Ncol*sized);
  cudaMalloc( (void**)&d_sr2, Nsr*Ncol*sized);
  cudaMalloc( (void**)&d_sn1, Nsn*Ncol*sized);
  cudaMalloc( (void**)&d_sn2, Nsn*Ncol*sized);
  cudaMalloc( (void**)&d_srCln, Ncol*sizel);
  cudaMalloc( (void**)&d_idx, Nsr*sizel);
  cudaMalloc( (void**)&d_ttl, Nsn*Ncol*sizel);
  cudaMalloc( (void**)&d_snCln, Ncol*sizel);
  cudaMalloc( (void**)&d_snTar, Ntar*sizel);
  cudaMalloc( (void**)&d_Nsr,  sizel);
  cudaMalloc( (void**)&d_Ncol, sizel);
  cudaMalloc( (void**)&d_Ntar, sizel);


  cudaMemcpy( d_sr1,   sr1,   Nsr*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_sr2,   sr2,   Nsr*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_sn1,   sn1,   Nsn*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_sn2,   sn2,   Nsn*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_srCln, srCln, Ncol*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_idx,   idx,   Nsr*sizel,      cudaMemcpyHostToDevice);
  cudaMemcpy( d_ttl,   ttl,   Nsn*Ncol*sizel, cudaMemcpyHostToDevice);
  cudaMemcpy( d_snCln, snCln, Ncol*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_snTar, snTar, Ntar*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_Nsr,   &Nsr,   sizel,          cudaMemcpyHostToDevice);
  cudaMemcpy( d_Ncol,  &Ncol,  sizel,          cudaMemcpyHostToDevice);
  cudaMemcpy( d_Ntar,  &Ntar,  sizel,          cudaMemcpyHostToDevice);

  rpg2ker<<<Nblocks,Nthreads>>>(d_sr1,d_sr2,d_srCln,d_Nsr,d_Ncol,d_idx,d_ttl,d_sn1,d_sn2,d_snCln,d_snTar,d_Ntar);

  cudaDeviceSynchronize();

  cudaMemcpy( sn1, d_sn1, Nsn*Ncol*sized, cudaMemcpyDeviceToHost);
  cudaMemcpy( sn2, d_sn2, Nsn*Ncol*sized, cudaMemcpyDeviceToHost);
  cudaMemcpy( ttl, d_ttl, Nsn*Ncol*sized, cudaMemcpyDeviceToHost);

  cudaFree( d_Ntar );
  cudaFree( d_Ncol );
  cudaFree( d_Nsr  );
  cudaFree( d_snTar);
  cudaFree( d_snCln);
  cudaFree( d_ttl  );
  cudaFree( d_idx  );
  cudaFree( d_srCln);
  cudaFree( d_sn2  );
  cudaFree( d_sn1  );
  cudaFree( d_sr2  );
  cudaFree( d_sr1  );
}
