#include <stddef.h>
//#include <stdio.h>
//#include <cuda.h>

#define NUM_THREADS 512
// cuda (gpu) reduce_pixels! for MeasuredSymmetric source and sink
__global__ void rpg1ker(double *sr1, long long *srCln, long long *Nsr, long long *Ncol,
                        long long *idx, long long *ttl,double *sn1, long long *snCln,
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
          ttl[l]+=1LL; // add one to the number of contributing source pixels
        }
      }
    }
  }
}

extern "C" void rpg1(double *sr1, long long *srCln, long long Nsr, long long Ncol,
          long long *idx, long long *ttl, double *sn1, long long *snCln, long long Nsn,
          long long *snTar, long long Ntar)
{
  double *d_sr1, *d_sn1;
  long long *d_srCln, *d_idx, *d_ttl, *d_snCln, *d_snTar;
  long long *d_Nsr, *d_Ncol, *d_Ntar;

  int sized = sizeof(double);
  int sizel = sizeof(long long);
  int Nthreads = NUM_THREADS;
  // determine how many blocks are necessary to visit all Ntar target
  // indexes using Nthreads threads.
  int Nblocks = (Ntar%Nthreads > 0) ? 1+Ntar/Nthreads : Ntar/Nthreads;

  cudaMalloc( (void**)&d_sr1, Nsr*Ncol*sized);
  cudaMalloc( (void**)&d_sn1, Nsn*Ncol*sized);
  cudaMalloc( (void**)&d_srCln, Ncol*sizel);
  cudaMalloc( (void**)&d_idx, Nsr*sizel);
  cudaMalloc( (void**)&d_ttl, Nsn*Ncol*sizel);
  cudaMalloc( (void**)&d_snCln, Ncol*sizel);
  cudaMalloc( (void**)&d_snTar, Ntar*sizel);
  cudaMalloc( (void**)&d_Nsr,  sizel);
  cudaMalloc( (void**)&d_Ncol, sizel);
  cudaMalloc( (void**)&d_Ntar, sizel);


  cudaMemcpy( d_sr1,   sr1,   Nsr*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_sn1,   sn1,   Nsn*Ncol*sized, cudaMemcpyHostToDevice);
  cudaMemcpy( d_srCln, srCln, Ncol*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_idx,   idx,   Nsr*sizel,      cudaMemcpyHostToDevice);
  cudaMemcpy( d_ttl,   ttl,   Nsn*Ncol*sizel, cudaMemcpyHostToDevice);
  cudaMemcpy( d_snCln, snCln, Ncol*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_snTar, snTar, Ntar*sizel,     cudaMemcpyHostToDevice);
  cudaMemcpy( d_Nsr,   &Nsr,   sizel,          cudaMemcpyHostToDevice);
  cudaMemcpy( d_Ncol,  &Ncol,  sizel,          cudaMemcpyHostToDevice);
  cudaMemcpy( d_Ntar,  &Ntar,  sizel,          cudaMemcpyHostToDevice);

  rpg1ker<<<Nblocks,Nthreads>>>(d_sr1,d_srCln,d_Nsr,d_Ncol,d_idx,d_ttl,d_sn1,d_snCln,d_snTar,d_Ntar);

  cudaDeviceSynchronize();

  cudaMemcpy( sn1, d_sn1, Nsn*Ncol*sized, cudaMemcpyDeviceToHost);
  cudaMemcpy( ttl, d_ttl, Nsn*Ncol*sized, cudaMemcpyDeviceToHost);

  cudaFree( d_Ntar );
  cudaFree( d_Ncol );
  cudaFree( d_Nsr  );
  cudaFree( d_snTar);
  cudaFree( d_snCln);
  cudaFree( d_ttl  );
  cudaFree( d_idx  );
  cudaFree( d_srCln);
  cudaFree( d_sn1  );
  cudaFree( d_sr1  );
}
