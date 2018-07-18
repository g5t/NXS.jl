#include <stddef.h>
#include <omp.h>
// multi-processor reduce_pixels! for MeasuredSymmetric source and sink
// (not used at this point)

extern "C" void reduce_pixels_2_float_int(float *sr1, float *sr2, int *srCln, int Nsr, int Ncol, int *idx, int *ttl, float *sn1, float *sn2, int *snCln, int *snTar, int Ntar);
extern "C" void reduce_pixels_2_double_int(double *sr1, double *sr2, int *srCln, int Nsr, int Ncol, int *idx, int *ttl, double *sn1, double *sn2, int *snCln, int *snTar, int Ntar);
extern "C" void reduce_pixels_2_float_longlong(float *sr1, float *sr2, long long *srCln, long long Nsr, long long Ncol, long long *idx, long long *ttl, float *sn1, float *sn2, long long *snCln, long long *snTar, long long Ntar);
extern "C" void reduce_pixels_2_double_longlong(double *sr1, double *sr2, long long *srCln, long long Nsr, long long Ncol, long long *idx, long long *ttl, double *sn1, double *sn2, long long *snCln, long long *snTar, long long Ntar);

template <typename T, typename I> void reduce_pixels_2(T *sr1, T *sr2, I *srCln, I Nsr, I Ncol, I *idx, I *ttl, T *sn1, T *sn2, I *snCln, I *snTar, I Ntar)
{
  I i,j,k,l,m,mytarget;
# pragma omp parallel shared(sr1,sr2,srCln,Nsr,Ncol,idx,ttl,sn1,sn2,snCln,snTar,Ntar) private(i,j,k,l,m,mytarget)
# pragma omp for
  for (i=0;i<Ntar;i++) // loop over all target sink pixels that have contributing source pixels
  {
    mytarget=snTar[i]; // C-indexed target sink pixel (-1 since julia has 1-based indexing)
    //printf("---- looking for matches to sink pixel %d\n",mytarget+1);
    for (j=0;j<Nsr;j++) // loop through all source pixel indicies
    {
      if (idx[j]==mytarget) // if the source pixel index (-1 since julia has 1-based indexing) matches our target
      {
        for (k=0;k<Ncol;k++) //loop over the columns for that pixel and
        {
          l=snCln[k]+mytarget; // the index of the target's kth column
          m=srCln[k]+j;        // the index of the source's kth column
          sn1[l]+=sr1[m]; // add to the target sink1 pixel columns
          sn2[l]+=sr2[m]; // add to the target sink2 pixel columns
          ttl[l]+=1; // add one to the number of contributing source pixels
        }
      }
    }
  }
}

void reduce_pixels_2_float_int(float *sr1, float *sr2, int *srCln, int Nsr, int Ncol, int *idx, int *ttl, float *sn1, float *sn2, int *snCln, int *snTar, int Ntar){
  reduce_pixels_2<float,int>(sr1,sr2,srCln,Nsr,Ncol,idx,ttl,sn1,sn2,snCln,snTar,Ntar);
}
void reduce_pixels_2_double_int(double *sr1, double *sr2, int *srCln, int Nsr, int Ncol, int *idx, int *ttl, double *sn1, double *sn2, int *snCln, int *snTar, int Ntar){
  reduce_pixels_2<double,int>(sr1,sr2,srCln,Nsr,Ncol,idx,ttl,sn1,sn2,snCln,snTar,Ntar);
}
void reduce_pixels_2_float_longlong(float *sr1, float *sr2, long long *srCln, long long Nsr, long long Ncol, long long *idx, long long *ttl, float *sn1, float *sn2, long long *snCln, long long *snTar, long long Ntar){
  reduce_pixels_2<float,long long>(sr1,sr2,srCln,Nsr,Ncol,idx,ttl,sn1,sn2,snCln,snTar,Ntar);
}
void reduce_pixels_2_double_longlong(double *sr1, double *sr2, long long *srCln, long long Nsr, long long Ncol, long long *idx, long long *ttl, double *sn1, double *sn2, long long *snCln, long long *snTar, long long Ntar){
  reduce_pixels_2<double,long long>(sr1,sr2,srCln,Nsr,Ncol,idx,ttl,sn1,sn2,snCln,snTar,Ntar);
}
