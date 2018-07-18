#include <stddef.h>
#include <omp.h>
// multi-processor reduce_pixels! for Real source and sink
void rpp1(double *sr1, long long *srCln, long long Nsr, long long Ncol,
          long long *idx, long long *ttl,
          double *sn1, long long *snCln,
          long long *snTar, long long Ntar
          )
{
  long long i,j,k,l,m,mytarget;
# pragma omp parallel shared(sr1,srCln,Nsr,Ncol,idx,ttl,sn1,snCln,snTar,Ntar) private(i,j,k,l,m,mytarget)
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
          sn1[l]+=sr1[m]; // add to the target sink pixel columns
          ttl[l]+=1LL; // add one to the number of contributing source pixels
        }
      }
    }
  }
}
