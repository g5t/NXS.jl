#include <stddef.h>
#include <omp.h>
// multi-processor distribute_pixels! for Real source and sink
void dpp1(double *src, long long nosr, long long *srcln,
          long long *idx, long long maxcnt, long long *ttl,
          double *snk, long long nosn, long long *sncln,
          long long nocol)
{
  long long i,j,k,matches;
  long long thisidx[maxcnt]; // the number of indicies that match any one sink pixel can not exceed the maximum calculated in julia
  omp_set_num_threads(8);
# pragma omp parallel shared ( src, nosr, srcln, idx, ttl, snk, nosn, sncln, nocol ) private ( i, j, k, matches, thisidx )
# pragma omp for
  for (i=0;i<nosn;i++) // loop over all sink pixels and find the contributing source pixels
  {
    // loop through all source pixel indicies and find which ones match this sink pixel index
    matches=0LL; // (re)initialize the matches counter
    for (j=0;j<nosr;j++)
    {
      if (idx[j]==i+1) // i+1 due to julia using 1-based indicies
      {
        thisidx[matches]=j;
        matches+=1;
      }
    }
    /* We now have linear indexing into the first N=D-1 dimesions of src for
       the `matches` number of source pixels that go into the sink pixel `i` */
    for (j=0;j<matches;j++) // loop over the matched source pixels
    {
      for (k=0;k<nocol;k++) // and the columns
      {
        snk[i+sncln[k]]+=src[thisidx[j]+srcln[k]]; // and add their values to the sink pixel columns'
        ttl[i+sncln[k]]+=1;
      }
    }
  }
}
