#include <cmath>
#include <stdlib.h>
#include <stdio.h>

extern "C" {
 void givensc_( double *a, int* lda_in, int* n_in, double* tol_in )
{
  double tol = *tol_in;
  int n = *n_in;
  int lda = *lda_in;
  double t1,t2,cx,sx;
  double r,rs;

  //  #pragma omp target map(tofrom: a) map(to:n,lda,tol)
  {
    for( int j=0; j< n-1; j++ )
    {
      for( int i=n-1; i> j; i-- )
	{
	  cx = a[(i-1)*lda+j];
	  sx = a[i*lda+j];

	  r = cx*cx + sx*sx;
	  rs = 1.0/sqrt(r);
	  if  (   fabs( r ) > tol  ) {
	    cx = cx*rs;
	    sx = sx*rs;
	  }
	  else
	  {
	      cx = 1.0;
	      sx = 0.0;
	  }
	  double tmp_sx = sx;
	  const int tmp_i = i;
	  double tmp_cx = cx;
	  const int tmp_n = n;
	  const int tmp_lda = lda;
	  for( int k=0; k<tmp_n; k++)
	    {
	      double t1 = a[tmp_lda*(tmp_i-1)+k];
	      double t2 = a[tmp_i*tmp_lda+k];
	      a[(tmp_i-1)*tmp_lda+k] = t1*tmp_cx + t2*tmp_sx;
	      a[tmp_i*tmp_lda+k] = t2*tmp_cx - t1*tmp_sx;
	    }
	}
    }
  }
  return;
}
}
