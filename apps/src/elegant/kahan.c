/* file: Kahan.c
 * This is the implementation of Kahan's algorithm. 
 * It has been tested with GNU compiler without optimization flag. 
 * Other compiler or flags could get unexpected results.
 * Yusong Wang, 2005
 */

#include "track.h"

double Kahan (long length, double a[], double *error)
{
  
  double sum = 0.0, C = 0.0, Y, T;
  long i; 

  /* Kahan's summation formula */
  for (i=0; i<length; i++) {
     Y = a[i] - C;
     T = sum + Y;
     C = (T-sum) - Y;
     sum = T;
  }

  *error = -C;
  return sum;
}


double KahanPlus (double oldSum, double b, double *error)
{
  double sum, Y;

  Y = b + *error;
  sum = oldSum + Y;
  *error = Y - (sum-oldSum);
 
  return sum; 
}

#if USE_MPI
/* calculate sum of one dimensional array on several processors with Kahan's algorithm */
double KahanParallel (double sum,  double error, MPI_Comm comm)
{
  long i, offset=0;
  double *sumArray, *errorArray;
  double total = 0.0, error_sum = 0.0;

  sumArray = malloc(sizeof(double) * n_processors);
  errorArray = malloc(sizeof(double) * n_processors);

  MPI_Allgather(&sum,1,MPI_DOUBLE,sumArray,1,MPI_DOUBLE,comm);
  MPI_Allgather(&error,1,MPI_DOUBLE,errorArray,1,MPI_DOUBLE,comm);
  if (comm==workers)
    offset=1;
  for (i=1-offset; i<n_processors-offset; i++) {
    total = KahanPlus(total, sumArray[i], &error_sum);
    total = KahanPlus(total, errorArray[i], &error_sum);
  }
  tfree(sumArray); sumArray = NULL;
  tfree(errorArray); errorArray = NULL;
  return total;
}
#endif
