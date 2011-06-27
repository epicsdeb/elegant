/* file: median_oag.c 
  * contents: compute_median(), compute_percentile(), compute_average(), compute_middle()
  * see also: find_XX() routines in rowmedian.c  These return the position of the
  *     median, etc. in parallel elegant
  * Yusong Wang, 2007
  */

#if SDDS_MPI_IO

#include "track.h"

 long approximate_percentiles_p(double *position, double *percent, long positions, double *x, long n, 
                              long bins)
 {
   double *hist, *cdf, xMin, xMax, xCenter, xRange;
   long i, j, k;

   if (!notSinglePart)
     return approximate_percentiles(position, percent, positions, x, n, bins);

   if (bins<2 || positions<=0) /* In the parallel version, n=0 will not be a condition to return, as no particle on Master */
     return 0;

   if (!(hist = calloc(sizeof(*hist), bins)))
     return 0;
   if (isMaster && notSinglePart)
     n = 0;
   find_min_max(&xMin, &xMax, x, n);
   if(isMaster) {
     xMin = DBL_MAX;
     xMax = -DBL_MAX;
   }
   find_global_min_max(&xMin, &xMax, n, MPI_COMM_WORLD);    
    
   xCenter = (xMax+xMin)/2;
   xRange = (xMax-xMin)*(1+1./bins)/2;
   xMin = xCenter-xRange;
   xMax = xCenter+xRange;
   if (isSlave)
     make_histogram(hist, bins, xMin, xMax, x, n, 1);
   else {
     if (!(cdf = calloc(sizeof(*cdf), bins)))
       return 0;
   }

   MPI_Reduce(hist, cdf, bins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   
   if (isMaster) {
     for (i=1; i<bins; i++)
       cdf[i] += cdf[i-1];
     for (i=0; i<bins; i++)
       cdf[i] /= cdf[bins-1];
     
     for (j=0; j<positions; j++) {
       for (i=k=0; i<bins; i++) {
	 if (cdf[i]<percent[j]/100.0)
	   k = i;
	 else
	   break;
       }
       /* printf ("p%ld,cdf[%ld] = %lf %", j, k, cdf[k]*100); */
       position[j] = xMin + (k*(xMax-xMin))/bins;
     }
     free (cdf);
   }
   MPI_Bcast(position, positions, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   free(hist);
   return 1;               
 }

long find_median_of_row_p(double *best_particle, double **x, long index, long n, long n_total)
{
    static double **data = NULL;
    static long last_n = 0;
    long i, j, global_median_position = n_total/2, last_index, first_index, median_found, best_i;
    double median_i, local_median;
    static double *median_array = NULL;
    static long *median_positions = NULL;
    static long *global_median_positions = NULL;

    if (index<0 && n_total<=0)
        return(-1);

    if (n_total < n_processors-1) {
      long *n_vector = (long*)tmalloc(n_processors*sizeof(*n_vector));
      long *offset = (long*)tmalloc(n_processors*sizeof(*offset));
      double median;

     
	static long last_n_total = 0;

	if (!last_n_total && data)
	  free_zarray_2d((void**)data, last_n_total, 7);
	data = (double**)zarray_2d(sizeof(**data), n_total, 7);
	last_n_total = n_total;
    
      
	MPI_Allgather(&n, 1, MPI_LONG, n_vector, 1, MPI_LONG, MPI_COMM_WORLD);
      offset[0] = 0;
      for (i=0; i<n_processors-1; i++) {
	n_vector[i] *= 7;
	offset[i+1] = offset[i] + n_vector[i];
      }	

      if (!n) /* We need allocate dummy memory to make MPI_Gatherv work */
	x = (double**)zarray_2d(sizeof(**data), 1, 7);
      MPI_Allgatherv (&x[0][0], 7*n, MPI_DOUBLE, &data[0][0], (int*)n_vector, (int*)offset, MPI_DOUBLE, MPI_COMM_WORLD);

      if (isMaster)
	best_i = find_median_of_row(&median, data, index, n_total)+1;
      /* We assume no more than one particle per CPU and it is filled from lower rank to higher rank */
      MPI_Bcast(&best_i, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      if (best_i == myid) {
	memcpy(best_particle, x[0], sizeof(*best_particle)*COORDINATES_PER_PARTICLE);
      }
      MPI_Bcast(best_particle, COORDINATES_PER_PARTICLE, MPI_DOUBLE, best_i, MPI_COMM_WORLD);

      return 1;
    }

    if (n>last_n) {
        if (data)
            free_zarray_2d((void**)data, last_n, 2);
        data = (double**)zarray_2d(sizeof(**data), n, 2);
        last_n = n;
    }

    if (!median_array)
      median_array = tmalloc(sizeof(*median_array)*n_processors);
    if (!global_median_positions)
	global_median_positions = tmalloc(sizeof(*global_median_positions)*n_processors);
    
    for (i=0; i<n; i++) {
        data[i][0] = x[i][index];
        data[i][1] = i;
        }

    /* Sort locally only once on each of the processors */
    set_up_row_sort(0, 2, sizeof(**data), double_cmpasc);
    qsort((void*)data, n, sizeof(*data), row_compare);

    first_index = 0; 
    last_index = n-1;

    median_found = 0;
    
    while (!median_found){   
      local_median = data[(last_index-first_index)/2][0];
      
      /* Gather all the medians */
      MPI_Allgather(&local_median, 1, MPI_DOUBLE, median_array, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      qsort((void*)(median_array+1), n_processors-1, sizeof(*median_array), double_cmpasc);
      
    
      /* find out the positions of each local median in the global array */
      j = 0;
      for (i = 0; i<n_processors; i++) {
	median_i = median_array[i];	
	while ( data[j][0] < median_i)
	  j++;
	median_positions[i] = j;
      }
      MPI_Allreduce (median_positions, global_median_positions, n_processors, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      for (i=0; i<n_processors; i++) {
#if MPI_DEBUG
	fprintf (stderr, "global_median_positions[%ld]=%ld\t", i, global_median_positions[i]);
#endif
	if (global_median_positions[i] == global_median_position) {
	  median_found =1;
	  /* The median should be on the i-th processor */
	  if (i == myid) {
	    best_i = (long) data[median_positions[i]][1];
	    memcpy(best_particle, x[best_i], sizeof(*best_particle)*COORDINATES_PER_PARTICLE);
	  }
	  MPI_Bcast(best_particle, COORDINATES_PER_PARTICLE, MPI_DOUBLE, i, MPI_COMM_WORLD);
	  break;
	}
	else if (global_median_positions[i] > global_median_position) {
	  break;
	}
      if (global_median_positions[myid] > global_median_position) /* The global median could be in the first part on this processor */
	last_index = median_positions[myid];
      else
	first_index = median_positions[myid];
      }
    
    }
    
    return 1;
    
}
#endif
