/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: compute_centroids()
 * contents: compute_centroids(), compute_sigmas(), zero_beam_sums(),
 *           accumulate_beam_sums()
 *
 * Michael Borland, 1989
 */


/* routine: compute_centroids()
 * purpose: compute centroids for (x, x', y, y', s, deltap/p)
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

void compute_centroids(
		       double *centroid,
		       double **coordinates,
		       long n_part
		       )
{
  long i_part, i_coord;
  double sum[6], *part;
  long active = 1; 
#ifdef USE_KAHAN
  double error[6];
#endif

#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
  long n_total;
#ifdef USE_KAHAN
  long j;
  double error_sum=0.0, error_total=0.0,
    **sumMatrix, **errorMatrix,
    *sumArray, *errorArray;
  sumMatrix = (double**)czarray_2d(sizeof(**sumMatrix), n_processors, 6); 
  errorMatrix = (double**)czarray_2d(sizeof(**errorMatrix), n_processors, 6);
  sumArray = malloc(sizeof(double) * n_processors);
  errorArray = malloc(sizeof(double) * n_processors);
#endif
  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }
#endif
 
  for (i_coord=0; i_coord<6; i_coord++) {
    sum[i_coord] = centroid[i_coord] = 0;
#ifdef USE_KAHAN
    error[i_coord] = 0.0;
#endif
  }

  if (active) {
    for (i_part=0; i_part<n_part; i_part++) {
      part = coordinates[i_part];
      for (i_coord=0; i_coord<6; i_coord++) 
#ifndef USE_KAHAN
	sum[i_coord] += part[i_coord];
#else
	sum[i_coord] = KahanPlus(sum[i_coord], part[i_coord], &error[i_coord]); 
#endif	
    }
  }
  if (!USE_MPI || !notSinglePart) {
    if (n_part)
      for (i_coord=0; i_coord<6; i_coord++)
	centroid[i_coord] = sum[i_coord]/n_part;

  }
#if USE_MPI
  if (notSinglePart) {
    if (parallelStatus!=trueParallel) {
      if (isMaster && n_part)
	for (i_coord=0; i_coord<6; i_coord++)
	  centroid[i_coord] = sum[i_coord]/n_part;
    }
    else {
      if (isMaster) {
	n_part = 0;
      }
      /* compute centroid sum over processors */
#ifndef USE_KAHAN  
      MPI_Allreduce(sum,centroid,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
      MPI_Allgather(sum,6,MPI_DOUBLE,&sumMatrix[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD);
      /* compute error sum over processors */
      MPI_Allgather(error,6,MPI_DOUBLE,&errorMatrix[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD);

      for (i_coord=0; i_coord<6; i_coord++) {
        error_sum = 0.0;
	/* extract the columnwise array from the matrix */
	for (j=0; j<n_processors; j++) {         
	  sumArray[j] = sumMatrix[j][i_coord];             
	  errorArray[j] = errorMatrix[j][i_coord];
	}
	centroid[i_coord] = Kahan(n_processors-1,&sumArray[1],&error_sum);
	error_total = Kahan(n_processors-1,&errorArray[1],&error_sum);
	centroid[i_coord] += error_total; 
      }
#endif 
      /* compute total number of particles over processors */
      MPI_Allreduce(&n_part, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);  
      if (n_total)
	for (i_coord=0; i_coord<6; i_coord++)
	  centroid[i_coord] /= n_total;

    }
  }
#endif  
   
#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
#ifdef USE_KAHAN
  free_czarray_2d((void**)sumMatrix, n_processors, 6); 
  free_czarray_2d((void**)errorMatrix, n_processors, 6);
  free(sumArray);
  free(errorArray);
#endif
#endif

}

/* routine: compute_sigmas()
 * purpose: compute sigmas for (x, x', y, y', s, deltap/p)
 *
 * Michael Borland, 1989
 */

void compute_sigmas(
    double *emit,
    double *sigma,
    double *centroid,
    double **coordinates,
    long n_part
    )
{
    long i_part, i_coord;
    double sum2[6], *part, value;
    long active = 1;
#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
  long n_total;
  double sum2_total[6];

  if (notSinglePart) {
    if (isMaster)
      n_part = 0;
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }
#endif

  if (active) {
    for (i_coord=0; i_coord<6; i_coord++)
        sum2[i_coord] = 0;
    if (emit)
      for (i_coord=0; i_coord<3;  i_coord++)
        emit[i_coord] = 0;
    
    for (i_part=0; i_part<n_part; i_part++) {
        part = coordinates[i_part];
        for (i_coord=0; i_coord<6; i_coord++) 
            sum2[i_coord] += sqr(part[i_coord]-centroid[i_coord]);
        }
#if USE_MPI
    if (notSinglePart) {
      /* compute total number of particles over processors */
      MPI_Allreduce(sum2, sum2_total, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&n_part, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);      
    }
    if (n_total) {
        for (i_coord=0; i_coord<6; i_coord++)
            if ((value=sum2_total[i_coord])>0)
                sigma[i_coord] = sqrt(value/n_total);
            else
                sigma[i_coord] = 0;
        if (emit) {
          for (i_coord=0; i_coord<6; i_coord+=2) {
            double sum12, sum12_total;
            sum12 = 0;
            for (i_part=0; i_part<n_part; i_part++) {
              part = coordinates[i_part];
              sum12 += (part[i_coord]-centroid[i_coord])*(part[i_coord+1]-centroid[i_coord+1]);
            }
	    MPI_Allreduce(&sum12, &sum12_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if ((emit[i_coord/2] = sqr(sigma[i_coord]*sigma[i_coord+1]) - sqr(sum12_total/n_total))>0)
              emit[i_coord/2] = sqrt(emit[i_coord/2]);
            else
              emit[i_coord/2] = 0;
          }
        }
#else
      if (n_part) {
        for (i_coord=0; i_coord<6; i_coord++)
            if ((value=sum2[i_coord])>0)
                sigma[i_coord] = sqrt(value/n_part);
            else
                sigma[i_coord] = 0;
        if (emit) {
          for (i_coord=0; i_coord<6; i_coord+=2) {
            double sum12;
            sum12 = 0;
            for (i_part=0; i_part<n_part; i_part++) {
              part = coordinates[i_part];
              sum12 += (part[i_coord]-centroid[i_coord])*(part[i_coord+1]-centroid[i_coord+1]);
            }
            if ((emit[i_coord/2] = sqr(sigma[i_coord]*sigma[i_coord+1]) - sqr(sum12/n_part))>0)
              emit[i_coord/2] = sqrt(emit[i_coord/2]);
            else
              emit[i_coord/2] = 0;
          }
        }
#endif
    }
  }
}

void zero_beam_sums(
                    BEAM_SUMS *sums,
                    long n
                    )
{
  long i, j, k;
  for (i=0; i<n; i++) {
    for (j=0; j<6; j++)
      sums[i].maxabs[j] = 0;
    for (j=0; j<6; j++)
      sums[i].min[j] = -(sums[i].max[j] = -DBL_MAX);
    for (j=0; j<6; j++)
      sums[i].centroid[j] = 0;
    for (j=0; j<6; j++)
      for (k=j; k<6; k++)
        sums[i].sigma[j][k] = 0;
    sums[i].n_part = sums[i].z = sums[i].p0 = 0;
  }
}

void accumulate_beam_sums(
                          BEAM_SUMS *sums,
                          double **coord,
                          long n_part,
                          double p_central
                          )
{
  long i_part, i, j;
  double centroid[6];
  double value;
  long active = 1;
  double Sij;
#ifdef USE_KAHAN
  double errorCen[6], errorSig[21];
#endif

#if USE_MPI  /* In the non-parallel mode, it will be same with the serial version */ 
  double buffer[6], Sij_p[21], Sij_total[21];
  long n_total, offset=0, index;  
#ifdef USE_KAHAN
  double error_sum=0.0, error_total=0.0,
    **sumMatrixCen, **errorMatrixCen,
    **sumMatrixSig, **errorMatrixSig,
    *sumArray, *errorArray;
  long k;
  sumMatrixCen = (double**)czarray_2d(sizeof(**sumMatrixCen), n_processors, 6); 
  errorMatrixCen = (double**)czarray_2d(sizeof(**errorMatrixCen), n_processors, 6);
  sumMatrixSig = (double**)czarray_2d(sizeof(**sumMatrixSig), n_processors, 21);
  errorMatrixSig = (double**)czarray_2d(sizeof(**errorMatrixSig), n_processors, 21);
  sumArray = malloc(sizeof(double) * n_processors);
  errorArray = malloc(sizeof(double) * n_processors);
#endif
  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }
#endif
    
  if (!sums->n_part)
    sums->p0 = p_central;

  if ((!USE_MPI && n_part) || USE_MPI) {
    if (active) {
      /* maximum amplitudes */
      for (i=0; i<6; i++) {
	if (i==4)
	  continue;  /* done below */
	for (i_part=0; i_part<n_part; i_part++) {
	  if ((value=fabs(coord[i_part][i]))>sums->maxabs[i])
	    sums->maxabs[i] = value;
	  if ((value=coord[i_part][i]) > sums->max[i])
	    sums->max[i] = value;
	  if ((value=coord[i_part][i]) < sums->min[i])
	    sums->min[i] = value;
	}
      }
      /* compute centroids for present beam and add in to existing centroid data */
      for (i=0; i<6; i++) {
#ifdef USE_KAHAN
        errorCen[i] = 0.0;
#endif
	for (centroid[i]=i_part=0; i_part<n_part; i_part++) {
#ifndef USE_KAHAN
	  centroid[i] += coord[i_part][i];
#else
	centroid[i] = KahanPlus(centroid[i], coord[i_part][i], &errorCen[i]); 
#endif
	}
	if (!USE_MPI || !notSinglePart) {	
	  sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+n_part);
	  centroid[i] /= n_part;
	}

#if USE_MPI
        if (notSinglePart) {
	  if ((parallelStatus!=trueParallel) && isMaster) {
	    if (n_part) {
	      sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+n_part);
	      centroid[i] /= n_part;
	    }
	  }
	}
#endif
      }
    }
#if USE_MPI
    if (notSinglePart) {
      if(parallelStatus==trueParallel) {
	if (isMaster) { 
	  n_part = 0;  /* All the particles have been distributed to the slave processors */
	  memset(centroid, 0.0,  sizeof(double)*6);
	}
	/* compute centroid sum over processors */
#ifndef USE_KAHAN 
	MPI_Allreduce(centroid, buffer, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	memcpy(centroid, buffer, sizeof(double)*6);
#else
	MPI_Allgather(centroid,6,MPI_DOUBLE,&sumMatrixCen[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD);
	/* compute error sum over processors */
      	MPI_Allgather(errorCen,6,MPI_DOUBLE,&errorMatrixCen[0][0],6,MPI_DOUBLE,MPI_COMM_WORLD); 
	for (i=0; i<6; i++) {
	  error_sum = 0.0;
          centroid[i] = 0.0;
          /* extract the columnwise array from the matrix */
          for (j=1; j<n_processors; j++) {
            centroid[i] = KahanPlus(centroid[i], sumMatrixCen[j][i], &error_sum);
            centroid[i] = KahanPlus(centroid[i], errorMatrixCen[j][i], &error_sum);
	  }
	}
#endif 
	/* compute total number of particles over processors */
	MPI_Allreduce(&n_part, &n_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	if (n_total)
	  for (i=0; i<6; i++) {
	    sums->centroid[i] = (sums->centroid[i]*sums->n_part+centroid[i])/(sums->n_part+n_total);
	    centroid[i] /= n_total;
	  }     
      }
      else
	if (isMaster) {
	  n_total = n_part;
	}
    }
#endif
    if (active) {	        
      i = 4;
      for (i_part=0; i_part<n_part; i_part++) {
        if ((value=fabs(coord[i_part][i]-centroid[i]))>sums->maxabs[i])
          sums->maxabs[i] = value;
        if ((value=coord[i_part][i]-centroid[i])>sums->max[i])
  	  sums->max[i] = value;
        if ((value=coord[i_part][i]-centroid[i])<sums->min[i])
	  sums->min[i] = value;
      }
    }
#if USE_MPI
    if (notSinglePart) {
      if (parallelStatus==trueParallel) {
	/* compute sums->maxabs over processors*/
	MPI_Allreduce(sums->maxabs, buffer, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
	memcpy(sums->maxabs, buffer, sizeof(double)*6);     
	/* compute sums->max over processors */
	MPI_Allreduce(sums->max, buffer, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
	memcpy(sums->max, buffer, sizeof(double)*6);      
	/* compute sums->min over processors */
	MPI_Allreduce(sums->min, buffer, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
	memcpy(sums->min, buffer, sizeof(double)*6);  
      }    
    }
#endif    
    /* compute Sigma[i][j] for present beam and add to existing data */
    if (active) {
      for (i=0; i<6; i++) {
#if USE_MPI
	if ((parallelStatus==trueParallel) && notSinglePart)
	  if (i>=1)        
	    offset += i-1;
#endif
	for (j=i; j<6; j++) {
#ifdef USE_KAHAN
          errorSig[j]=0.0;
#endif

#if (!USE_MPI) 
	  for (Sij=i_part=0; i_part<n_part; i_part++) {
#ifndef USE_KAHAN
	    Sij += (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]);
#else
	    Sij = KahanPlus(Sij, (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]), &errorSig[j]); 
#endif
	  }
	  sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sij)/(sums->n_part+n_part);
#else 
	  if (notSinglePart) {
	    if (parallelStatus==trueParallel) {
	      index = 5*i+j-offset;
	      Sij_p[index] = 0;
#ifdef USE_KAHAN
              errorSig[index] = 0.0;
#endif
	      for (i_part=0; i_part<n_part; i_part++) {
#ifndef USE_KAHAN
		Sij_p[index] += (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]);
#else
	        Sij_p[index] = KahanPlus(Sij_p[index], (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]), &errorSig[index]); 
#endif
	      }
	    }
	    else if (isMaster) {
	      for (Sij=i_part=0; i_part<n_part; i_part++) {
#ifndef USE_KAHAN
		Sij += (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]);
#else
		Sij = KahanPlus(Sij, (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]), &errorSig[j]); 
#endif
	      }
	      if (n_part)
		sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sij)/(sums->n_part+n_part);
	    }
	  }
	  else { /* Single particle case */
	    for (Sij=i_part=0; i_part<n_part; i_part++) {
#ifndef USE_KAHAN
	      Sij += (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]);
#else
	      Sij = KahanPlus(Sij, (coord[i_part][i]-centroid[i])*(coord[i_part][j]-centroid[j]), &errorSig[j]); 
#endif
	    }
	    if (n_part)
	      sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sij)/(sums->n_part+n_part);
	  }
#endif
	}
      }
    }
#if USE_MPI
    if (notSinglePart) {
      if (parallelStatus==trueParallel) {
	if (isMaster) {
	  memset(Sij_p, 0.0,  sizeof(double)*21);
#ifdef USE_KAHAN
          memset(errorSig, 0.0,  sizeof(double)*21);
#endif
	}
	/* compute Sij sum over processors */
#ifndef USE_KAHAN
	MPI_Allreduce(Sij_p, Sij_total, 21, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        MPI_Allgather(Sij_p, 21, MPI_DOUBLE, &sumMatrixSig[0][0], 21, MPI_DOUBLE, MPI_COMM_WORLD);
	/* compute error sum over processors */
	MPI_Allgather(errorSig, 21, MPI_DOUBLE, &errorMatrixSig[0][0], 21, MPI_DOUBLE,MPI_COMM_WORLD);
        offset = 0;
	for (i=0; i<6; i++) {
	  if (i>=1)        
	    offset += i-1;
	  for (j=i; j<6; j++) {
	    index = 5*i+j-offset;
	    error_sum = 0.0;
	    /* extract the columnwise array from the matrix */
	    for (k=0; k<n_processors; k++) {         
	      sumArray[k] = sumMatrixSig[k][index];             
	      errorArray[k] = errorMatrixSig[k][index];
	    }
	    Sij_total[index] = Kahan(n_processors-1,&sumArray[1],&error_sum);
	    error_total = Kahan(n_processors-1,&errorArray[1],&error_sum);
	    Sij_total[index] += error_total;
	  } 
	}
      
#endif
	if (n_total) {
	  offset = 0; 
	  for (i=0; i<6; i++) {
	    if (i>=1)
	      offset += i-1;
	    for (j=i; j<6; j++) {
	      index = 5*i+j-offset;
	      sums->sigma[i][j] = (sums->sigma[i][j]*sums->n_part+Sij_total[index])/(sums->n_part+n_total);
	    }
	  }
	}
      }
    }
#endif
  }

#if (!USE_MPI)    
  sums->n_part += n_part;
#else
  if (!notSinglePart)
    sums->n_part += n_part; 
  else if (!SDDS_MPI_IO) {
    if (parallelStatus==trueParallel) 
      sums->n_part += n_total;
    else if (isMaster)
      sums->n_part += n_part;
  } else {
    if (isMaster)
      sums->n_part += n_total;
    else
      sums->n_part += n_part;
  
  }
   
 
    
#endif    


#if USE_MPI
#ifdef USE_KAHAN
  free_czarray_2d((void**)sumMatrixCen, n_processors, 6); 
  free_czarray_2d((void**)errorMatrixCen, n_processors, 6);
  free_czarray_2d((void**)sumMatrixSig, n_processors, 21);
  free_czarray_2d((void**)errorMatrixSig, n_processors, 21);
  free(sumArray);
  free(errorArray);
#endif
#endif

}

void copy_beam_sums(
    BEAM_SUMS *target,
    BEAM_SUMS *source
    )
{
  long i, j;

  for (i=0; i<6; i++) {
    target->maxabs[i] = source->maxabs[i];
  }
  for (i=0; i<6; i++) {
    target->centroid[i] = source->centroid[i];
    for (j=i; j<6; j++)
      target->sigma[i][j] = source->sigma[i][j];
  }
  target->n_part = source->n_part;
  target->z      = source->z;
  target->p0     = source->p0;
}

long computeSliceMoments(double C[6], double S[6][6], 
			 double **part, long np, 
			 double sMin, double sMax)
{
  long i, j, k, count = 0;
  if (!part)
    bombElegant("NULL pointer passed to computeSliceMoments", NULL);
  for (j=0; j<6; j++) {
    C[j] = 0;
    for (k=0; k<6; k++)
      S[j][k] = 0;
  }
  if (!np)
    return 0;

  for (i=0; i<np; i++) {
    if (sMin!=sMax && (part[i][4]<sMin || part[i][4]>sMax))
      continue;
    count++;
    for (j=0; j<6; j++) {
      C[j] += part[i][j];
    }
  }
  if (!count)
    return 0;

  for (j=0; j<6; j++)
    C[j] /= count;

  for (i=0; i<np; i++) {
    if (sMin!=sMax && (part[i][4]<sMin || part[i][4]>sMax))
      continue;
    for (j=0; j<6; j++)
      for (k=0; k<=j; k++)
	S[j][k] += (part[i][j]-C[j])*(part[i][k]-C[k]);
  }

  for (j=0; j<6; j++)
    for (k=0; k<=j; k++) {
      S[j][k] /= count; 
      S[k][j] = S[j][k];
    }
  return count;
}

	 
