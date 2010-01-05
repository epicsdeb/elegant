/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge);
void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2);


void track_through_wake(double **part, long np, WAKE *wakeData, double *PoInput,
                        RUN *run, long i_pass, CHARGE *charge
                        )
{
  static double *Itime = NULL;           /* array for histogram of particle density */
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static long max_np = 0;
  static short shortBunchWarning = 0;
  long ib, nb, n_binned;
  double factor, tmin, tmax, tmean, dt, Po, rampFactor;
#if USE_MPI
  double *buffer;
#endif

  set_up_wake(wakeData, run, i_pass, np, charge);
  rampFactor = 0;
  if (i_pass>=(wakeData->rampPasses-1))
    rampFactor = 1;
  else
    rampFactor = (i_pass+1.0)/wakeData->rampPasses;
  Po = *PoInput;

  if (!USE_MPI || !notSinglePart) {
    if (np>max_np) {
      pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
      time = trealloc(time, sizeof(*time)*max_np);
    }
  }
#if USE_MPI
  else if (USE_MPI) {
      long np_total;
      if (isSlave) {
	MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	if (np_total>max_np) { 
	  /* if the total number of particles is increased, we do reallocation for every CPU */
	  pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
	  time = trealloc(time, sizeof(*time)*max_np);
	  max_np = np_total; /* max_np should be the sum across all the processors */
	}
      }
    }
#endif 

  /* Compute time coordinate of each particle */
  if (isSlave ||  !notSinglePart)
    tmean = computeTimeCoordinates(time, Po, part, np);

  find_min_max(&tmin, &tmax, time, np);
#if USE_MPI
  if (isSlave && notSinglePart)
    find_global_min_max(&tmin, &tmax, np, workers);      
#endif
  if (isSlave ||  !notSinglePart) {
    if ((tmax-tmin) > (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0])) {
      if (!wakeData->allowLongBeam) {
	fprintf(stderr, "Error: The beam is longer than the longitudinal wake function.\nThis may produce unphysical results.\n");
	fprintf(stderr, "The beam length is %le s, while the wake length is %le s\n",
		tmax-tmin, wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
	exit(1);
      }
      fprintf(stdout, "Warning: The beam is longer than the longitudinal wake function.\nThis may produce unphysical results.\n");
      fprintf(stdout, "The beam length is %le s, while the wake length is %le s\n",
	      tmax-tmin, wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
      if (abs(tmax-tmean)<abs(tmin-tmean)) 
	tmin = tmax - (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
      else
	tmax = tmin + (wakeData->t[wakeData->wakePoints-1]-wakeData->t[0]);
    }

    dt = wakeData->dt;
    if (np>1 && (tmax-tmin)<20*dt && !shortBunchWarning) {
      fprintf(stdout, "Warning: The beam is shorter than 20*DT, where DT is the spacing of the wake points.\n");
      fprintf(stdout, "         Depending on the longitudinal distribution and shape of the wake, this may produce poor results.\n");
      fprintf(stdout, "         Consider using a wake with finer time spacing in WAKE elements.\n");
      fflush(stdout);
      shortBunchWarning = 1;
    }
    
    if (wakeData->n_bins) {
      nb = wakeData->n_bins;
      tmin = tmean-dt*nb/2.0;
    }
    else {
      nb = (tmax-tmin)/dt+3;
      tmin -= dt;
      tmax += dt;
    }
    if (nb<=0) {
      fprintf(stdout, "Warning: Number of wake bins is 0 or negative\n");
      fprintf(stderr, "probably indicating an extremely long bunch\n");
      fprintf(stderr, "Wake ignored!\n");
      return;
    }

    if (nb>max_n_bins) {
      Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins=nb));
      Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
    }
 
    n_binned = binTimeDistribution(Itime, pbin, tmin, dt, nb, time, part, Po, np);
  }

  if (!USE_MPI || !notSinglePart) {
    if (n_binned!=np) {
      fprintf(stdout, "warning: only %ld of %ld particles were binned (WAKE)\n", n_binned, np);
      fprintf(stdout, "consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
      fflush(stdout);
      return;
    }
  }
#if USE_MPI
  else if (isSlave) {
    int all_binned, result = 1;

    result = ((n_binned==np) ? 1 : 0);		             
    MPI_Allreduce(&result, &all_binned, 1, MPI_INT, MPI_LAND, workers);
    if (!all_binned) {
      if (myid==1) {  
	/* This warning will be given only if the flag MPI_DEBUG is defined for the Pelegant */ 
	fprintf(stdout, "warning: Not all of %ld particles were binned (WAKE)\n", np);
	fprintf(stdout, "consider setting n_bins=0 in WAKE definition to invoke autoscaling\n");
	fflush(stdout); 
      }
    }
  }
  
  if (isSlave && notSinglePart) {
    buffer = malloc(sizeof(double) * nb);
    MPI_Allreduce(Itime, buffer, nb, MPI_DOUBLE, MPI_SUM, workers);
    memcpy(Itime, buffer, sizeof(double)*nb);
    free(buffer);
  }
#endif
  if (isSlave || !notSinglePart) {
    if (wakeData->smoothing && nb>=(2*wakeData->SGHalfWidth+1)) {
      if (!SavitzyGolaySmooth(Itime, nb, wakeData->SGOrder, wakeData->SGHalfWidth, wakeData->SGHalfWidth, 0)) {
	fprintf(stderr, "Problem with smoothing for WAKE element (file %s)\n",
		wakeData->inputFile);
	fprintf(stderr, "Parameters: nbins=%ld, order=%ld, half-width=%ld\n",
		nb, wakeData->SGOrder, wakeData->SGHalfWidth);
	exit(1);
      }
    }
  }

  /* Do the convolution of the particle density and the wake function,
     V(T) = Integral[W(T-t)*I(t)dt, t={-infinity, T}]
     Note that T<0 is the head of the bunch.
     For the wake, the argument is the normal convention wherein larger
     arguments are later times.
     */
  if (isSlave || !notSinglePart) {
    Vtime[nb] = 0;
    convolveArrays(Vtime, nb,
		   Itime, nb, 
		   wakeData->W, wakeData->wakePoints);

    factor = wakeData->macroParticleCharge*particleRelSign*wakeData->factor*rampFactor;
    for (ib=0; ib<nb; ib++)
      Vtime[ib] *= factor;
  
    applyLongitudinalWakeKicks(part, time, pbin, np, Po, 
			       Vtime, nb, tmin, dt, wakeData->interpolate);
  }

  if (wakeData->change_p0)
    do_match_energy(part, np, PoInput, 0);
  
#if defined(MINIMIZE_MEMORY)
  free(Itime);
  free(Vtime);
  free(pbin);
  free(time);
  Itime = Vtime = time = NULL;
  pbin = NULL;
  max_n_bins =  max_np= 0;
#endif
}

void applyLongitudinalWakeKicks(double **part, double *time, long *pbin, long np, double Po,
                                double *Vtime, long nb, double tmin, double dt,
                                long interpolate)
{
  long ip, ib;
  double dt1, dgam;
  
  /* change particle momentum offsets to reflect voltage in relevant bin */
  for (ip=0; ip<np; ip++) {
    if ((ib=pbin[ip])>=0 && ib<=nb-1) {
      if (interpolate) {
        dt1 = time[ip]-(tmin+dt*ib);  /* distance to bin center */
        if ((dt1<0 && ib) || ib==nb-1) {
          ib--;
          dt1 += dt;
        }
        dgam = (Vtime[ib]+(Vtime[ib+1]-Vtime[ib])/dt*dt1)/(1e6*particleMassMV*particleRelSign);
      }
      else
        dgam = Vtime[ib]/(1e6*particleMassMV*particleRelSign);
      if (dgam) {
        /* Put in minus sign here as the voltage decelerates the beam */
	add_to_particle_energy(part[ip], time[ip], Po, -dgam); 
      }
    }
  }
}

typedef struct {
  char *filename;
  long points;
  double *t, *W;
} WAKE_DATA;

static WAKE_DATA *storedWake = NULL;
static long storedWakes = 0;

void set_up_wake(WAKE *wakeData, RUN *run, long pass, long particles, CHARGE *charge)
{
  SDDS_DATASET SDDSin;
  double tmin, tmax;
  long iw;
#if SDDS_MPI_IO 
/* All the processes will read the wake file, but not in parallel.
   Zero the Memory when call  SDDS_InitializeInput */
  SDDSin.parallel_io = 0; 
#endif
  
  if (charge) {
    wakeData->macroParticleCharge = charge->macroParticleCharge;
  } else if (pass==0) {
    wakeData->macroParticleCharge = 0;
    if (wakeData->charge<0)
      bomb("WAKE charge parameter should be non-negative.  Use change_particle to set particle charge state.", NULL);
#if (!USE_MPI)
    if (particles)
      wakeData->macroParticleCharge = wakeData->charge/particles;
#else
    if (USE_MPI) {
      long particles_total;
      if (isSlave) {
	MPI_Allreduce(&particles, &particles_total, 1, MPI_LONG, MPI_SUM, workers);
	if (particles_total)
	  wakeData->macroParticleCharge = wakeData->charge/particles_total;  
      }
    } 
#endif
  }
  
  if (wakeData->initialized)
    return;
  wakeData->initialized = 1;
  wakeData->W = wakeData->t = NULL;

  if (wakeData->n_bins<2 && wakeData->n_bins!=0)
    bomb("n_bins must be >=2 or else 0 (autoscale) for WAKE element", NULL);

  if (!wakeData->inputFile || !strlen(wakeData->inputFile))
    bomb("supply inputFile for WAKE element", NULL);
  if (!wakeData->tColumn || !strlen(wakeData->tColumn))
    bomb("supply tColumn for WAKE element", NULL);
  if (!wakeData->WColumn || !strlen(wakeData->WColumn))
    bomb("supply WColumn for WAKE element", NULL);
  
  for (iw=0; iw<storedWakes; iw++) {
    if (strcmp(storedWake[iw].filename, wakeData->inputFile)==0)
      break;
  }
  
  if (iw==storedWakes) {
    /* read in a new wake */
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, wakeData->inputFile) || SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stderr, "Error: unable to open or read WAKE file %s\n", wakeData->inputFile);
      exit(1);
    }
    if ((wakeData->wakePoints=SDDS_RowCount(&SDDSin))<0) {
      fprintf(stderr, "Error: no data in WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (wakeData->wakePoints<2) {
      fprintf(stderr, "Error: too little data in WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDSin, wakeData->tColumn, "s", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK) {
      fprintf(stderr, "Error: problem with time column in WAKE file %s.  Check existence, type, and units.\n",  wakeData->inputFile);
      exit(1);
    }
    if (!(wakeData->t=SDDS_GetColumnInDoubles(&SDDSin, wakeData->tColumn))) {
      fprintf(stderr, "Error: problem retrieving time data from WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    if (SDDS_CheckColumn(&SDDSin, wakeData->WColumn, "V/C", SDDS_ANY_FLOATING_TYPE, 
                         stdout)!=SDDS_CHECK_OK) {
      fprintf(stderr, "Error: problem with wake column in WAKE file %s.  Check existence, type, and units.\n",  wakeData->inputFile);
      exit(1);
    }
    if (!(wakeData->W=SDDS_GetColumnInDoubles(&SDDSin, wakeData->WColumn))) {
      fprintf(stderr, "Error: problem retrieving wake data from WAKE file %s\n",  wakeData->inputFile);
      exit(1);
    }
    SDDS_Terminate(&SDDSin);

    /* record in wake storage */
    if (!(storedWake=SDDS_Realloc(storedWake, sizeof(*storedWake)*(storedWakes+1))) || 
        !SDDS_CopyString(&storedWake[storedWakes].filename, wakeData->inputFile))
      SDDS_Bomb("Memory allocation failure (WAKE)");
    storedWake[storedWakes].t = wakeData->t;
    storedWake[storedWakes].W = wakeData->W;
    storedWake[storedWakes].points = wakeData->wakePoints;
    wakeData->isCopy = 0;
    storedWakes++;
  }
  else {
    /* point to an existing wake */
    wakeData->t = storedWake[iw].t;
    wakeData->W = storedWake[iw].W;
    wakeData->wakePoints = storedWake[iw].points;
    wakeData->isCopy = 1;
  }
  find_min_max(&tmin, &tmax, wakeData->t, wakeData->wakePoints);
#if USE_MPI
  if (isSlave)
    find_global_min_max(&tmin, &tmax, wakeData->wakePoints, workers);      
#endif
  if (tmin>=tmax) {
    fprintf(stderr, "Error: zero or negative time span in WAKE file %s\n",  wakeData->inputFile);
    exit(1);
  }
  if (tmin!=0) {
    fprintf(stderr, "Error: WAKE function does not start at t=0 for file %s\n",  wakeData->inputFile);
    exit(1);
  }
  wakeData->dt = (tmax-tmin)/(wakeData->wakePoints-1);
}

void convolveArrays(double *output, long outputs, 
                    double *a1, long n1,
                    double *a2, long n2)
{
  long ib, ib1, ib2;
  for (ib=0; ib<outputs; ib++) {
    output[ib] = 0;
    ib2 = ib;
    ib1 = 0;
    if (ib2>=n2) {
      ib2 = n2-1;
      ib1 = ib-ib2;
      if (ib1>=n1)
        continue;
    }
    for (; ib1<=ib; ib1++, ib2--)
      output[ib] += a1[ib1]*a2[ib2];
  }
}

long binTimeDistribution(double *Itime, long *pbin, double tmin,
                         double dt, long nb, double *time, double **part, double Po, long np)
{
  long ib, ip, n_binned;
  
  for (ib=0; ib<nb; ib++)
    Itime[ib] = 0;

  for (ip=n_binned=0; ip<np; ip++) {
    pbin[ip] = -1;
    /* Bin CENTERS are at tmin+ib*dt */
    ib = (time[ip]-tmin)/dt+0.5;
    if (ib<0)
      continue;
    if (ib>nb - 1)
      continue;
    Itime[ib] += 1;
    pbin[ip] = ib;
    n_binned++;
  }
  return n_binned;
}
