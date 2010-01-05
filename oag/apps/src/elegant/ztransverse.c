/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: zlongit.c
 * contents: track_through_zlongit()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

#define WAKE_COLUMNS 5
static SDDS_DEFINITION wake_column[WAKE_COLUMNS] = {
    {"Deltat", "&column name=Deltat, symbol=\"$gD$rt\", units=s, type=double, description=\"Time after head of bunch\" &end"},
    {"Ix", "&column name=Ix, symbol=\"<I*x>\", units=C*m/s, type=double, description=\"Transverse horizontal moment\" &end"},
    {"Wx", "&column name=Wx, symbol=\"W$bx$n\", units=V/m, type=double, description=\"Transverse horizontal wake\" &end"},
    {"Iy", "&column name=Iy, symbol=\"<I*y>\", units=C*m/s, type=double, description=\"Transverse vertical moment\" &end"},
    {"Wy", "&column name=Wy, symbol=\"W$by$n\", units=V/m, type=double, description=\"Transverse vertical wake\" &end"},
    };

#define WAKE_PARAMETERS 5
#define BB_WAKE_PARAMETERS 5
#define NBB_WAKE_PARAMETERS 2
static SDDS_DEFINITION wake_parameter[WAKE_PARAMETERS] = {
    {"Pass", "&parameter name=Pass, type=long &end"},
    {"q", "&parameter name=q, units=C, type=double, description=\"Total charge\" &end"},
    {"Rs", "&parameter name=Rs, symbol=\"R$bs$n\", units=\"$gW$r\", type=double, description=\"Broad-band impedance\" &end"},
    {"fo", "&parameter name=fo, symbol=\"f$bo$n\", units=Hz, type=double, description=\"Frequency of BB resonator\" &end"},
    {"Deltaf", "&parameter name=Deltaf, symbol=\"$gD$rf\", units=Hz, type=double, description=\"Frequency sampling interval\" &end"},
    } ;


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge,
                        double timeSpan);
double *getTransverseImpedance(SDDS_DATASET *SDDSin, char *ZName);

void track_through_ztransverse(double **part, long np, ZTRANSVERSE *ztransverse, double Po,
                               RUN *run, long i_pass, CHARGE *charge
                               )
{
  static double *posItime[2] = {NULL, NULL};     /* array for particle density times x, y*/
  static double *posIfreq;                       /* array for FFT of particle density times x or y*/
  static double *Vtime = NULL;           /* array for voltage acting on each bin */
  static long max_n_bins = 0;
  static long *pbin = NULL;              /* array to record which bin each particle is in */
  static double *time = NULL;            /* array to record arrival time of each particle */
  static double *pz = NULL;
  static long max_np = 0;
  double *Vfreq, *iZ;
  long i, ib, nb, n_binned, nfreq, iReal, iImag, plane, first;
  double factor, tmin, tmax, tmean, dt, userFactor[2], rampFactor=1;
  static long not_first_call = -1;
  long ip, ip1, ip2, bunches, bunch, npb;
  long bucketEnd[MAX_BUCKETS];
#if USE_MPI
  float *buffer_send, *buffer_recv;
#endif
#if defined(DEBUG)
  FILE *fp;
#endif
  if ((i_pass -= ztransverse->startOnPass)<0)
    return;

  if (i_pass>=(ztransverse->rampPasses-1))
    rampFactor = 1;
  else
    rampFactor = (i_pass+1.0)/ztransverse->rampPasses;
  
  not_first_call += 1;

#if (!USE_MPI)
  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
    pz = trealloc(pz, sizeof(*pz)*max_np);
  }
#else
    if (USE_MPI) {
      long np_total;
      if (isSlave) {
	MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	if (np_total>max_np) { 
	  /* if the total number of particles is increased, we do reallocation for every CPU */
	  pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
	  time = trealloc(time, sizeof(*time)*max_np);
	  pz = trealloc(pz, sizeof(*pz)*max_np);
	  max_np = np_total; /* max_np should be the sum across all the processors */
	}
      }
    } 
#endif
  
  if ((part[0][6]-floor(part[0][6]))!=0) {
    /* This is a kludgey way to determine that particles have been assigned to buckets */
    printf("Bunched beam detected\n"); fflush(stdout);
#if USE_MPI
    if (n_processors!=1) {
      printf("Error (ZTRANSVERSE): must have bunch_frequency=0 for parallel mode.\n");
      MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
#endif
    /* Start by sorting in bucket order */
    qsort(part[0], np, COORDINATES_PER_PARTICLE*sizeof(double), comp_BucketNumbers);
    /* Find the end of the buckets in the particle list */
    bunches = 0;
    for (ip=1; ip<np; ip++) {
      if ((part[ip-1][6]-floor(part[ip-1][6]))!=(part[ip][6]-floor(part[ip][6]))) {
        /* printf("Bucket %ld ends with ip=%ld\n", bunches, ip-1); fflush(stdout); */
        bucketEnd[bunches++] = ip-1;
        if (bunches>=MAX_BUCKETS) {
            bomb("Error (wake): maximum number of buckets was exceeded", NULL);
          }
      }
    }
    bucketEnd[bunches++] = np-1;
  } else {
    bunches = 1;
    bucketEnd[0] = np-1;
  }
  /* printf("Bucket %ld ends with ip=%ld\n", bunches, bucketEnd[bunches-1]); fflush(stdout); */

  ip2 = -1;
  for (bunch=0; bunch<bunches; bunch++) {
    ip1 = ip2+1;
    ip2 = bucketEnd[bunch];
    npb = ip2-ip1+1;
    /* printf("Processing bunch %ld with %ld particles\n", bunch, npb); fflush(stdout); */
    /* Compute time coordinate of each particle */
    tmean = computeTimeCoordinates(time+ip1, Po, part+ip1, npb);
    find_min_max(&tmin, &tmax, time+ip1, npb);
#if USE_MPI
    find_global_min_max(&tmin, &tmax, np, workers);      
#endif
    if (bunch==0) {
      /* use np here since we need to compute the macroparticle charge */
      set_up_ztransverse(ztransverse, run, i_pass, np, charge, tmax-tmin);
    }
    nb = ztransverse->n_bins;
    dt = ztransverse->bin_size;
    tmin -= dt;
    tmax -= dt;
    if ((tmax-tmin)*2>nb*dt) {
      TRACKING_CONTEXT tcontext;
      getTrackingContext(&tcontext);
      fprintf(stderr, "%s %s: Time span of bunch (%le s) is more than half the total time span (%le s).\n",
              entity_name[tcontext.elementType],
              tcontext.elementName,
              tmax-tmin, nb*dt);
      fprintf(stderr, "If using broad-band impedance, you should increase the number of bins and rerun.\n");
      fprintf(stderr, "If using file-based impedance, you should increase the number of data points or decrease the frequency resolution.\n");
      exit(1);
    }

    if (nb>max_n_bins) {
      posItime[0] = trealloc(posItime[0], 2*sizeof(**posItime)*(max_n_bins=nb));
      posItime[1] = trealloc(posItime[1], 2*sizeof(**posItime)*(max_n_bins=nb));
      posIfreq = trealloc(posIfreq, 2*sizeof(*posIfreq)*(max_n_bins=nb));
      Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
    }

    for (ib=0; ib<nb; ib++)
      posItime[0][2*ib] = posItime[0][2*ib+1] = 
        posItime[1][2*ib] = posItime[1][2*ib+1] = 0;

    /* make arrays of I(t)*x and I(t)*y */
    n_binned = binTransverseTimeDistribution(posItime, pz, pbin, tmin, dt, nb, 
                                             time+ip1, part+ip1, Po, npb,
                                             ztransverse->dx, ztransverse->dy, 1, 1);
#if (!USE_MPI)
    if (n_binned!=npb) {
      fprintf(stdout, "Warning: only %ld of %ld particles were binned (ZTRANSVERSE)!\n", n_binned, npb);
      if (!not_first_call) {
        fprintf(stdout, "*** This may produce unphysical results.  Your wake needs smaller frequency\n");
        fprintf(stdout, "    spacing to cover a longer time span.\n");
      }
      fflush(stdout);
    }  
#else
    if (USE_MPI) {
      int all_binned, result = 1;
      if (isSlave)
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
#endif

    userFactor[0] = ztransverse->factor*ztransverse->xfactor*rampFactor;
    userFactor[1] = ztransverse->factor*ztransverse->yfactor*rampFactor;

    first = 1;
    for (plane=0; plane<2; plane++) {
#if USE_MPI
      /* This could be good for some advanced network structures */
      /*
        if (isSlave) {
        double buffer[nb];
        MPI_Allreduce(posItime[plane], buffer, nb, MPI_DOUBLE, MPI_SUM, workers);
        memcpy(posItime[plane], buffer, sizeof(double)*nb);
        }
        */
      if (isSlave) {
        buffer_send = malloc(sizeof(float) * nb);
        buffer_recv = malloc(sizeof(float) * nb);
        for (i=0; i<nb; i++)
          buffer_send[i] = posItime[plane][i];
        MPI_Reduce(buffer_send, buffer_recv, nb, MPI_FLOAT, MPI_SUM, 1, workers);
        MPI_Bcast(buffer_recv, nb, MPI_FLOAT, 1, workers);
        for (i=0; i<nb; i++)
          posItime[plane][i] = buffer_recv[i];
        free(buffer_send);
        free(buffer_recv);
      }

#endif
      if (userFactor[plane]==0) {
        for (ib=0; ib<nb; ib++)
          Vtime[ib] = 0;
      } else {
        if (ztransverse->smoothing)
          SavitzyGolaySmooth(posItime[plane], nb, ztransverse->SGOrder,
                             ztransverse->SGHalfWidth, ztransverse->SGHalfWidth, 0);
        
        /* Take the FFT of (x*I)(t) to get (x*I)(f) */
        memcpy(posIfreq, posItime[plane], 2*nb*sizeof(*posIfreq));
        realFFT(posIfreq, nb, 0);
        
        /* Compute V(f) = i*Z(f)*(x*I)(f), putting in a factor 
         * to normalize the current waveform
         */
        Vfreq = Vtime;
        factor = ztransverse->macroParticleCharge*particleRelSign/dt*userFactor[plane];
        iZ = ztransverse->iZ[plane];
        Vfreq[0] = posIfreq[0]*iZ[0]*factor;
        nfreq = nb/2 + 1;
        if (nb%2==0)
          /* Nyquist term */
          Vfreq[nb-1] = posIfreq[nb-1]*iZ[nb-1]*factor;
        for (ib=1; ib<nfreq-1; ib++) {
          iImag = (iReal = 2*ib-1)+1;
          /* The signs are chosen here to get agreement with TRFMODE.
             In particular, test particles following closely behind the 
             drive particle get defocused.
             */
          Vfreq[iReal] =  (posIfreq[iReal]*iZ[iImag] + posIfreq[iImag]*iZ[iReal])*factor; 
          Vfreq[iImag] = -(posIfreq[iReal]*iZ[iReal] - posIfreq[iImag]*iZ[iImag])*factor;
        }
        
        /* Compute inverse FFT of V(f) to get V(t) */
        realFFT(Vfreq, nb, INVERSE_FFT);
        Vtime = Vfreq;
        
        /* change particle transverse momenta to reflect voltage in relevant bin */
        applyTransverseWakeKicks(part+ip1, time+ip1, pz, pbin, npb, 
                                 Po, plane, 
                                 Vtime, nb, tmin, dt, ztransverse->interpolate);
      }

      if (ztransverse->SDDS_wake_initialized && ztransverse->wakes) {
        /* wake potential output */
        factor = ztransverse->macroParticleCharge*particleRelSign/dt;
        if (ztransverse->wake_interval<=0 || (i_pass%ztransverse->wake_interval)==0) {
          if (first && !SDDS_StartTable(&ztransverse->SDDS_wake, nb)) {
            SDDS_SetError("Problem starting SDDS table for wake output (track_through_ztransverse)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          for (ib=0; ib<nb; ib++) {
            if (!SDDS_SetRowValues(&ztransverse->SDDS_wake, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ib,
                                   0, ib*dt, 
                                   1+plane*2, posItime[plane][ib]*factor,  
                                   2+plane*2, Vtime[ib], -1)) {
              SDDS_SetError("Problem setting rows of SDDS table for wake output (track_through_ztransverse)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (!SDDS_SetParameters(&ztransverse->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                  "Pass", i_pass, "q", ztransverse->macroParticleCharge*particleRelSign*np, NULL)) {
            SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_ztransverse)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          if (ztransverse->broad_band) {
            if (!SDDS_SetParameters(&ztransverse->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                    "Rs", ztransverse->Rs, "fo", ztransverse->freq, 
                                    "Deltaf", ztransverse->bin_size, NULL)) {
              SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_ztransverse)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (!first) {
            if (!SDDS_WriteTable(&ztransverse->SDDS_wake)) {
              SDDS_SetError("Problem writing SDDS table for wake output (track_through_ztransverse)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
            if (!inhibitFileSync)
              SDDS_DoFSync(&ztransverse->SDDS_wake);
          }
        }
      }
      first = 0;
    }
  }
  
#if defined(MIMIMIZE_MEMORY)
  free(posItime[0]);
  free(posItime[1]);
  free(posIfreq);
  free(Vtime);
  free(pbin);
  free(time);
  free(pz);
  posItime[0] = posItime[1] = Vtime = time = pz = posIfreq = NULL;
  pbin = NULL;
  max_n_bins = max_np = 0 ;
#endif

}


void set_up_ztransverse(ZTRANSVERSE *ztransverse, RUN *run, long pass, long particles, CHARGE *charge,
                        double timeSpan)
{
  long i, nfreq;
  double df;

  if (charge) {
    ztransverse->macroParticleCharge = charge->macroParticleCharge;
  } else if (pass==0) {
    ztransverse->macroParticleCharge = 0;
    if (ztransverse->charge<0)
      bomb("ZTRANSVERSE charge parameter should be non-negative. Use change_particle to set particle charge state.", NULL);
#if (!USE_MPI)
    if (particles)
      ztransverse->macroParticleCharge = ztransverse->charge/particles;
#else
      if (USE_MPI) {
	long particles_total;

	MPI_Allreduce(&particles, &particles_total, 1, MPI_LONG, MPI_SUM, workers);
	if (particles_total)
	  ztransverse->macroParticleCharge = ztransverse->charge/particles_total;  
      } 
#endif
  }

  if (ztransverse->initialized)
    return;

  ztransverse->SDDS_wake_initialized = 0;

  if (ztransverse->broad_band) {
    /* Use impedance Z = -i*wr/w*Rs/(1 + i*Q(w/wr-wr/w))
       */
    double term;
    if (ztransverse->bin_size<=0)
      bomb("bin_size must be positive for ZTRANSVERSE element", NULL);
    if (ztransverse->n_bins%2!=0)
      bomb("ZTRANSVERSE element must have n_bins divisible by 2", NULL);
    if (ztransverse->ZxReal || ztransverse->ZxImag ||
        ztransverse->ZyReal || ztransverse->ZyImag )
      bomb("can't specify both broad_band impedance and Z(f) files for ZTRANSVERSE element", NULL);

    optimizeBinSettingsForImpedance(timeSpan, ztransverse->freq, ztransverse->Q,
                                    &(ztransverse->bin_size), &(ztransverse->n_bins),
                                    ztransverse->max_n_bins);
    
    nfreq = ztransverse->n_bins/2 + 1;
    ztransverse->iZ[0] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    ztransverse->iZ[1] = tmalloc(sizeof(**(ztransverse->iZ))*ztransverse->n_bins);
    /* df is the frequency spacing normalized to the resonant frequency */
    df = 1/(ztransverse->n_bins*ztransverse->bin_size)/(ztransverse->freq);
    /* DC term of iZ is 0  */
    ztransverse->iZ[0][0] = ztransverse->iZ[1][0] = 0;
    for (i=1; i<nfreq-1; i++) {
      term = ztransverse->Q*(i*df-1.0/(i*df));
      /* real part of i*Z */
      ztransverse->iZ[0][2*i-1] =  
        ztransverse->iZ[1][2*i-1] =  
	ztransverse->Rs/(i*df)/(1+term*term);
      /* imaginary part of i*Z is -Real[i*Z]*term */
      ztransverse->iZ[0][2*i] = 
        ztransverse->iZ[1][2*i] = 
	-term*ztransverse->iZ[0][2*i-1];
    }
    /* Nyquist term--real part of iZ only */
    term = ztransverse->Q*(1.0/(nfreq*df)-nfreq*df);
    ztransverse->iZ[0][ztransverse->n_bins-1] = 
      ztransverse->iZ[1][ztransverse->n_bins-1] = 
      ztransverse->Rs/(nfreq*df)/(1+term*term);
    df *= ztransverse->freq;
  } else {
    double *ZReal[2], *ZImag[2], *freqData;
    double df_spect;
    long n_spect;
    SDDS_DATASET SDDSin;
    if (!ztransverse->freqColumn || !ztransverse->inputFile)
      bomb("you must give an inputFile and freqColumn, or use a broad band model (ZTRANSVERSE)", NULL);
    if (!ztransverse->ZxReal && !ztransverse->ZxImag &&
        !ztransverse->ZyReal && !ztransverse->ZxImag)
      bomb("you must either give broad_band=1, or Z[xy]Real and/or Z[xy]Imag (ZTRANSVERSE)", NULL);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, ztransverse->inputFile) || !SDDS_ReadPage(&SDDSin)) {
      fprintf(stdout, "unable to read file %s\n", ztransverse->inputFile);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors); 
      exit(1);
    }
    if ((n_spect=SDDS_RowCount(&SDDSin))<4) {
      fprintf(stdout, "too little data in %s\n", ztransverse->inputFile);
      fflush(stdout);
      exit(1);
    }
    if (!power_of_2(n_spect-1))
      bomb("number of spectrum points must be 2^n+1, n>1 (ZTRANSVERSE)", NULL);
    ZReal[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxReal);
    ZImag[0] = getTransverseImpedance(&SDDSin, ztransverse->ZxImag);
    ZReal[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyReal);
    ZImag[1] = getTransverseImpedance(&SDDSin, ztransverse->ZyImag);
    if (!(freqData=SDDS_GetColumnInDoubles(&SDDSin, ztransverse->freqColumn))) {
      fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n", ztransverse->freqColumn);
      fflush(stdout);
      exit(1);
    }
    if (!checkPointSpacing(freqData, n_spect, 1e-6)) {
      fprintf(stdout, "Frequency values are not equispaced (ZTRANSVERSE)\n");
      fflush(stdout);
      exit(1);
    }
    if ((df_spect = (freqData[n_spect-1]-freqData[0])/(n_spect-1))<=0) {
      fprintf(stdout, "Zero or negative frequency spacing in %s (ZTRANSVERSE)\n",
              ztransverse->inputFile);
      fflush(stdout);
      exit(1);
    }
    df = df_spect;
    nfreq = n_spect;
    ztransverse->n_bins = 2*(n_spect-1);
    ztransverse->bin_size = 1.0/(ztransverse->n_bins*df_spect);
    if (!SDDS_Terminate(&SDDSin)) {
      fprintf(stdout, "Error closing data set %s\n",
              ztransverse->inputFile);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!(ztransverse->iZ[0] =
          calloc(sizeof(*ztransverse->iZ[0]), n_spect*2)) ||
        !(ztransverse->iZ[1] =
          calloc(sizeof(*ztransverse->iZ[1]), n_spect*2)))
      bomb("memory allocation failure (ZTRANSVERSE)", NULL);
    for (i=0; i<n_spect; i++) {
      if (i==0) {
        /* DC term */
        ztransverse->iZ[0][i] = -ZImag[0][i];
        ztransverse->iZ[1][i] = -ZImag[1][i];
      } else if (i==n_spect-1 && ztransverse->n_bins%2==0) {
        /* Nyquist term */
        ztransverse->iZ[0][2*i-1] = -ZImag[0][i];
        ztransverse->iZ[1][2*i-1] = -ZImag[1][i];
      } else {
        /* real part of iZ */
        ztransverse->iZ[0][2*i-1] = -ZImag[0][i];
        ztransverse->iZ[1][2*i-1] = -ZImag[1][i];
        /* imaginary part of iZ */
        ztransverse->iZ[0][2*i  ] = ZReal[0][i];
        ztransverse->iZ[1][2*i  ] = ZReal[1][i];
      }
    }
    free(ZReal[0]);
    free(ZReal[1]);
    free(ZImag[0]);
    free(ZImag[1]);
  }

#if (!USE_MPI)
  /* Only the serial version will dump this part of output */
  if (ztransverse->wakes) {
    ztransverse->wakes = compose_filename(ztransverse->wakes, run->rootname);
    if (ztransverse->broad_band) 
      SDDS_ElegantOutputSetup(&ztransverse->SDDS_wake, ztransverse->wakes, SDDS_BINARY, 
			      1, "transverse wake",
			      run->runfile, run->lattice, wake_parameter, BB_WAKE_PARAMETERS,
			      wake_column, WAKE_COLUMNS, "set_up_ztransverse", 
			      SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    else {
      SDDS_ElegantOutputSetup(&ztransverse->SDDS_wake, ztransverse->wakes, SDDS_BINARY, 
			      1, "transverse wake",
			      run->runfile, run->lattice, wake_parameter, NBB_WAKE_PARAMETERS,
			      wake_column, WAKE_COLUMNS, "set_up_ztransverse", 
			      SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    }
    ztransverse->SDDS_wake_initialized = 1;
  }
#endif

  if (ztransverse->highFrequencyCutoff0>0) {
    applyLowPassFilterToImpedance(ztransverse->iZ[0], nfreq,
                                  ztransverse->highFrequencyCutoff0, 
                                  ztransverse->highFrequencyCutoff1);
    applyLowPassFilterToImpedance(ztransverse->iZ[1], nfreq,
                                  ztransverse->highFrequencyCutoff0, 
                                  ztransverse->highFrequencyCutoff1);
  }

#if 0
  if (!ztransverse->initialized) {
    FILE *fp;
    fp = fopen_e("ztransverse.sdds", "w", 0);
    fprintf(fp, "SDDS1\n&column name=f units=Hz type=double &end\n");
    fprintf(fp, "&column name=ZReal type=double &end\n");
    fprintf(fp, "&column name=ZImag type=double &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
    for (i=0; i<nfreq; i++) 
      fprintf(fp, "%21.15e %21.15e %21.15e\n",
              i*df, ztransverse->iZ[0][2*i], i>0?-ztransverse->iZ[0][2*i-1]:0);
    fclose(fp);
  }
#endif
  
  ztransverse->initialized = 1;
}

double *getTransverseImpedance(SDDS_DATASET *SDDSin, 
                               char *ZName)
{
  long rows;
  double *Z;
  rows = SDDS_RowCount(SDDSin);
  if (!ZName || !strlen(ZName)) {
    Z = calloc(sizeof(*Z), rows);
    return Z;
  }
  if (!(Z=SDDS_GetColumnInDoubles(SDDSin, ZName))) {
    fprintf(stdout, "Unable to read column %s (ZTRANSVERSE)\n",
            ZName);
    fflush(stdout);
    exit(1);
  }
  return Z;
}

void optimizeBinSettingsForImpedance(double timeSpan, double freq, double Q,
                                     double *binSize, long *nBins, long maxBins)
{
  long n_bins, maxBins2;
  double bin_size, factor, factor1, factor2;
  TRACKING_CONTEXT tcontext;
  
  n_bins = *nBins;
  bin_size = *binSize;
  getTrackingContext(&tcontext);

  if (maxBins<=0)
    maxBins2 = pow(2, 20);
  else
    maxBins2 = pow(2, (long)(log(maxBins)/log(2)+1));
  if (maxBins>0 && maxBins!=maxBins2)
    fprintf(stdout, "Adjusted maximum number of bins for %s %s to %ld\n",
            entity_name[tcontext.elementType],
            tcontext.elementName, maxBins2);
  
  if (1/(2*freq*bin_size)<10) {
    /* want maximum frequency in Z > 10*fResonance */
    fprintf(stdout, "%s %s has excessive bin size for given resonance frequency\n",
            entity_name[tcontext.elementType],
            tcontext.elementName);
    bin_size = 1./freq/20;
    fprintf(stdout, "  Bin size adjusted to %e\n", bin_size);
    fflush(stdout);
  }
  if (2*timeSpan>bin_size*n_bins) {
    fprintf(stdout, "%s %s has insufficient time span for initial bunch\n",
            entity_name[tcontext.elementType],
            tcontext.elementName);
    n_bins = pow(2,
                 (long)(log(2*timeSpan*1.05/bin_size)/log(2)+1));
    if (maxBins2<n_bins) {
      fprintf(stderr, "  Maximum number of bins does not allow sufficient time span!\n");
      exit(1);
    }
    fprintf(stdout, "  Number of bins adjusted to %ld\n",
            n_bins);
    fflush(stdout);
  }
  if (Q<1) 
    /* Ideally, want frequency resolution < fResonance/200 and < fResonanceWidth/200 */
    factor = 200/(n_bins*bin_size*freq/Q);
  else
    factor = 200/(n_bins*bin_size*freq);
  if (factor>1) {
    fprintf(stdout, "%s %s has too few bins or excessively small bin size for given frequency\n",
            entity_name[tcontext.elementType],
            tcontext.elementName);
    if (n_bins*factor>maxBins2) {
      if ((n_bins*factor)/maxBins2>50) {
        if (Q<1)
          fprintf(stdout, " With %ld bins, the frequency resolution is only %.1g times the resonant frequency\n",
                  maxBins2, freq*maxBins2*bin_size);
        else
          fprintf(stdout, " With %ld bins, the frequency resolution is only %.1g times the resonance width\n",
                  maxBins2, freq/Q*maxBins2*bin_size);
        fprintf(stdout, " It isn't possible to model this situation accurately with %ld bins.  Consider the RFMODE or TRFMODE element.\n",
                maxBins2);
        fprintf(stdout, " Alternatively, consider increasing your bin size or maximum number of bins\n");
        exit(1);
      }
      n_bins = maxBins2;
      fprintf(stdout, "  Number of bins adjusted to %ld\n", n_bins);
      
    } else {
      n_bins = pow(2, (long)(log(n_bins*factor)/log(2)+1));
      fprintf(stdout, "  Number of bins adjusted to %ld\n",
              n_bins);
    }
    fflush(stdout);
  }
  
  *nBins = n_bins;
  *binSize = bin_size;
}
