/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: ftrfmode.c
 * contents: track_through_ftrfmode(), set_up_ftrfmode(),
 *
 * Michael Borland, 1993, 2003
 */
#include "mdb.h"
#include "track.h"

#if USE_MPI
typedef struct {
  long index;  /* Records the real index in the whole histogram */
  long count;
  double xsum, ysum;
} HISTOGRAM_DATA;

#endif

void track_through_ftrfmode(
                            double **part, long np, FTRFMODE *trfmode, double Po,
                            char *element_name, double element_z, long pass, long n_passes,
                            CHARGE *charge
                            )
{
  static unsigned long *count = NULL;
  static double *xsum = NULL;              /* sum of x coordinate in each bin = N*<x> */
  static double *ysum = NULL;              /* sum of y coordinate in each bin = N*<y> */
  static double *Vxbin = NULL;             /* array for voltage acting on each bin MV */
  static double *Vybin = NULL;             /* array for voltage acting on each bin MV */
  static double *Vzbin = NULL;             /* array for voltage acting on each bin MV */
  static long *pbin = NULL;                /* array to record which bin each particle is in */
  static double *time = NULL;              /* array to record arrival time of each particle */
  static long max_np = 0;
  long ip, ib;
  double tmin, tmax, tmean, dt, P;
  double Vxb, Vyb, V, omega, phase, t, k, omegaOverC, damping_factor, tau;
  double Px, Py, Pz;
  double Q, Qrp;
  long lastBin, imode;
  double rampFactor;
  long nonEmptyBins = 0;
#if USE_MPI
  static long *nonEmptyArray = NULL;
  double *buffer = NULL;
  long  *reverseMap = NULL;
  long np_total, i, j, map_j=0;
  long firstBin = trfmode->n_bins, firstBin_global, lastBin_global;
  long nonEmptyBins_total = 0, offset = 0;
  MPI_Status status;
  HISTOGRAM_DATA *subHis; /* a compressed histogram with non-zero bins only */ 
#ifdef  USE_MPE /* use the MPE library */
  int event1a, event1b, event2a, event2b;
  event1a = MPE_Log_get_event_number();
  event1b = MPE_Log_get_event_number();
  event2a = MPE_Log_get_event_number();
  event2b = MPE_Log_get_event_number();
  MPE_Describe_state(event1a, event1b, "Count nonEmptyBins", "red");
  MPE_Describe_state(event2a, event2b, "Compress-unCompress", "yellow");
#endif
#endif

  if (charge)
    trfmode->mp_charge = charge->macroParticleCharge;
  else
    bombElegant("CHARGE element required to use FTRFMODE", NULL);
  if (trfmode->mp_charge==0 || (trfmode->xfactor==0 && trfmode->yfactor==0))
    return;
  
  if (!trfmode->initialized)
    bombElegant("track_through_ftrfmode called with uninitialized element", NULL);

  if (trfmode->outputFile && pass==0 && !SDDS_StartPage(&trfmode->SDDSout, n_passes))
    SDDS_Bomb("Problem starting page for FTRFMODE output file");

  
  xsum = calloc(trfmode->n_bins, sizeof(*xsum));
  ysum = calloc(trfmode->n_bins, sizeof(*ysum));
  count = calloc(trfmode->n_bins, sizeof(*count));
#if (!USE_MPI)
  Vxbin = calloc(trfmode->n_bins, sizeof(*Vxbin));
  Vybin = calloc(trfmode->n_bins, sizeof(*Vybin));
  Vzbin = calloc(trfmode->n_bins, sizeof(*Vzbin));
  if (!(xsum && ysum && count && Vxbin && Vybin && Vzbin))
    bomb ("Memory allocation failure in track_through_ftrfmod", NULL);
#endif

  if (!(xsum && ysum && count))
    bomb ("Memory allocation failure in track_through_ftrfmod", NULL);

  if (np>max_np) {
    pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
    time = trealloc(time, sizeof(*time)*max_np);
  }

  tmean = 0;
  if (isSlave) {
    for (ip=0; ip<np; ip++) {
      P = Po*(part[ip][5]+1);
      time[ip] = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
      tmean += time[ip];
    }
  }
#if USE_MPI
  if (isSlave) {
    double t_total;
    MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
    MPI_Allreduce(&tmean, &t_total, 1, MPI_DOUBLE, MPI_SUM, workers);
    tmean = t_total;
  }
  tmean /= np_total;      
#else
  tmean /= np;
#endif

    tmin = tmean - trfmode->bin_size*trfmode->n_bins/2.;
    tmax = tmean + trfmode->bin_size*trfmode->n_bins/2.;

    dt = (tmax - tmin)/trfmode->n_bins;
    lastBin = -1;

  if (isSlave) {  
    for (ip=0; ip<np; ip++) {
      pbin[ip] = -1;
      ib = (time[ip]-tmin)/dt;
      if (ib<0)
	continue;
      if (ib>trfmode->n_bins - 1)
	continue;
      xsum[ib] += part[ip][0]-trfmode->dx;
      ysum[ib] += part[ip][2]-trfmode->dy;
      count[ib] += 1;
      if (count[ib]==1) /* Count nonEmptyBins for the first time only */
        nonEmptyBins++;
      pbin[ip] = ib;
      if (ib>lastBin)
	lastBin = ib;
#if USE_MPI
      if (ib<firstBin)    
	firstBin = ib;
#endif
    }
  }

#if USE_MPI
  MPI_Reduce(&lastBin, &lastBin_global, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&firstBin, &firstBin_global, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    if (isMaster) {
      lastBin = lastBin_global;
      firstBin = firstBin_global; 
    }

    if (!nonEmptyArray)
      nonEmptyArray = malloc(sizeof(*nonEmptyArray)*n_processors);

    MPI_Gather(&nonEmptyBins,1,MPI_LONG,nonEmptyArray,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (isMaster){
      for (i=1; i<n_processors; i++) 
	nonEmptyBins_total += nonEmptyArray[i];
    }
    MPI_Bcast(&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  
    if (!(subHis = calloc(nonEmptyBins_total, sizeof(*subHis))))
      bomb ("Memory allocation failure in track_through_ftrfmod", NULL);
 
    Vxbin = calloc(nonEmptyBins_total, sizeof(*Vxbin));
    Vybin = calloc(nonEmptyBins_total, sizeof(*Vybin));
    Vzbin = calloc(nonEmptyBins_total, sizeof(*Vzbin));
    reverseMap = calloc(trfmode->n_bins, sizeof(*reverseMap));

    if (isSlave) {
      /* figure out the overlapped part of histogram between neighbor processors */
      /* MPI_Sendrecv(&lastBin, 1, MPI_LONG, myid+1, 109, &ovlap_f, myid-1, 109, workers, &status);
      MPI_Sendrecv(&firstBin, 1, MPI_LONG, myid-1, 110, &ovlap_e, myid+1, 110, workers, &status); 
      */
      for (i=0,ib=firstBin; ib<=lastBin; ib++) {
	if (count[ib]){        
          subHis[i].index = ib;
	  subHis[i].count = count[ib];
	  subHis[i].xsum = xsum[ib];
       	  subHis[i].ysum = ysum[ib]; 
	  i++;
	}
      }
      MPI_Send(subHis, nonEmptyBins*sizeof(*subHis), MPI_BYTE, 0, 108, MPI_COMM_WORLD);
    }
    else {
      for (i=1; i<n_processors; i++) {
        if (i>1)
	  offset += nonEmptyArray[i-1];
	MPI_Recv (&subHis[offset], nonEmptyArray[i]*sizeof(*subHis), MPI_BYTE, i, 108, MPI_COMM_WORLD, &status); 
	for (j=offset; j<nonEmptyArray[i]+offset; j++) {
          #define current subHis[j]
          map_j = current.index;
          count[map_j] += current.count;
          xsum[map_j] += current.xsum;
	  ysum[map_j] += current.ysum; 
	}
      }
#ifdef  USE_MPE
      MPE_Log_event(event1a, 0, "start initialize");
#endif       
      for (i=0, ib=firstBin; ib<=lastBin; ib++) { 
        if (count[ib]) {
          subHis[i].index = ib;
	  subHis[i].count = count[ib];
	  subHis[i].xsum = xsum[ib];
       	  subHis[i].ysum = ysum[ib];
          i++;  
        }
      } 
#ifdef  USE_MPE
      MPE_Log_event(event1b, 0, "end initialize");
#endif 
      /* If there are overlapped bins between different processors, the number should be less than the original
	 nonEmptyBins_total */
      nonEmptyBins_total = i;  
    }
    MPI_Bcast (&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast (subHis, nonEmptyBins_total*sizeof(*subHis), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    if (isSlave) {
      for (i=0; i<nonEmptyBins_total; i++) {      
	count[i] = subHis[i].count;
	xsum[i] = subHis[i].xsum;
	ysum[i] = subHis[i].ysum;
        reverseMap[subHis[i].index] = i;  
      }
    }
#endif
    
    rampFactor = 0;
    if (pass > (trfmode->rampPasses-1)) 
      rampFactor = 1;
    else
      rampFactor = (pass+1.0)/trfmode->rampPasses;

    if (isSlave){
      double last_t = trfmode->last_t; /* Save it and use it later for different modes */

      for (imode=0; imode<trfmode->modes; imode++) 
#if USE_MPI
	if (myid == (imode%(n_processors-1)+1))    /* Decompse the modes on different processors */ 
#endif
      {
	if (trfmode->cutoffFrequency>0 && (trfmode->omega[imode] > PIx2*trfmode->cutoffFrequency))
	  continue;
	if (!trfmode->doX[imode] && !trfmode->doY[imode])
          continue;

	trfmode->last_t = last_t;
#if USE_MPI
	for (i=0; i<nonEmptyBins_total; i++) {  /* For parallel version, we use non-empty part of array only */
	  ib = i;
	  t = tmin+(subHis[i].index+0.5)*dt;           /* middle arrival time for this bin */
#else
	for (ib=0; ib<=lastBin; ib++) {
	  if (!count[ib] || (!xsum[ib] && !ysum[ib]))
	    continue;
	  t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */
#endif
      if (t<trfmode->last_t) {
        trfmode->last_t = t;
	/*  fprintf(stdout, "*** Warning: reference time reset for FTRFMODE.  Should only happen once per step.\n");
	    fflush(stdout); */
      }            
	omega = trfmode->omega[imode];
	Q = trfmode->Q[imode]/(1+trfmode->beta[imode]);
	tau = 2*Q/omega;
	Qrp = sqrt(Q*Q - 0.25);
	k = omega/2*trfmode->Rs[imode]/trfmode->Q[imode];
	/* These adjustments per Zotter and Kheifets, 3.2.4, 3.3.2 */
	k *= Q/Qrp;
	omega *= Qrp/Q;
        omegaOverC = omega/c_mks;
        
	damping_factor = exp(-(t-trfmode->last_t)/tau);
	if (trfmode->doX[imode]) {
	  /* -- x plane */
	  /* advance the phasor */
	  phase = trfmode->lastPhasex[imode] + omega*(t - trfmode->last_t);
	  V = trfmode->Vx[imode]*damping_factor;
	  trfmode->Vxr[imode] = V*cos(phase);
	  trfmode->Vxi[imode] = V*sin(phase);
	  trfmode->lastPhasex[imode] = phase;
	  /* add this cavity's contribution to this bin */
	  Vxbin[ib] += trfmode->Vxr[imode];

	  /* compute beam-induced voltage for this bin */
	  Vxb = 2*k*trfmode->mp_charge*particleRelSign*xsum[ib]*trfmode->xfactor*rampFactor; 
          Vzbin[ib] += omegaOverC*(xsum[ib]/count[ib])*(trfmode->Vxi[imode] - Vxb/2);

	  /* add beam-induced voltage to cavity voltage---it is imaginary as
	   * the voltage is 90deg out of phase 
	   */
	  trfmode->Vxi[imode] -= Vxb;

          /* update the phasor */
	  if (trfmode->Vxi[imode]==0 && trfmode->Vxr[imode]==0)
	    trfmode->lastPhasex[imode] = 0;
	  else
	    trfmode->lastPhasex[imode] = atan2(trfmode->Vxi[imode], trfmode->Vxr[imode]);
	  trfmode->Vx[imode] = sqrt(sqr(trfmode->Vxr[imode])+sqr(trfmode->Vxi[imode]));
	}
	if (trfmode->doY[imode]) {
	  /* -- y plane */
	  /* advance the phasor */
	  phase = trfmode->lastPhasey[imode] + omega*(t - trfmode->last_t);
	  V = trfmode->Vy[imode]*damping_factor;
	  trfmode->Vyr[imode] = V*cos(phase);
	  trfmode->Vyi[imode] = V*sin(phase);
	  trfmode->lastPhasey[imode] = phase;
	  /* add this cavity's contribution to this bin */
	  Vybin[ib] += trfmode->Vyr[imode];

	  /* compute beam-induced voltage for this bin */
	  Vyb = 2*k*trfmode->mp_charge*particleRelSign*ysum[ib]*trfmode->yfactor*rampFactor;
          Vzbin[ib] += omegaOverC*(ysum[ib]/count[ib])*(trfmode->Vyi[imode] - Vyb/2);

	  /* add beam-induced voltage to cavity voltage---it is imaginary as
	   * the voltage is 90deg out of phase 
	   */
	  trfmode->Vyi[imode] -= Vyb;

          /* update the phasor */
	  if (trfmode->Vyi[imode]==0 && trfmode->Vyr[imode]==0)
	    trfmode->lastPhasey[imode] = 0;
	  else
            trfmode->lastPhasey[imode] = atan2(trfmode->Vyi[imode], trfmode->Vyr[imode]);
	  trfmode->Vy[imode] = sqrt(sqr(trfmode->Vyr[imode])+sqr(trfmode->Vyi[imode]));
	}
     trfmode->last_t = t;
      } /* loop over bins */

    } /* loop over modes */
    }

#if (USE_MPI)
      if (isSlave) {
	/* Sum up the result for different modes on different CPUs */ 
       	buffer = malloc(sizeof(double)*nonEmptyBins_total);
	MPI_Allreduce(Vxbin, buffer, nonEmptyBins_total, MPI_DOUBLE, MPI_SUM, workers);
        memcpy(Vxbin, buffer, sizeof(*Vxbin)*nonEmptyBins_total);
	MPI_Allreduce(Vybin, buffer, nonEmptyBins_total, MPI_DOUBLE, MPI_SUM, workers);
        memcpy(Vybin, buffer, sizeof(*Vybin)*nonEmptyBins_total);
	MPI_Allreduce(Vzbin, buffer, nonEmptyBins_total, MPI_DOUBLE, MPI_SUM, workers);
        memcpy(Vzbin, buffer, sizeof(*Vzbin)*nonEmptyBins_total); 
      }
      if (trfmode->outputFile) {
	printf (" the trfmode->output will not work for the current version of SDDS in Pelegant\n");
        exit (1);
      }
#endif   
    if (trfmode->outputFile) {
      for (imode=0; imode<trfmode->modes; imode++) {
        if (!SDDS_SetRowValues(&trfmode->SDDSout, 
                               SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, pass,
                               trfmode->xModeIndex[imode], trfmode->Vx[imode],
                               trfmode->yModeIndex[imode], trfmode->Vy[imode], -1))
          SDDS_Bomb("Problem writing data to FTRFMODE output file");
      }
    }

    if (0) {
      FILE *fp;
      fp = fopen("modeData.sdds", "w");
      fprintf(fp, "SDDS1\n&column name=bin type=long &end\n");
      fprintf(fp, "&column name=Voltage, type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (ib=0; ib<trfmode->n_bins; ib++)
        fprintf(fp, "%ld %e\n", ib, Vybin[ib]);
      fclose(fp);
    }
     
    /* change particle slopes to reflect voltage in relevant bin */
#if  (!USE_MPI)
    for (ip=0; ip<np; ip++) {
      if (pbin[ip]>=0) {
	P = Po*(1+part[ip][5]);
	Pz = P/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3])) + Vzbin[pbin[ip]]/(1e6*particleMassMV*particleRelSign);
	Px = part[ip][1]*Pz + Vxbin[pbin[ip]]/(1e6*particleMassMV*particleRelSign);
	Py = part[ip][3]*Pz + Vybin[pbin[ip]]/(1e6*particleMassMV*particleRelSign);
	P  = sqrt(Pz*Pz+Px*Px+Py*Py);
	part[ip][1] = Px/Pz;
	part[ip][3] = Py/Pz;
	part[ip][5] = (P-Po)/Po;
	part[ip][4] = time[ip]*c_mks*P/sqrt(sqr(P)+1);
      }
    }
#else
    if (isSlave) {
      for (ip=0; ip<np; ip++) {
	if (pbin[ip]>=0) {
          i = reverseMap[pbin[ip]];
	  P = Po*(1+part[ip][5]);
	  Pz = P/sqrt(1+sqr(part[ip][1])+sqr(part[ip][3])) + Vzbin[i]/(1e6*particleMassMV*particleRelSign);
	  Px = part[ip][1]*Pz + Vxbin[i]/(1e6*particleMassMV*particleRelSign);
	  Py = part[ip][3]*Pz + Vybin[i]/(1e6*particleMassMV*particleRelSign);
	  P  = sqrt(Pz*Pz+Px*Px+Py*Py);
	  part[ip][1] = Px/Pz;
	  part[ip][3] = Py/Pz;
	  part[ip][5] = (P-Po)/Po;
	  part[ip][4] = time[ip]*c_mks*P/sqrt(sqr(P)+1);
	}
      }
    } 

    if (myid == 1) /* We let the first slave to dump the parameter */
#endif
  if (trfmode->outputFile) {
    if (!SDDS_SetRowValues(&trfmode->SDDSout, 
                           SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, pass,
                           "Pass", pass, NULL))
      SDDS_Bomb("Problem writing data to FTRFMODE output file");
    if ((trfmode->flushInterval<1 || pass%trfmode->flushInterval==0 || pass==(n_passes-1)) &&
        !SDDS_UpdatePage(&trfmode->SDDSout, 0))
      SDDS_Bomb("Problem writing data to FTRFMODE output file");
  }
     
#ifdef  USE_MPE
      MPE_Log_event(event2a, 0, "start histogram"); 
#endif
  free(xsum);
  free(ysum);
  free(count);
  free(Vxbin);
  free(Vybin);
  free(Vzbin);
#if USE_MPI
  #if defined(MINIMIZE_MEMORY)
    free(nonEmptyArray); nonEmptyArray = NULL;
  #endif
  free(reverseMap);
  free(buffer);
  free(subHis);
#endif

#if defined(MINIMIZE_MEMORY)
  free(pbin);
  free(time);
  time = NULL;
  pbin = NULL;
  max_np = 0;
#endif
#ifdef  USE_MPE
    MPE_Log_event(event2b, 0, "end histogram"); 
#endif
}

void set_up_ftrfmode(FTRFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                     RUN *run, long n_particles,
                     double Po, double total_length)
{
  long imode;
  double T;
  SDDS_DATASET SDDSin;
  
  if (rfmode->initialized)
    return;
#if SDDS_MPI_IO
  if (isSlave)
#endif
  if (n_particles<1)
    bombElegant("too few particles in set_up_ftrfmode()", NULL);
  if (rfmode->n_bins<2)
    bombElegant("too few bins for FTRFMODE", NULL);
  if (!rfmode->filename ||
      !SDDS_InitializeInput(&SDDSin, rfmode->filename))
    bombElegant("unable to open file for FTRFMODE element", NULL);
  /* check existence and properties of required columns  */
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Frequency", "Hz", SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Frequency column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "Q", NULL, SDDS_ANY_FLOATING_TYPE,
                       stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Error: problem with Q column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
    exitElegant(1);
  }
  if (rfmode->useSymmData) {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "$gW$r/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedanceSymm", "Ohms/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedanceSymm column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exitElegant(1);
    }
  }
  else {
    if (SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "$gW$r/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin, "ShuntImpedance", "Ohms/m", SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with ShuntImpedance column for FTRFMODE file %s.  Check existence, type, and units.\n", rfmode->filename);
      exitElegant(1);
    }
  }

  if (!SDDS_ReadPage(&SDDSin))
    SDDS_Bomb("unable to read page from file for FTRFMODE element");
  if ((rfmode->modes = SDDS_RowCount(&SDDSin))<1) {
    fprintf(stdout, "Error: no data in FTRFMODE file %s\n", rfmode->filename);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "beta", NULL, SDDS_ANY_FLOATING_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "beta", NULL, SDDS_ANY_FLOATING_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"beta\" column for FRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exitElegant(1);
    }
  }
  if (SDDS_CheckColumn(&SDDSin, "xMode", NULL, SDDS_ANY_INTEGER_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "xMode", NULL, SDDS_ANY_INTEGER_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"doX\" column for FTRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exitElegant(1);
    }
  }
  if (SDDS_CheckColumn(&SDDSin, "yMode", NULL, SDDS_ANY_INTEGER_TYPE,
                       NULL)!=SDDS_CHECK_NONEXISTENT) {
    if (SDDS_CheckColumn(&SDDSin, "yMode", NULL, SDDS_ANY_INTEGER_TYPE,
                         NULL)!=SDDS_CHECK_OK) {
      fprintf(stdout, "Error: problem with \"doY\" column for FTRFMODE file %s.  Check type and units.\n", rfmode->filename);
      exitElegant(1);
    }
  }
  if (!(rfmode->omega = SDDS_GetColumnInDoubles(&SDDSin, "Frequency")) ||
      !(rfmode->Q = SDDS_GetColumnInDoubles(&SDDSin, "Q")) ||
      (rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedanceSymm"))) ||
      (!rfmode->useSymmData &&
       !(rfmode->Rs = SDDS_GetColumnInDoubles(&SDDSin, "ShuntImpedance")))) 
    SDDS_Bomb("Problem getting data from FTRFMODE file");
  if (!(rfmode->beta = SDDS_GetColumnInDoubles(&SDDSin, "beta"))) {
    if (!(rfmode->beta = malloc(sizeof(*(rfmode->beta))*rfmode->modes)))
      bombElegant("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->beta[imode] = 0;
  }
  if (!(rfmode->doX = SDDS_GetColumnInLong(&SDDSin, "xMode"))) {
    if (!(rfmode->doX = malloc(sizeof(*(rfmode->doX))*rfmode->modes)))
      bombElegant("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->doX[imode] = 1;
  }
  if (!(rfmode->doY = SDDS_GetColumnInLong(&SDDSin, "yMode"))) {
    if (!(rfmode->doY = malloc(sizeof(*(rfmode->doY))*rfmode->modes)))
      bombElegant("memory allocation failure (FTRFMODE)", NULL);
    for (imode=0; imode<rfmode->modes; imode++)
      rfmode->doY[imode] = 1;
  }

  if (!(rfmode->Vx  = malloc(sizeof(*(rfmode->Vx ))*rfmode->modes)) ||
      !(rfmode->Vxr = malloc(sizeof(*(rfmode->Vxr))*rfmode->modes)) ||
      !(rfmode->Vxi = malloc(sizeof(*(rfmode->Vxi))*rfmode->modes)) ||
      !(rfmode->Vy  = malloc(sizeof(*(rfmode->Vy ))*rfmode->modes)) ||
      !(rfmode->Vyr = malloc(sizeof(*(rfmode->Vyr))*rfmode->modes)) ||
      !(rfmode->Vyi = malloc(sizeof(*(rfmode->Vyi))*rfmode->modes)) ||
      !(rfmode->lastPhasex = malloc(sizeof(*(rfmode->lastPhasex))*rfmode->modes)) ||
      !(rfmode->lastPhasey = malloc(sizeof(*(rfmode->lastPhasey))*rfmode->modes)))
    bombElegant("memory allocation failure (FTRFMODE)", NULL);
  
  for (imode=0; imode<rfmode->modes; imode++) {
    rfmode->omega[imode] *= PIx2;
    rfmode->Vx[imode] = rfmode->Vxr[imode] = rfmode->Vxi[imode] = 0;
    rfmode->lastPhasex[imode] = 0;
    rfmode->Vy[imode] = rfmode->Vyr[imode] = rfmode->Vyi[imode] = 0;
    rfmode->lastPhasey[imode] = 0;
  }

  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      T = rfmode->bin_size*rfmode->n_bins;
      rfmode->bin_size = 0.1/(rfmode->omega[imode]/PIx2);
      rfmode->n_bins = T/rfmode->bin_size+1;
      rfmode->bin_size = T/rfmode->n_bins;
      fprintf(stdout, "The FTRFMODE %s bin size is too large for mode %ld--setting to %e and increasing to %ld bins\n",
              element_name, imode, rfmode->bin_size, rfmode->n_bins);
      fflush(stdout);
    }
  }

  for (imode=0; imode<rfmode->modes; imode++) {
    if (rfmode->bin_size*rfmode->omega[imode]/PIx2>0.1) {
      fprintf(stdout, "Error: FTRFMODE bin size adjustment failed\n");
      exitElegant(1);
    }
  }

  rfmode->last_t = element_z/c_mks;

#if (USE_MPI)
  if (myid == 1) {/* We let the first slave to dump the parameter */
    dup2(fd,fileno(stdout));
#endif
  if (rfmode->outputFile) {
    TRACKING_CONTEXT context;
    char *filename;
    getTrackingContext(&context);
    filename = compose_filename(rfmode->outputFile, context.rootname);
    if (!SDDS_InitializeOutput(&rfmode->SDDSout, SDDS_BINARY, 0, NULL, NULL, 
                               filename) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSout, "Pass", NULL, SDDS_LONG)) {
      fprintf(stderr, "Problem initializing file %s for FTRFMODE element\n", filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    if (!(rfmode->xModeIndex = malloc(sizeof(*(rfmode->xModeIndex))*rfmode->modes)) ||
        !(rfmode->yModeIndex = malloc(sizeof(*(rfmode->yModeIndex))*rfmode->modes))) {
      fprintf(stderr, "Memory allocation failure for TFRFMODE element\n");
      exitElegant(1);
    }
    for (imode=0; imode<rfmode->modes; imode++) {
      char sx[100], sy[100];
      sprintf(sx, "VxMode%03ld", imode);
      sprintf(sy, "VyMode%03ld", imode);
      if ((rfmode->xModeIndex[imode]
           =SDDS_DefineColumn(&rfmode->SDDSout, sx, NULL, "V", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          ((rfmode->yModeIndex[imode]
            =SDDS_DefineColumn(&rfmode->SDDSout, sy, NULL, "V", NULL, NULL, SDDS_DOUBLE, 0))<0)) {
        fprintf(stderr, "Problem initializing file %s for FTRFMODE element\n", filename);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exitElegant(1);
      }
    }
    if (!SDDS_WriteLayout(&rfmode->SDDSout)) {
      fprintf(stderr, "Problem initializing file %s for FTRFMODE element\n", filename);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    if (filename!=rfmode->outputFile)
      free(filename);
  }
#if (USE_MPI)
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout);  
#endif
  } /* We let the first slave to dump the parameter */
#endif

  rfmode->initialized = 1;
}

