/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: rfmode.c
 * contents: track_through_rfmode()
 *
 * Michael Borland, 1993
 */
#if USE_MPI 
#include "mpi.h" /* Defines the interface to MPI allowing the use of all MPI functions. */
#if USE_MPE
#include "mpe.h" /* Defines the MPE library */ 
#endif
#endif
#include <complex>
#include "mdb.h"
#include "track.h"

#ifdef __cplusplus
extern "C" {
#endif
long find_nearby_array_entry(double *entry, long n, double key);
double linear_interpolation(double *y, double *t, long n, double t0, long i);
#ifdef __cplusplus
}
#endif

#if USE_MPI
void histogram_sums(long nonEmptyBins, long firstBin, long *lastBin, long *his);
#endif

void runBinlessRfMode(double **part, long np, RFMODE *rfmode, double Po,
		      char *element_name, double element_z, long pass, long n_passes,
		      CHARGE *charge);

void track_through_rfmode(
                          double **part, long np, RFMODE *rfmode, double Po,
                          char *element_name, double element_z, long pass, long n_passes,
                          CHARGE *charge
    )
{
    static long *Ihist = NULL;               /* array for histogram of particle density */
    static double *Vbin = NULL;              /* array for voltage acting on each bin */
    static long max_n_bins = 0;
    static long *pbin = NULL;                /* array to record which bin each particle is in */
    static double *time = NULL;              /* array to record arrival time of each particle */
    static long max_np = 0;
    long ip, ib, nb2, lastBin=0, n_binned=0;
    double tmin=0, tmax, tmean, dt=0, P;
    double Vb, V, omega=0, phase, t, k, damping_factor, tau;
    double V_sum, Vr_sum, phase_sum;
    double Q_sum, dgamma;
    long n_summed, max_hist, n_occupied;
    static long been_warned = 0;
    double Qrp, VbImagFactor, Q=0;
    long deltaPass;
#if USE_MPI
    long nonEmptyBins = 0;
    long np_total;
    long firstBin = rfmode->n_bins;
#endif

    if (rfmode->binless) { /* This can't be done in parallel mode */
#if USE_MPI
    fprintf(stdout, (char*)"binless in rfmode is not supported in the current parallel version.\n");
    fprintf(stdout, (char*)"Please use serial version.\n");
    fflush(stdout);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 9);
#endif
      runBinlessRfMode(part, np, rfmode, Po, element_name, element_z, pass, n_passes, charge);
      return;
    }
    if (charge) {
      rfmode->mp_charge = charge->macroParticleCharge;
    } else if (pass==0) {
      rfmode->mp_charge = 0;
      if (rfmode->charge<0)
        bombElegant((char*)"RFMODE charge parameter should be non-negative. Use change_particle to set particle charge state.", NULL);
#if (!USE_MPI) 
      if (np)
        rfmode->mp_charge = rfmode->charge/np;
#else
      if (notSinglePart) {
	MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	if (np_total)
	  rfmode->mp_charge = rfmode->charge/np_total; 
      } else {
        if (np)
          rfmode->mp_charge = rfmode->charge/np;
      }
#endif
    }

    if (pass%rfmode->pass_interval)
      return;
    
    if (isMaster && (!been_warned)) {        
        if (rfmode->freq<1e3 && rfmode->freq)  {
            fprintf(stdout, (char*)"\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
            fflush(stdout);
            been_warned = 1;
            }
        if (been_warned) {
            fprintf(stdout, (char*)"units of parameters for RFMODE are as follows:\n");
            fflush(stdout);
            print_dictionary_entry(stdout, T_RFMODE, 0, 0);
            }
        }

    if (rfmode->mp_charge==0) {
      return ;
    }
    if (rfmode->detuned_until_pass>pass) {
      return ;
    }

    if (!rfmode->initialized)
        bombElegant((char*)"track_through_rfmode called with uninitialized element", NULL);

    if (rfmode->n_bins>max_n_bins) {
      Ihist = (long*)trealloc(Ihist, sizeof(*Ihist)*(max_n_bins=rfmode->n_bins));
      Vbin = (double*)trealloc(Vbin, sizeof(*Vbin)*max_n_bins);
    }

    if (np>max_np) {
      pbin = (long*)trealloc(pbin, sizeof(*pbin)*(max_np=np));
      time = (double*)trealloc(time, sizeof(*time)*max_np);
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
    if (notSinglePart) {
      if (isSlave) {
	double t_total;
	MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	MPI_Allreduce(&tmean, &t_total, 1, MPI_DOUBLE, MPI_SUM, workers);
	tmean = t_total;
      }
      tmean /= np_total;
    } else
      tmean /= np;
#else
    tmean /= np;
#endif

   if (isSlave) {
      tmin = tmean - rfmode->bin_size*rfmode->n_bins/2.;
      tmax = tmean + rfmode->bin_size*rfmode->n_bins/2.;

      for (ib=0; ib<rfmode->n_bins; ib++)
        Ihist[ib] = 0;
      
      dt = (tmax - tmin)/rfmode->n_bins;
      n_binned = lastBin = 0;
      for (ip=0; ip<np; ip++) {
        pbin[ip] = -1;
        ib = (long)((time[ip]-tmin)/dt);
        if (ib<0)
            continue;
        if (ib>rfmode->n_bins - 1)
            continue;
        Ihist[ib] += 1;
#if USE_MPI
	if (Ihist[ib]==1)
	  nonEmptyBins++;
#endif
        pbin[ip] = ib;
        if (ib>lastBin)
          lastBin = ib;
#if USE_MPI
      if (ib<firstBin)    
	firstBin = ib;
#endif
        n_binned++;
      }
      V_sum = Vr_sum = phase_sum = Q_sum = 0;
      n_summed = max_hist = n_occupied = 0;
      nb2 = rfmode->n_bins/2;
    
      /* find frequency and Q at this time */
      omega = PIx2*rfmode->freq;
      if (rfmode->nFreq) {
	double omegaFactor;
	ib = find_nearby_array_entry(rfmode->tFreq, rfmode->nFreq, tmean);
	omega *= (omegaFactor=linear_interpolation(rfmode->fFreq, rfmode->tFreq, rfmode->nFreq, tmean, ib));
	/* keeps stored energy constant for constant R/Q */
	rfmode->V *= sqrt(omegaFactor);
      }
      Q = rfmode->Q/(1+rfmode->beta);
      if (rfmode->nQ) {
	ib = find_nearby_array_entry(rfmode->tQ, rfmode->nQ, tmean);
	Q *= linear_interpolation(rfmode->fQ, rfmode->tQ, rfmode->nQ, tmean, ib);
      }
    }
#if (!USE_MPI)
    if (Q<0.5) {
      fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
      fflush(stdout);
      exitElegant(1);
    }
#else
    if (myid == 1) { /* Let the first slave processor write the output */
      if (Q<0.5) {
	dup2(fd,fileno(stdout)); 
	fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
	fflush(stdout);
	close(fd);
	MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);
      }
    }
#endif
      tau = 2*Q/omega;
      k = omega/4*(rfmode->RaInternal)/rfmode->Q;
      if ((deltaPass = (pass-rfmode->detuned_until_pass)) <= (rfmode->rampPasses-1)) 
	k *= (deltaPass+1.0)/rfmode->rampPasses;

      /* These adjustments per Zotter and Kheifets, 3.2.4 */
      Qrp = sqrt(Q*Q - 0.25);
      VbImagFactor = 1/(2*Qrp);
      omega *= Qrp/Q;

      if (rfmode->single_pass) {
	rfmode->V = rfmode->last_phase = 0;
	rfmode->last_t = tmin + 0.5*dt;
      }
#if USE_MPI
      histogram_sums(nonEmptyBins, firstBin, &lastBin, Ihist);

      /*      if (isSlave) {
        long lastBin_global;         
        MPI_Allreduce(&lastBin, &lastBin_global, 1, MPI_LONG, MPI_MAX, workers);
        lastBin = lastBin_global;
      }
      if(isSlave) {
        buffer = (double*)tmalloc(sizeof(double) * (lastBin+1)); 
	MPI_Allreduce(Ihist, buffer, lastBin+1, MPI_LONG, MPI_SUM, workers);
        memcpy(Ihist, buffer, sizeof(long)*(lastBin+1));
	free(buffer);
      }
      */
     
#endif 
      for (ib=0; ib<=lastBin; ib++) {
	if (!Ihist[ib])
          continue;
        if (Ihist[ib]>max_hist)
          max_hist = Ihist[ib];
        n_occupied++;

        t = tmin+(ib+0.5)*dt;           /* middle arrival time for this bin */

        /* advance cavity to this time */
        phase = rfmode->last_phase + omega*(t - rfmode->last_t);
        damping_factor = exp(-(t-rfmode->last_t)/tau);
        rfmode->last_t = t;
        rfmode->last_phase = phase;
        V = rfmode->V*damping_factor;
        rfmode->Vr = V*cos(phase);
        rfmode->Vi = V*sin(phase);

        /* compute beam-induced voltage for this bin */
        Vb = 2*k*rfmode->mp_charge*particleRelSign*rfmode->pass_interval*Ihist[ib];
        Vbin[ib] = rfmode->Vr - Vb/2;
        
        /* add beam-induced voltage to cavity voltage */
        rfmode->Vr -= Vb;
        rfmode->Vi -= Vb*VbImagFactor;
        rfmode->last_phase = atan2(rfmode->Vi, rfmode->Vr);
        rfmode->V = sqrt(sqr(rfmode->Vr)+sqr(rfmode->Vi));
        
        V_sum  += Ihist[ib]*rfmode->V;
        Vr_sum += Ihist[ib]*rfmode->Vr;
        phase_sum += Ihist[ib]*rfmode->last_phase;
        Q_sum += Ihist[ib]*rfmode->mp_charge*particleRelSign;
        n_summed  += Ihist[ib];
      }

      if (rfmode->rigid_until_pass<=pass) {
	/* change particle momentum offsets to reflect voltage in relevant bin */
	/* also recompute slopes for new momentum to conserve transverse momentum */
	for (ip=0; ip<np; ip++) {
	  if (pbin[ip]>=0) {
	    /* compute new momentum and momentum offset for this particle */
	    dgamma = Vbin[pbin[ip]]/(1e6*particleMassMV*particleRelSign);
	    add_to_particle_energy(part[ip], time[ip], Po, dgamma);
	  }
	}
      }
    
    if (rfmode->record) {
#if (USE_MPI)
      if (myid == 1) /* We let the first slave to dump the parameter */
#endif
      if ((pass%rfmode->sample_interval)==0 && 
          (!SDDS_SetRowValues(&rfmode->SDDSrec, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
				(pass/rfmode->sample_interval),
			      (char*)"Pass", pass, (char*)"NumberOccupied", n_occupied,
			      (char*)"FractionBinned", np?(1.0*n_binned)/np:0.0,
			      (char*)"VPostBeam", rfmode->V, (char*)"PhasePostBeam", rfmode->last_phase,
			      (char*)"tPostBeam", rfmode->last_t,
			      (char*)"V", n_summed?V_sum/n_summed:0.0,
			      (char*)"VReal", n_summed?Vr_sum/n_summed:0.0,
			      (char*)"Phase", n_summed?phase_sum/n_summed:0.0, NULL) ||
	   !SDDS_UpdatePage(&rfmode->SDDSrec, 0))) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb((char*)"problem setting up data for RFMODE record file");
      }
      if (pass==n_passes-1 && !SDDS_Terminate(&rfmode->SDDSrec)) {
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb((char*)"problem writing data for RFMODE record file");
      }
    }
 
#if defined(MINIMIZE_MEMORY)
    free(Ihist);
    free(Vbin);
    free(pbin);
    free(time);
    Ihist = pbin = NULL;
    Vbin = time = NULL;
    max_n_bins = max_np = 0;
#endif

  }

void set_up_rfmode(RFMODE *rfmode, char *element_name, double element_z, long n_passes, 
                   RUN *run, long n_particles,
                   double Po, double total_length)
{
  long i, n;
  double T;
  TABLE data;
  
  if (rfmode->initialized)
    return;
  
  rfmode->initialized = 1;
  if (rfmode->pass_interval<=0)
    bombElegant((char*)"pass_interval <= 0 for RFMODE", NULL);
#if !SDDS_MPI_IO
  if (n_particles<1)
    bombElegant((char*)"too few particles in set_up_rfmode()", NULL);
#endif
  if (rfmode->n_bins<2)
    bombElegant((char*)"too few bins for RFMODE", NULL);
  if (rfmode->bin_size<=0 && !rfmode->binless)
    bombElegant((char*)"bin_size must be positive for RFMODE", NULL);
  if (rfmode->Ra && rfmode->Rs) 
    bombElegant((char*)"RFMODE element may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
  if (rfmode->Ra)
    rfmode->RaInternal = rfmode->Ra;
  else
    rfmode->RaInternal = 2*rfmode->Rs;
  if (rfmode->bin_size*rfmode->freq>0.1) {
    T = rfmode->bin_size*rfmode->n_bins;
    rfmode->bin_size = 0.1/rfmode->freq;
    rfmode->n_bins = ((long)(T/rfmode->bin_size+1));
    rfmode->bin_size = T/rfmode->n_bins;
    fprintf(stdout, (char*)"The RFMODE %s bin size is too large--setting to %e and increasing to %ld bins\n",
            element_name, rfmode->bin_size, rfmode->n_bins);
    fprintf(stdout, (char*)"Total span changed from %le to %le\n",
            T, rfmode->n_bins*rfmode->bin_size);
    fflush(stdout);
  }
  if (rfmode->sample_interval<=0)
    rfmode->sample_interval = 1;
  if (rfmode->record && !(rfmode->fileInitialized)) {
    rfmode->record = compose_filename(rfmode->record, run->rootname);
    n = n_passes/rfmode->sample_interval;
#if (USE_MPI)
    if (myid == 1) /* We let the first slave to dump the parameter */
#endif
    if (!SDDS_InitializeOutput(&rfmode->SDDSrec, SDDS_BINARY, 1, NULL, NULL, rfmode->record) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"NumberOccupied", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"FractionBinned", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"V", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VReal", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"Phase", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"VPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"PhasePostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&rfmode->SDDSrec, (char*)"tPostBeam", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&rfmode->SDDSrec) ||
        !SDDS_StartPage(&rfmode->SDDSrec, n+1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem setting up RFMODE record file");
    }
    rfmode->fileInitialized = 0;
  }
  if (rfmode->preload && rfmode->charge) {
    double Vb, omega, To, tau;
    std::complex <double> Vc;
    if (rfmode->fwaveform || rfmode->Qwaveform) {
      printf((char*)"Warning: preloading of RFMODE doesn't work properly with frequency or Q waveforms\n");
      printf((char*)"unless the initial values of the frequency and Q factors are 1.\n");
    }
    To = total_length/(Po*c_mks/sqrt(sqr(Po)+1));
    omega = rfmode->freq*PIx2;
    Vb = 2 * omega/4*(rfmode->Ra)/rfmode->Q * rfmode->charge * rfmode->preload_factor * particleRelSign;
    tau = 2*rfmode->Q/(omega*(1+rfmode->beta));

    Vc = -Vb/(1.0-cexpi(omega*To)*exp(-To/tau));
    rfmode->V = std::abs<double>(Vc);
    rfmode->last_phase = atan2(Vc.imag(), Vc.real());
    fprintf(stdout, (char*)"RFMODE %s at z=%fm preloaded:  V = (%e, %e) V  =  %eV at %fdeg \n",
            element_name, element_z, Vc.real(), Vc.imag(),
            rfmode->V, rfmode->last_phase*180/PI);
    fflush(stdout);
    fprintf(stdout, (char*)"omega=%21.15e To=%21.15es, Vb = %21.15eV, tau = %21.15es\n", omega, To, Vb, tau);
    fflush(stdout);
    rfmode->last_t = element_z/c_mks - To;
  }
  else if (rfmode->initial_V) {
    rfmode->last_phase = rfmode->initial_phase;
    rfmode->V  = rfmode->initial_V;
    rfmode->last_t = rfmode->initial_t;
  }
  else 
    rfmode->last_t = element_z/c_mks;
  rfmode->tFreq = rfmode->fFreq = NULL;
  rfmode->nFreq = 0;
  if (rfmode->fwaveform) {
    if (!getTableFromSearchPath(&data, rfmode->fwaveform, 1, 0)) 
      bombElegant((char*)"unable to read frequency waveform for RFMODE", NULL);
    if (data.n_data<=1)
      bombElegant((char*)"RFMODE frequency table contains less than 2 points", NULL);
    rfmode->tFreq = data.c1;
    rfmode->fFreq = data.c2;
    rfmode->nFreq = data.n_data;
    for (i=1; i<rfmode->nFreq; i++)
      if (rfmode->tFreq[i-1]>rfmode->tFreq[i])
        bombElegant((char*)"time values are not monotonically increasing in RFMODE frequency waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    data.c1 = data.c2 = NULL;
  }
  rfmode->tQ = rfmode->fQ = NULL;
  rfmode->nQ = 0;
  if (rfmode->Qwaveform) {
    if (!getTableFromSearchPath(&data, rfmode->Qwaveform, 1, 0)) 
      bombElegant((char*)"unable to read Q waveform for RFMODE", NULL);
    if (data.n_data<=1)
      bombElegant((char*)"RFMODE Q table contains less than 2 points", NULL);
    rfmode->tQ = data.c1;
    rfmode->fQ = data.c2;
    rfmode->nQ = data.n_data;
    for (i=1; i<rfmode->nQ; i++)
      if (rfmode->tQ[i-1]>rfmode->tQ[i])
        bombElegant((char*)"time values are not monotonically increasing in RFMODE Q waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    data.c1 = data.c2 = NULL;
  }
}

void runBinlessRfMode(
		      double **part, long np, RFMODE *rfmode, double Po,
		      char *element_name, double element_z, long pass, long n_passes,
		      CHARGE *charge
		      )
{
  static TIMEDATA *tData;
  static long max_np = 0;
  long ip, ip0, ib;
  double P;
  double Vb, V, omega, phase, t, k, damping_factor, tau;
  double V_sum, Vr_sum, phase_sum;
  double Vc, Vcr, Q_sum, dgamma, Vb_sum;
  long n_summed, max_hist, n_occupied;
  static long been_warned = 0;
  double Qrp, VbImagFactor, Q;
  double tmean;
  long deltaPass;
#if USE_MPI
  long np_total;

  if (isSlave)
    MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
#endif

  if (charge) {
    rfmode->mp_charge = charge->macroParticleCharge;
  } else if (pass==0) {
    rfmode->mp_charge = 0;
    if (rfmode->charge<0)
      bombElegant((char*)"RFMODE charge parameter should be non-negative. Use change_particle to set particle charge state.", NULL);
#if !USE_MPI
    if (np)
      rfmode->mp_charge = rfmode->charge/np;
#else
    if (notSinglePart) {
      if (np_total)
	rfmode->mp_charge = rfmode->charge/np_total;
    } else {
      if (np)
	rfmode->mp_charge = rfmode->charge/np;  
    }
#endif
  }

  if (pass%rfmode->pass_interval)
    return;
    
  if (!been_warned) {        
    if (rfmode->freq<1e3 && rfmode->freq)  {
      fprintf(stdout, (char*)"\7\7\7warning: your RFMODE frequency is less than 1kHz--this may be an error\n");
      fflush(stdout);
      been_warned = 1;
    }
    if (been_warned) {
      fprintf(stdout, (char*)"units of parameters for RFMODE are as follows:\n");
      fflush(stdout);
      print_dictionary_entry(stdout, T_RFMODE, 0, 0);
    }
  }

  if (rfmode->mp_charge==0) {
    return ;
  }
  if (rfmode->detuned_until_pass>pass) {
    return ;
  }

  if (!rfmode->initialized)
    bombElegant((char*)"track_through_rfmode called with uninitialized element", NULL);

  if (np>max_np) 
    tData = (TIMEDATA*)trealloc(tData, sizeof(*tData)*(max_np=np));

  tmean = 0;
  for (ip=0; ip<np; ip++) {
    P = Po*(part[ip][5]+1);
    tData[ip].t = part[ip][4]*sqrt(sqr(P)+1)/(c_mks*P);
    tData[ip].ip = ip;
    if (isSlave) {
      tmean += tData[ip].t;
    }
  }
  qsort(tData, np, sizeof(*tData), compTimeData);

#if USE_MPI
  if (notSinglePart) {
    if (isSlave) {
      double t_total;

      MPI_Allreduce(&tmean, &t_total, 1, MPI_DOUBLE, MPI_SUM, workers);
      tmean = t_total;
    }
    tmean /= np_total;
  } else
    tmean /= np;
#else
  tmean /= np;
#endif

  V_sum = Vr_sum = phase_sum = Vc = Vcr = Q_sum = Vb_sum = 0;
  n_summed = max_hist = n_occupied = 0;
    
  /* find frequency and Q at this time */
  omega = PIx2*rfmode->freq;
  if (rfmode->nFreq) {
    double omegaFactor;
    ib = find_nearby_array_entry(rfmode->tFreq, rfmode->nFreq, tmean);
    omega *= (omegaFactor=linear_interpolation(rfmode->fFreq, rfmode->tFreq, rfmode->nFreq, tmean, ib));
    /* keeps stored energy constant for constant R/Q */
    rfmode->V *= sqrt(omegaFactor);
  }
  Q = rfmode->Q/(1+rfmode->beta);
  if (rfmode->nQ) {
    ib = find_nearby_array_entry(rfmode->tQ, rfmode->nQ, tmean);
    Q *= linear_interpolation(rfmode->fQ, rfmode->tQ, rfmode->nQ, tmean, ib);
  }
  if (Q<0.5) {
    fprintf(stdout, (char*)"The effective Q<=0.5 for RFMODE.  Use the ZLONGIT element.\n");
    fflush(stdout);
    exitElegant(1);
  }
  tau = 2*Q/omega;
  k = omega/4*(rfmode->RaInternal)/rfmode->Q;
  if ((deltaPass = (pass-rfmode->detuned_until_pass)) <= (rfmode->rampPasses-1)) 
    k *= (deltaPass+1.0)/rfmode->rampPasses;

  /* These adjustments per Zotter and Kheifets, 3.2.4 */
  Qrp = sqrt(Q*Q - 0.25);
  VbImagFactor = 1/(2*Qrp);
  omega *= Qrp/Q;

  if (rfmode->single_pass) {
    rfmode->V = rfmode->last_phase = 0;
    rfmode->last_t = tData[0].t;
  }
    
  for (ip0=0; ip0<np; ip0++) {
    ip = tData[ip0].ip;
    t = tData[ip0].t;
        
    /* advance cavity to this time */
    phase = rfmode->last_phase + omega*(t - rfmode->last_t);
    damping_factor = exp(-(t-rfmode->last_t)/tau);
    rfmode->last_t = t;
    rfmode->last_phase = phase;
    V = rfmode->V*damping_factor;
    rfmode->Vr = V*cos(phase);
    rfmode->Vi = V*sin(phase);
        
    /* compute beam-induced voltage from this particle */
    Vb = 2*k*rfmode->mp_charge*particleRelSign*rfmode->pass_interval;
    /* compute voltage seen by this particle */
    V = rfmode->Vr - Vb/2;
        
    /* add beam-induced voltage to cavity voltage */
    rfmode->Vr -= Vb;
    rfmode->Vi -= Vb*VbImagFactor;
    rfmode->last_phase = atan2(rfmode->Vi, rfmode->Vr);
    rfmode->V = sqrt(sqr(rfmode->Vr)+sqr(rfmode->Vi));
        
    V_sum  += rfmode->V;
    Vb_sum += Vb;
    Vr_sum += rfmode->Vr;
    phase_sum += rfmode->last_phase;
    Q_sum += rfmode->mp_charge*particleRelSign;
    n_summed  += 1;
        
    if (rfmode->rigid_until_pass<=pass) {
      /* change the particle energy */
      dgamma = V/(1e6*particleMassMV*particleRelSign);
      add_to_particle_energy(part[ip], tData[ip0].t, Po, dgamma);
    }
  }

  if (rfmode->record) {
    if ((pass%rfmode->sample_interval)==0 && 
	(!SDDS_SetRowValues(&rfmode->SDDSrec, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    (pass/rfmode->sample_interval),
			    (char*)"Pass", pass, (char*)"NumberOccupied", n_occupied,
			    (char*)"FractionBinned", 1.0,
			    (char*)"VPostBeam", rfmode->V, (char*)"PhasePostBeam", rfmode->last_phase,
			    (char*)"tPostBeam", rfmode->last_t,
			    (char*)"V", n_summed?V_sum/n_summed:0.0,
			    (char*)"VReal", n_summed?Vr_sum/n_summed:0.0,
			    (char*)"Phase", n_summed?phase_sum/n_summed:0.0, 
			    NULL) ||
	 !SDDS_UpdatePage(&rfmode->SDDSrec, 0))) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem setting up data for RFMODE record file");
    }
    if (pass==n_passes-1 && !SDDS_Terminate(&rfmode->SDDSrec)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb((char*)"problem writing data for RFMODE record file");
    }
  }
 
#if defined(MINIMIZE_MEMORY)
  free(tData);
  tData = NULL;
  max_np = 0;
#endif

}


#if USE_MPI
typedef struct {
  long index;  /* Records the real index in the whole histogram */
  long sum;
} SUB_HISTOGRAM;

 void histogram_sums(long nonEmptyBins, long firstBin, long *lastBin, long *his)
{
  static long *nonEmptyArray = NULL;
  long lastBin_global, firstBin_global;
  long nonEmptyBins_total = 0, offset = 0;
  long i, j, map_j, ib;
  static SUB_HISTOGRAM *subHis; /* a compressed histogram with non-zero bins only */
  static long max_nonEmptyBins_total = 0;
  MPI_Status status;

  MPI_Reduce(lastBin, &lastBin_global, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&firstBin, &firstBin_global, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
  if (isMaster) {
    *lastBin = lastBin_global;
    firstBin = firstBin_global; 
  }

  if (!nonEmptyArray)
    nonEmptyArray =(long*) tmalloc(sizeof(*nonEmptyArray)*n_processors);

  MPI_Gather(&nonEmptyBins,1,MPI_LONG,nonEmptyArray,1,MPI_LONG,0,MPI_COMM_WORLD);
  if (isMaster){
    for (i=1; i<n_processors; i++) 
      nonEmptyBins_total += nonEmptyArray[i];
  }
  MPI_Bcast(&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
 
  if (nonEmptyBins_total>max_nonEmptyBins_total) {
    if (!(subHis = (SUB_HISTOGRAM*) trealloc(subHis, sizeof(*subHis)*(max_nonEmptyBins_total=nonEmptyBins_total))))
      bomb ("Memory allocation failure in track_through_ftrfmod", NULL);
  }

  if (isSlave) {
    for (i=0,ib=firstBin; ib<=*lastBin; ib++) {
      if (his[ib]){        
	subHis[i].index = ib;
	subHis[i].sum = his[ib]; 
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
	his[map_j] += current.sum; 
      }
    }
      
    for (i=0, ib=firstBin; ib<=*lastBin; ib++) { 
      if (his[ib]) {
	subHis[i].index = ib;
	subHis[i].sum = his[ib];
	i++;  
      }
    } 
    /* If there are overlapped bins between different processors, the number should be less than the original
       nonEmptyBins_total */
    nonEmptyBins_total = i;  
  }
  MPI_Bcast (&nonEmptyBins_total, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast (subHis, nonEmptyBins_total*sizeof(*subHis), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast (lastBin, 1, MPI_LONG, 0, MPI_COMM_WORLD);

  if (isSlave) {
    for (i=0; i<nonEmptyBins_total; i++) {      
      his[subHis[i].index] = subHis[i].sum; 
    }
  } 
}
#endif
