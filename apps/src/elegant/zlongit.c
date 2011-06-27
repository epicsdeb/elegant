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

/* Sign conventions in this module:
 * V = I*Z 
 * V(t) = V(w)*exp(i*w*t)  I(t) = I(w)*exp(i*w*t)
 * For inductor:   Z = i*w*L
 * For capacitor:  Z = -i/(w*C)
 * For resistor:   Z = R
 * For resonator:  Z = R/(1+i*Q*(w/wr-wr/w))
 *
 * The minus sign to get energy loss instead of energy gain is 
 * used internally and must not be included in the impedance.
 *
 */

#include "mdb.h"
#include "track.h"
#include "table.h"
#include "fftpackC.h"

#define WAKE_COLUMNS 3
static SDDS_DEFINITION wake_column[WAKE_COLUMNS] = {
    {"Deltat", "&column name=Deltat, symbol=\"$gD$rt\", units=s, type=double, description=\"Time after head of bunch\" &end"},
    {"Wz", "&column name=Wz, symbol=\"W$bz$n\", units=V, type=double, description=\"Longitudinal wake\" &end"},
    {"LinearDensity", "&column name=LinearDensity, units=C/s, type=double &end"},
    };

#define WAKE_PARAMETERS 5
#define BB_WAKE_PARAMETERS 5
#define NBB_WAKE_PARAMETERS 2
static SDDS_DEFINITION wake_parameter[WAKE_PARAMETERS] = {
    {"Pass", "&parameter name=Pass, type=long &end"},
    {"q", "&parameter name=q, units=C, type=double, description=\"Total charge\" &end"},
    {"Ra", "&parameter name=Ra, symbol=\"R$ba$n\", units=\"$gW$r\", type=double, description=\"Broad-band impedance\" &end"},
    {"fo", "&parameter name=fo, symbol=\"f$bo$n\", units=Hz, type=double, description=\"Frequency of BB resonator\" &end"},
    {"Deltaf", "&parameter name=Deltaf, symbol=\"$gD$rf\", units=Hz, type=double, description=\"Frequency sampling interval\" &end"},
    } ;

#define DEBUG 0

void set_up_zlongit(ZLONGIT *zlongit, RUN *run, long pass, long particles, CHARGE *charge,
                    double timeSpan);

void track_through_zlongit(double **part, long np, ZLONGIT *zlongit, double Po,
    RUN *run, long i_pass, CHARGE *charge
    )
{
    static double *Itime = NULL;           /* array for histogram of particle density */
    static double *Ifreq = NULL;           /* array for FFT of histogram of particle density */
    static double *Vtime = NULL;           /* array for voltage acting on each bin */
    static long max_n_bins = 0;
    static long *pbin = NULL;              /* array to record which bin each particle is in */
    static double *time = NULL;           /* array to record arrival time of each particle */
    static long max_np = 0;
    double *Vfreq, *Z;
    long ip, ib, nb, n_binned, nfreq, iReal, iImag;
    double factor, tmin, tmax, tmean, dt, dt1, dgam, rampFactor;
    long ip1, ip2, bunches, bunch, npb, i_pass0;
    long bucketEnd[MAX_BUCKETS];
    static long not_first_call = -1;
#if USE_MPI
    double *buffer;
    double tmin_part, tmax_part;           /* record the actual tmin and tmax for particles to reduce communications */
    long offset, length;
#endif
   
#ifdef  USE_MPE /* use the MPE library */
  int event1a, event1b, event2a, event2b;
  event1a = MPE_Log_get_event_number();
  event1b = MPE_Log_get_event_number();
  event2a = MPE_Log_get_event_number();
  event2b = MPE_Log_get_event_number();
  MPE_Describe_state(event1a, event1b, "SavitzyGolaySmooth", "red");
  MPE_Describe_state(event2a, event2b, "fft_inverse", "yellow");
#endif

    i_pass0 = i_pass;
    if ((i_pass -= zlongit->startOnPass)<0)
      return;

    rampFactor = 0;
    if (i_pass>=(zlongit->rampPasses-1))
      rampFactor = 1;
    else
      rampFactor = (i_pass+1.0)/zlongit->rampPasses;
    
    not_first_call += 1;

#if (!USE_MPI)    
    if (np>max_np) {
      pbin = trealloc(pbin, sizeof(*pbin)*(max_np=np));
      time = trealloc(time, sizeof(*time)*max_np);
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
        printf("Error (ZLONGIT): must have bunch_frequency=0 for parallel mode.\n");
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
            bombElegant("Error (wake): maximum number of buckets was exceeded", NULL);
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
      tmean = computeTimeCoordinates(time+ip1, Po, part+ip1, npb);
      find_min_max(&tmin, &tmax, time+ip1, npb);
#if USE_MPI
      find_global_min_max(&tmin, &tmax, np, workers); 
      tmin_part = tmin;
      tmax_part = tmax;     
#endif
      if (bunch==0) {
        /* use np here since we need to compute the macroparticle charge */
        set_up_zlongit(zlongit, run, i_pass, np, charge, tmax-tmin);
      }
  
      nb = zlongit->n_bins;
      dt = zlongit->bin_size;
      if ((tmax-tmin)*2>nb*dt) {
        TRACKING_CONTEXT tcontext;
        getTrackingContext(&tcontext);
        fprintf(stderr, "%s %s: Time span of bunch (%le s) is more than half the total time span (%le s).\n",
                entity_name[tcontext.elementType],
                tcontext.elementName, tmax-tmin, nb*dt);
        fprintf(stderr, "If using broad-band impedance, you should increase the number of bins and rerun.\n");
        fprintf(stderr, "If using file-based impedance, you should increase the number of data points or decrease the frequency resolution.\n");
        exitElegant(1);
      }
      
      if (zlongit->n_bins>max_n_bins) {
        Itime = trealloc(Itime, 2*sizeof(*Itime)*(max_n_bins=zlongit->n_bins));
        Ifreq = trealloc(Ifreq, 2*sizeof(*Ifreq)*(max_n_bins=zlongit->n_bins));
        Vtime = trealloc(Vtime, 2*sizeof(*Vtime)*(max_n_bins+1));
      }
      
      if (zlongit->reverseTimeOrder) {
        for (ip=ip1; ip<=ip2; ip++)
          time[ip] = 2*tmean-time[ip];
      }
      tmin = tmean - dt*zlongit->n_bins/2.0;
      
      for (ib=0; ib<zlongit->n_bins; ib++)
        Itime[2*ib] = Itime[2*ib+1] = 0;
      
      n_binned=0; 
      for (ip=ip1; ip<=ip2; ip++) {
        pbin[ip] = -1;
        ib = (time[ip]-tmin)/dt;
        if (ib<0)
          continue;
        if (ib>nb - 1)
          continue;
        if (zlongit->area_weight && ib>1 && ib<(nb-1)) {
          double dist;
          dist = (time[ip]-((ib+0.5)*dt+tmin))/dt;
          Itime[ib] += 0.5;
          Itime[ib-1] += 0.25-0.5*dist;
          Itime[ib+1] += 0.25+0.5*dist;
        }
        else 
          Itime[ib] += 1;
        pbin[ip] = ib;
        n_binned++;
      }
#if (!USE_MPI)
      if (n_binned!=npb) {
        fprintf(stdout, "Warning: only %ld of %ld particles were binned (ZLONGIT)!\n", n_binned, ip2-ip1+1);
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

#if USE_MPI 
      offset = ((long)((tmin_part-tmin)/dt)-1 ? (long)((tmin_part-tmin)/dt)-1:0);
      length = ((long)((tmax_part-tmin_part)/dt)+2 < nb ? (long)((tmax_part-tmin_part)/dt)+2:nb);
      if (isSlave) {
        buffer = malloc(sizeof(double) * length);;
        MPI_Allreduce(&Itime[offset], buffer, length, MPI_DOUBLE, MPI_SUM, workers);
        memcpy(&Itime[offset], buffer, sizeof(double)*length);
        free(buffer);
      }
#endif

#ifdef  USE_MPE
      MPE_Log_event(event1a, 0, "start zlongit"); /* record time spent on I/O operations */
#endif
      if (zlongit->smoothing)
        SavitzyGolaySmooth(Itime, nb, zlongit->SGOrder, 
                           zlongit->SGHalfWidth, zlongit->SGHalfWidth, 0);
#ifdef  USE_MPE
      MPE_Log_event(event1b, 0, "end zlongit"); /* record time spent on I/O operations */
#endif

#if DEBUG
      /* Output the time-binned data */
      if (1) {
        FILE *fp;
        fp = fopen("zlongit.tbin", "w");
        fprintf(fp, "SDDS1\n&column name=t type=double units=s &end\n&column name=I type=double &end\n&data mode=ascii &end\n");
        fprintf(fp, "%ld\n", nb);
        for (ib=0; ib<nb; ib++) 
          fprintf(fp, "%e %e\n",
                  ib*dt+tmin, Itime[ib]*zlongit->macroParticleCharge*particleRelSign/dt);
        fclose(fp);
      }
#endif

      /* Take the FFT of I(t) to get I(f) */
      memcpy(Ifreq, Itime, 2*zlongit->n_bins*sizeof(*Ifreq));
      realFFT(Ifreq, nb, 0);

      /* Compute V(f) = Z(f)*I(f), putting in a factor 
       * to normalize the current waveform.
       */
      Vfreq = Vtime;
      factor = zlongit->macroParticleCharge*particleRelSign/dt*zlongit->factor*rampFactor;
      Z = zlongit->Z;
      Vfreq[0] = Ifreq[0]*Z[0]*factor;
      nfreq = nb/2 + 1;
      if (nb%2==0)
        /* Nyquist term */
        Vfreq[nb-1] = Ifreq[nb-1]*Z[nb-1]*factor;
      for (ib=1; ib<nfreq-1; ib++) {
        iImag = (iReal = 2*ib-1)+1;
        Vfreq[iReal] = (Ifreq[iReal]*Z[iReal] - Ifreq[iImag]*Z[iImag])*factor;
        Vfreq[iImag] = (Ifreq[iReal]*Z[iImag] + Ifreq[iImag]*Z[iReal])*factor; 
      }
      
      /* Compute inverse FFT of V(f) to get V(t) */
#ifdef  USE_MPE
      MPE_Log_event(event2a, 0, "start zlongit"); /* record time spent on I/O operations */
#endif
      realFFT(Vfreq, nb, INVERSE_FFT);
#ifdef  USE_MPE
      MPE_Log_event(event2b, 0, "start zlongit"); /* record time spent on I/O operations */
#endif
      Vtime = Vfreq;
      
      if (zlongit->SDDS_wake_initialized && zlongit->wakes) {
        /* wake potential output */
        factor = zlongit->macroParticleCharge*particleRelSign/dt;
        if ((zlongit->wake_interval<=0 || ((i_pass0-zlongit->wake_start)%zlongit->wake_interval)==0) &&
            i_pass0>=zlongit->wake_start && i_pass0<=zlongit->wake_end) {
          if (!SDDS_StartTable(&zlongit->SDDS_wake, nb)) {
            SDDS_SetError("Problem starting SDDS table for wake output (track_through_zlongit)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          for (ib=0; ib<nb; ib++) {
            if (!SDDS_SetRowValues(&zlongit->SDDS_wake, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ib,
                                   0, ib*dt, 1, Vtime[ib], 2, Itime[ib]*factor, -1)) {
              SDDS_SetError("Problem setting rows of SDDS table for wake output (track_through_zlongit)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (!SDDS_SetParameters(&zlongit->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                  "Pass", i_pass0, "q", zlongit->macroParticleCharge*particleRelSign*np, NULL)) {
            SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_zlongit)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          if (zlongit->broad_band) {
            if (!SDDS_SetParameters(&zlongit->SDDS_wake, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                    "Ra", zlongit->Ra, "fo", zlongit->freq, "Deltaf", zlongit->bin_size, NULL)) {
              SDDS_SetError("Problem setting parameters of SDDS table for wake output (track_through_zlongit)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (!SDDS_WriteTable(&zlongit->SDDS_wake)) {
            SDDS_SetError("Problem writing SDDS table for wake output (track_through_zlongit)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
          if (!inhibitFileSync)
            SDDS_DoFSync(&zlongit->SDDS_wake);
        }
      }

      /* put zero voltage in Vtime[nb] for use in interpolation */
      Vtime[nb] = 0;
      /* change particle momentum offsets to reflect voltage in relevant bin */
      for (ip=ip1; ip<=ip2; ip++) {
        if ((ib=pbin[ip])>=0 && ib<=nb-1) {
          if (zlongit->interpolate && ib>0) {
            /* dt/2 offset is so that center of bin is location where
             * particle sees voltage for that bin only
             */
            dt1 = time[ip]-(tmin+dt/2.0+dt*ib);
            if (dt1<0)
              dt1 += dt;
            else
              ib += 1;
            if (ib<nb)
              dgam = (Vtime[ib-1]+(Vtime[ib]-Vtime[ib-1])/dt*dt1)/(1e6*particleMassMV*particleRelSign);
            else
              continue;
          }
          else
            dgam = Vtime[ib]/(1e6*particleMassMV*particleRelSign);
          if (dgam) {
            if (zlongit->reverseTimeOrder)
              time[ip] = 2*tmean - time[ip];
            /* Put in minus sign here as the voltage decelerates the beam */
            add_to_particle_energy(part[ip], time[ip], Po, -dgam);
          }
        }
      }
    }
    

#if defined(MINIMIZE_MEMORY)
    free(Itime);
    free(Vtime);
    free(pbin);
    free(time);
    Itime = Vtime = time = NULL;
    pbin = NULL;
    max_np = max_n_bins = 0;
#endif
  }

void set_up_zlongit(ZLONGIT *zlongit, RUN *run, long pass, long particles, CHARGE *charge,
                    double timeSpan)
{
    long i, nfreq;
    double df;

    if (charge) {
      zlongit->macroParticleCharge = charge->macroParticleCharge;
    } else if (pass==0) {
      zlongit->macroParticleCharge = 0;
      if (zlongit->charge<0)
        bombElegant("ZLONGIT charge parameter should be non-negative.  Use change_particle to set particle charge state.", NULL);
#if (!USE_MPI)
      if (particles)
        zlongit->macroParticleCharge = zlongit->charge/particles;
#else
      if (USE_MPI) {
	long particles_total;

	MPI_Allreduce(&particles, &particles_total, 1, MPI_LONG, MPI_SUM, workers);
	if (particles_total)
	  zlongit->macroParticleCharge = zlongit->charge/particles_total;  
      } 
#endif
    }

    if (zlongit->initialized)
      return ;

    if (zlongit->broad_band) {
      /* compute impedance for a resonator.  Recall that I use V(t) = Vo*exp(i*w*t) convention,
       * so the impedance is Z(w) = (Ra/2)*(1 + i*T)/(1+T^2), where T=Q*(wo/w-w/wo).
       * The imaginary and real parts are positive for small w.
       */
        double term;
        if (zlongit->bin_size<=0)
          bombElegant("bin_size must be positive for ZLONGIT element", NULL);
        if (zlongit->Ra && zlongit->Rs) 
          bombElegant("ZLONGIT element broad-band resonator may have only one of Ra or Rs nonzero.  Ra is just 2*Rs", NULL);
        if (!zlongit->Ra)
          zlongit->Ra = 2*zlongit->Rs;
        if (zlongit->n_bins%2!=0)
            bombElegant("ZLONGIT element must have n_bins divisible by 2", NULL);
        if (zlongit->Zreal  || zlongit->Zimag) 
            bombElegant("can't specify both broad_band impedance and Z(f) files for ZLONGIT element", NULL);
        optimizeBinSettingsForImpedance(timeSpan, zlongit->freq, zlongit->Q,
                                        &(zlongit->bin_size), &(zlongit->n_bins), zlongit->max_n_bins);
        df = 1/(zlongit->n_bins*zlongit->bin_size)/(zlongit->freq);
        nfreq = zlongit->n_bins/2 + 1;
        fprintf(stdout, "ZLONGIT has %ld frequency points with df=%e\n",
                nfreq, df);
        fflush(stdout);
        zlongit->Z = tmalloc(sizeof(*(zlongit->Z))*zlongit->n_bins);
        zlongit->Z[0] = 0;
        zlongit->Z[zlongit->n_bins-1] = 0;    /* Nyquist term */
        for (i=1; i<nfreq-1; i++) {
            term = zlongit->Q*(1.0/(i*df)-i*df);
            zlongit->Z[2*i-1] =  zlongit->Ra/2/(1+sqr(term));
            zlongit->Z[2*i  ] =  zlongit->Z[2*i-1]*term;
            }
        if (0) {
            FILE *fp;
            fp = fopen("zbb.debug", "w");
            fputs("SDDS1\n&column name=Index, type=long &end\n", fp);
            fputs("&column name=zReal, type=double &end\n", fp);
            fputs("&column name=zImag, type=double &end\n", fp);
            fputs("&data mode=ascii, no_row_counts=1 &end\n", fp);
            fprintf(fp, "0 %e 0\n", zlongit->Z[0]);
            for (i=1; i<nfreq-1; i++)
                fprintf(fp, "%ld %e %e\n",
                        i, zlongit->Z[2*i-1], zlongit->Z[2*i]);
            fprintf(fp, "%ld %e 0\n", i, zlongit->Z[zlongit->n_bins-1]);
            fclose(fp);
            }
        df *= zlongit->freq;
        }
    else {
        TABLE Zr_data, Zi_data;
        double *Zr=NULL, *Zi=NULL;
        double df_spect=0;
        long n_spect=0;
        if (!zlongit->Zreal && !zlongit->Zimag)
            bombElegant("you must either give broad_band=1, or Zreal and/or Zimag (ZLONGIT)", NULL);
        if (zlongit->Zreal && !getTableFromSearchPath(&Zr_data, zlongit->Zreal, 1, 0))
            bombElegant("unable to read real impedance function (ZLONGIT)", NULL);
        if (zlongit->Zimag && !getTableFromSearchPath(&Zi_data, zlongit->Zimag, 1, 0))
            bombElegant("unable to read imaginary impedance function (ZLONGIT)", NULL);
        if (zlongit->Zreal && !zlongit->Zimag) {
            if (!checkPointSpacing(Zr_data.c1, Zr_data.n_data, 1e-6))
                bombElegant("frequency values not equally spaced for real data (ZLONGIT)",  NULL);
            Zr = Zr_data.c2;
            if ((n_spect = Zr_data.n_data)<2)
                bombElegant("too little data in real impedance input file (ZLONGIT)", NULL);
            df_spect = Zr_data.c1[1]-Zr_data.c1[0];
            Zi = tmalloc(sizeof(*Zi)*n_spect);
            for (i=0; i<n_spect; i++)
                Zi[i] = 0;
            }
        else if (zlongit->Zimag && !zlongit->Zreal) {
            if (!checkPointSpacing(Zi_data.c1, Zi_data.n_data, 1e-6))
                bombElegant("frequency values not equally spaced for real data (ZLONGIT)",  NULL);
            Zi = Zi_data.c2;
            if ((n_spect = Zi_data.n_data)<2)
                bombElegant("too little data in imaginary impedance input file (ZLONGIT)", NULL);
            df_spect = Zi_data.c1[1]-Zi_data.c1[0];
            Zr = tmalloc(sizeof(*Zr)*n_spect);
            for (i=0; i<n_spect; i++)
                Zr[i] = 0;
            }
        else if (zlongit->Zimag && zlongit->Zreal) {
            if (!checkPointSpacing(Zr_data.c1, Zr_data.n_data, 1e-6))
                bombElegant("frequency values not equally spaced for real data (ZLONGIT)",  NULL);
            if (!checkPointSpacing(Zi_data.c1, Zi_data.n_data, 1e-6))
                bombElegant("frequency values not equally spaced for real data (ZLONGIT)",  NULL);
            if (Zi_data.n_data!=Zr_data.n_data)
                bombElegant("real and imaginary impedance files have different amounts of data (ZLONGIT)", NULL);
            n_spect = Zi_data.n_data;
            df_spect = Zi_data.c1[1]-Zi_data.c1[0];
            if (df_spect!=(Zi_data.c1[1]-Zi_data.c1[0]))
                bombElegant("real and imaginary impedance files have different frequency spacing (ZLONGIT)", NULL);
            Zi = Zi_data.c2;
            Zr = Zr_data.c2;
            }
        if (Zi[0])
            bombElegant("impedance spectrum has non-zero imaginary DC term (ZLONGIT)", NULL);
        if (!power_of_2(n_spect-1))
            bombElegant("number of spectrum points must be 2^n+1, n>1 (ZLONGIT)", NULL);
        zlongit->n_bins = 2*(n_spect-1);
        zlongit->bin_size = 1.0/(zlongit->n_bins*df_spect);
        nfreq = n_spect;
        fprintf(stdout, "Using Nb=%ld and dt=%e s (span of %e s) in ZLONGIT\n",
                zlongit->n_bins, zlongit->bin_size, zlongit->n_bins*zlongit->bin_size);
        fflush(stdout);
        zlongit->Z = tmalloc(sizeof(*zlongit->Z)*2*zlongit->n_bins);
        for (i=0; i<n_spect; i++) {
            if (i==0)
                /* DC term */
                zlongit->Z[0] = Zr[0];
            else if (i==n_spect-1 && zlongit->n_bins%2==0)
                /* Nyquist term */
                zlongit->Z[2*i-1] = Zr[i];
            else {
                zlongit->Z[2*i-1] = Zr[i];
                zlongit->Z[2*i  ] = Zi[i];
                }
            }
        df = df_spect;
      }

    if (zlongit->SDDS_wake_initialized && !SDDS_Terminate(&zlongit->SDDS_wake)) {
        SDDS_SetError("Problem terminating SDDS output (set_up_zlongit)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    zlongit->SDDS_wake_initialized = 0;

#if (!USE_MPI)  
    /* Only the serial version will dump this part of output */
    if (zlongit->wakes) {
        zlongit->wakes = compose_filename(zlongit->wakes, run->rootname);
        if (zlongit->broad_band) 
          SDDS_ElegantOutputSetup(&zlongit->SDDS_wake, zlongit->wakes, SDDS_BINARY, 1, "longitudinal wake",
                                  run->runfile, run->lattice, wake_parameter, BB_WAKE_PARAMETERS,
                                  wake_column, WAKE_COLUMNS, "set_up_zlongit", SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
        else {
          SDDS_ElegantOutputSetup(&zlongit->SDDS_wake, zlongit->wakes, SDDS_BINARY, 1, "longitudinal wake",
                                  run->runfile, run->lattice, wake_parameter, NBB_WAKE_PARAMETERS,
                                  wake_column, WAKE_COLUMNS, "set_up_zlongit", SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
        }
        zlongit->SDDS_wake_initialized = 1;
      }
#endif

    if (zlongit->highFrequencyCutoff0>0)
      applyLowPassFilterToImpedance(zlongit->Z, nfreq,
                                    zlongit->highFrequencyCutoff0, 
                                    zlongit->highFrequencyCutoff1);

#if 0
    if (!zlongit->initialized) {
      FILE *fp;
      fp = fopen_e("zlongit.sdds", "w", 0);
      fprintf(fp, "SDDS1\n&column name=f units=Hz type=double &end\n");
      fprintf(fp, "&column name=ZReal type=double &end\n");
      fprintf(fp, "&column name=ZImag type=double &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (i=0; i<nfreq; i++) 
        if (i==0)
          fprintf(fp, "%21.15e %21.15e %21.15e\n",
                  i*df, zlongit->Z[0], 0.0);
        else 
          fprintf(fp, "%21.15e %21.15e %21.15e\n",
                  i*df, zlongit->Z[2*i-1], zlongit->Z[2*i]);
      fclose(fp);
    }
#endif
  
    zlongit->initialized = 1;
  }

void applyLowPassFilterToImpedance(double *Z, long nfreq, double cutoff0, double cutoff1)
{
  long i;
  double f;

  for (i=1; i<nfreq-1; i++) {
    f = (i*1.0)/nfreq;
    if (f<cutoff0)
      continue;
    else if (f>cutoff1 || cutoff1<=cutoff0)
      Z[2*i-1] = Z[2*i] = 0;
    else {
      Z[2*i-1] *= 1-(f-cutoff0)/(cutoff1-cutoff0);
      Z[2*i  ] *= 1-(f-cutoff0)/(cutoff1-cutoff0);
    }
  }
}


long checkPointSpacing(double *x, long n, double tolerance)
{
  double dx, dx0, range;
  long i;
  
  if (n<3)
    return 1;
  if ((range = x[n-1] - x[0])<=0)
    return 0;
  dx0 = (x[1] - x[0])/range;
  for (i=1; i<n-1; i++) {
    dx = (x[i+1]-x[i])/range;
    if (fabs(dx-dx0)>tolerance)
      return 0;
  }
  return 1;
}

