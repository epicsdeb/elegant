/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: do_tracking()
 * purpose: track a collection of particles through a beamline
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#ifdef USE_GSL
#include "gsl/gsl_poly.h"
#endif
/* #include "smath.h" */
void flushTransverseFeedbackDriverFiles(TFBDRIVER *tfbd);
void set_up_frfmode(FRFMODE *rfmode, char *element_name, double element_z, long n_passes,  RUN *run, long n_particles, double Po, double total_length);
void track_through_frfmode(double **part, long np, FRFMODE *rfmode, double Po,char *element_name, double element_z, long pass, long n_passes,CHARGE *charge);
void set_up_ftrfmode(FTRFMODE *rfmode, char *element_name, double element_z, long n_passes,RUN *run, long n_particles,double Po, double total_length);
void track_through_ftrfmode(double **part, long np, FTRFMODE *trfmode, double Po,char *element_name, double element_z, long pass, long n_passes,CHARGE *charge);
void transformEmittances(double **coord, long np, double pCentral, EMITTANCEELEMENT *ee);

ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr);
void trackLongitudinalOnlyRing(double **part, long np, VMATRIX *M, double *alpha);
void store_fitpoint_matrix_values(MARK *fpt, char *name, long occurence, VMATRIX *M);
long trackWithIndividualizedLinearMatrix(double **particle, long particles,
                                    double **accepted, double Po, double z,
                                    ELEMENT_LIST *eptr,
                                    TWISS *twiss0, double *tune0,
                                    double *chrom, double *chrom2, double *chrom3,
                                    double *dbeta_dPoP, double *dalpha_dPoP,
                                    double *alphac, double *eta2, 
                                    ILMATRIX *ilmat);
void matr_element_tracking(double **coord, VMATRIX *M, MATR *matr,
                           long np, double z);
void ematrix_element_tracking(double **coord, VMATRIX *M, EMATRIX *matr,
                              long np, double z, double *Pcentral);
void distributionScatter(double **part, long np, double Po, DSCATTER *scat, long iPass);
void recordLostParticles(double **lossBuffer, double **coord, long *nLost, long nLeft,long  nToTrack, long pass);
void storeMonitorOrbitValues(ELEMENT_LIST *eptr, double **part, long np);
void mhist_table(ELEMENT_LIST *eptr0, ELEMENT_LIST *eptr, long step, long pass, double **coord, long np, 
                double Po, double length, double charge, double z);
void set_up_mhist(MHISTOGRAM *mhist, RUN *run, long occurence);
void findMinMax (double **coord, long np, double *min, double *max, double *c0, double Po);

void interpolateFTable(double *B, double *xyz, FTABLE *ftable);
void rotate_coordinate(double **A, double *x, long inverse);
void ftable_frame_converter(double **coord, long np, FTABLE *ftable, long entrance_exit);
double choose_theta(double rho, double x0, double x1, double x2);
void track_through_multipole_deflector(
                                double **final, 
                                MRFDF *rf_param,
                                double **initial,
                                long n_particles,
                                double pc_central
                                );

#if USE_MPI
typedef enum balanceMode {badBalance, startMode, goodBalance} balance;
void scatterParticles(double **coord, long *nToTrack, double **accepted,
                      long n_processors, int myid, balance balanceStatus, 
                      double my_rate, double nParPerElements, double round,
                      int lostSinceSeqMod,int *distributed, 
                      long *reAllocate, double *P_central);
void gatherParticles(double ***coord, long **lostOnPass, long *nToTrack, 
                     long *nLost, double ***accepted, long n_processors, 
                     int myid, double *round);
/* Avoid unnecessary communications by checking if an operation will be executed in advance*/
int usefulOperation (ELEMENT_LIST *eptr, unsigned long flags, long i_pass);
balance checkBalance(double my_wtime, int myid, long n_processors);
#endif

#ifdef SORT   
int comp_IDs(const void *coord1, const void *coord2);
#endif
static TRACKING_CONTEXT trackingContext;


double beta_from_delta(double p, double delta)
{
  p *= 1+delta;
  return( p/sqrt(p*p+1));
}


/* This is used if one needs to wedge a function into the lattice at a specific
 * location
 */

static void (*trackingWedgeFunction)(double **part, long np, long pass, double *pCentral) = NULL;
static ELEMENT_LIST *trackingWedgeElement = NULL;
void setTrackingWedgeFunction(void (*wedgeFunc)(double **part, long np, long pass, double *pCentral),
                              ELEMENT_LIST *eptr)
{
  trackingWedgeFunction = wedgeFunc;
  trackingWedgeElement = eptr;
}

long do_tracking(
                 /* Either the beam pointer or the coord pointer must be supplied, but not both */
                 BEAM *beam,  
                 double **coord,
                 long nOriginal,   /* Used only if coord is supplied */
                 long *effort,
                 LINE_LIST *beamline,
                 double *P_central,    /* beta*gamma for central particle */
                 double **accepted,
                 BEAM_SUMS **sums_vs_z,
                 long *n_z_points,
                 TRAJECTORY *traj_vs_z,
                 RUN *run,
                 long step,
                 unsigned long flags,
                 long n_passes,
                 long passOffset,
                 SASEFEL_OUTPUT *sasefel,
		 SLICE_OUTPUT *sliceAnalysis,
                 double *finalCharge,
		 double **lostParticles,
		 ELEMENT_LIST *startElem
                 )
{
  RFMODE *rfmode; TRFMODE *trfmode;
  FRFMODE *frfmode; FTRFMODE *ftrfmode;
  WATCH *watch;
  STRAY *stray;
  HISTOGRAM *histogram;
  MHISTOGRAM *mhist;
  FTABLE *ftable;
  ENERGY *energy;
  MAXAMP *maxamp;
  MALIGN *malign;
  ELEMENT_LIST *eptr, *eptrPred, *eptrCLMatrix=NULL;
  long nToTrack;  /* number of particles being tracked */
  long nLeft;     /* number of those that are left after a tracking routine returns */
  long nLost=0;     /* accumulated number lost */
  long nMaximum=0;  /* maximum number of particles seen */
  long show_dE, maxampOpenCode=0, maxampExponent=0, maxampYExponent=0;
  double dgamma, dP[3], z, z_recirc, last_z;
  long i, j, i_traj=0, i_sums, i_pass, isConcat;
  long i_sums_recirc, saveISR=0;
  long watch_pt_seen, feedbackDriverSeen;
  double sum, x_max, y_max;
  long elliptical;
  double et1, et2=0;
  long is_batch = 0, last_type;
  static long is_ansi_term = -1;
  char s[100], *name;
  long check_nan, sums_allocated = 0;
  long elementsTracked, sliceAnDone = 0;
  CHARGE *charge;
  static long warnedAboutChargePosition = 0;
  unsigned long classFlags = 0;
  long nParticlesStartPass = 0;
  int myid = 0, active = 1, lostSinceSeqMode = 0, needSort = 0;
#ifdef SORT
  int nToTrackAtLastSort;
#endif
#if USE_MPI 
  long old_nToTrack = 0, nParElements, nElements; 
  int checkFlags;
  double my_wtime, start_wtime, end_wtime, nParPerElements, my_rate;
  double round = 0.5;
  balance balanceStatus;
#if SDDS_MPI_IO
  int distributed = 1;
  long total_nOriginal;
  long total_nToTrack;
#else
  int distributed = 0; /* indicate if the particles have been scattered */
#endif
  long reAllocate = 0; /* indicate if new memory needs to be allocated */
#ifdef  USE_MPE /* use the MPE library */
  int event1a, event1b, event2a, event2b;
  MPE_LOG_BYTES  bytebuf;
  int            bytebuf_pos = 0;
  event1a = MPE_Log_get_event_number();
  event1b = MPE_Log_get_event_number();
  event2a = MPE_Log_get_event_number(); 
  event2b = MPE_Log_get_event_number();
  if(isMaster) {
    MPE_Describe_state(event1a, event1b, "Watch", "red");
    MPE_Describe_info_state( event2a, event2b, "Tracking_element", "orange",
			     "Element: %s" );
  }
#endif
  balanceStatus = startMode;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* get ID number for each processor */
  trackingContext.myid = myid;
  if (myid==0) 
    my_rate = 0.0;
  else
    my_rate = 1.0;
  if (notSinglePart && partOnMaster) /* This is a special case when the first beam is fiducial. We need scatter the beam in the second step. */
    distributed = 0;
#endif 
  strncpy(trackingContext.rootname, run->rootname, CONTEXT_BUFSIZE);
  if (!coord && !beam)
    bombElegant("Null particle coordinate array and null beam pointer! (do_tracking)", NULL);
  if (coord && beam)
    bombElegant("Particle coordinate array and beam pointer both supplied!  (do_tracking)", NULL);
  if (beam) {
    coord = beam->particle;
    nOriginal = beam->n_to_track;  /* used only for computing macroparticle charge */
  }

#if SDDS_MPI_IO
  if (notSinglePart && !partOnMaster) {
    if (isMaster )
      nOriginal = 0; 
    MPI_Allreduce(&nOriginal, &total_nOriginal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  else /* single partticle case where all the processors track the same particle(s), or particles are on master when the first beam is fiducial */
    total_nOriginal = nOriginal;
#endif
  
#ifdef WATCH_MEMORY
  fprintf(stdout, "start do_tracking():  CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
          cpu_time()/100.0, page_faults(), memory_count());
  fflush(stdout);
#endif
  
#if defined(UNIX) || defined(_WIN32)
  if (is_ansi_term==-1) {
    char *ptr;
    is_ansi_term = 1;
    if (!(ptr=getenv("TERM")))
      is_ansi_term = 0;
    else if (strcmp(ptr, "emacs")==0)
      is_ansi_term = 0;
  }
#endif

#if SDDS_MPI_IO
  if (isSlave || (!notSinglePart))
#else 
  if (isMaster)
#endif
    if (accepted)
      copy_particles(accepted, coord, nOriginal);

#ifdef VAX_VMS
  is_batch = job_mode(getpid())==2?1:0;
#endif
  
  z = z_recirc = last_z =  0;
  i_sums = i_sums_recirc = 0;
  x_max = y_max = 0;
  nToTrack = nLeft = nMaximum = nOriginal;
#if USE_MPI
  if (!partOnMaster && notSinglePart) {
    if (isMaster) nToTrack = 0; 
    MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  } else { /* singlePart tracking or partOnMaster */
    if (beam)	
      beam->n_to_track_total = nToTrack;
  }
#endif
  et1 = -2.0;
  elliptical = isConcat = 0;
  watch_pt_seen = feedbackDriverSeen = 0;

#ifdef SORT
  nToTrackAtLastSort = nToTrack;
#endif
  
  check_nan = 1;
  eptr = &(beamline->elem);

  if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN)) {
    /* this is required just in case rf elements etc. have previously
     * been activated by computation of correction matrices or trajectory
     * correction.
     */
    delete_phase_references();
    reset_special_elements(beamline, 1);
  }
  reset_driftCSR();

  while (eptr) {
    if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN))
      eptr->Pref_output_fiducial = 0;
    eptr = eptr->succ;
  }
  if (!(flags&FIDUCIAL_BEAM_SEEN) && flags&PRECORRECTION_BEAM)
    flags &= ~FIRST_BEAM_IS_FIDUCIAL; 
  if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN) && !(flags&SILENT_RUNNING)) {
    fprintf(stdout, "This step establishes energy profile vs s (fiducial beam).\n");
    fflush(stdout);
  }
  
  log_exit("do_tracking.1");
  log_entry("do_tracking.2");
  name = "_BEG_";
  last_type = sums_allocated = 0;
  charge = NULL;
  if (finalCharge)
    *finalCharge = 0;  

#if SDDS_MPI_IO
  if (isSlave || (!notSinglePart) || partOnMaster) {
#else
  if (isMaster) {   /* As the particles have not been distributed, only master needs to do these computation */
#endif 
    if (check_nan) {
      nLeft = limit_amplitudes(coord, DBL_MAX, DBL_MAX, nToTrack, accepted, z, *P_central, 0,
					  0);
      if (nLeft!=nToTrack)
	recordLostParticles(lostParticles, coord, &nLost, nLeft, nToTrack, 0);
      nToTrack = nLeft;
    }
    if (run->apertureData.initialized)  {
      nLeft = imposeApertureData(coord, nToTrack, accepted, 0.0, *P_central, &(run->apertureData));
      if (nLeft!=nToTrack)
	recordLostParticles(lostParticles, coord, &nLost, nLeft, nToTrack, 0);
      nToTrack = nLeft;
    }
  }
  
  for (i_pass=passOffset; i_pass<n_passes+passOffset; i_pass++) {
    log_entry("do_tracking.2.1");
    if (run->stopTrackingParticleLimit>0) {
#if !USE_MPI
      if (nToTrack<run->stopTrackingParticleLimit) 
#else
      MPI_Allreduce(&nToTrack, &total_nToTrack, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (total_nToTrack<run->stopTrackingParticleLimit)
#endif
      {
	/* force loss of all the particles */
	for (i=0; i<nToTrack; i++) {
	  coord[i][4] = z;
	  coord[i][5] = *P_central*(1+coord[i][5]);
	}
	nLeft = 0;
	if (nLeft!=nToTrack)
	  recordLostParticles(lostParticles, coord, &nLost, nLeft, nToTrack, i_pass);
	nToTrack = 0;
      }
    }

    ResetNoiseGroupValues();
    if (applyElementModulations(&(run->modulationData), *P_central, coord, nToTrack, run, i_pass)) {
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
    }
    if (applyElementRamps(&(run->rampData), *P_central, run, i_pass)) {
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
    }
    
    if (beamline->links) {
      sprintf(s, "%.15e sto p_central  %ld sto turn", *P_central, i_pass);
      rpn(s);
      rpn_clear();
      if (rpn_check_error()) exitElegant(1);
      if (assert_element_links(beamline->links, run, beamline, TURN_BY_TURN_LINK)) {
        beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
        beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
      }
    }
    if (beamline->flags&BEAMLINE_TWISS_WANTED && !(beamline->flags&BEAMLINE_TWISS_CURRENT)
        && !(flags&TEST_PARTICLES)) {
      update_twiss_parameters(run, beamline, NULL);
    }
    if (run->concat_order && !(flags&TEST_PARTICLES) && 
        !(beamline->flags&BEAMLINE_CONCAT_CURRENT) ) {
      /* form concatenated beamline for tracking */
      if (getSCMULTSpecCount())
        bombElegant("space charge calculation can not work together with matrix concatenation tracking. \n Please remove concat_order from run_setup", NULL);
      concatenate_beamline(beamline, run);
    }

    if (run->concat_order && beamline->flags&BEAMLINE_CONCAT_DONE &&
        !(flags&TEST_PARTICLES)) {
      if (beamline->ecat_recirc && (i_pass || flags&BEGIN_AT_RECIRC))
        eptr = beamline->ecat_recirc;
      else
        eptr = &(beamline->ecat);
      isConcat = 1;
    }
    else if (beamline->elem_recirc && (i_pass || flags&BEGIN_AT_RECIRC))
      eptr = beamline->elem_recirc;
    else
      eptr = &(beamline->elem);

    if (i_pass==0) {
      if (flags&LINEAR_CHROMATIC_MATRIX) {
        if (!isConcat) {
          fprintf(stdout, "Error: in order to use the \"linear chromatic matrix\" for\n");
          fflush(stdout);
          fprintf(stdout, "tracking, you must ask for matrix concatenation in the run_setup.\n");
          fflush(stdout);
          exitElegant(1);
        }
        eptrCLMatrix = findBeamlineMatrixElement(eptr);
      }
      if (flags&LONGITUDINAL_RING_ONLY) {
        if (!isConcat) {
          fprintf(stdout, "Error: in order to use the \"longitudinal ring\" mode of\n");
          fflush(stdout);
          fprintf(stdout, "tracking, you must ask for matrix concatenation in the run_setup.\n");
          fflush(stdout);
          exitElegant(1);
        }
        eptrCLMatrix = findBeamlineMatrixElement(eptr);
      }
    }
      
    if (sums_vs_z && n_z_points) {
      if (!sums_allocated && !*sums_vs_z) {
        /* allocate storage for beam sums */
        if (!isConcat)
          *n_z_points = beamline->n_elems + 1 + 
            (run->wrap_around?0:(n_passes-1)*(beamline->n_elems-beamline->i_recirc));
        else
          *n_z_points = beamline->ncat_elems + 1 +
            (run->wrap_around?0:(n_passes-1)*(beamline->ncat_elems-beamline->i_recirc));
        if (flags&FINAL_SUMS_ONLY)
          *n_z_points = 0;
        *sums_vs_z = tmalloc(sizeof(**sums_vs_z)*(*n_z_points+1));
        zero_beam_sums(*sums_vs_z, *n_z_points+1);
        sums_allocated = 1;
      }
      else if (!run->combine_bunch_statistics && i_pass==0)
        zero_beam_sums(*sums_vs_z, *n_z_points+1);
    }

    if (run->wrap_around) {
      i_sums = i_sums_recirc;  /* ==0 for i_pass==0 */
      z = z_recirc;            /* ditto */
      last_z = z;
    }
    if (run->final_pass && sums_vs_z && n_z_points)
      zero_beam_sums(*sums_vs_z, *n_z_points+1);

    log_exit("do_tracking.2.1");
    log_entry("do_tracking.2.2");
    if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)
        && !(run->tracking_updates==0)) {
#if defined(VAX_VMS)
      sprintf(s, "%ld particles present after pass %ld        ",
              nToTrack, i_pass);
      fputs(s, stdout);
      if (is_ansi_term)
        backspace(strlen(s));
      else
        fputc('\n', stdout);
      fflush(stdout);
      et1 = et2;
#endif
#if defined(UNIX) || defined(_WIN32)
#if !SDDS_MPI_IO
      if ((et2=delapsed_time())-et1>2.0) {
        sprintf(s, "%ld particles present after pass %ld        ", 
                nToTrack, i_pass);
#else
      if (i_pass%20==0) {
	if (!partOnMaster && notSinglePart) {
	  /* We have to collect information from all the processors to print correct info during tracking */
	  if (isMaster) nToTrack = 0; 
	  MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	  sprintf(s, "%ld particles present after pass %ld        ", 
		  beam->n_to_track_total, i_pass);
	}
	else { /* singlePart tracking or partOnMaster */
	  sprintf(s, "%ld particles present after pass %ld        ", 
		  nToTrack, i_pass);
	}
#endif
        fputs(s, stdout);
        if (is_ansi_term)
          backspace(strlen(s));
        else
          fputc('\n', stdout);
        fflush(stdout);
        et1 = et2;
      }
#else
      sprintf(s, "%ld particles present after pass %ld        ", 
              nToTrack, i_pass);
      fputs(s, stdout);
      if (is_ansi_term)
        backspace(strlen(s));
      else
        fputc('\n', stdout);
      fflush(stdout);
#endif 
    }
    elementsTracked = -1;
    eptrPred = eptr;
#if USE_MPI
    if (notSinglePart) {
      my_wtime = 0.0;
      nParElements = 0;
      nElements = 0;
    }
#endif
    nParticlesStartPass = nToTrack;

    if (getSCMULTSpecCount()) 
      /* prepare space charge effects calculation  */
      initializeSCMULT(eptr, coord, nToTrack, *P_central, i_pass);

    if (i_pass==0 && startElem) {
      /* start tracking from an interior point in the beamline */
      while (eptr && eptr!=startElem) {
        if (eptr->type==T_MAXAMP) {
          maxamp = (MAXAMP*) eptr->p_elem;
          x_max = maxamp->x_max;
          y_max = maxamp->y_max;
          elliptical = maxamp->elliptical;
          maxampOpenCode = determineOpenSideCode(maxamp->openSide);
          maxampExponent = maxamp->exponent;
          maxampYExponent = maxamp->yExponent;
        }
        eptr = eptr->succ;
      }
      z = startElem->end_pos;
      startElem = NULL; 
    }

    while (eptr && (nToTrack || (USE_MPI && notSinglePart))) {
      if (trackingWedgeFunction && eptr==trackingWedgeElement)
        (*trackingWedgeFunction)(coord, nToTrack, i_pass, P_central);

      classFlags = entity_description[eptr->type].flags;
      elementsTracked++;
      log_entry("do_tracking.2.2.0");
      if (!eptr->name) {
        fprintf(stdout, "error: element ending at %em has NULL name pointer\n", eptr->end_pos);
        fflush(stdout);
        if (eptr->pred && eptr->pred->name) {
          fprintf(stdout, "previous element is %s\n", eptr->pred->name);
          fflush(stdout);
        }
        else if (eptr->succ && eptr->succ->name) {
          fprintf(stdout, "next element is %s\n", eptr->succ->name);
          fflush(stdout);
        }
        abort();
      }
      if (!eptr->p_elem && !run->concat_order) {
        fprintf(stdout, "element %s has NULL p_elem pointer", eptr->name);
        fflush(stdout);
        exitElegant(1);
      }
      if (eptr->type<=0 || eptr->type>=N_TYPES) {
        fprintf(stdout, "element %s has type %ld--not recognized/not allowed\n", eptr->name, eptr->type);
        fflush(stdout);
        exitElegant(1);
      }
      log_exit("do_tracking.2.2.0");

      log_entry("do_tracking.2.2.1");

#ifdef SORT
      if (!USE_MPI || needSort)
	if (nToTrackAtLastSort > nToTrack) {/* indicates more particles are lost, need sort */
          if (beam && beam->bunchFrequency!=0)
            fprintf(stdout, "*** Warning: particle ID sort not being performed because bunch frequency is nonzero\n");
          else {
            qsort(coord[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
            if (accepted!=NULL)
              qsort(accepted[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
            nToTrackAtLastSort = nToTrack;
            needSort = 0;
          }   
        }
#endif   
      if (sums_vs_z && *sums_vs_z && (!run->final_pass || i_pass==n_passes-1) && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES)) {
        if (i_sums<0)
          bombElegant("attempt to accumulate beam sums with negative index!", NULL);
        accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
        (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
        fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
                i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
        fflush(stdout);
#endif
        i_sums++;
      }
#if USE_MPI
      if (notSinglePart) {
	active = 0;
	if (classFlags&UNIPROCESSOR) {
	  /* This element cannot be done in parallel. Only the master CPU will work. */
	  if (myid == 0)
	    active = 1;
	  else if (eptr->type==T_SCRIPT) /* The slave processors will be notified if they need allocate new memory */
	    active = 1; 
	  else 
	    active = 0;
	  if (parallelStatus==trueParallel) {
	    if (!partOnMaster) {
              if(usefulOperation(eptr, flags, i_pass)) {
		printf("Warning: %s (%s) is a serial element. It is not recommended for the simulation with a large number of particles because of memory issue.\n", eptr->name, entity_name[eptr->type]);
		gatherParticles(&coord, NULL, &nToTrack, &nLost, &accepted, 
				n_processors, myid, &round);
		if (isMaster)
		  nMaximum = nToTrack;
                partOnMaster = 1;
	      }
	      /* update the nMaximum for recording the nLost on all the slave processors */
	      if (myid!=0)
		nMaximum = nToTrack;
	    }

	    /* The element will change the state of particles. Scatter is required
	       for parallel computation */
	    if (!(classFlags&(UNIDIAGNOSTIC&(~UNIPROCESSOR)))) {
	      parallelStatus = notParallel;           
	    }
	  }     
	} 
	else {
	  /* This element can be done in parallel. Only the slave CPUS will work. */
	  if (myid != 0)
	    active = 1; 
	  else {  /* myid == 0 */
	    if (!(classFlags&MPALGORITHM)) {
	      active = 0;
	    }
	    else /* The master CPU needs to participate communications */
	      active = 1;
	  }
	  if ((balanceStatus==badBalance) && (parallelStatus==trueParallel)) {
	    gatherParticles(&coord, NULL, &nToTrack, &nLost, &accepted, n_processors, myid, &round);
	    nMaximum = nToTrack;
	  } 
	  /* Particles will be scattered in startMode, bad balancing status or notParallel state */  
	  if ((balanceStatus==badBalance) || (parallelStatus==notParallel)) {
	    scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
			     balanceStatus, my_rate, nParPerElements, round, lostSinceSeqMode, &distributed, &reAllocate, P_central);
	    if (myid != 0) {
	      /* update the nMaximum for recording the nLost on all the slave processors */
	      nMaximum = nToTrack;  
	    }
	    if (balanceStatus!=startMode)
	      balanceStatus = goodBalance;
	  }
	  else if (balanceStatus==startMode) { 
	    /* For the first pass, scatter when it is not in parallel mode */
	    if (parallelStatus!=trueParallel) {
	      scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
			       balanceStatus, my_rate, nParPerElements, round, lostSinceSeqMode, &distributed, &reAllocate, P_central);
	      if (myid != 0) {
		/* update the nMaximum for recording the nLost on all the slave processors */
		nMaximum = nToTrack; 
	      }
	    }
	  }
	  parallelStatus = trueParallel;         
	  lostSinceSeqMode = 0;
	  if (myid != 0) {
	    if (!(classFlags&MPALGORITHM)) {
	      /* We do not count time spent on those elements which need collective communications,
		 as it will be the time spent on synchronization instead of computation */ 
	      nElements++;
	      /* count the total number of particles tracked by all of the elements for each pass */
	      nParElements = nParElements+nToTrack;
	      start_wtime = MPI_Wtime();  
	    }
	  }
	  partOnMaster = 0;
	}
      } 
#endif

      name = eptr->name;
      last_z = z;
      if (entity_description[eptr->type].flags&HAS_LENGTH && eptr->p_elem)
        z += ((DRIFT*)eptr->p_elem)->length;
      else {
        if (eptr->pred)
          z += eptr->end_pos - eptr->pred->end_pos;
        else
          z += eptr->end_pos;
      }
      /* fill a structure that can be used to pass to other routines 
       * information on the tracking context 
       */
      strncpy(trackingContext.elementName, eptr->name, CONTEXT_BUFSIZE);
      trackingContext.elementOccurrence = eptr->occurence;
      trackingContext.sliceAnalysis = sliceAnalysis?
	(sliceAnalysis->finalValuesOnly?NULL:sliceAnalysis):NULL;
      trackingContext.zStart = last_z;
      trackingContext.zEnd = z;
      trackingContext.step = step;
      trackingContext.elementType = eptr->type;

      log_exit("do_tracking.2.2.1");
      if (eptr->p_elem || eptr->matrix) {
#ifdef VAX_VMS
        if (run->print_statistics && !(flags&TEST_PARTICLES))
          fprintf(stdout, "tracking through %s%c", eptr->name, ' ');
	fflush(stdout);
#else
        if (run->print_statistics && !(flags&TEST_PARTICLES))
          fprintf(stdout, "tracking through %s%c", eptr->name, '\n');
	fflush(stdout);
#endif
        show_dE = 0;
        nLeft = nToTrack;  /* in case it isn't set by the element tracking */
        if (eptr==eptrCLMatrix) {
          /* This element is the place-holder for the chromatic linear matrix or
           * the longitudinal-only matrix 
           */
          if ((!USE_MPI || !notSinglePart) || (USE_MPI && (myid!=0))) {
            /* Only the slave CPUs will work on this part */ 
            if (flags&LINEAR_CHROMATIC_MATRIX) 
              nLeft
	        = trackWithIndividualizedLinearMatrix(coord, nToTrack, accepted,
		       			         *P_central, z, eptrCLMatrix,
					         beamline->twiss0, beamline->tune,
					         beamline->chromaticity,
					         beamline->chrom2, beamline->chrom3,
					         beamline->dbeta_dPoP, beamline->dalpha_dPoP,
					         beamline->alpha, beamline->eta2, NULL);
            else 
              trackLongitudinalOnlyRing(coord, nToTrack, 
                                        eptrCLMatrix->matrix,
                                        beamline->alpha);
          }
	}
        else if (entity_description[eptr->type].flags&MATRIX_TRACKING &&
		 !(flags&IBS_ONLY_TRACKING)) {
          if (!(entity_description[eptr->type].flags&HAS_MATRIX))
            bombElegant("attempt to matrix-multiply for element with no matrix!",  NULL);
          if (!eptr->matrix) {
            if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
              bombElegant("no matrix for element that must have matrix", NULL);
          }
          if (eptr->matrix->C[5]!=0) {
            fprintf(stdout, "Warning: matrix with C5!=0 detected in matrix multiplier--this shouldn't happen!\nAll particles considered lost!\n");
            fprintf(stdout, "Element in question is %s, C5=%le\n", eptr->name, eptr->matrix->C[5]);
            fflush(stdout);
            nLeft = 0;
          } else {
            if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
              fprintf(stdout, "Tracking matrix for %s\n", eptr->name);
              fflush(stdout);
              print_elem(stdout, eptr);
              print_matrices(stdout, "", eptr->matrix);
            }
            if (flags&CLOSED_ORBIT_TRACKING) {
              switch (eptr->type) {
              case T_MONI:
              case T_HMON:
              case T_VMON:
                storeMonitorOrbitValues(eptr, coord, nToTrack);
                break;
              default:
                break;
              }
            }
            /* Only the slave CPUs will track */ 
            if ((!USE_MPI || !notSinglePart) || (USE_MPI && (myid!=0))) 
              track_particles(coord, eptr->matrix, coord, nToTrack);
          }
        }
        else {
          long type;
          if (run->print_statistics>1 && !(flags&TEST_PARTICLES)) {
            fprintf(stdout, "Tracking element: ");
            fflush(stdout);
            print_elem(stdout, eptr);
          }
	  type = eptr->type;
	  if (flags&IBS_ONLY_TRACKING) {
	    switch (type) {
	    case T_IBSCATTER:
	    case T_WATCH:
	    case T_CLEAN:
	    case T_RCOL:
	    case T_CHARGE:
	      break;
	    default:
	      type = -1;
	      break;
	    }
	  }
#ifdef  USE_MPE
	  bytebuf_pos = 0;
	  MPE_Log_event( event2a, 0, NULL );
	  MPE_Log_pack( bytebuf, &bytebuf_pos, 's', strlen(entity_name[eptr->type]), entity_name[eptr->type]); 
#endif
	  if (active && (((!USE_MPI || !notSinglePart) && nParticlesStartPass) || nToTrack || 
	     (USE_MPI && (classFlags&RUN_ZERO_PARTICLES)))) {
	    switch (type) {
	    case -1:
	      break;
	    case T_CHARGE:
	      if (i_pass==0) {
		if (elementsTracked!=0 && !warnedAboutChargePosition) {
		  warnedAboutChargePosition = 1;
		  fprintf(stdout, "Warning: the CHARGE element is not at the start of the beamline.\n");
		  fflush(stdout);
		}
		if (charge!=NULL) {
		  fprintf(stdout, "Fatal error: multipole CHARGE elements in one beamline.\n");
		  fflush(stdout);
		  exitElegant(1);
		}
		charge = (CHARGE*)eptr->p_elem;
		charge->macroParticleCharge = 0;
#if !SDDS_MPI_IO
		if (nOriginal)
		  charge->macroParticleCharge = charge->charge/(nOriginal);
#else
		if (notSinglePart) {
		  if (total_nOriginal)
		    charge->macroParticleCharge = charge->charge/(total_nOriginal);
		} else {
		  if (nOriginal)
		    charge->macroParticleCharge = charge->charge/(nOriginal);
		}
#endif
		if (charge->chargePerParticle)
		  charge->macroParticleCharge = charge->chargePerParticle;
                if (charge->macroParticleCharge<0) 
                  bombElegant("Error: CHARGE element should specify the quantity of charge (in Coulombs) without the sign", NULL);
	      }
	      break;
	    case T_MARK:
	      if (((MARK*)eptr->p_elem)->fitpoint && i_pass==n_passes-1) {
		/*
		  if (beamline->flags&BEAMLINE_TWISS_WANTED) {
		  if (!(beamline->flags&BEAMLINE_TWISS_DONE))
                  update_twiss_parameters(run, beamline, NULL);
		  store_fitpoint_twiss_parameters((MARK*)eptr->p_elem, eptr->name, 
		  eptr->occurence, eptr->twiss);
		  }
		*/
                if (isMaster || !notSinglePart)
       		  store_fitpoint_matrix_values((MARK*)eptr->p_elem, eptr->name, 
					       eptr->occurence, eptr->accumMatrix);
		store_fitpoint_beam_parameters((MARK*)eptr->p_elem, eptr->name,eptr->occurence, 
					       coord, nToTrack, *P_central); 
		if (flags&CLOSED_ORBIT_TRACKING)
		  storeMonitorOrbitValues(eptr, coord, nToTrack);
	      }
	      break;
	    case T_RECIRC:
	      /* Recognize and record recirculation point.  */
	      if (i_pass==0) {
		i_sums_recirc = i_sums-1;
		z_recirc = last_z;
	      }
	      break;
	    case T_RFDF:
	      if (!(flags&TIME_DEPENDENCE_OFF) || (flags&CLOSED_ORBIT_TRACKING))
		track_through_rf_deflector(coord, (RFDF*)eptr->p_elem,
                                           coord, nToTrack, *P_central,
                                           beamline->revolution_length, eptr->end_pos,
                                           i_pass);
	      else
		exactDrift(coord, nToTrack, ((RFDF*)eptr->p_elem)->length);
	      break;
	    case T_MRFDF:
	      if (!(flags&TIME_DEPENDENCE_OFF))
		track_through_multipole_deflector(coord, (MRFDF*)eptr->p_elem,
                                           coord, nToTrack, *P_central);
	      break;
	    case T_RFTM110:
	      if (!(flags&TIME_DEPENDENCE_OFF))
		track_through_rftm110_deflector(coord, (RFTM110*)eptr->p_elem,
						coord, nToTrack, *P_central,
						beamline->revolution_length, eptr->end_pos,
						i_pass);
	      break;
	    case T_RMDF:
	      if (!(flags&TIME_DEPENDENCE_OFF))
		track_through_ramped_deflector(coord, (RMDF*)eptr->p_elem,
					       coord, nToTrack, *P_central);
	      else
		drift_beam(coord, nToTrack, ((RMDF*)eptr->p_elem)->length, run->default_order);
	      break;
	    case T_RFTMEZ0:
              nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
              show_dE = 1;
	      break;
	    case T_TMCF:
	    case T_CEPL:
	    case T_TWPL:
	      if (!(flags&TIME_DEPENDENCE_OFF)) {
		nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			       &dgamma, dP, accepted, last_z);
		show_dE = 1;
	      }
	      else
		drift_beam(coord, nToTrack, ((TW_LINAC*)eptr->p_elem)->length, run->default_order);
	      break;
	    case T_MAPSOLENOID:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      break;
	    case T_TWLA:
	    case T_TWMTA:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      show_dE = 1;
	      break;
	    case T_RCOL:
	      if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
		drift_beam(coord, nToTrack, ((RCOL*)eptr->p_elem)->length, run->default_order);
	      else {
		nLeft = rectangular_collimator(coord, (RCOL*)eptr->p_elem, nToTrack, accepted, last_z, *P_central);
	      }
	      break;
	    case T_ECOL:
	      if (flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))
		drift_beam(coord, nToTrack, ((ECOL*)eptr->p_elem)->length, run->default_order);
	      else
		nLeft = elliptical_collimator(coord, (ECOL*)eptr->p_elem, nToTrack, accepted, last_z, *P_central);
	      break;
	    case T_CLEAN:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
		nLeft = remove_outlier_particles(coord, (CLEAN*)eptr->p_elem, 
						 nToTrack, accepted, z, *P_central);
	      break;
	    case T_SCRAPER:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))) {
		nLeft = beam_scraper(coord, (SCRAPER*)eptr->p_elem, nToTrack, accepted, last_z, *P_central);
	      }
	      break;
	    case T_PFILTER:
	      if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES)))
		nLeft = track_through_pfilter(coord, (PFILTER*)eptr->p_elem, nToTrack, 
					      accepted, z, *P_central);
	      break;
	    case T_CENTER:
	      center_beam(coord, (CENTER*)eptr->p_elem, nToTrack, i_pass, *P_central);
	      break;
	    case T_REMCOR:
	      remove_correlations(coord, (REMCOR*)eptr->p_elem, nToTrack);
	      break;
	    case T_RFCA:
	      nLeft = simple_rf_cavity(coord, nToTrack, (RFCA*)eptr->p_elem, accepted, P_central, z);
	      break;
	    case T_RFCW:
	      nLeft = track_through_rfcw(coord, nToTrack, (RFCW*)eptr->p_elem, accepted, P_central, z,
					 run, i_pass, charge);
	      break;
	    case T_MODRF:
	      modulated_rf_cavity(coord, nToTrack, (MODRF*)eptr->p_elem, *P_central, z);
	      break;
	    case T_WATCH:
#ifdef  USE_MPE
	      MPE_Log_event(event1a, 0, "start watch"); /* record time spent on I/O operations */
#endif
#if USE_MPI
	      if (!notSinglePart) /* When each processor tracks the beam independently, the watch point will be disabled in Pelegant */
		break;
	      if (!partOnMaster && notSinglePart) { /* Update the total particle number to get the correct charge */
		if (isMaster) nToTrack = 0;
		MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	      } else { /* singlePart tracking or partOnMaster */
		beam->n_to_track_total = nToTrack;
	      }
#endif
	      if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
	        watch = (WATCH*)eptr->p_elem;
	        if (!watch->disable) {
	          watch_pt_seen = 1;
	          if (!watch->initialized) 
	            set_up_watch_point(watch, run, eptr->occurence, eptr->pred?eptr->pred->name:NULL);
	          if (i_pass==0 && (n_passes/watch->interval)==0)
	            fprintf(stdout, "warning: n_passes = %ld and WATCH interval = %ld--no output will be generated!\n",
	     	     n_passes, watch->interval);
#if SDDS_MPI_IO
		  if(watch_not_allowed) {
		    dup2(fd,fileno(stdout));
		    printf("/****************************************************************************/\n");
		    printf ("Watch point can not be used for dynamic aperture searching with Pelegant\n");
		    printf("/****************************************************************************/\n");
		    fflush(stdout);
		    MPI_Abort(MPI_COMM_WORLD, 1);
		  }
#endif
		  fflush(stdout);
		  if (i_pass>=watch->start_pass && (i_pass-watch->start_pass)%watch->interval==0 &&
                      (watch->end_pass<0 || i_pass<=watch->end_pass)) {
	            switch (watch->mode_code) {
	            case WATCH_COORDINATES:
#if !SDDS_MPI_IO
		      dump_watch_particles(watch, step, i_pass, coord, nToTrack, *P_central,
			        	   beamline->revolution_length, 
					   charge?charge->macroParticleCharge*nToTrack:0.0, z);
#else
		      dump_watch_particles(watch, step, i_pass, coord, nToTrack, *P_central,
			        	   beamline->revolution_length, 
					   charge?charge->macroParticleCharge*beam->n_to_track_total:0.0, z);
#endif
		      break;
		    case WATCH_PARAMETERS:
		    case WATCH_CENTROIDS:
#if SDDS_MPI_IO
		      dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, total_nOriginal, *P_central,
					    beamline->revolution_length, z);
#else
		      dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central,
					    beamline->revolution_length, z);
#endif
		      break;
		    case WATCH_FFT:
#if SDDS_MPI_IO
		      /* This part will be done in serial for now. A parallel version of FFT could be used here */
		      if (!partOnMaster && notSinglePart) {
			printf("Warning: %s (%s FFT) is a serial element. It is not recommended for the simulation with a large number of particles because of memory issue.\n", eptr->name, entity_name[eptr->type]);
		      }
			gatherParticles(&coord, NULL, &nToTrack, &nLost, &accepted, n_processors, myid, &round);
		      if (isMaster)
#endif
		      dump_watch_FFT(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central);
#if SDDS_MPI_IO
		      if (!partOnMaster && notSinglePart) 
			scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
					 balanceStatus, my_rate, nParPerElements, round, 
					 lostSinceSeqMode, &distributed, &reAllocate, P_central);
#endif
		      break;
		    }
		  }
		}
	      }
#ifdef  USE_MPE
	      MPE_Log_event(event1b, 0, "end watch");
#endif
	      break;
	    case T_HISTOGRAM:
	      if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
		histogram = (HISTOGRAM*)eptr->p_elem;
		if (!histogram->disable) {
		  watch_pt_seen = 1;   /* yes, this should be here */
		  if (!histogram->initialized) 
		    set_up_histogram(histogram, run, eptr->occurence);
		  if (i_pass==0 && (n_passes/histogram->interval)==0)
		    fprintf(stdout, "warning: n_passes = %ld and HISTOGRAM interval = %ld--no output will be generated!\n",
			    n_passes, histogram->interval);
		  fflush(stdout);
		  if (i_pass>=histogram->startPass && (i_pass-histogram->startPass)%histogram->interval==0) {
#if !SDDS_MPI_IO   
		    dump_particle_histogram(histogram, step, i_pass, coord, nToTrack, *P_central,
					    beamline->revolution_length, 
					    charge?charge->macroParticleCharge*nToTrack:0.0, z);
#else
		    dump_particle_histogram(histogram, step, i_pass, coord, nToTrack, *P_central,
					    beamline->revolution_length, 
					    charge?charge->macroParticleCharge*beam->n_to_track_total:0.0, z);

#endif
		  }
		}
	      }
	      break;
            case T_MHISTOGRAM:
              if (!eptr->pred)
                bombElegant("MHISTOGRAM should not be the first element of the beamline.", NULL);
	      if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
                mhist = (MHISTOGRAM*)eptr->p_elem;
                if (!mhist->disable) {
                  watch_pt_seen = 1;   /* yes, this should be here */
                  if (i_pass==0 && (n_passes/mhist->interval)==0)
                    fprintf(stdout, "warning: n_passes = %ld and MHISTOGRAM interval = %ld--no output will be generated!\n",
                            n_passes, mhist->interval);
                  fflush(stdout);
                  if (i_pass>=mhist->startPass && (i_pass-mhist->startPass)%mhist->interval==0) {

                    ELEMENT_LIST *eptr0;
                    if (mhist->lumped) {
                      if (!mhist->initialized && eptr->occurence==1)
                        set_up_mhist(mhist, run, 0); 
                      eptr0 = &(beamline->elem);
                      while (eptr0) {
                        if (eptr0->type == eptr->type && strcmp(eptr0->name, eptr->name)==0)
                          break;
                        eptr0 = eptr0->succ;
                      }
                    }
                    else {
                      if (!mhist->initialized) 
                        set_up_mhist(mhist, run, eptr->occurence);
                      eptr0 = NULL;
                    }

                    mhist_table(eptr0, eptr, step, i_pass, coord, nToTrack, *P_central,
                                beamline->revolution_length, 
                                charge?charge->macroParticleCharge*nToTrack:0.0, z);
                  }
                }
	      }
	      break;
            case T_FTABLE:
              ftable = (FTABLE*)eptr->p_elem;
              field_table_tracking(coord, nToTrack, ftable, *P_central, run);
              break;       
	    case T_MALIGN:
	      malign = (MALIGN*)eptr->p_elem;
	      if (malign->on_pass==-1 || malign->on_pass==i_pass)
		offset_beam(coord, nToTrack, (MALIGN*)eptr->p_elem, *P_central);
	      break;
	    case T_PEPPOT:
	      nLeft = pepper_pot_plate(coord, (PEPPOT*)eptr->p_elem, nToTrack, accepted);
	      break;
	    case T_ENERGY:
	      energy = (ENERGY*)eptr->p_elem;
	      if (energy->match_beamline) {
		if ((flags&FIDUCIAL_BEAM_SEEN) && eptr->Pref_output_fiducial>0)
		  /* Beamline momentum is defined.  Change particle reference momentum to match. */
		  set_central_momentum(coord, nToTrack, eptr->Pref_output_fiducial, P_central);
		else
		  /* Compute new central momentum to match the average momentum of the particles. */
		  do_match_energy(coord, nToTrack, P_central, 0);
		if (energy->match_particles)
		  bombElegant("can't match_beamline AND match_particles for ENERGY element", NULL);
	      }
	      else if (energy->match_particles) {
		/* change the particle momenta so that the centroid is the central momentum */
		do_match_energy(coord, nToTrack, P_central, 1);
	      }
	      else if (energy->central_energy)
		/* Change particle reference momentum to match the given energy */
		set_central_momentum(coord, nToTrack, sqrt(sqr(energy->central_energy+1)-1), 
				     P_central);
	      else if (energy->central_momentum)
		/* Change particle reference momentum to match the given value */
		set_central_momentum(coord, nToTrack, energy->central_momentum, P_central);
	      break;
	    case T_MAXAMP:
	      maxamp = (MAXAMP*) eptr->p_elem;
	      x_max = maxamp->x_max;
	      y_max = maxamp->y_max;
	      elliptical = maxamp->elliptical;
	      maxampOpenCode = determineOpenSideCode(maxamp->openSide);
	      maxampExponent = maxamp->exponent;
	      maxampYExponent = maxamp->yExponent;
	      break;
	    case T_TRCOUNT:
	      /* >>>>> Needs to be updated */
	      /* 
	       *n_original = nLeft; 
	       if (accepted && i_pass==0)
               copy_particles(accepted, coord, *n_original);
              */
	      break;
	    case T_ALPH:
	      if (!eptr->matrix && !(eptr->matrix=compute_matrix(eptr, run, NULL)))
		bombElegant("no matrix for alpha magnet", NULL);
	      nLeft = alpha_magnet_tracking(coord, eptr->matrix, (ALPH*)eptr->p_elem, nToTrack,
					    accepted, *P_central, z);
	      break;
	    case T_MATR:
	      if (!eptr->matrix)
		eptr->matrix = compute_matrix(eptr, run, NULL);
	      matr_element_tracking(coord, eptr->matrix, (MATR*)eptr->p_elem, nToTrack,
				    z);
	      break;
	    case T_EMATRIX:
	      if (!eptr->matrix)
		eptr->matrix = compute_matrix(eptr, run, NULL);
	      ematrix_element_tracking(coord, eptr->matrix, (EMATRIX*)eptr->p_elem, nToTrack,
				       z, P_central);
	      break;
            case T_ILMATRIX:
              nLeft = trackWithIndividualizedLinearMatrix(coord, nToTrack, accepted, *P_central,
                                                     z, eptr, 
                                                     NULL,
                                                     ((ILMATRIX*)eptr->p_elem)->tune,
                                                     ((ILMATRIX*)eptr->p_elem)->chrom,
                                                     ((ILMATRIX*)eptr->p_elem)->chrom2,
                                                     ((ILMATRIX*)eptr->p_elem)->chrom3,
                                                     ((ILMATRIX*)eptr->p_elem)->beta1,
                                                     ((ILMATRIX*)eptr->p_elem)->alpha1,
                                                     ((ILMATRIX*)eptr->p_elem)->alphac,
                                                     ((ILMATRIX*)eptr->p_elem)->eta1,
                                                     ((ILMATRIX*)eptr->p_elem));
              break;
	    case T_MULT:
	      nLeft = multipole_tracking(coord, nToTrack, (MULT*)eptr->p_elem, 0.0,
					 *P_central, accepted, last_z);
	      break;
	    case T_FMULT:
	      nLeft = fmultipole_tracking(coord, nToTrack, (FMULT*)eptr->p_elem, 0.0,
					  *P_central, accepted, last_z);
	      break;
	    case T_TAYLORSERIES:
	      nLeft = taylorSeries_tracking(coord, nToTrack, (TAYLORSERIES*)eptr->p_elem, 0.0,
					    *P_central, accepted, z);
	      break;
	    case T_KICKER:
	      if (flags&TIME_DEPENDENCE_OFF)
		drift_beam(coord, nToTrack, ((KICKER*)eptr->p_elem)->length, run->default_order);
	      else
		track_through_kicker(coord, nToTrack, (KICKER*)eptr->p_elem, *P_central, i_pass, run->default_order);
	      break;
	    case T_MKICKER:
	      if (flags&TIME_DEPENDENCE_OFF)
		drift_beam(coord, nToTrack, ((MKICKER*)eptr->p_elem)->length, run->default_order);
	      else
		track_through_mkicker(coord, nToTrack, (MKICKER*)eptr->p_elem, *P_central, i_pass, run->default_order);
	      break;
	    case T_KSBEND:
	      nLeft = track_through_kick_sbend(coord, nToTrack, (KSBEND*)eptr->p_elem, 0.0,
					       *P_central, accepted, z);
	      break;
	    case T_CSBEND:
              ((CSBEND*)eptr->p_elem)->edgeFlags = determine_bend_flags(eptr, ((CSBEND*)eptr->p_elem)->edge1_effects,
                                                                        ((CSBEND*)eptr->p_elem)->edge2_effects);
	      if (flags&TEST_PARTICLES) {
		saveISR = ((CSBEND*)eptr->p_elem)->isr;
		((CSBEND*)eptr->p_elem)->isr = 0;
	      }
	      nLeft = track_through_csbend(coord, nToTrack, (CSBEND*)eptr->p_elem, 0.0,
					   *P_central, accepted, last_z, NULL);
	      if (flags&TEST_PARTICLES)
		((CSBEND*)eptr->p_elem)->isr = saveISR;	  
	      break;
	    case T_CSRCSBEND:
              ((CSRCSBEND*)eptr->p_elem)->edgeFlags = determine_bend_flags(eptr, ((CSRCSBEND*)eptr->p_elem)->edge1_effects,
                                                                        ((CSRCSBEND*)eptr->p_elem)->edge2_effects);
	      if (flags&TEST_PARTICLES) {
		saveISR = ((CSRCSBEND*)eptr->p_elem)->isr;
		((CSRCSBEND*)eptr->p_elem)->isr = 0;
	      }
	      nLeft = track_through_csbendCSR(coord, nToTrack, (CSRCSBEND*)eptr->p_elem, 0.0,
					      *P_central, accepted, last_z, z, charge, run->rootname);
	      if (flags&TEST_PARTICLES)
		((CSRCSBEND*)eptr->p_elem)->isr = saveISR;
	      break;
	    case T_CSRDRIFT:
	      nLeft = track_through_driftCSR(coord, nToTrack, (CSRDRIFT*)eptr->p_elem,
					     *P_central, accepted, last_z, 
					     beamline->revolution_length,
					     run->rootname);
	      break;
	    case T_LSCDRIFT:
	      track_through_lscdrift(coord, nToTrack, (LSCDRIFT*)eptr->p_elem, *P_central, charge);
	      break;
	    case T_SCMULT:
	      if (getSCMULTSpecCount()) trackThroughSCMULT(coord, nToTrack, eptr);
	      break;
	    case T_EDRIFT:
	      exactDrift(coord, nToTrack, ((EDRIFT*)eptr->p_elem)->length);
	      break;
	    case T_TUBEND:
	      nLeft = track_through_tubend(coord, nToTrack, 
					   (TUBEND*)eptr->p_elem, 0.0,
					   *P_central, accepted, z);
	      break;
	    case T_KQUAD:
              if (flags&TEST_PARTICLES) {
                saveISR = ((KQUAD*)eptr->p_elem)->isr;
                ((KQUAD*)eptr->p_elem)->isr = 0;
              }
	      nLeft = multipole_tracking2(coord, nToTrack, eptr, 0.0,
                                          *P_central, accepted, last_z,
                                          x_max, y_max, elliptical,
                                          &(run->apertureData), NULL);
              if (flags&TEST_PARTICLES)
                ((KQUAD*)eptr->p_elem)->isr = saveISR;
              break;
	    case T_KSEXT:
              if (flags&TEST_PARTICLES) {
                saveISR = ((KSEXT*)eptr->p_elem)->isr;
                ((KSEXT*)eptr->p_elem)->isr = 0;
              }
	      nLeft = multipole_tracking2(coord, nToTrack, eptr, 0.0,
                                          *P_central, accepted, last_z,
                                          x_max, y_max, elliptical,
                                          &(run->apertureData), NULL);
              if (flags&TEST_PARTICLES)
                ((KSEXT*)eptr->p_elem)->isr = saveISR;
	      break;
	    case T_KOCT:
              if (flags&TEST_PARTICLES) {
                saveISR = ((KOCT*)eptr->p_elem)->isr;
                ((KOCT*)eptr->p_elem)->isr = 0;
              }
	      nLeft = multipole_tracking2(coord, nToTrack, eptr, 0.0,
                                          *P_central, accepted, last_z,
                                          x_max, y_max, elliptical,
                                          &(run->apertureData), NULL);
              if (flags&TEST_PARTICLES)
                ((KOCT*)eptr->p_elem)->isr = saveISR;
	      break;
	    case T_KQUSE:
              if (((KQUSE*)eptr->p_elem)->matrixTracking) {
                if (!eptr->matrix)
                  eptr->matrix = compute_matrix(eptr, run, NULL);
                track_particles(coord, eptr->matrix, coord, nToTrack);
              } else {
                if (flags&TEST_PARTICLES) {
                  saveISR = ((KQUSE*)eptr->p_elem)->isr;
                  ((KQUSE*)eptr->p_elem)->isr = 0;
                }
                nLeft = multipole_tracking2(coord, nToTrack, eptr, 0.0,
                                            *P_central, accepted, last_z,
                                            x_max, y_max, elliptical,
                                            &(run->apertureData), NULL);
                if (flags&TEST_PARTICLES)
                  ((KQUSE*)eptr->p_elem)->isr = saveISR;
              }
	      break;
	    case T_SAMPLE:
	      if (!(flags&TEST_PARTICLES))
		nLeft = sample_particles(coord, (SAMPLE*)eptr->p_elem, nToTrack, accepted, z, *P_central);
	      break;
	    case T_SCATTER:
	      if (!(flags&TEST_PARTICLES))
		scatter(coord, nToTrack, *P_central, (SCATTER*)eptr->p_elem);
	      break;
	    case T_DSCATTER:
	      if (!(flags&TEST_PARTICLES))
		distributionScatter(coord, nToTrack, *P_central, (DSCATTER*)eptr->p_elem, i_pass);
	      break;
	    case T_TSCATTER:
              break;
	    case T_NIBEND:
	      nLeft = lorentz(coord, nToTrack, (NIBEND*)eptr->p_elem, T_NIBEND, *P_central, accepted);
	      break;
	    case T_NISEPT:
	      nLeft = lorentz(coord, nToTrack, (NISEPT*)eptr->p_elem, T_NISEPT, *P_central, accepted);
	      break;
	    case T_BMAPXY:
	      nLeft = lorentz(coord, nToTrack, (BMAPXY*)eptr->p_elem, T_BMAPXY, *P_central, accepted);
	      break;
	    case T_KPOLY:
	      nLeft = polynomial_kicks(coord, nToTrack, (KPOLY*)eptr->p_elem, 0.0,
				       *P_central, accepted, z);
	      break;
	    case T_RAMPRF:
	      ramped_rf_cavity(coord, nToTrack, (RAMPRF*)eptr->p_elem, *P_central, beamline->revolution_length,
			       eptr->end_pos, i_pass);
	      break;
	    case T_RAMPP:
	      ramp_momentum(coord, nToTrack, (RAMPP*)eptr->p_elem, P_central, i_pass);
	      break;
	    case T_SOLE:
	      if (((SOLE*)eptr->p_elem)->B) {
		SOLE *sptr;
		double ks;
		sptr = (SOLE*)eptr->p_elem;
		if ((ks = -sptr->B/(*P_central*particleMass*c_mks/particleCharge))!=sptr->ks) {
		  sptr->ks = ks;
		  if (eptr->matrix) {
		    free_matrices(eptr->matrix);
                    free(eptr->matrix);
                    eptr->matrix = NULL;
                  }
		  if (!(eptr->matrix = compute_matrix(eptr, run, NULL)))
		    bombElegant("no matrix for element that must have matrix", NULL);
		}
		sptr->ks = 0;  /* reset so it is clear that B is fundamental quantity */
	      }
	      if (!eptr->matrix) {
		if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
		  bombElegant("no matrix for element that must have matrix", NULL);
	      }
	      track_particles(coord, eptr->matrix, coord, nToTrack);
	      break;
	    case T_MATTER:
	      track_through_matter(coord, nToTrack, (MATTER*)eptr->p_elem, *P_central);
	      break;
	    case T_RFMODE:
	      rfmode = (RFMODE*)eptr->p_elem;
	      if (!rfmode->initialized)
		set_up_rfmode(rfmode, eptr->name, eptr->end_pos, n_passes, run, 
			      nOriginal, *P_central,
			      beamline->revolution_length);
	      track_through_rfmode(coord, nToTrack, (RFMODE*)eptr->p_elem, *P_central,
				   eptr->name, eptr->end_pos, i_pass, n_passes,
				   charge);
	      break;
	    case T_FRFMODE:
	      frfmode = (FRFMODE*)eptr->p_elem;
	      if (!frfmode->initialized)
		set_up_frfmode(frfmode, eptr->name, eptr->end_pos, n_passes, run, 
			       nOriginal, *P_central,
			       beamline->revolution_length);
	      track_through_frfmode(coord, nToTrack, frfmode, *P_central,
				    eptr->name, eptr->end_pos, i_pass, n_passes,
				    charge);
	      break;
	    case T_TRFMODE:
	      trfmode = (TRFMODE*)eptr->p_elem;
	      if (!trfmode->initialized)
		set_up_trfmode(trfmode, eptr->name, eptr->end_pos, n_passes, run, nOriginal);
	      track_through_trfmode(coord, nToTrack, (TRFMODE*)eptr->p_elem, *P_central,
				    eptr->name, eptr->end_pos, i_pass, n_passes,
				    charge);
	      break;
	    case T_FTRFMODE:
	      ftrfmode = (FTRFMODE*)eptr->p_elem;
	      if (!ftrfmode->initialized)
		set_up_ftrfmode(ftrfmode, eptr->name, eptr->end_pos, n_passes, run, 
				nOriginal, *P_central,
				beamline->revolution_length);
	      track_through_ftrfmode(coord, nToTrack, ftrfmode, *P_central,
				     eptr->name, eptr->end_pos, i_pass, n_passes,
				     charge);
	      break;
	    case T_ZLONGIT:
	      track_through_zlongit(coord, nToTrack, (ZLONGIT*)eptr->p_elem, *P_central, run, i_pass,
				    charge);
	      break;
	    case T_ZTRANSVERSE:
	      track_through_ztransverse(coord, nToTrack, (ZTRANSVERSE*)eptr->p_elem, *P_central, run, i_pass,
					charge);
	      break;
	    case T_CORGPIPE:
              nLeft = elimit_amplitudes(coord, ((CORGPIPE*)eptr->p_elem)->radius, ((CORGPIPE*)eptr->p_elem)->radius, 
                                        nToTrack, accepted, z-((CORGPIPE*)eptr->p_elem)->length, *P_central, 0, 0, 2, 2);
	      track_through_corgpipe(coord, nLeft, (CORGPIPE*)eptr->p_elem, P_central, run, i_pass,
				 charge);
              nLeft = elimit_amplitudes(coord, ((CORGPIPE*)eptr->p_elem)->radius, ((CORGPIPE*)eptr->p_elem)->radius, 
                                        nLeft, accepted, z, *P_central, 1, 0, 2, 2);
	      break;
	    case T_WAKE:
	      track_through_wake(coord, nToTrack, (WAKE*)eptr->p_elem, P_central, run, i_pass,
				 charge);
	      break;
	    case T_TRWAKE:
	      track_through_trwake(coord, nToTrack, (TRWAKE*)eptr->p_elem, *P_central, run, i_pass, 
				   charge);
	      break;
	    case T_SREFFECTS:
              track_SReffects(coord, nToTrack, (SREFFECTS*)eptr->p_elem, *P_central, eptr->twiss, &(beamline->radIntegrals),
                              flags&TEST_PARTICLES);
	      break;
	    case T_IBSCATTER:
	      if (!(flags&TEST_PARTICLES))
		track_IBS(coord, nToTrack, (IBSCATTER*)eptr->p_elem,
			  *P_central, &(beamline->elem), charge, i_pass, n_passes, run);
	      break;
	    case T_SCRIPT:
#if !USE_MPI
	if (nLeft<nMaximum)
	      if (((SCRIPT*)eptr->p_elem)->verbosity>1)
		fprintf(stdout, "nLost=%ld, beam->n_particle=%ld, beam->n_to_track=%ld, nLeft=%ld, nToTrack=%ld, nMaximum=%ld\n",
			nLost, beam->n_particle, beam->n_to_track, nLeft, nToTrack, nMaximum);
#endif
	      nLeft = transformBeamWithScript((SCRIPT*)eptr->p_elem, *P_central, charge, 
					      beam, coord, nToTrack, &nLost, 
					      run->rootname, i_pass, run->default_order);
#if USE_MPI
	      nToTrack = nLeft;
	      /* As the particles could be redistributed across processors, we need adjust the beam->n_to_track to dump lost particle coordinate at the end */ 
	      beam->n_to_track = nLeft+nLost;
#endif
	      if (((SCRIPT*)eptr->p_elem)->verbosity>2)
		fprintf(stdout, "nLost=%ld, beam->n_particle=%ld, beam->n_to_track=%ld, nLeft=%ld, nToTrack=%ld, nMaximum=%ld\n",
			nLost, beam->n_particle, beam->n_to_track, nLeft, nToTrack, nMaximum);
	      if (beam && coord!=beam->particle) {
		/* particles were created and so the particle array was changed */
		coord = beam->particle;
#if !USE_MPI
		if (nLost != (beam->n_to_track - nLeft)) {
		  fprintf(stderr, "Particle accounting problem after return from script.\n");
		  fprintf(stderr, "nLost=%ld, beam->n_particle=%ld, nLeft=%ld\n",
			  nLost, beam->n_particle, nLeft);
		  exitElegant(1);         	                        
		}
#endif
	      }
	      if (beam && lostParticles!=beam->lost)
		lostParticles = beam->lost;

	      if (nMaximum<beam->n_to_track)
		nMaximum = beam->n_to_track;
#if !USE_MPI
	      if (((SCRIPT*)eptr->p_elem)->verbosity>1)
		fprintf(stdout, "nLost=%ld, beam->n_particle=%ld, beam->n_to_track=%ld, nLeft=%ld, nMaximum=%ld\n\n",
			nLost, beam->n_particle, beam->n_to_track, nLeft, nMaximum);
#endif
	      break;
	    case T_FLOORELEMENT:
	      break;
	    case T_TFBPICKUP:
	      if (!(flags&TEST_PARTICLES))
		transverseFeedbackPickup((TFBPICKUP*)eptr->p_elem, coord, nToTrack, i_pass);
	      break;
	    case T_STRAY:
	      if (eptr->matrix) {
		free_matrices(eptr->matrix);
                free(eptr->matrix);
                eptr->matrix = NULL;
              }
	      stray = (STRAY*)eptr->p_elem;
	      eptr->matrix = stray_field_matrix(stray->length, &stray->lBx, &stray->gBx, 
						eptr->end_theta, stray->order?stray->order:run->default_order,
						*P_central, 
						stray->Wi);
	      track_particles(coord, eptr->matrix, coord, nToTrack);
	      break;
	    case T_TFBDRIVER:
	      if (!(flags&TEST_PARTICLES))
		transverseFeedbackDriver((TFBDRIVER*)eptr->p_elem, coord, nToTrack, beamline, i_pass, n_passes, run->rootname);
	      feedbackDriverSeen = 1;
	      break;
	    case T_LSRMDLTR:
	      nLeft = motion(coord, nToTrack, eptr->p_elem, eptr->type, P_central, 
			     &dgamma, dP, accepted, last_z);
	      show_dE = 1;
	      break;
	    case T_CWIGGLER:
              if (flags&TEST_PARTICLES) {
                saveISR = ((CWIGGLER*)eptr->p_elem)->isr;
                ((CWIGGLER*)eptr->p_elem)->isr = 0;
              }
	      GWigSymplecticPass(coord, nToTrack, *P_central, (CWIGGLER*)eptr->p_elem);
              if (flags&TEST_PARTICLES)
                ((CWIGGLER*)eptr->p_elem)->isr = saveISR;
	      break;
	    case T_APPLE:
              if (flags&TEST_PARTICLES) {
                saveISR = ((APPLE*)eptr->p_elem)->isr;
                ((APPLE*)eptr->p_elem)->isr = 0;
              }
	      APPLE_Track(coord, nToTrack, *P_central, (APPLE*)eptr->p_elem);
              if (flags&TEST_PARTICLES)
                ((APPLE*)eptr->p_elem)->isr = saveISR;
	      break;
            case T_UKICKMAP:
              nLeft = trackUndulatorKickMap(coord, accepted, nToTrack, *P_central, (UKICKMAP*)eptr->p_elem, 
                                            last_z);
              break;
            case T_TWISSELEMENT:
	      if (((TWISSELEMENT*)eptr->p_elem)->disable)
		break;
              if ( ((TWISSELEMENT*)eptr->p_elem)->applyOnce==0 || i_pass==passOffset) {
                /* If applying once, do so on the first pass through only */
                if ( ((TWISSELEMENT*)eptr->p_elem)->fromBeam ) {
                  /* Compute the transformation from the beam, rather than the lattice twiss parameters */
                  if ( ((TWISSELEMENT*)eptr->p_elem)->computeOnce==0 || 
                      ((TWISSELEMENT*)eptr->p_elem)->transformComputed==0) {
                    TWISS beamTwiss;
                    if (((TWISSELEMENT*)eptr->p_elem)->verbose)
                      printf("* Computing beam-based twiss transformation matrix for %s at z=%e m\n",
                             eptr->name, eptr->end_pos);
#if SDDS_MPI_IO
                    if (!partOnMaster && notSinglePart) {
                      if (isMaster) nToTrack = 0;
                      MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
                    } else { /* singlePart tracking or partOnMaster */
                      beam->n_to_track_total = nToTrack;
                    }
		    if (isMaster && (beam->n_to_track_total<10)) {
#else
                    if (nToTrack<10) {
#endif		     
                      printf("*** Error: too few particles (<10) for computation of twiss parameters from beam\n");
                      exitElegant(1);
                    }
                    computeBeamTwissParameters(&beamTwiss, coord, nToTrack);
                    if (eptr->matrix) {
                      free_matrices(eptr->matrix);
                      free(eptr->matrix);
                      eptr->matrix = NULL;
                    }
                    eptr->matrix = twissTransformMatrix((TWISSELEMENT*)eptr->p_elem, &beamTwiss);
                    ((TWISSELEMENT*)eptr->p_elem)->transformComputed = 1;
                  }
		}
                if (((TWISSELEMENT*)eptr->p_elem)->transformComputed==0) {
                  if (((TWISSELEMENT*)eptr->p_elem)->from0Values) {
		    if (eptr->matrix) {
		      free_matrices(eptr->matrix);
                      free(eptr->matrix);
                      eptr->matrix = NULL;
                    }
		    eptr->matrix = twissTransformMatrix1(&(((TWISSELEMENT*)eptr->p_elem)->twiss), &(((TWISSELEMENT*)eptr->p_elem)->twiss0));
		  } else {
		    printf("Error: The twiss parameter transformation matrix was not computed for element %s at z=%e m\n",
			   eptr->name, eptr->end_pos);
		    printf("This means you set FROM_BEAM=0 but didn't issue a twiss_output command.\n");
		    exitElegant(1);
		  }
                }
                if (eptr->matrix==NULL) {
                  printf("Error: twiss parameter transformation matrix was not computed for element %s at z=%e m\n",
                         eptr->name, eptr->end_pos);
                  printf("and this wasn't properly detected.  Please send your input files to borland@aps.anl.gov.\n");
                  exitElegant(1);
                }
                if (((TWISSELEMENT*)eptr->p_elem)->verbose) {
                  TWISS beamTwiss;
                  printf("* Applying twiss parameter transformation matrix (%s at z=%e m) to beam.\n", eptr->name, eptr->end_pos);
                  computeBeamTwissParameters(&beamTwiss, coord, nToTrack);
                  printf("  * Initial twiss parameters:\n");
                  printf("  betax = %le  alphax = %le  etax = %le, etaxp = %le\n",
                         beamTwiss.betax, beamTwiss.alphax, beamTwiss.etax, beamTwiss.etapx);
                  printf("  betay = %le  alphay = %le  etay = %le, etayp = %le\n",
                         beamTwiss.betay, beamTwiss.alphay, beamTwiss.etay, beamTwiss.etapy);
                  fflush(stdout);
                }
                track_particles(coord, eptr->matrix, coord, nToTrack);
                if (((TWISSELEMENT*)eptr->p_elem)->verbose) {
                  TWISS beamTwiss;
                  computeBeamTwissParameters(&beamTwiss, coord, nToTrack);
                  printf("  * Final twiss parameters:\n");
                  printf("  betax = %le  alphax = %le  etax = %le, etaxp = %le\n",
                         beamTwiss.betax, beamTwiss.alphax, beamTwiss.etax, beamTwiss.etapx);
                  printf("  betay = %le  alphay = %le  etay = %le, etayp = %le\n",
                         beamTwiss.betay, beamTwiss.alphay, beamTwiss.etay, beamTwiss.etapy);
                  fflush(stdout);
                }
              }
              break;
            case T_EMITTANCE:
              transformEmittances(coord, nToTrack, *P_central, (EMITTANCEELEMENT*)eptr->p_elem);
              break;
            case T_MRADINTEGRALS:
              break;
            case T_HCOR:
            case T_VCOR:
            case T_HVCOR:
              if (!(entity_description[eptr->type].flags&HAS_MATRIX))
		bombElegant("attempt to matrix-multiply for element with no matrix!",  NULL);
	      if (!eptr->matrix) {
		if (!(eptr->matrix=compute_matrix(eptr, run, NULL)))
		  bombElegant("no matrix for element that must have matrix", NULL);
	      }
	      /* Only the slave CPUs will track */ 
	      if ((!USE_MPI || !notSinglePart) || (USE_MPI && (myid!=0))) 
		track_particles(coord, eptr->matrix, coord, nToTrack);
	      switch (type) {
	      case T_HCOR:
		if (((HCOR*)(eptr->p_elem))->synchRad)
		  addCorrectorRadiationKick(coord, nToTrack, eptr, type, *P_central, NULL, 
					    flags&(CLOSED_ORBIT_TRACKING+TEST_PARTICLES));
		break;
	      case T_VCOR:
		if (((VCOR*)(eptr->p_elem))->synchRad)
		  addCorrectorRadiationKick(coord, nToTrack, eptr, type, *P_central, NULL, 
					    flags&(CLOSED_ORBIT_TRACKING+TEST_PARTICLES));
		break;
	      case T_HVCOR:
		if (((HVCOR*)(eptr->p_elem))->synchRad)
		  addCorrectorRadiationKick(coord, nToTrack, eptr, type, *P_central, NULL,
					    flags&(CLOSED_ORBIT_TRACKING+TEST_PARTICLES));
		break;
	      }
	      break;
	    default:
	      fprintf(stdout, "programming error: no tracking statements for element %s (type %s)\n",
		      eptr->name, entity_name[eptr->type]);
	      fflush(stdout);
	      exitElegant(1);
	      break;
	    }
	  }
#ifdef USE_MPE
	      MPE_Log_event( event2b, 0, bytebuf );
#endif
	}
#if USE_MPI
	if ((myid==0) && notSinglePart && (!usefulOperation(eptr, flags, i_pass)))
	  active = 0;
#endif
	if ((!USE_MPI || !notSinglePart ) || (USE_MPI && active)) {
	  if (!(flags&TEST_PARTICLES && !(flags&TEST_PARTICLE_LOSSES))) {
            if (x_max || y_max) {
              if (!elliptical) 
                nLeft = limit_amplitudes(coord, x_max, y_max, nLeft, accepted, z, *P_central, 
                                         eptr->type==T_DRIF || eptr->type==T_STRAY,
                                         maxampOpenCode);
              else
                nLeft = elimit_amplitudes(coord, x_max, y_max, nLeft, accepted, z, *P_central, 
                                          eptr->type==T_DRIF || eptr->type==T_STRAY,
                                          maxampOpenCode, maxampExponent, maxampYExponent);
            }
            if (run->apertureData.initialized) 
              nLeft = imposeApertureData(coord, nLeft, accepted, z, *P_central, 
                                 &(run->apertureData));
          }
	}
        if (run->print_statistics && !(flags&TEST_PARTICLES)) {
          report_stats(stdout, ": ");
          fprintf(stdout, "central momentum is %e    zstart = %em  zend = %em\n", *P_central, last_z, z);
          fflush(stdout);
          if (nLeft!=nToTrack)
            fprintf(stdout, "%ld particles left\n", nLeft);
	  fflush(stdout);
        }
      }
      else if (!(flags&TEST_PARTICLES)) {
        fprintf(stdout, "element %s was ignored in tracking.\n",
                eptr->name);
        fflush(stdout);
      }
      if (flags&FIRST_BEAM_IS_FIDUCIAL && !(flags&FIDUCIAL_BEAM_SEEN)) {
        if (!(flags&RESTRICT_FIDUCIALIZATION) ||
            (entity_description[eptr->type].flags&MAY_CHANGE_ENERGY)) {
	  if (!(classFlags&(UNIDIAGNOSTIC&(~UNIPROCESSOR))))
	    /* If it is a Diagnostic element, nothing needs to be done */	    
	    do_match_energy(coord, nLeft, P_central, 0);
        }
        eptr->Pref_output_fiducial = *P_central;
      } else if (flags&FIDUCIAL_BEAM_SEEN) {
        if (*P_central!=eptr->Pref_output_fiducial)
          set_central_momentum(coord, nLeft, eptr->Pref_output_fiducial, P_central);
      }
      else if (run->always_change_p0)
	if (!(classFlags&(UNIDIAGNOSTIC&(~UNIPROCESSOR))))
	  /* If it is a Diagnostic element, nothing needs to be done */
	  do_match_energy(coord, nLeft, P_central, 0);
      if (i_pass==0 && traj_vs_z) {
        /* collect trajectory data--used mostly by trajectory correction routines */
        if (!traj_vs_z[i_traj].centroid) {
          fprintf(stdout, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                  eptr->name);
          fflush(stdout);
          exitElegant(1);
        }
        traj_vs_z[i_traj].elem = eptr;
        if (!(traj_vs_z[i_traj].n_part=nLeft)) {
          for (i=0; i<6; i++)
            traj_vs_z[i_traj].centroid[i] = 0;
        }
        else {
          for (i=0; i<6; i++) {
            for (j=sum=0; j<nToTrack; j++)
              sum += coord[j][i];
            traj_vs_z[i_traj].centroid[i] = sum/nLeft;
          }
        }
        i_traj++;
      }
      if (!(flags&TEST_PARTICLES) && sliceAnalysis && sliceAnalysis->active && !sliceAnalysis->finalValuesOnly) {
#if USE_MPI
	if (!(classFlags&UNIPROCESSOR)) { /* This function will be parallelized in the future */
	  fprintf(stdout, "performSliceAnalysisOutput is not supported in parallel mode currently.\n");
	  MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
	  MPI_Abort(MPI_COMM_WORLD, 1); 
	}
#endif
	performSliceAnalysisOutput(sliceAnalysis, coord, nToTrack, 
				   !sliceAnDone, step, 
				   *P_central, 
				   charge?charge->macroParticleCharge*nToTrack:0.0, 
				   eptr->name, eptr->end_pos, 0); 
	sliceAnDone = 1;
      }
#if USE_MPI
      if (notSinglePart) {
	if (!(classFlags&(UNIPROCESSOR|MPALGORITHM))) {
	  end_wtime = MPI_Wtime();
	  my_wtime = my_wtime+end_wtime-start_wtime; 
	}
	else if (!(classFlags&((UNIDIAGNOSTIC&(~UNIPROCESSOR))|MPALGORITHM))) { 
	  /* a non-diagnostic uniprocessor element */
	  if ((myid == 0) && (nMaximum!=(nLeft+nLost)))         
	    /* there are additional losses occurred */
	    lostSinceSeqMode = needSort= 1;
	  MPI_Bcast (&lostSinceSeqMode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	if (classFlags&MPALGORITHM && isMaster) {
	  /* Master does not need to do limit_amplitudes for MPALGORITHM elements */
	  active = 0; 
	}
      }
#endif
      if ((!USE_MPI || !notSinglePart) || (USE_MPI && active)) {
        nLeft = limit_amplitudes(coord, DBL_MAX, DBL_MAX, nLeft, accepted, z, *P_central, 0, 0);

	if (eptr->type!=T_SCRIPT) { /* For the SCRIPT element, the lost particle coordinate will be recorded inside the element */
	  if (nLeft!=nToTrack)
	    recordLostParticles(lostParticles, coord, &nLost, nLeft, nToTrack, i_pass);
	}
      }  

      if (getSCMULTSpecCount() && entity_description[eptr->type].flags&HAS_LENGTH) {
	/* calcaulate beam size at exit of element for use in space  charge calculation with SCMULT */
	/* need special care for element with 0 length but phase space rotation */
      	if (((DRIFT*)eptr->p_elem)->length > 0.0) 
        	accumulateSCMULT(coord, nToTrack, eptr);
      }
      last_type = eptr->type;
      eptrPred = eptr;
      eptr = eptr->succ;
      nToTrack = nLeft;
    } /* end of the while loop */
    if (!(flags&TEST_PARTICLES) && sliceAnalysis && sliceAnalysis->active && !sliceAnalysis->finalValuesOnly) {
#if USE_MPI
      if (notSinglePart) {
	if (!(classFlags&UNIPROCESSOR)) { /* This function will be parallelized in the future */
	  fprintf(stdout, "performSliceAnalysisOutput is not supported in parallel mode currently.\n");
	  MPI_Abort(MPI_COMM_WORLD, 1); 
	}
      }
#endif
      performSliceAnalysisOutput(sliceAnalysis, coord, nToTrack, 
				 !sliceAnDone, step, 
				 *P_central, 
				 charge?charge->macroParticleCharge*nToTrack:0.0, 
				 eptrPred->name, eptrPred->end_pos, 0);
      
      sliceAnDone = 1;
    }
   
    if (effort) {
#if !SDDS_MPI_IO 
      *effort += nLeft;
#else
      if ((isMaster&&(partOnMaster||!notSinglePart)) || (isSlave&&!partOnMaster))
	*effort += nLeft;
#endif
    }


    log_exit("do_tracking.2.2");
#ifdef WATCH_MEMORY
    fprintf(stdout, "main tracking loop done: CPU: %6.2lf  PF: %6ld  MEM: %6ld\n",
            cpu_time()/100.0, page_faults(), memory_count());
    fflush(stdout);
#endif
    if ((!USE_MPI || !notSinglePart) && (i_pass==0 || watch_pt_seen || feedbackDriverSeen)) {
      /* if eptr is not NULL, then all particles have been lost */
      /* some work still has to be done, however. */
      while (eptr) {
        if (sums_vs_z && *sums_vs_z && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES)) {
          if (i_sums<0)
            bombElegant("attempt to accumulate beam sums with negative index!", NULL);
          accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
          (*sums_vs_z)[i_sums].z = z;
          i_sums++;
        }
        if (entity_description[eptr->type].flags&HAS_LENGTH && eptr->p_elem)
          z += ((DRIFT*)eptr->p_elem)->length;
        else {
          if (eptr->pred)
            z += eptr->end_pos - eptr->pred->end_pos;
          else
            z += eptr->end_pos;
        }
        switch (eptr->type) {
        case T_TFBDRIVER:
          flushTransverseFeedbackDriverFiles((TFBDRIVER *)(eptr->p_elem));
          break;
        case T_WATCH:
#if USE_MPI
	      if (!notSinglePart) /* When each processor tracks the beam independently, the watch point will be disabled in Pelegant */
		break;
#endif
          if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
            watch = (WATCH*)eptr->p_elem;
            if (!watch->initialized) 
              set_up_watch_point(watch, run, eptr->occurence, eptr->pred?eptr->pred->name:NULL);
	    if (!watch->disable) {
	      if (i_pass%watch->interval==0) {
		switch (watch->mode_code) {
		case WATCH_COORDINATES:
		  break;
		case WATCH_PARAMETERS:
		case WATCH_CENTROIDS:
#if SDDS_MPI_IO
		  dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, total_nOriginal, *P_central,
					beamline->revolution_length, z);
#else
		  dump_watch_parameters(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central,
					beamline->revolution_length, z);
#endif
		  break;
		case WATCH_FFT:
		  dump_watch_FFT(watch, step, i_pass, n_passes, coord, nToTrack, nOriginal, *P_central);
		  break;
		}
	      }
	    }
	  }
          break;
        default:
          break;
        }
        if (i_pass==0 && traj_vs_z) {
          /* collect trajectory data--used mostly by trajectory correction routines */
          if (!traj_vs_z[i_traj].centroid) {
            fprintf(stdout, "error: the trajectory centroid array for %s is NULL (do_tracking)",
                    eptr->name);
            fflush(stdout);
            exitElegant(1);
          }
          traj_vs_z[i_traj].elem = eptr;
          traj_vs_z[i_traj].n_part = 0;
          for (i=0; i<6; i++)
            traj_vs_z[i_traj].centroid[i] = 0;
          i_traj++;
        }
        eptr = eptr->succ;
      }
    }
 
    if (sums_vs_z && (*sums_vs_z) && !(flags&FINAL_SUMS_ONLY) && !(flags&TEST_PARTICLES) &&
        (run->wrap_around || i_pass==n_passes-1)) {
      if (i_sums<0)
        bombElegant("attempt to accumulate beam sums with negative index!", NULL);
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for %s at z=%em, sx=%e\n", 
              i_sums, name, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
      fflush(stdout);
#endif
      i_sums++;
    }
 
#if USE_MPI
    if (notSinglePart) {
      if (run->load_balancing_on) {  /* User can choose if load balancing needs to be done */
	if (balanceStatus==startMode) { 
	  balanceStatus = checkBalance (my_wtime, myid, n_processors);  
	  /* calculate the rate for all of the slave processors */
	  if (myid==0) {
	    my_rate = 0.0;
	    nParPerElements = 0.0;
	  }
	  else {
	    nParPerElements = (double)nParElements/(double)nElements;
	    if (my_wtime!=0.0) 
	      my_rate = nParPerElements/my_wtime;
	    else 
	      my_rate = 1.; 
	  } 
	  /*  lostSinceSeqMode = 1; */ /* set flag to distribute jobs according to  the speed.
                                          The default redistribution for the first turn is disabled
                                          as it might cause some problems for random number generator */
	}
	else { /* The workload balancing will be checked for every pass by default.
		  If user defined the CHECKFLAGS, the balance will be checked only 
		  when the nToTrack is changed. */
#ifdef CHECKFLAGS 
	  if (myid==0) {
	    if (old_nToTrack!=nToTrack) {
	      checkFlags = 1;
	    }
	    else
	      checkFlags = 0;
	  }
	  MPI_Bcast(&checkFlags, 1, MPI_INT, 0, MPI_COMM_WORLD);
#else
	  checkFlags = 1; /* the default option */
#endif
	  if (checkFlags)
	    balanceStatus = checkBalance (my_wtime, myid, n_processors);   
	  if (balanceStatus == badBalance) {
	    if (myid==0) {
	      my_rate = 0.0;
	      nParPerElements = 0.0;
	    }
	    else {
	      nParPerElements = (double)nParElements/(double)nElements;
	      if (my_wtime!=0.0)
		my_rate = nParPerElements/my_wtime;
	      else  
		/* set the speed to be euqal for the special case where all the elements are UNIPROCESSOR or MPALGORITHM */
		my_rate = 1.;
	    } 
	  }
#ifdef MPI_DEBUG  
	  fprintf(stdout, "\n\nmyid=%d, nParPerElements=%e, my_time=%lf, my_rate=%lf\n",
		  myid, nParPerElements, my_wtime, nParPerElements/my_wtime);
	  fprintf(stdout, "nParElements=%ld, nElements=%ld\n",nParElements, nElements);
#endif
	}
      }
      else 
        balanceStatus = goodBalance;
      if (myid==0)
	old_nToTrack = nToTrack;
    }
#endif
  } /* end of the for loop for n_passes*/


#ifdef SORT   /* Sort the particles when the particles are lost at the very last element */
      if (!USE_MPI || needSort)
	if (nToTrackAtLastSort > nToTrack)  {/* indicates more particles are lost, need sort */
          if (beam && beam->bunchFrequency!=0)
            fprintf(stdout, "*** Warning: particle ID sort not being performed because bunch frequency is nonzero\n");
          else { 
            qsort(coord[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
            if (accepted!=NULL)
              qsort(accepted[0], nToTrack, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs);
            nToTrackAtLastSort = nToTrack;
          }   
        }
#endif 

#if USE_MPI 
  #if  !SDDS_MPI_IO  
   if (notSinglePart)
      /* change back to sequential mode before leaving the do_tracking function */
      if (parallelStatus==trueParallel && notSinglePart) {
	gatherParticles(&coord, NULL, &nToTrack, &nLost, &accepted, n_processors, myid, &round);
      	MPI_Bcast(&nToTrack, 1, MPI_LONG, 0, MPI_COMM_WORLD); 
	parallelStatus = notParallel ;
	partOnMaster = 1;
      }
  #else
      /* Make sure that the particles are distributed to the slave processors for parallel IO */ 
      if (partOnMaster && notSinglePart) {
	scatterParticles(coord, &nToTrack, accepted, n_processors, myid,
			       balanceStatus, my_rate, nParPerElements, round, lostSinceSeqMode, &distributed, &reAllocate, P_central);
        parallelStatus = trueParallel;
      }
  #endif

#endif

  /* do this here to get report of CSR drift normalization */
  reset_driftCSR();

  log_exit("do_tracking.2");
  log_entry("do_tracking.3");

  if (((!USE_MPI && nLeft) || USE_MPI) && sums_vs_z && *sums_vs_z && !(flags&TEST_PARTICLES)) {
    if (flags&FINAL_SUMS_ONLY) {
      log_entry("do_tracking.3.1");
      i_sums = 0;
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
      (*sums_vs_z)[i_sums].z = z;
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
      fflush(stdout);
#endif
      log_exit("do_tracking.3.1");
    }
    else if (run->wrap_around) {
      log_entry("do_tracking.3.2");
      if (i_sums<0)
        bombElegant("attempt to accumulate beam sums with negative index!", NULL);
      /* accumulate sums for final output */
      accumulate_beam_sums(*sums_vs_z+i_sums, coord, nToTrack, *P_central);
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums accumulated in slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, z, sqrt((*sums_vs_z)[i_sums].sum2[0]/nLeft));
      fflush(stdout);
#endif
      log_exit("do_tracking.3.2");
    }
    else {
      log_entry("do_tracking.3.3");
      if (i_sums<0)
        bombElegant("attempt to accumulate beam sums with negative index!", NULL);
      copy_beam_sums(*sums_vs_z+i_sums, *sums_vs_z+i_sums-1);
#if defined(BEAM_SUMS_DEBUG)
      fprintf(stdout, "beam sums copied to slot %ld from slot %ld for final sums at z=%em, sx=%e\n", 
              i_sums, i_sums-1, z, (*sums_vs_z)[i_sums].sum2[0]);
      fflush(stdout);
#endif
      log_exit("do_tracking.3.3");
    }
  }

  if (sasefel && sasefel->active) {
    if (!charge) {
      fprintf(stdout, "Can't compute SASE FEL---no CHARGE element seen");
      fflush(stdout);
      exitElegant(1);
    }
#if SDDS_MPI_IO
  if (!partOnMaster && notSinglePart) {
    if (isMaster) nToTrack = 0;
    MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  } else { /* singlePart tracking or partOnMaster */
    beam->n_to_track_total = nToTrack;
  }
    /* The charge is correct on master only, but it should not affect the output */
    computeSASEFELAtEnd(sasefel, coord, nToTrack, *P_central, charge->macroParticleCharge*beam->n_to_track_total);
#else
    computeSASEFELAtEnd(sasefel, coord, nToTrack, *P_central, charge->macroParticleCharge*nToTrack);
#endif
  }
  
  log_exit("do_tracking.3");
  log_entry("do_tracking.4");
  if (!(flags&SILENT_RUNNING) && !is_batch && n_passes!=1 && !(flags&TEST_PARTICLES)) {
#if !SDDS_MPI_IO
    fprintf(stdout, "%ld particles present after pass %ld        \n", 
            nToTrack, i_pass);
#else
    if (!partOnMaster && notSinglePart) {
      /* We have to collect information from all the processors to print correct info during tracking */
      if (isMaster) nToTrack = 0; 
      MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      fprintf(stdout, "%ld particles present after pass %ld        \n", 
	      beam->n_to_track_total, i_pass);
    }
    else {
      beam->n_to_track_total = nToTrack;
      fprintf(stdout, "%ld particles present after pass %ld        \n", 
	      nToTrack, i_pass);
    }
#endif
    fflush(stdout);
  }

#ifdef MPI_DEBUG  
  #ifdef CHECKFLAGS 
    printf("Balance is checked for the first pass and when particles are lost only.\n"); 
    fflush(stdout);
  #else
    if (run->load_balancing_on) {
      printf("Balance is checked for every pass.\n"); 
      fflush(stdout);
    }
  #endif
#endif

  log_exit("do_tracking.4");

  log_exit("do_tracking");
 
  if (charge && finalCharge) {
#if !SDDS_MPI_IO
    *finalCharge = nToTrack*charge->macroParticleCharge;
#else
    if (!partOnMaster) {
      /* We have to collect information from all the processors to print correct info after tracking */
      if (isMaster) nToTrack = 0; 
      MPI_Reduce (&nToTrack, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
    else {
      beam->n_to_track_total = nToTrack;
    }
    /* Only Master will have the correct information */
    *finalCharge = beam->n_to_track_total*charge->macroParticleCharge;
#endif
  }

  return(nToTrack);
}

void offset_beam(
                 double **coord,
                 long nToTrack, 
                 MALIGN *offset,
                 double P_central
                 )
{
  long i_part;
  double *part, pc, beta, gamma, t;
  double ds;
  
  log_entry("offset_beam");
  
  for (i_part=nToTrack-1; i_part>=0; i_part--) {
    part = coord[i_part];
    if (offset->dz)
      ds = offset->dz*sqrt(1+sqr(part[1])+sqr(part[3]));
    else
      ds = 0;
    part[0] += offset->dx + offset->dz*part[1];
    part[1] += offset->dxp;
    part[2] += offset->dy + offset->dz*part[3];
    part[3] += offset->dyp;
    part[4] += ds;
    if (offset->dt || offset->dp || offset->de) {
      pc = P_central*(1+part[5]);
      beta = pc/(gamma=sqrt(1+pc*pc));
      t = part[4]/(beta*c_mks) + offset->dt;
      if (offset->dp) {
        part[5] += offset->dp;
        pc = P_central*(1+part[5]);
        beta = pc/sqrt(1+pc*pc);
      }
      if (offset->de) {
        gamma += offset->de*gamma;
        pc = sqrt(gamma*gamma-1);
        beta = pc/gamma;
        part[5] = (pc-P_central)/P_central;
      }
      part[4] = t*beta*c_mks;
    }
  }
  log_exit("offset_beam");
}

void do_match_energy(
                     double **coord, 
                     long np,
                     double *P_central,
                     long change_beam
                     )
{
  long ip;
  double P_average, dP_centroid, P, t;
  long active = 1;
#ifdef USE_KAHAN
  double error = 0.0;
#endif
#if USE_MPI
  long np_total;
  double P_total = 0.0;
  if (notSinglePart) {
    if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster))
      active = 1;
    else 
      active = 0;
  }  
#endif

  log_entry("do_match_energy");

#if (!USE_MPI)  
  if (!np) {
    log_exit("do_match_energy");
    return;
  }
#else
  if (notSinglePart) {
    if (parallelStatus!=trueParallel) {
      if (!np) {
	log_exit("do_match_energy");   
	return;   
      }
    }
    else {
      if (isMaster) 
	np = 0; /* All the particles have been distributed to the slave processors */    
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);        
    }    
  }
  else if (!np) {
    log_exit("do_match_energy");
    return;
  }
#endif

  if (!change_beam) {
    /* change the central momentum so that it matches the beam's centroid */
    P_average = 0;
    if (active) {
      for (ip=0; ip<np; ip++) {
#ifndef USE_KAHAN	     
	P_average += (*P_central)*(1+coord[ip][5]);
#else
	P_average = KahanPlus(P_average,(*P_central)*(1+coord[ip][5]), &error); 
#endif	
      }
    }
#if (!USE_MPI)
    P_average /= np;
#else
    if (notSinglePart) {
      if (parallelStatus!=trueParallel) {
	if (isMaster)            
	  P_average /= np; 
      }
      else {
#ifndef USE_KAHAN    
	MPI_Allreduce(&P_average, &P_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  
#else
        P_total = KahanParallel (P_average, error, MPI_COMM_WORLD);

#endif 
	P_average = P_total/np_total;
      }
    }
    else /* Single particle case, all the processors will do the same as in serial version */
      P_average /= np;
#endif 
#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e in %s at %e to match beam\n",
              *P_central, P_average, trackingContext.elementName, trackingContext.zEnd);
#endif
    if (fabs(P_average-(*P_central))/(*P_central)>1e-14){ 
     /* if (P_average!= *P_central) { */
      if (active) {
	for (ip=0; ip<np; ip++)
	  coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_average)/ P_average;
      }
      *P_central =  P_average;
    }
  }
  else {
    /* change the particle momenta so that the centroid is the central momentum */
    /* the path length is adjusted so that the time-of-flight at the current
       velocity is fixed */
    P_average = 0;
    if (active) {
      for (ip=0; ip<np; ip++) {
#ifndef USE_KAHAN	     
        P_average += (*P_central*(1+coord[ip][5]));
#else
	P_average = KahanPlus(P_average, (*P_central*(1+coord[ip][5])), &error); 
#endif	
      }
    }
#if (!USE_MPI)
    P_average /= np;
#else
    if (notSinglePart) { 
      if (parallelStatus!=trueParallel) {
	if (isMaster)
	  P_average /= np; 
      }
      else {
#ifndef USE_KAHAN    
	MPI_Allreduce(&P_average, &P_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  
#else
        P_total = KahanParallel (P_average, error, MPI_COMM_WORLD);
#endif
	P_average = P_total/np_total;
      }
    }
    else
      P_average /= np;
#endif       
    if (active) {
      dP_centroid =  *P_central - P_average;
      for (ip=0; ip<np; ip++) {
	P = (1+coord[ip][5])*(*P_central);
	t = coord[ip][4]/(P/sqrt(P*P+1));
	P += dP_centroid;
	coord[ip][5] = (P - *P_central)/ (*P_central);
	coord[ip][4] = t*(P/sqrt(P*P+1));
#if defined(IEEE_MATH)
	if (isnan(coord[ip][4]) || isinf(coord[ip][4])) {
	  long i;
	  fprintf(stdout, "error: bad time coordinate for particle %ld\n", ip);
	  fflush(stdout);
	  for (i=0; i<6; i++)
	    fprintf(stdout, "%15.8e ", coord[ip][i]);
	  fflush(stdout);
	  fputc('\n', stdout);
	  fprintf(stdout, "P_average = %e  P_central = %e  t = %e  dP_centroid = %e\n",
		  P_average, *P_central, t, dP_centroid);
	  fflush(stdout);
#if (USE_MPI)
	  if (active)
	    MPI_Abort(MPI_COMM_WORLD, 1);
#endif    
	  abort();
	}
#endif
      }
    }
  }

  log_exit("do_match_energy");
}

void set_central_energy(
                        double **coord, 
                        long np,
                        double new_energy,  /* new central gamma - 1*/
                        double *P_central
                        )
{
  
  log_entry("set_central_energy");
  set_central_momentum(coord, np, sqrt(sqr(new_energy+1)-1), P_central);
  log_exit("set_central_energy");
}

void set_central_momentum(
                          double **coord, 
                          long np,
                          double  P_new,  /* new central beta*gamma */
                          double *P_central
                          )
{
  long ip;

#if (!USE_MPI)  
  if (!np) {
    *P_central =  P_new;
    return;
  }
  if (*P_central != P_new) {
    for (ip=0; ip<np; ip++)
      coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_new)/P_new;
    *P_central =  P_new;
  }
#else
  if (notSinglePart) {
    if (!np) 
      *P_central =  P_new;

    if (*P_central != P_new) {
#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e in %s at %e to match beam\n",
              *P_central, P_new, trackingContext.elementName, trackingContext.zEnd);
#endif
      if (((parallelStatus==trueParallel) && isSlave) || ((parallelStatus!=trueParallel) && isMaster)) {
	for (ip=0; ip<np; ip++)
	  coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_new)/P_new;
      }
      *P_central =  P_new;
    }
  }
  else {
    if (!np) {
      *P_central =  P_new;
      return;
    }
    if (*P_central != P_new) {
#ifdef DEBUG_FIDUCIALIZATION
      fprintf(stdout, "Changing reference momentum from %e to %e in %s\n",
              *P_central, P_new, trackingContext.elementName);
#endif
      for (ip=0; ip<np; ip++)
	coord[ip][5] = ((1+coord[ip][5])*(*P_central) - P_new)/P_new;
      *P_central =  P_new;
    }
  }
#endif
}

void remove_correlations(double **part, REMCOR *remcor, long np)
{
  double sumxy, sumy2, ratio;
  long ip, ic, wc;
  long removeFrom[4];
  
  if (!np) 
    return;
  
  removeFrom[0] = remcor->x;
  removeFrom[1] = remcor->xp;
  removeFrom[2] = remcor->y;
  removeFrom[3] = remcor->yp;
  wc = remcor->with-1;

  for (ip=sumy2=0; ip<np; ip++)
    sumy2 += part[ip][wc]*part[ip][wc];
  if (!sumy2)
    return;
  
  for (ic=0; ic<4; ic++) {
    if (!removeFrom[ic] || ic==wc)
      continue;
    if (!remcor->ratioSet[ic]) {
      for (ip=sumxy=0; ip<np; ip++)
        sumxy += part[ip][ic]*part[ip][wc];
      ratio = sumxy/sumy2;
      if (remcor->onceOnly) {
        remcor->ratio[ic] = ratio;
        remcor->ratioSet[ic] = 1;
      }
    }
    else 
      ratio = remcor->ratio[ic];
    for (ip=0; ip<np; ip++)
      part[ip][ic] -= part[ip][wc]*ratio;
  }
}


void center_beam(double **part, CENTER *center, long np, long iPass, double p0)
{
  double sum, offset;
  long i, ic;

  if (!np) {
    return;
  }

  if (center->onPass>=0 && iPass!=center->onPass)
    return;

  for (ic=0; ic<6; ic++) {
    if (center->doCoord[ic]) {
      if (!center->deltaSet[ic]) {
        for (i=sum=0; i<np; i++)
          sum += part[i][ic];
#if USE_MPI
	if (notSinglePart) {
	  double sum_total;
          long np_total;

          MPI_Allreduce (&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);
          sum = sum_total;
	  MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
          center->delta[ic] = offset = sum/np_total;
	} else
	  center->delta[ic] = offset = sum/np;
#else
        center->delta[ic] = offset = sum/np;
#endif
        if (center->onceOnly)
          center->deltaSet[ic] = 1;
      } else 
        offset = center->delta[ic];
/*      printf("centering coordinate %ld by subtracting %le\n", ic, offset); */
      for (i=0; i<np; i++)
        part[i][ic] -= offset;
    }
  }

  if (center->doCoord[ic=6]) {
    /* Special treatment for time coordinate */
    double *timeCoord;
    timeCoord = tmalloc(sizeof(*timeCoord)*np);
    offset = computeTimeCoordinates(timeCoord, p0, part, np);
    if (center->deltaSet[ic])
      offset = center->delta[ic];
    for (i=0; i<np; i++)
      timeCoord[i] -= offset;
    computeDistanceCoordinates(timeCoord, p0, part, np);
    if (center->onceOnly && !center->deltaSet[ic]) {
      center->delta[ic] = offset;
      center->deltaSet[ic] = 1;
    }
  }
}


void drift_beam(double **part, long np, double length, long order)
{
  VMATRIX *M;
  
  log_entry("drift_beam");
  
  if (length) {
    M = drift_matrix(length, order);
    track_particles(part, M, part, np);
    free_matrices(M);
    tfree(M);
    M = NULL;
  }
  log_exit("drift_beam");
}

void scatter(double **part, long np, double Po, SCATTER *scat)
{
  long i, ip;
  double t, P, beta;
  double sigma[4];

  if (!np)
    return;
  
  log_entry("scatter");
  sigma[0] = scat->x;
  sigma[1] = scat->xp;
  sigma[2] = scat->y;
  sigma[3] = scat->yp;
  for (ip=0; ip<np; ip++) {
    if (scat->probability<1 && random_2(1)>scat->probability)
      continue;
    for (i=0; i<4; i++) {
      if (!sigma[i])
        continue;
      part[ip][i] += gauss_rn(0, random_2)*sigma[i];
    }
    if (scat->dp) {
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      t = part[ip][4]/beta;
      part[ip][5] += scat->dp*gauss_rn(0, random_2);
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      part[ip][4] = t*beta;
    }
  }
  
  log_exit("scatter");
}

void store_fitpoint_matrix_values(MARK *fpt, char *name, long occurence, VMATRIX *M)
{
  char buffer[1000];
  long i, j, k, l, count;

  if (!M) 
    return;
  if (!M->R)
    bombElegant("NULL R matrix passed to store_fitpoint_matrix_values", NULL);

  if (!(fpt->init_flags&8)) {
    if (M->order==1) {
      if (!(fpt->matrix_mem = malloc(sizeof(*(fpt->matrix_mem))*(6+36))))
        bombElegant("memory allocation failure (store_fitpoint_matrix_values)", NULL);
    } else if (M->order==2) {
      if (!(fpt->matrix_mem = malloc(sizeof(*(fpt->matrix_mem))*(6+36+126))))
        bombElegant("memory allocation failure (store_fitpoint_matrix_values)", NULL);
    } else
      if (!(fpt->matrix_mem = malloc(sizeof(*(fpt->matrix_mem))*(6+36+126+336))))
        bombElegant("memory allocation failure (store_fitpoint_matrix_values)", NULL);
    for (i=count=0; i<6; i++) {
      sprintf(buffer, "%s#%ld.C%ld", name, occurence, i+1);
      fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
    }
    for (i=0; i<6; i++) {
      for (j=0; j<6; j++) {
        sprintf(buffer, "%s#%ld.R%ld%ld", name, occurence, i+1, j+1);
        fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
      }
    }
    if (M->order>1) {
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          for (k=0; k<=j; k++) {
            sprintf(buffer, "%s#%ld.T%ld%ld%ld", name, occurence, i+1, j+1, k+1);
            fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
          }
        }
      }
    }
    if (M->order>2) {
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          for (k=0; k<=j; k++) {
            for (l=0; l<=k; l++) {
              sprintf(buffer, "%s#%ld.U%ld%ld%ld%ld", name, occurence, i+1, j+1, k+1, l+1);
              fpt->matrix_mem[count++] = rpn_create_mem(buffer, 0);
            }
          }
        }
      }
    }    
    fpt->init_flags |= 8;
  }
  
  for (i=count=0; i<6; i++)
    rpn_store(M->C[i], NULL, fpt->matrix_mem[count++]);
  for (i=0; i<6; i++)
    for (j=0; j<6; j++)
      rpn_store(M->R[i][j], NULL, fpt->matrix_mem[count++]);
  if (M->order>1)
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++) 
          rpn_store(M->T[i][j][k], NULL, fpt->matrix_mem[count++]);
  if (M->order>2)
    for (i=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++) 
            for (l=0; l<=k; l++) 
              rpn_store(M->Q[i][j][k][l], NULL, fpt->matrix_mem[count++]);
}

void store_fitpoint_beam_parameters(MARK *fpt, char *name, long occurence, double **coord, long np, double Po)
{
  long i, j, k;
  static double emit[3], sigma[6], centroid[6], beta[3], alpha[3], emitc[3];
  static BEAM_SUMS sums;
  static char *centroid_name_suffix[8] = {
    "Cx", "Cxp", "Cy", "Cyp", "Cs", "Cdelta", "pCentral", "Particles" };
  static char *sigma_name_suffix[6] = {
    "Sx", "Sxp", "Sy", "Syp", "Ss", "Sdelta" };
  static char *emit_name_suffix[5] = {
    "ex", "ey", "es", "ecx", "ecy"};
  static char *beta_name_suffix[2] = {
    "betaxBeam", "betayBeam", 
  };
  static char *alpha_name_suffix[2] = {
    "alphaxBeam", "alphayBeam",
  };
  static char s[1000];

  zero_beam_sums(&sums, 1);
  accumulate_beam_sums(&sums, coord, np, Po);
  if (isMaster || !notSinglePart) {
    for (i=0; i<6; i++) {
      centroid[i] = sums.centroid[i];
      sigma[i] = sqrt(sums.sigma[i][i]);
      if (i%2==0) {
        beta[i/2] = alpha[i/2] = emitc[i/2] = emit[i/2] = 0;
        computeEmitTwissFromSigmaMatrix(emit+i/2, emitc+i/2, beta+i/2, alpha+i/2, sums.sigma, i);
/*        printf("%s#%ld : emit = %e, beta = %e, alpha = %e\n",
               name, occurence, emit[i/2], beta[i/2], alpha[i/2]); */
      }
    }
    
    if (!(fpt->init_flags&2)) {
      fpt->centroid_mem = tmalloc(sizeof(*fpt->centroid_mem)*8);
      fpt->sigma_mem = tmalloc(sizeof(*fpt->sigma_mem)*6);
      fpt->emit_mem = tmalloc(sizeof(*fpt->emit_mem)*5);
      fpt->betaBeam_mem = tmalloc(sizeof(*fpt->betaBeam_mem)*2);
      fpt->alphaBeam_mem = tmalloc(sizeof(*fpt->alphaBeam_mem)*2);
      fpt->sij_mem = tmalloc(sizeof(*fpt->sigma_mem)*15);
      for (i=0; i<8; i++) {
	sprintf(s, "%s#%ld.%s", name, occurence, centroid_name_suffix[i]);
	fpt->centroid_mem[i] = rpn_create_mem(s, 0);
      }
      for (i=0; i<6; i++) {
	sprintf(s, "%s#%ld.%s", name, occurence, sigma_name_suffix[i]);
	fpt->sigma_mem[i] = rpn_create_mem(s, 0);
      }
      for (i=0; i<5; i++) {
	sprintf(s, "%s#%ld.%s", name, occurence, emit_name_suffix[i]);
	fpt->emit_mem[i] = rpn_create_mem(s, 0);
      }
      for (i=0; i<2; i++) {
	sprintf(s, "%s#%ld.%s", name, occurence, beta_name_suffix[i]);
	fpt->betaBeam_mem[i] = rpn_create_mem(s, 0);
	sprintf(s, "%s#%ld.%s", name, occurence, alpha_name_suffix[i]);
	fpt->alphaBeam_mem[i] = rpn_create_mem(s, 0);
      }
      for (i=k=0; i<6; i++) {
	for (j=i+1; j<6; j++, k++) {
	  sprintf(s, "%s#%ld.s%ld%ld", name, occurence, i+1, j+1);
	  fpt->sij_mem[k] = rpn_create_mem(s, 0);
	}
      }
      fpt->init_flags |= 2;
    }
    for (i=0; i<6; i++) {
      rpn_store(centroid[i], NULL, fpt->centroid_mem[i]);
      rpn_store(sigma[i], NULL, fpt->sigma_mem[i]);
    }
    for (i=0; i<3; i++)
      rpn_store(emit[i], NULL, fpt->emit_mem[i]);
    for (i=0; i<2; i++) {
      rpn_store(emitc[i], NULL, fpt->emit_mem[i+3]);
      rpn_store(beta[i], NULL, fpt->betaBeam_mem[i]);
      /* printf("%s#%ld.%s = %e\n", name, occurence, beta_name_suffix[i], beta[i]); */
      rpn_store(alpha[i], NULL, fpt->alphaBeam_mem[i]);
    }
    for (i=k=0; i<6; i++)
      for (j=i+1; j<6; j++, k++)
	rpn_store(sums.sigma[i][j], NULL, fpt->sij_mem[k]);
    rpn_store(Po, NULL, fpt->centroid_mem[6]);
    rpn_store((double)np, NULL, fpt->centroid_mem[7]);
  }
}

ELEMENT_LIST *findBeamlineMatrixElement(ELEMENT_LIST *eptr)
{
  ELEMENT_LIST *eptr0=NULL, *eptrPassed;
  long matrixSeen = 0;
  eptrPassed = eptr;
  while (eptr) {
    if (eptr->type==T_MATR) {
      eptr0 = eptr;
      matrixSeen = 1;
      eptr = eptr->succ;
      break;
    }
    eptr = eptr->succ;
  }
  if (!matrixSeen)
    bombElegant("Can't do \"linear chromatic\" or \"longitudinal-only\" matrix tracking---no matrices!", NULL);
  while (eptr) {
    if ((eptr->p_elem || eptr->matrix) && eptr->type==T_MATR) {
      fprintf(stderr, "***** WARNING ****\n");
      fprintf(stderr, "Possible problem with \"linear chromatic\" or \"longitudinal-only\" matrix tracking\n");
      fprintf(stderr, "Concatenation resulted in more than one matrix.  Make sure the additional matrices do\n");
      fprintf(stderr, "not affect the revolution matrix!\n");
      print_elem_list(stderr, eptrPassed);
      fprintf(stderr, "***** WARNING ****\n");
      break;
    }
    eptr = eptr->succ;
  }
  return eptr0;
}

long trackWithIndividualizedLinearMatrix(double **particle, long particles, double **accepted,
                                    double Po, double z, ELEMENT_LIST *eptr,
                                    TWISS *twiss,
                                    double *tune0,
                                    double *chrom,    /* d   nu /ddelta   */
                                    double *chrom2,   /* d^2 nu /ddelta^2 */
                                    double *chrom3,   /* d^3 nu /ddelta^3 */
                                    double *dbeta_dPoP, 
                                    double *dalpha_dPoP,
                                    double *alphac,   /* Cs = Cs(0) + delta*alphac[0] + delta^2*alphac[1] */
                                    double *eta2,     /* x = x(0) + eta*delta + eta2*delta^2 */
                                    ILMATRIX *ilmat   /* used only if twiss==NULL */
                                    )
{
  long ip, plane, offset, i, j, itop, is_lost;
  double *coord, deltaPoP, tune2pi, sin_phi, cos_phi;
  double alpha[2], beta[2], eta[4], beta1, alpha1, A[2];
  double R11, R22, R12;
  long allowResonanceCrossing = 0;
  static VMATRIX *M1 = NULL;
  double det;
  
  if (ilmat)
    allowResonanceCrossing = ilmat->allowResonanceCrossing;
  
  if (!M1) {
    M1 = tmalloc(sizeof(*M1));
    initialize_matrices(M1, 1);
  }
  if (twiss) {
    beta[0] = twiss->betax;
    beta[1] = twiss->betay;
    alpha[0] = twiss->alphax;
    alpha[1] = twiss->alphay;
    eta[0] = twiss->etax;
    eta[1] = twiss->etapx;
    eta[2] = twiss->etay;
    eta[3] = twiss->etapy;
  } else {
    for (i=0; i<2; i++) {
      beta[i] = ilmat->beta[i];
      alpha[i] = ilmat->alpha[i];
    }
    for (i=0; i<4; i++)
      eta[i] = ilmat->eta[i];
  }
    
  for (i=0; i<6; i++) {
    M1->C[i] = eptr->matrix->C[i];
    for (j=0; j<6; j++)
      M1->R[i][j] = i==j?1:0;
  }

  if (ilmat && ilmat->tilt) 
    rotateBeamCoordinates(particle, particles, ilmat->tilt);
  
  itop = particles-1;
  for (ip=0; ip<particles; ip++) {
    coord = particle[ip];
    deltaPoP = coord[5];
    /* remove the dispersive orbit from the particle coordinates */
    coord[5] -= deltaPoP;
    for (plane=0; plane<2; plane++) {
      coord[2*plane]   -= deltaPoP*(eta[2*plane]   + deltaPoP*eta2[2*plane]);
      coord[2*plane+1] -= deltaPoP*(eta[2*plane+1] + deltaPoP*eta2[2*plane+1]);
      /* compute the betatron amplitude, if needed */
      A[plane] = 0;
      if (ilmat)
        A[plane] = (sqr(coord[2*plane]) + ipow(alpha[plane]*coord[2*plane]+beta[plane]*coord[2*plane+1], 2))/beta[plane];
    }
    is_lost = 0;
    for (plane=0; !is_lost && plane<2; plane++) {
      tune2pi = PIx2*(tune0[plane] + 
                      deltaPoP*(chrom[plane] +
                                deltaPoP/2*(chrom2[plane] + 
                                            deltaPoP/3*chrom3[plane])));
      if (ilmat)
        tune2pi += PIx2*(A[0]*ilmat->tswax[plane] + A[1]*ilmat->tsway[plane]);
      offset = 2*plane;
      if ((beta1 = beta[plane]+dbeta_dPoP[plane]*deltaPoP)<=0) {
        fprintf(stdout, "nonpositive beta function for particle with delta=%le\n",
                deltaPoP);
        fprintf(stdout, "particle is lost\n");
        is_lost = 1;
        continue;
      }
      if (!allowResonanceCrossing && fabs( ((long)(2*tune2pi/PIx2)) - ((long)(2*tune0[plane]))) != 0) {
        fprintf(stdout, "particle with delta=%le crossed integer or half-integer resonance\n",
                deltaPoP);
        fprintf(stdout, "particle is lost\n");
        is_lost = 1;
        continue;
      }
      /* R11=R22 or R33=R44 */
      sin_phi = sin(tune2pi);
      cos_phi = cos(tune2pi);
      alpha1 = alpha[plane]+dalpha_dPoP[plane]*deltaPoP;
      /* R11 or R33 */
      R11 = M1->R[0+offset][0+offset] = cos_phi + alpha1*sin_phi;
      /* R22 or R44 */
      R22 = M1->R[1+offset][1+offset] = cos_phi - alpha1*sin_phi;
      /* R12 or R34 */
      if ((R12 = M1->R[0+offset][1+offset] = beta1*sin_phi)) {
        /* R21 or R43 */
        M1->R[1+offset][0+offset] = (R11*R22-1)/R12;
      }
      else {
        bombElegant("divided by zero in trackWithChromaticLinearMatrix", NULL);
      }
      det = M1->R[0+offset][0+offset]*M1->R[1+offset][1+offset] -
        M1->R[0+offset][1+offset]*M1->R[1+offset][0+offset];
      if (fabs(det-1)>1e-6) {
        fprintf(stdout, "Determinant is suspect for particle with delta=%e\n", deltaPoP);
        fprintf(stdout, "particle is lost\n");
        is_lost = 1;
        continue;
      }
    }
    if (is_lost) {
      swapParticles(particle[ip], particle[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      particle[itop][4] = z;
      particle[itop][5] = Po*(1+deltaPoP);
      --itop;
      --ip;
      --particles;
    } else {
      /* momentum-dependent pathlength --- note that other path-length terms are ignored ! */
      M1->C[4] = eptr->matrix->C[4]*(1 + deltaPoP*(alphac[0] + deltaPoP*alphac[1]));
      track_particles(&coord, M1, &coord, 1);
      /* add back the dispersive orbit to the particle coordinates */
      coord[5] += deltaPoP;
      for (plane=0; plane<2; plane++) {
        coord[2*plane]   += deltaPoP*(eta[2*plane]   + deltaPoP*eta2[2*plane]);
        coord[2*plane+1] += deltaPoP*(eta[2*plane+1] + deltaPoP*eta2[2*plane+1]);
      }
    }
  }

  if (ilmat && ilmat->tilt) 
    rotateBeamCoordinates(particle, particles, -ilmat->tilt);

  return particles;
}
  

void trackLongitudinalOnlyRing(double **part, long np, VMATRIX *M, double *alpha)
{
  long ip;
  double *coord, length, alpha1, alpha2;
  
  alpha1 = alpha[0];
  alpha2 = alpha[1];
  length = M->C[4];
  for (ip=0; ip<np; ip++) {
    coord = part[ip];
    coord[0] = coord[1] = coord[2] = coord[3] = 0;
    coord[4] += length*(1+(alpha1+alpha2*coord[5])*coord[5]);
  }
}

void getTrackingContext(TRACKING_CONTEXT *trackingContext0) 
{
  memcpy(trackingContext0, &trackingContext, sizeof(trackingContext));
}

void matr_element_tracking(double **coord, VMATRIX *M, MATR *matr,
                           long np, double z)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  long i;
#if !USE_MPI
  if (!np)
    return;
#endif
  if (!matr) {
    track_particles(coord, M, coord, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      for (i=0; i<np; i++)
        sum += coord[i][4];
#if !USE_MPI
      matr->sReference = sum/np;
#else
      if (notSinglePart) {
	if (isSlave) {
	  double sum_total;
	  long np_total;

	  MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	  MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);	      
	  matr->sReference = sum_total/np_total;
	}
      } else
	matr->sReference = sum/np;
#endif
      matr->fiducialSeen = 1;
    }
    for (i=0; i<np; i++)
      coord[i][4] -= matr->sReference;
    track_particles(coord, M, coord, np);
    for (i=0; i<np; i++)
      coord[i][4] += matr->sReference;
  }
}

void ematrix_element_tracking(double **coord, VMATRIX *M, EMATRIX *matr,
			      long np, double z, double *P_central)
/* subtract off <s> prior to using a user-supplied matrix to avoid possible
 * problems with R5? and T?5? elements
 */
{
  long i;

#if !USE_MPI
  if (!np)
    return;
#endif
  if (!matr) {
    fprintf(stderr, "ematrix_element_tracking: matr=NULL, tracking with M (%ld order)\n",
            M->order);
    track_particles(coord, M, coord, np);
  } else {
    if (!matr->fiducialSeen) {
      double sum = 0;
      for (i=0; i<np; i++)
        sum += coord[i][4];
#if !USE_MPI
      matr->sReference = sum/np;
#else
      if (notSinglePart) {
	if (isSlave) {
	  double sum_total;
	  long np_total;
	
	  MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
	  MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, workers);	      
	  matr->sReference = sum_total/np_total;
	}
      } else
	matr->sReference = sum/np;
#endif
      matr->fiducialSeen = 1;
    }
    for (i=0; i<np; i++)
      coord[i][4] -= matr->sReference;
    track_particles(coord, M, coord, np);
    for (i=0; i<np; i++)
      coord[i][4] += matr->sReference;
  }
  if (matr->deltaP)
    *P_central += matr->deltaP;
}

long transformBeamWithScript(SCRIPT *script, double pCentral, CHARGE *charge, 
                             BEAM *beam, double **part, long np, long *nLost,
                             char *mainRootname, long iPass, long driftOrder)
{
  char *rootname, *input, *output=NULL;
  char *cmdBuffer0, *cmdBuffer1=NULL;
  SDDS_DATASET SDDSout, SDDSin;
  double *data;
  char *dataname[6] = {"x","xp","y","yp","t","p"};
  long i, j, npNew, nameLength, doDrift;
  char passString[20];
#if !USE_MPI
  long k, lostIndex;
#else
  long npTotal=0, rootnameLength;

  if (notSinglePart) { 
    MPI_Allreduce (&np, &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (!beam->n_to_track_total)
       return 0;
  }	
#endif

  if (!script->rootname || !strlen(script->rootname)) {
    /* generate random rootname */
    if (isMaster)
    if (!(rootname = tmpname(NULL)))
      bombElegant("problem generating temporary filename for script", NULL);
#if SDDS_MPI_IO
    if (isMaster)
      rootnameLength = strlen(rootname)+1;
    MPI_Bcast(&rootnameLength, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (isSlave) /* Master and slave could have different rootname length if use C-shell, which calls tmpname one more time on master for execution.*/
      rootname = malloc(rootnameLength); 
  /* As different processors will have different process names, we need
     make sure they have the same file name for parallel I/O */
    MPI_Bcast(rootname, rootnameLength, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
  } else 
    rootname = compose_filename(script->rootname, mainRootname);
  nameLength = (script->directory?strlen(script->directory):0) + \
    strlen(rootname) + strlen(script->inputExtension) +
    strlen(script->outputExtension) + 4;
  if (!(input = malloc(sizeof(*input)*nameLength)) ||
      !(output = malloc(sizeof(*output)*nameLength)))
    bombElegant("problem generating temporary filename for script", NULL);

  doDrift = 0;
  if (script->onPass>=0) {
    if (script->onPass!=iPass)
      doDrift = 1;
  } else if (script->startPass>=0) {
    if (script->startPass>iPass)
      doDrift = 1;
  }

  if (doDrift) {
    drift_beam(part, np, script->length, driftOrder);
    return np;
  }
  
  /* prepare command */
  if (script->directory && strlen(script->directory)) {
#if defined(_WIN32)
    sprintf(input, "%s\\%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s\\%s.%s", script->directory, rootname, script->outputExtension);
#else
    sprintf(input, "%s/%s.%s", script->directory, rootname, script->inputExtension);
    sprintf(output, "%s/%s.%s", script->directory, rootname, script->outputExtension);
#endif
  } else {
    sprintf(input, "%s.%s", rootname, script->inputExtension);
    sprintf(output, "%s.%s", rootname, script->outputExtension);
  }
  if (rootname!=script->rootname)
    free(rootname);

 sprintf(passString, "%ld", iPass);

  if (!(cmdBuffer0=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))) ||
      !(cmdBuffer1=malloc(sizeof(char)*(strlen(script->command)+10*strlen(input)+10*strlen(output)+strlen(passString)))))
    bombElegant("memory allocation failure making command buffer for script", NULL);
  replaceString(cmdBuffer0, script->command, "%i", input, 9, 0);
  replaceString(cmdBuffer1, cmdBuffer0, "%o", output, 9, 0);
 
  replaceString(cmdBuffer0, cmdBuffer1, "%p", passString, 9, 0);
  strcpy_ss(cmdBuffer1, cmdBuffer0);
  
  /* substitute numerical parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], value[25], *ptr;
    sprintf(tag, "%%np%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    sprintf(value, "%21.15e", script->NP[i]);
    if (!(cmdBuffer0 = SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))) ||
        !(cmdBuffer1 = SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*25+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, value, count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  /* substitute string parameters */
  for (i=0; i<10; i++) {
    long count = 0;
    char tag[10], *ptr;
    if (!script->SP[i] || strlen(script->SP[i])==0)
      continue;
    sprintf(tag, "%%sp%ld", i);
    ptr = cmdBuffer1;
    while ((ptr=strstr(ptr, tag))) {
      count ++;
      ptr += 3;
    }
    if (!count) continue;
    if (!(cmdBuffer0 = 
          SDDS_Realloc(cmdBuffer0, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))) ||
        !(cmdBuffer1 = 
          SDDS_Realloc(cmdBuffer1, sizeof(*cmdBuffer1)*(strlen(cmdBuffer1)+count*strlen(script->SP[i])+1))))
      SDDS_Bomb("memory allocation failure");
    replaceString(cmdBuffer0, cmdBuffer1, tag, script->SP[i], count, 0);
    strcpy_ss(cmdBuffer1, cmdBuffer0);
  }
  interpret_escaped_quotes(cmdBuffer1);
  
  if (script->verbosity>0) {
    fprintf(stdout, "%s\n", cmdBuffer1);
    fflush(stdout);
  }
 
  /* dump the data to script input file */
#if USE_MPI
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  {
    SDDS_ForceInactive(&SDDSout);
    SDDS_PhaseSpaceSetup(&SDDSout, input, SDDS_BINARY, 1, "script input", 
			 "unknown", "unknown",
			 "transformBeamWithScript");
    dump_phase_space(&SDDSout, part, np, 0, pCentral, charge?charge->macroParticleCharge*np:0.0);
    
    if (!SDDS_Terminate(&SDDSout))
      SDDS_Bomb("problem terminating script input file");
  }
#if defined(CONDOR_COMPILE)
  _condor_ckpt_disable();
#endif

  /* run the script */
  if (isMaster) /* This will be done on the master */
  {  
    if (script->useCsh)
      executeCshCommand(cmdBuffer1);
    else 
      system(cmdBuffer1);
  }
#if defined(CONDOR_COMPILE)
  _condor_ckpt_enable();
#endif

  if (script->verbosity>0) {
    fprintf(stdout, "Command completed\n");
    fflush(stdout);
  }

  /* read the data from script output file */
#if SDDS_MPI_IO
  MPI_Barrier(MPI_COMM_WORLD);

  if (notSinglePart) {
    if (!fexists(output)) 
      SDDS_Bomb("unable to find script output file");
    SDDSin.parallel_io = 1;
    /* set up parallel IO information */      
    SDDS_MPI_Setup(&SDDSin, 1, n_processors, myid, MPI_COMM_WORLD, 0);
    if (!SDDS_MPI_InitializeInput(&SDDSin, output)) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  } else if (isMaster){
    if (!fexists(output)) 
      SDDS_Bomb("unable to find script output file");
    if (!SDDS_InitializeInput(&SDDSin, output)) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
#else  
  if (!fexists(output)) 
    SDDS_Bomb("unable to find script output file");
  if (!SDDS_InitializeInput(&SDDSin, output)) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
#endif

#if USE_MPI
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m") ||
      !check_sdds_column(&SDDSin, "xp", NULL) ||
      !check_sdds_column(&SDDSin, "yp", NULL) ||
      !check_sdds_column(&SDDSin, "p", "m$be$nc") ||
      !check_sdds_column(&SDDSin, "t", "s")) {
    if (!check_sdds_column(&SDDSin, "p", "m$be$nc") &&
        check_sdds_column(&SDDSin, "p", NULL)) {
      fprintf(stdout, "Warning: p has no units in script output file.  Expected m$be$nc\n");
      fflush(stdout);
    } else {
      fprintf(stdout, 
              "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in script output");
      fflush(stdout);
      exitElegant(1);
    }
  }
 
#if !SDDS_MPI_IO  
  if (SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("Unable to read script output file");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
#else
  if (notSinglePart)
    SDDS_MPI_ReadPage(&SDDSin);
  else if (isMaster) 
    SDDS_ReadPage(&SDDSin);
  if (notSinglePart || (!notSinglePart&&isMaster))
#endif
  npNew = SDDS_RowCount(&SDDSin);

#if !USE_MPI
  if (script->verbosity>0) {
    fprintf(stdout, "%ld particles in script output file (was %ld)\n", npNew, np);
    fflush(stdout);
  }

  if (!npNew) {
    return 0;
  }
  if (npNew>np)
#else
  if (!notSinglePart) 
    MPI_Bcast(&npNew, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  else {
    if (beam) {
      npTotal = beam->n_to_track_total;
      beam->n_to_track_total = SDDS_MPI_TotalRowCount(&SDDSin);
      if (script->verbosity>0) {
        fprintf(stdout, "%ld particles in script output file (was %ld)\n", beam->n_to_track_total, npTotal);
        fflush(stdout);
      }	
    }
  }
  if (beam && script->verbosity>0) {
    fprintf(stdout, "%ld particles in script output file (was %ld)\n", npTotal, beam->n_to_track_total);
    fflush(stdout);
  }
  if ((!notSinglePart && (npNew>np))||(notSinglePart && (beam->n_to_track_total>npTotal)))
#endif
  if (script->noNewParticles)
    bombElegant("The number of particles increased after the SCRIPT element without seting the correct flag.\n Please set the correct NO_NEW_PARTICLES flag for the SCRIPT element!", NULL);

  if (!script->useParticleID) {
    if ((!script->noNewParticles || (npNew!=np)) && (iPass==0)){
      fprintf (stdout, "Warning: There is no particle ID available to find which particles are lost. The particles lost in the SCRIPT element will not be recorded!\n");
    }
    if ((!script->noNewParticles) && (npNew < np)) {
      /* Do particle accounting only */
      beam->n_to_track = npNew+*nLost;
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=0; i<np-npNew; i++) {
	  swapParticles(beam->lost[npNew+i], beam->lost[beam->n_to_track+i]);
	}
    }
  }    
      
  if (npNew>np) {
    /* We may have to resize the arrays in the BEAM structure */
    
    fprintf(stdout, "Increasing number of particles from %ld (%ld active) to %ld (%ld active)\n",
            np+*nLost, np, npNew+*nLost, npNew);
    
    if (!beam) {
      fprintf(stderr, "Error: script element increased the number of particles from %ld to %ld\n.",
              np, npNew);
      fprintf(stderr, "This happened (apparently) during a pre-tracking stage, which isn't allowed\n");
      exitElegant(1);
    }
    /* Check that previous particle counts are correct */
    if ((np+*nLost)!=beam->n_to_track) {
      fprintf(stderr, "Particle accounting problem in SCRIPT element:\n");
      fprintf(stderr, "np = %ld, *nLost = %ld, beam->n_to_track = %ld, beam->n_particle=%ld\n",
              np, *nLost, beam->n_to_track, beam->n_particle);
      fprintf(stderr, "This could happen if the particleID is not unique.\n");
      exitElegant(1);
    }

    if ((npNew+*nLost) > beam->n_particle) {
      if (beam->original==beam->particle) {
        /* This means, oddly enough, that the particle array and original array are the same because the
         * separate original array wasn't needed.  n_original gives the size of both arrays (including
         * live and lost particles).  To avoid confusion, we'll copy the data to a new array before
         * doing anything else, even though it means the original array is not used for anything and
         * contains a useless frozen copy of the present beam.
         * Use n_original since that's the size of the array, including lost particles. 
         */
        beam->particle = (double**)czarray_2d(sizeof(double), beam->n_original, 7);
        copy_particles(beam->particle, beam->original, beam->n_original);
      }
      /* resize the particle array, leaving space for the lost particle data at the top */
      if (!(beam->particle = (double**)resize_czarray_2d((void**)beam->particle,sizeof(double), npNew+*nLost, 7)) ||
          !(beam->lost = (double**)resize_czarray_2d((void**)beam->lost,sizeof(double), npNew+*nLost, 8))) {
        fprintf(stderr, "Memory allocation failure increasing particle array size to %ld\n",
                npNew+*nLost);
        exitElegant(1);
      }
      beam->n_particle = npNew+*nLost;
#if !USE_MPI
      /* move lost particles into the upper part of the arrays */
      if (beam->lost)
	for (i=*nLost-1; i>=0; i--) {
	  swapParticles(beam->lost[np+i], beam->lost[npNew+i]);
      }
#endif
    }
    
    if (beam->accepted)  {
      /* this data is invalid when particles are added */
      free_czarray_2d((void**)beam->accepted, np+*nLost, 7);
      beam->accepted = NULL;
    }
    beam->n_to_track = npNew+*nLost;
    fprintf(stdout, "beam->n_particle = %ld, beam->n_to_track = %ld\n",
	    beam->n_particle, beam->n_to_track);
   
    part = beam->particle;
  }

#if USE_MPI
  /* Particles could be redistributed, move lost particles into the upper part of the arrays */
  if ((np != npNew) && beam && nLost && beam->lost)
    for (i=0; i<=*nLost-1; i++) {
      swapParticles(beam->lost[np+i], beam->lost[npNew+i]);
    }
  if ((isSlave && notSinglePart) || (isMaster && !notSinglePart)) 
#endif
  for (i=0; i<6; i++) {
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, dataname[i]))) {
      SDDS_SetError("Unable to read script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    for (j=0; j<npNew; j++) {
      if (part[j][i] != data[j])
      part[j][i] = data[j];
    }
    free(data);
  }

#if USE_MPI
  if (!notSinglePart)
    MPI_Bcast(part[0], 7*npNew, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (script->useParticleID ) {
#if USE_MPI
    if ((isSlave && notSinglePart) || (isMaster && !notSinglePart))
#endif
    if (!(data = SDDS_GetColumnInDoubles(&SDDSin, "particleID"))) {
      SDDS_SetError("Unable to read particleID from script output file. Please set USE_PARTICLE_ID=0 if this is desired.\n");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    if (!script->noNewParticles && (iPass==0))
	fprintf (stdout, "Warning: New particles are added in the SCRIPT element. The particles lost in the SCRIPT element will not be recorded!\n");
    else {
#if !USE_MPI
      if (npNew<np) {
	/* Find out which particles are lost */
	lostIndex = npNew;
	for (i=0; i<np; i++) {
	  for (j=0; j<npNew; j++) {
	    if (part[i][6]==data[j])
	      break; /* this particle survived */
	  }
	  if (j==npNew) { /* record the lost particle */
	    /* Copy the lost particles in the SCRIPT element into the upper part of the arrays. */
	    (*nLost)++;
	    if (beam->lost) {
	      for (k=0; k<7; k++) {
		beam->lost[lostIndex][k] = part[i][k];
	      }
	      beam->lost[lostIndex][7] = (double) iPass; 
	      lostIndex++;
	    }
	  }
	}
      }
#else
      /* Even though particle ID is available, we will not record lost particle coordinates due to particle
	 redistribution in Pelegant. The particles move from one processor to another make it very complicated 
	 to find out which particles are lost. While it should be not hard to find the lost particles within 
	 the script if there is such an interest. */
      if (beam && (beam->n_to_track_total<npTotal) && (iPass==0))
	fprintf (stdout, "Warning: Lost particle coordinates in the SCRIPT element will not be recorded.\n");
#endif
    } 
  }

  /* assign new particle IDs if there are new particles */
#if !USE_MPI
  if (npNew>np && script->useParticleID) {
    for (j=0; j<npNew; j++)
      part[j][6] = j+1;
  }
#else
  if (isSlave && notSinglePart) {
    if (beam && (beam->n_to_track_total>npTotal) && script->useParticleID) {
      long sum=0, tmp, my_offset=0, *offset = tmalloc(n_processors*sizeof(*offset));
      MPI_Allgather (&npNew, 1, MPI_LONG, offset, 1, MPI_LONG, workers);
      tmp = offset[0];
      for (i=1; i<n_processors; i++) {
        sum += tmp;
        tmp = offset[i];
        offset[i] = sum;
      }
      offset[0] = 0;
      my_offset = offset[myid-1];
      tfree(offset);
      for (j=0; j<npNew; j++){
	part[j][6] = j+1+my_offset;
      }
    }
  }
#endif

#if USE_MPI
  if ((!notSinglePart&&isMaster) || notSinglePart)
#endif
  if (charge) {
    double totalCharge;
    if (!SDDS_GetParameterAsDouble(&SDDSin, "Charge", &totalCharge)) {
      SDDS_SetError("Unable to read Charge parameter from script output file");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    charge->charge = totalCharge;
    charge->macroParticleCharge = 0;
#if USE_MPI
    if (!notSinglePart) {
      if (npNew)
	charge->macroParticleCharge = totalCharge/npNew;
    } else
      if (beam->n_to_track_total)
	charge->macroParticleCharge = totalCharge/beam->n_to_track_total;
#else
    if (npNew)
      charge->macroParticleCharge = totalCharge/npNew;
#endif 
 }

#if USE_MPI
  if (charge && notSinglePart) {
    MPI_Bcast(&(charge->charge), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&charge->macroParticleCharge, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if (notSinglePart || (isMaster && !notSinglePart)) 
#endif
  {
  if (SDDS_ReadPage(&SDDSin)!=-1)
    SDDS_Bomb("Script output file has multiple pages");
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (script->verbosity) {
    fprintf(stdout, "done with file\n");
    fflush(stdout);
  }

  /* convert (t, p) data to (s, delta) */
  for (j=0; j<npNew; j++) {
    double p, beta;
    p = part[j][5];
    part[j][5] = (p-pCentral)/pCentral;
    beta = p/sqrt(sqr(p)+1);
    part[j][4] *= beta*c_mks;
  }

#if USE_MPI
  if (isMaster)
#endif
  if (!script->keepFiles) {
    /* delete the input and output files */
    remove(input);
    remove(output);
  }

  /* clean up */
  free(cmdBuffer0);
  free(cmdBuffer1);
  
  return npNew;
}

void distributionScatter(double **part, long np, double Po, DSCATTER *scat, long iPass)
{
  static DSCATTER_GROUP *dscatterGroup = NULL;
  static long dscatterGroups = 0;
  long i, ip, interpCode, nScattered, nLeftThisPass;
  double t, P, beta, amplitude, cdf;
  TRACKING_CONTEXT context;
    
  if (!np)
    return;
  getTrackingContext(&context);

  if (!scat->initialized) {
    SDDS_DATASET SDDSin;
    static char *planeName[3] = {"xp", "yp", "dp"};
    static short planeIndex[3] = {1, 3, 5};
    scat->initialized = 1;
    scat->indepData = scat->cdfData = NULL;
    scat->groupIndex = -1;
    if ((i=match_string(scat->plane, planeName, 3, MATCH_WHOLE_STRING))<0) {
      fprintf(stderr, "Error for %s: plane is not valid.  Give xp, yp, or dp\n",
              context.elementName);
      exitElegant(1);
    }
    scat->iPlane = planeIndex[i];
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, scat->fileName) ||
        SDDS_ReadPage(&SDDSin)!=1) {
      fprintf(stderr, "Error for %s: file is not valid.\n", context.elementName);
      exitElegant(1);
    }
    if ((scat->nData=SDDS_RowCount(&SDDSin))<2) {
      fprintf(stderr, "Error for %s: file contains insufficient data.\n", context.elementName);
      exitElegant(1);
    }
    /* Get independent data */
    if (!(scat->indepData=SDDS_GetColumnInDoubles(&SDDSin, scat->valueName))) {
      fprintf(stderr, "Error for %s: independent variable data is invalid.\n",
              context.elementName);
      exitElegant(1);
    }
    /* Check that independent data is monotonically increasing */
    for (i=1; i<scat->nData; i++)
      if (scat->indepData[i]<=scat->indepData[i-1]) {
        fprintf(stderr, "Error for %s: independent variable data is not monotonically increasing.\n",
                context.elementName);
        exitElegant(1);
      }
    /* Get CDF or PDF data */
    if (!(scat->cdfData=SDDS_GetColumnInDoubles(&SDDSin, 
                                                scat->cdfName?scat->cdfName:scat->pdfName))) {
      fprintf(stderr, "Error for %s: CDF/PDF data is invalid.\n",
              context.elementName);
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    SDDS_Terminate(&SDDSin);
    if (!(scat->cdfName)) {
      /* must integrate to get the CDF */
      double *integral, *ptr;
      if (!(integral=malloc(sizeof(*integral)*scat->nData))) {
        fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
        exitElegant(1);
      }
      trapazoidIntegration1(scat->indepData, scat->cdfData, scat->nData, integral);
      ptr = scat->cdfData;
      scat->cdfData = integral;
      free(ptr);
    }
    /* Check that CDF is monotonically increasing */
    for (i=0; i<scat->nData; i++) {
      if (scat->cdfData[i]<0 || (i && scat->cdfData[i]<=scat->cdfData[i-1])) {
        long j;
        fprintf(stderr, "Error for %s: CDF not monotonically increasing at index %ld.\n", context.elementName,
                i);
        for (j=0; j<=i+1 && j<scat->nData; j++)
          fprintf(stderr, "%ld %21.15e\n", j, scat->cdfData[i]);
        exitElegant(1);
      }
    }
    /* Normalize CDF to 1 */
    for (i=0; i<scat->nData; i++) 
      scat->cdfData[i] /= scat->cdfData[scat->nData-1];

    if (scat->oncePerParticle) {
      /* Figure out scattering group assignment */
      if (scat->group>=0) {
        for (i=0; i<dscatterGroups; i++) {
          if (dscatterGroup[i].group==scat->group)
            break;
        }
      } else
        i = dscatterGroups;
      if (i==dscatterGroups) {
        /* make a new group */
        fprintf(stderr, "Making new scattering group structure (#%ld, group %ld) for %s #%ld\n",
                i, scat->group, context.elementName, context.elementOccurrence);
        if (!(dscatterGroup = SDDS_Realloc(dscatterGroup, sizeof(*dscatterGroup)*(dscatterGroups+1)))) {
          fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
          exitElegant(1);
        }
        dscatterGroup[i].particleIDScattered = NULL;
        dscatterGroup[i].nParticles = 0;
        dscatterGroup[i].group = scat->group;
        scat->firstInGroup = 1;
        dscatterGroups++;
      } else {
        fprintf(stderr, "Linking to scattering group structure (#%ld, group %ld) for %s #%ld\n",
                i, scat->group, context.elementName, context.elementOccurrence);
        scat->firstInGroup = 0;
      }
      scat->groupIndex = i;
    }
    if (scat->startOnPass<0)
      scat->startOnPass = 0;
  }


  if (iPass==0 && scat->oncePerParticle && scat->firstInGroup) {
    /* Initialize data for remembering which particles have been scattered */
    dscatterGroup[scat->groupIndex].nParticles = np;
    if (dscatterGroup[scat->groupIndex].particleIDScattered)
      free(dscatterGroup[scat->groupIndex].particleIDScattered);
    if (!(dscatterGroup[scat->groupIndex].particleIDScattered 
          = malloc(sizeof(*(dscatterGroup[scat->groupIndex].particleIDScattered))*np))) {
      fprintf(stderr, "Error for %s: memory allocation failure.\n", context.elementName);
      exitElegant(1);
    }
    dscatterGroup[scat->groupIndex].nScattered = 0;
    dscatterGroup[scat->groupIndex].allScattered = 0;
  }

  if (iPass<scat->startOnPass) 
    return;
  if (scat->endOnPass>=0 && iPass>scat->endOnPass)
    return;
  if (scat->oncePerParticle && dscatterGroup[scat->groupIndex].allScattered)
    return;

  if (iPass==scat->startOnPass) {
    scat->nLeft = scat->limitTotal>=0 ? scat->limitTotal : np;
    if (scat->nLeft>np)
      scat->nLeft = np;
  }
  nLeftThisPass = np;
  if (scat->limitTotal>=0)
    nLeftThisPass = scat->nLeft;
  if (scat->limitPerPass>=0) {
    if (scat->limitPerPass<nLeftThisPass)
      nLeftThisPass = scat->limitPerPass;
  }
  if (nLeftThisPass==0)
    return;

  if (scat->oncePerParticle && dscatterGroup[scat->groupIndex].nParticles<np) {
    fprintf(stderr, "Error for %s: number of particles is greater than the size of the particle ID array.\n",
            context.elementName);
    exitElegant(1);
  }

  nScattered = 0;
  for (ip=0; ip<np; ip++) {
    if (scat->probability<1 && random_2(1)>scat->probability)
      continue;
    if (nLeftThisPass==0)
      break;
    if (scat->oncePerParticle) {
      short found = 0;
      if (dscatterGroup[scat->groupIndex].nScattered>=dscatterGroup[scat->groupIndex].nParticles) {
        fprintf(stderr, "All particles scattered for group %ld (nscattered=%ld, np0=%ld)\n",
                scat->group, dscatterGroup[scat->groupIndex].nScattered, 
		dscatterGroup[scat->groupIndex].nParticles);
	dscatterGroup[scat->groupIndex].allScattered = 1;
	return ;
      }
      for (i=0; i<dscatterGroup[scat->groupIndex].nScattered; i++) {
        if (dscatterGroup[scat->groupIndex].particleIDScattered[i]==part[ip][6]) {
          found = 1;
          break;
        }
      }
      if (found)
        continue;
      dscatterGroup[scat->groupIndex].particleIDScattered[i] = part[ip][6];
      dscatterGroup[scat->groupIndex].nScattered += 1;
    }
    nScattered++;
    nLeftThisPass--;
    cdf = random_2(1);
    amplitude = scat->factor*interp(scat->indepData, scat->cdfData, scat->nData, cdf, 0, 1, &interpCode);
    if (scat->randomSign)
      amplitude *= random_2(1)>0.5 ? 1 : -1;
    if (!interpCode)
      fprintf(stderr, "Warning: interpolation error for %s.  cdf=%e\n",
              context.elementName, cdf);
    if (scat->iPlane==5) {
      /* momentum scattering */
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      t = part[ip][4]/beta;
      part[ip][5] += amplitude;
      P = (1+part[ip][5])*Po;
      beta = P/sqrt(sqr(P)+1);
      part[ip][4] = t*beta;
    } else
      part[ip][scat->iPlane] += amplitude;
  }
  scat->nLeft -= nScattered;
  /*
    fprintf(stderr, "%ld particles scattered by %s#%ld, group %ld, total=%ld\n", 
    nScattered, context.elementName, context.elementOccurrence, 
    scat->group, dscatterGroup[scat->groupIndex].nScattered);
  */
}

void recordLostParticles(double **lossBuffer, double **coord, long *nLost, long nLeft,long  nToTrack, long pass)
{
  long ip, j, nLost0;
  static FILE *fp = NULL;

/*  if (fp==NULL) {
    char s[1024];
    RUN run;
    getRunSetupContext(&run);
    sprintf(s, "%s.losdeb", run.rootname);
    fp = fopen(s, "w");
    fprintf(fp, "SDDS1\n&column name=s type=double units=m &end\n");
    fprintf(fp, "&column name=x type=double units=m &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }

  nLost0 = *nLost;
*/
  *nLost += (nToTrack - nLeft);

  if (!lossBuffer || !coord)
    return;

  if (nLeft==nToTrack) /* no additional losses occurred */
    return;
  
 /* printf("recording lost particles: nLost=%ld->%ld, nLeft=%ld, nToTrack=%ld\n", 
	 nLost0, *nLost, nLeft, nToTrack);
 */ 
  for (ip=nLeft; ip<nToTrack; ip++) {  /* copy the lost particle coordinates and pass information */
    for (j=0; j<7; j++)
      lossBuffer[ip][j] = coord[ip][j];
    lossBuffer[ip][7] = (double) pass;
  }  
}

void storeMonitorOrbitValues(ELEMENT_LIST *eptr, double **part, long np)
{
  MONI *moni;
  HMON *hmon;
  VMON *vmon;
  MARK *mark;
  char s[1000];
  double centroid[6];
  long i;
  
  if (!np)
    return;
  
  switch (eptr->type) {
  case T_MONI:
    moni = (MONI*)eptr->p_elem;
    if (!moni->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!moni->initialized) {
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[1] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[2] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      moni->coMemoryNumber[3] = rpn_create_mem(s, 0);
      moni->initialized = 1;
    }
    for (i=0; i<4; i++) 
      rpn_store(centroid[i], NULL, moni->coMemoryNumber[i]);
    break;
  case T_MARK:
    mark = (MARK*)eptr->p_elem;
    if (!mark->fitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (mark->co_mem==NULL) {
      mark->co_mem = tmalloc(4*sizeof(*(mark->co_mem)));
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      mark->co_mem[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      mark->co_mem[1] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      mark->co_mem[2] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      mark->co_mem[3] = rpn_create_mem(s, 0);
    }
    for (i=0; i<4; i++) 
      rpn_store(centroid[i], NULL, mark->co_mem[i]);
    break;
  case T_HMON:
    hmon = (HMON*)eptr->p_elem;
    if (!hmon->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!hmon->initialized) {
      sprintf(s, "%s#%ld.xco", eptr->name, eptr->occurence);
      hmon->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.xpco", eptr->name, eptr->occurence);
      hmon->coMemoryNumber[1] = rpn_create_mem(s, 0);
      hmon->initialized = 1;
    }
    for (i=0; i<2; i++) 
      rpn_store(centroid[i], NULL, hmon->coMemoryNumber[i]);
    break;
  case T_VMON:
    vmon = (VMON*)eptr->p_elem;
    if (!vmon->coFitpoint) 
      return;
    compute_centroids(centroid, part, np);
    if (!vmon->initialized) {
      sprintf(s, "%s#%ld.yco", eptr->name, eptr->occurence);
      vmon->coMemoryNumber[0] = rpn_create_mem(s, 0);
      sprintf(s, "%s#%ld.ypco", eptr->name, eptr->occurence);
      vmon->coMemoryNumber[1] = rpn_create_mem(s, 0);
      vmon->initialized = 1;
    }
    for (i=0; i<2; i++) 
      rpn_store(centroid[i+2], NULL, vmon->coMemoryNumber[i]);
    break;
  default:
    return ;
  }
}
 
#if USE_MPI
void scatterParticles(double **coord, long *nToTrack, double **accepted,
                      long n_processors, int myid, balance balanceStatus, 
                      double my_rate, double nParPerElements, double round, 
                      int lostSinceSeqMode, int *distributed, 
                      long *reAllocate, double *P_central)
{
  long work_processors = n_processors-1; 
  int root = 0, i, j;
  int my_nToTrack, nItems, *nToTrackCounts;
  double total_rate, constTime, *rateCounts;
  MPI_Status status;
  
  nToTrackCounts = malloc(sizeof(int) * n_processors);
  rateCounts = malloc(sizeof(double) * n_processors);

  /* The particles will be distributed to slave processors evenly for the first pass */
  if (((balanceStatus==startMode) && (!*distributed))) {
    /* || lostSinceSeqMode || *reAllocate )  we don't do redistribution for load balancing, as it could cause memory problem */ 
    MPI_Bcast(nToTrack, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (myid==0) 
      my_nToTrack = 0;
    else {
      my_nToTrack = *nToTrack/work_processors; 
      if (myid<=(*nToTrack%work_processors)) 
	my_nToTrack++;
    } 
    /* gather the number of particles to be sent to each processor */ 
    MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
    *distributed = 1; 
  }
  else if (balanceStatus == badBalance) {
    double minRate;
    if (myid==0)
      my_rate = 1.0;  /* set it to nonzero */
    MPI_Allreduce(&my_rate, &minRate, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
    if (myid==0)
      my_rate = 0.0;  /* set it back to zero */
    if (minRate==0.0) {  /* redistribute evenly when all particles are lost on any working processor */
      MPI_Bcast(nToTrack, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      if (myid==0) 
	my_nToTrack = 0;
      else {
	my_nToTrack = *nToTrack/work_processors; 
	if (myid<=(*nToTrack%work_processors)) 
	  my_nToTrack++;
      } 
      /* gather the number of particles to be sent to each processor */ 
      MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
    }
    else {
      /* calculating the number of jobs to be sent according to the speed of each processors */
      MPI_Gather(&my_rate, 1, MPI_DOUBLE, rateCounts, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Reduce(&my_rate, &total_rate, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (myid==0) {
	if (fabs(total_rate-0.0)>1e-12)
	  constTime = *nToTrack/total_rate;
	else
	  constTime = 0.0;
	my_nToTrack = 0;
	for (i=1; i<work_processors; i++) {
	  nToTrackCounts[i] =  rateCounts[i]*constTime+round; /* round to the nearest integer */
	  /* avoid sending more particles than the available ones */
	  if ((my_nToTrack+nToTrackCounts[i])>=*nToTrack) { 
	    nToTrackCounts[i] = *nToTrack-my_nToTrack; 
	    my_nToTrack = *nToTrack;
	    break;     
	  }
	  /* count the total number of particles that have been assigned */
	  my_nToTrack = my_nToTrack+nToTrackCounts[i];      
	}
	if (i<work_processors)  /* assign 0 particles to the remaining processors */
	  for (j=i; j<=work_processors; j++)
	    nToTrackCounts[j] = 0;
	else   /* The last processor will be responsible for all of the remaining particles */
	  nToTrackCounts[work_processors] = *nToTrack-my_nToTrack;
#ifdef MPI_DEBUG
	printf("total_rate=%lf, nToTrack=%ld, nToTrackForDistribution=%d\n",
	       total_rate, *nToTrack, my_nToTrack+ nToTrackCounts[work_processors] );
#endif
      }
      /* scatter the number of particles to be sent to each processor */ 
      MPI_Scatter(nToTrackCounts, 1, MPI_INT, &my_nToTrack, 1, MPI_INT, root, MPI_COMM_WORLD);
      *reAllocate = 0; /* set the flag back to 0 after scattering */ 
    }
  }
  else { /* keep the nToTrack unchanged */
    if (myid==0)
      my_nToTrack = 0;
    else       
      my_nToTrack =*nToTrack;
    /* gather the number of particles to be sent to each processor */ 
    MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
  }

  /* scatter particles to all of the slave processors */
  if (myid==0) {
    my_nToTrack = 0;
    for (i=1; i<=work_processors; i++) {
      /* calculate the number of elements that will be sent to each processor */
      nItems = nToTrackCounts[i]*COORDINATES_PER_PARTICLE;
      MPI_Send (&coord[my_nToTrack][0], nItems, MPI_DOUBLE, i, 104, MPI_COMM_WORLD); 
      if (accepted!=NULL)
        MPI_Send (&accepted[my_nToTrack][0], nItems, MPI_DOUBLE, i, 105, MPI_COMM_WORLD); 
      /* count the total number of particles that have been scattered */
      my_nToTrack = my_nToTrack+nToTrackCounts[i];      
    }
    *nToTrack = 0;
  } 
  else {
    MPI_Recv (&coord[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0,
              104, MPI_COMM_WORLD, &status); 
    if (accepted!=NULL)
      MPI_Recv (&accepted[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, 
                MPI_DOUBLE, 0, 105, MPI_COMM_WORLD, &status);
    *nToTrack = my_nToTrack;
  }
  /* broadcast the P_central to all the slave processors */
  MPI_Bcast(P_central, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
#ifdef MPI_DEBUG
  if (MPI_DEBUG) {
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    if (myid!=0)
      /*      if ((balanceStatus == badBalance) || lostSinceSeqMode) */ 
	{
	  MPI_Get_processor_name(processor_name,&namelen);
	  fprintf(stderr, "%d will be computed on %d (%s)\n",my_nToTrack,myid,processor_name);
	}
  }
#endif

  free(nToTrackCounts);
  free(rateCounts);
}

void gatherParticles(double ***coord, long **lostOnPass, long *nToTrack, long *nLost, double ***accepted, long n_processors, int myid, double *round)
{
  long work_processors = n_processors-1;
  int root = 0, i, nItems, displs ;
  int my_nToTrack, my_nLost, *nToTrackCounts, 
    *nLostCounts, current_nLost=0, nToTrack_total, nLost_total; 
 
  MPI_Status status;

  nToTrackCounts = malloc(sizeof(int) * n_processors);
  nLostCounts = malloc(sizeof(int) * n_processors);

  if (myid==0) {
    my_nToTrack = 0;  
    my_nLost = 0;
  }
  else {
    my_nToTrack = *nToTrack;
    my_nLost = *nLost;
  }

/* gather nToTrack and nLost from all of the slave processors to the master processors */ 
  MPI_Gather(&my_nToTrack, 1, MPI_INT, nToTrackCounts, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Gather(&my_nLost, 1, MPI_INT, nLostCounts, 1, MPI_INT, root, MPI_COMM_WORLD);

  MPI_Reduce (&my_nToTrack, &nToTrack_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (&my_nLost, &nLost_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (isMaster) {
    if(*coord==NULL)
      *coord = (double**)resize_czarray_2d((void**)(*coord), sizeof(double), nToTrack_total+nLost_total, 7);
    if(dumpAcceptance && (*accepted==NULL))
      *accepted = (double**)resize_czarray_2d((void**)(*accepted), sizeof(double), nToTrack_total+nLost_total, 7); 
  }

  if (myid==0) {
    for (i=1; i<=work_processors; i++) {
      /* the number of elements that are received from each processor (for root only) */
      nItems = nToTrackCounts[i]*COORDINATES_PER_PARTICLE;
      /* collect information for the left particles */
      MPI_Recv (&(*coord)[my_nToTrack][0], nItems, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status); 

      /* count the total number of particles to track and the total number of lost after the most recent scattering */
      my_nToTrack = my_nToTrack+nToTrackCounts[i];
      current_nLost = current_nLost+nLostCounts[i];      
    }
    *nLost = *nLost+current_nLost;
    *nToTrack = my_nToTrack;
  } 
  else {
    MPI_Send (&(*coord)[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD); 
  }

  /* collect information for the lost particles and gather the accepted array */
  
  MPI_Bcast(&current_nLost, 1, MPI_INT, root, MPI_COMM_WORLD);

    if (myid==0) {
      /* set up the displacement array and the number of elements that are received from each processor */ 
      nLostCounts[0] = 0;
      displs = my_nToTrack;
      my_nToTrack = 0;
      for (i=1; i<=work_processors; i++) {
        /* gather information for lost particles */  
  	displs = displs+nLostCounts[i-1];
        nItems = nLostCounts[i]*COORDINATES_PER_PARTICLE;
        MPI_Recv (&(*coord)[displs][0], nItems, MPI_DOUBLE, i, 102, MPI_COMM_WORLD, &status);
        if (*accepted!=NULL){
          MPI_Recv (&(*accepted)[my_nToTrack][0], nToTrackCounts[i]*COORDINATES_PER_PARTICLE, MPI_DOUBLE, i, 101, MPI_COMM_WORLD, &status); 
          MPI_Recv (&(*accepted)[displs][0], nItems, MPI_DOUBLE, i, 103, MPI_COMM_WORLD, &status);
        my_nToTrack = my_nToTrack+nToTrackCounts[i];
	}
      }
      /* update the round parameter to avoid more particles are distributed 
         than the available particles */
      if ((my_nToTrack/work_processors)<work_processors)
        *round = 0.0; 
    }
    else {
      /* send information for lost particles */
      MPI_Send (&(*coord)[my_nToTrack][0], my_nLost*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 102, MPI_COMM_WORLD);  
      if (*accepted!=NULL) {
        MPI_Send (&(*accepted)[0][0], my_nToTrack*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 101, MPI_COMM_WORLD);  
        MPI_Send (&(*accepted)[my_nToTrack][0], my_nLost*COORDINATES_PER_PARTICLE, MPI_DOUBLE, root, 103, MPI_COMM_WORLD); 
      }      
    }
    MPI_Bcast (round, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  free(nToTrackCounts);
  free(nLostCounts);
}

balance checkBalance (double my_wtime, int myid, long n_processors)
{
  double maxTime, minTime, *time;
  int i, balanceFlag = 1; 
  static int imbalanceCounter = 2; /* counter for the number of continuously imbalanced passes,
                                      the 1st pass is treated specially */

  time = malloc(sizeof(double) * n_processors);

  MPI_Gather (&my_wtime, 1, MPI_DOUBLE, time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  if (myid==0) {
    maxTime = minTime = time[1];
    for (i=2; i<n_processors; i++) {
      if (maxTime<time[i])
        maxTime = time[i];
      if (minTime>time[i])
        minTime = time[i]; 
    } 
    if ((maxTime-minTime)/minTime>0.10) {
      imbalanceCounter++;
      /* the workload balancing is treated as bad when there are 3 continuous imbalanced passes or 
         it spends more than 5 minutes for a pass */
      if ((imbalanceCounter>=3)||(maxTime>300)) {
         balanceFlag = 0;       
         imbalanceCounter = 0;  
      }
    }
    else
      imbalanceCounter = 0;
  }
  MPI_Bcast (&balanceFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);     

#ifdef MPI_DEBUG
  if (myid==0) {
    if ((maxTime-minTime)/minTime>0.10)
      printf("The balance is in bad status. ");
    else
      printf("The balance is in good status. ");
    if (balanceFlag==0)
      printf("We need redistribute the particles. ");
    printf("The difference is %4.2lf percent\n", (maxTime-minTime)/minTime*100);
    fflush(stdout);
  }
#endif

  free(time);

  if (balanceFlag==1) 
    return goodBalance;
  else
    return badBalance;
}

int usefulOperation (ELEMENT_LIST *eptr, unsigned long flags, long i_pass)  
{
  WATCH *watch;
  HISTOGRAM *histogram;

  if (eptr->type==T_WATCH) {
   if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
     watch = (WATCH*)eptr->p_elem;
     if (!watch->disable) { 
       if (!watch->initialized) {
         if (isSlave)   /* Slave processors will not go through the WATCH element, so we set the flag here */
	   watch->initialized = 1;
       }   
       if (i_pass>=watch->start_pass && (i_pass-watch->start_pass)%watch->interval==0 &&
	   (watch->end_pass<0 || i_pass<=watch->end_pass))
	 return 1;
     }
   }
   return 0;
  }
  else if (eptr->type==T_HISTOGRAM) {
    if (!(flags&TEST_PARTICLES) && !(flags&INHIBIT_FILE_OUTPUT)) {
      histogram = (HISTOGRAM*)eptr->p_elem;
      if (!histogram->disable) {
	if (!histogram->initialized) {
	  if (isSlave)   /* Slave processors will not go through the HISTOGRAM element */
	    histogram->initialized = 1;
	  return 1;
	}
	if (i_pass>=histogram->startPass && (i_pass-histogram->startPass)%histogram->interval==0)
          return 1;
      }
    }
    return 0 ;
  }
  return 1; /* This is default for all the other UNIPROCESSOR elements */
}
#endif

#ifdef SORT
  int comp_IDs(const void *coord1, const void *coord2)
  {
    if (((double*) coord1)[6] <((double*) coord2)[6]) 
      return -1;
    else if (((double*) coord1)[6] > ((double*) coord2)[6]) 
      return  1;
    else 
      return 0;
  }
#endif

void transformEmittances(double **coord, long np, double pCentral, EMITTANCEELEMENT *ee)
{
  double emit, emitc, factor, pAverage, eta, etap;
  long i, j;
  BEAM_SUMS sums;

#if SDDS_MPI_IO
  long npTotal;
  MPI_Reduce (&np, &npTotal, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (isMaster && npTotal<10) {
    printf("*** Error: too few particles (<10) for emittance modification\n");
    exitElegant(1);
  }
#endif

  zero_beam_sums(&sums, 1);
  accumulate_beam_sums(&sums, coord, np, pCentral);
  pAverage = pCentral*(1+sums.centroid[5]);
  
  for (i=0; i<2; i++) {
    computeEmitTwissFromSigmaMatrix(&emit, &emitc, NULL, NULL, sums.sigma, 2*i);
    eta = etap = 0;
    if (sums.sigma[5][5]) {
      eta  = sums.sigma[2*i+0][5]/sums.sigma[5][5];
      etap = sums.sigma[2*i+1][5]/sums.sigma[5][5];
    }
    if (ee->emit[i]>=0) {
      /* use geometric emittance */
      if (emitc>0)
        factor = sqrt(ee->emit[i]/emitc);
      else
        continue;
    } else if (ee->emitn[i]>=0) {
      /* use normalized emittance */
      emitc *= pAverage;
      if (emitc>0)
        factor = sqrt(ee->emitn[i]/emitc);
      else
        continue;
    } else {
      /* do nothing */
      continue;
    }
    for (j=0; j<np; j++) {
      coord[j][2*i+0] = factor*(coord[j][2*i+0]-eta*coord[j][5]) + eta*coord[j][5];
      coord[j][2*i+1] = factor*(coord[j][2*i+1]-eta*coord[j][5]) + eta*coord[j][5];
    }
  }
}

void set_up_mhist(MHISTOGRAM *mhist, RUN *run, long occurence)
{
  SDDS_DATASET SDDSin;

  if (mhist->disable)
    return;
  if (mhist->interval<=0)
    bombElegant("interval is non-positive for MHISTOGRAM element", NULL);
  if (!mhist->file1d && !mhist->file2dH && !mhist->file2dV && 
      !mhist->file2dL && !mhist->file4d && !mhist->file6d)
    bombElegant("No output request set to this mhistogram element. Use disable =1", NULL);  

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, mhist->inputBinFile) ||
      SDDS_ReadPage(&SDDSin)<=0 ||
      !(mhist->bins1d = (int32_t*)SDDS_GetColumnInLong(&SDDSin, "Bins_1D")) || 
      !(mhist->bins2d = (int32_t*)SDDS_GetColumnInLong(&SDDSin, "Bins_2D")) ||
      !(mhist->bins4d = (int32_t*)SDDS_GetColumnInLong(&SDDSin, "Bins_4D")) ||
      !(mhist->bins6d = (int32_t*)SDDS_GetColumnInLong(&SDDSin, "Bins_6D"))) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (mhist->file1d) {
    if (mhist->bins1d[0]<=2 && mhist->bins1d[1]<=2 && mhist->bins1d[2]<=2 &&
        mhist->bins1d[3]<=2 && mhist->bins1d[4]<=2 && mhist->bins1d[5]<=2)
      bombElegant("All 1D bin's value less than 2", NULL);  
    mhist->file1d = compose_filename_occurence(mhist->file1d, run->rootname, occurence);
  } 
  if (mhist->file2dH) {
    if (mhist->bins2d[0]<=2 || mhist->bins2d[1]<=2)
      bombElegant("2D x-x' bin's value less than 2", NULL);  
    mhist->file2dH = compose_filename_occurence(mhist->file2dH, run->rootname, occurence);    
  }
  if (mhist->file2dV) {
    if (mhist->bins2d[2]<=2 || mhist->bins2d[3]<=2)
      bombElegant("2D y-y' bin's value less than 2", NULL);  
    mhist->file2dV = compose_filename_occurence(mhist->file2dV, run->rootname, occurence);    
  }
  if (mhist->file2dL) {
    if (mhist->bins2d[4]<=2 || mhist->bins2d[5]<=2)
      bombElegant("2D dt-dp bin's value less than 2", NULL);  
    mhist->file2dL = compose_filename_occurence(mhist->file2dL, run->rootname, occurence);    
  }
  if (mhist->file4d) {
    if (mhist->bins4d[0]<=2 || mhist->bins4d[1]<=2 ||
        mhist->bins4d[2]<=2 || mhist->bins4d[3]<=2)
      bombElegant("4D x-x'-y-y' bin's value less than 2", NULL);  
    mhist->file4d = compose_filename_occurence(mhist->file4d, run->rootname, occurence);    
  }
  if (mhist->file6d) {
    if (mhist->bins6d[0]<=2 || mhist->bins6d[1]<=2 || mhist->bins6d[2]<=2 ||
        mhist->bins6d[3]<=2 || mhist->bins6d[4]<=2 || mhist->bins6d[5]<=2)
      bombElegant("6D x-x'-y-y'-dt-dp bin's value less than 2", NULL);  
    if (mhist->lumped)
      mhist->file6d = compose_filename_occurence(mhist->file6d, run->rootname, 0);    
    else
      mhist->file6d = compose_filename_occurence(mhist->file6d, run->rootname, occurence);    
  }

  mhist->initialized = 1;
  mhist->count = 0;

  return;
}

#define MHISTOGRAM_TABLE_PARAMETERS 12
static SDDS_DEFINITION mhist_table_para[MHISTOGRAM_TABLE_PARAMETERS] = {
  {"Step", "&parameter name=Step, symbol=\"m$be$nc\", description=\"Simulation step\", type=long &end"},
  {"pCentral", "&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double &end"},
  {"Charge", "&parameter name=Charge, units=\"C\", description=\"Beam charge\", type=double &end"},
  {"Particles", "&parameter name=Particles, description=\"Number of particles\",  type=long &end"},
  {"Pass", "&parameter name=Pass, type=long &end"},
  {"PassLength", "&parameter name=PassLength, units=\"m\", type=double &end"},
  {"PassCentralTime", "&parameter name=PassCentralTime, units=\"s\", type=double &end"},
  {"s", "&parameter name=s, units=\"m\", description=\"Location from beginning of beamline\", type=double &end"},
  {"ElementName", "&parameter name=ElementName, description=\"Previous element name\", type=string &end"},
  {"ElementOccurence", "&parameter name=ElementOccurence, description=\"Previous element's occurence\", format_string=%6ld, type=long, &end"},
  {"ElementType", "&parameter name=ElementType, description=\"Previous element-type name\", format_string=%10s, type=string, &end"},
  {"Normalize", "&parameter name=Normalize, description=\"If the table is normalized?\", type=long &end"},
  };
static void *mhist_table_paraValue[MHISTOGRAM_TABLE_PARAMETERS];

void mhist_table(ELEMENT_LIST *eptr0, ELEMENT_LIST *eptr, long step, long pass, double **coord, long np, 
                double Po, double length, double charge, double z)
{
  char *Name[6]={"x","xp","y","yp","dt","delta"};
  char *Units[6]={"m","","m","","s",""};
  double Min[6], Max[6], c0[6], t0;
  long i, j;
  double part[6], P, beta;
  MHISTOGRAM *mhist;

  if (eptr0)
    mhist = (MHISTOGRAM*)eptr0->p_elem;
  else
    mhist = (MHISTOGRAM*)eptr->p_elem;

  t0 = pass*length*sqrt(Po*Po+1)/(c_mks*(Po+1e-32));  
  mhist_table_paraValue[0] = (void*)(&step);
  mhist_table_paraValue[1] = (void*)(&Po);
  mhist_table_paraValue[2] = (void*)(&charge);
  mhist_table_paraValue[3] = (void*)(&np);
  mhist_table_paraValue[4] = (void*)(&pass);
  mhist_table_paraValue[5] = (void*)(&length);
  mhist_table_paraValue[6] = (void*)(&t0);
  mhist_table_paraValue[7] = (void*)(&z);
  mhist_table_paraValue[8] = (void*)(&eptr->pred->name);
  mhist_table_paraValue[9] = (void*)(&eptr->pred->occurence);
  mhist_table_paraValue[10] = (void*)(&entity_name[eptr->pred->type]);
  mhist_table_paraValue[11] = (void*)(&mhist->normalize);


  findMinMax (coord, np, Min, Max, c0, Po);
  Min[4] -= c0[4];
  Max[4] -= c0[4];

  if (mhist->file1d)
    mhist->x1d = chbook1m(Name, Units, Min, Max, mhist->bins1d, 6);

  if (mhist->file2dH)
    mhist->x2d = chbookn(Name, Units, 2, Min, Max, mhist->bins2d, 0);
  if (mhist->file2dV)
    mhist->y2d = chbookn(Name, Units, 2, Min, Max, mhist->bins2d, 2);
  if (mhist->file2dL)
    mhist->z2d = chbookn(Name, Units, 2, Min, Max, mhist->bins2d, 4);
  if (mhist->file4d)
    mhist->Tr4d = chbookn(Name, Units, 4, Min, Max, mhist->bins4d, 0);
  if (mhist->file6d)
    mhist->full6d = chbookn(Name, Units, 6, Min, Max, mhist->bins6d, 0);

  for (i=0; i<np; i++) {
    for (j=0; j<6; j++) {
      if (j==4) {
        P = Po*(1+coord[i][5]);
        beta = P/sqrt(P*P+1);
        part[4] = coord[i][4]/(beta*c_mks)-c0[4];
      }
      else
        part[j] = coord[i][j];
    }

    if (mhist->x1d)
      chfill1m(mhist->x1d, part, 1, mhist->bins1d, 6);
    if (mhist->x2d)
      chfilln(mhist->x2d, part, 1, 0);
    if (mhist->y2d)
      chfilln(mhist->y2d, part, 1, 2);
    if (mhist->z2d)
      chfilln(mhist->z2d, part, 1, 4);
    if (mhist->Tr4d)
      chfilln(mhist->Tr4d, part, 1, 0);
    if (mhist->full6d)
      chfilln(mhist->full6d, part, 1, 0);
  }

  if (mhist->x1d) {
    chprint1m(mhist->x1d, mhist->file1d, "One dimentional distribution", mhist_table_para, 
              mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbook1m(mhist->x1d);

  }
  if (mhist->x2d) {
    chprintn(mhist->x2d, mhist->file2dH, "x-xp distribution", mhist_table_para,
             mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbookn(mhist->x2d);
  }
  if (mhist->y2d) {
    chprintn(mhist->y2d, mhist->file2dV, "y-yp distribution", mhist_table_para,
             mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbookn(mhist->y2d);
  }
  if (mhist->z2d) {
    chprintn(mhist->z2d, mhist->file2dL, "dt-delta distribution", mhist_table_para,
             mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbookn(mhist->z2d);
  }
  if (mhist->Tr4d) {
    chprintn(mhist->Tr4d, mhist->file4d, "x-xp-y-yp distribution", mhist_table_para,
             mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbookn(mhist->Tr4d);
  }
  if (mhist->full6d) {
    chprintn(mhist->full6d, mhist->file6d, "full distribution", mhist_table_para,
             mhist_table_paraValue, MHISTOGRAM_TABLE_PARAMETERS, mhist->normalize, 0, mhist->count);
    free_hbookn(mhist->full6d);
  }

  mhist->count ++;

  return;
}

void findMinMax (double **coord, long np, double *min, double *max, double *c0, double Po)
{
  long i, j;
  double P, beta, time;

  for (j=0; j<6; j++) {
    max[j] = -(min[j] = DBL_MAX);
    c0[j] = 0;
    for (i=0; i<np; i++) {
      if (j==4) {
        P = Po*(1+coord[i][5]);
        beta = P/sqrt(P*P+1);
        time = coord[i][4]/(beta*c_mks);
        c0[j] += time;
        if (min[j] > time)
          min[j] = time;
        if (max[j] < time)
          max[j] = time;
      }
      else {
        c0[j] += coord[i][j];
        if (min[j] > coord[i][j])
          min[j] = coord[i][j];
        if (max[j] < coord[i][j])
          max[j] = coord[i][j];
      }
    }
    c0[j] /= (double)np; 
  }
  return;
}

void field_table_tracking(double **particle, long np, FTABLE *ftable, double Po, RUN *run)
{
  long ip, ik, nKicks, debug = ftable->verbose;
  double *coord, p0, factor;
  double xyz[3], p[3], B[3], BA, pA, **A, pz0;
  double rho, theta0, theta1, theta2, theta, tm_a, tm_b, tm_c;
  double step, eomc, s_location;
  char *rootname;

  static SDDS_TABLE test_output;
  static long first_time = 1;
  
  if (first_time && debug) {
    rootname = compose_filename("%s.phase", run->rootname);
    SDDS_PhaseSpaceSetup(&test_output, rootname, SDDS_BINARY, 1, "output phase space", run->runfile,
                         run->lattice, "setup_output");
    first_time =0 ;
  }
  
  if ((nKicks=ftable->nKicks)<1)
    bombElegant("N_KICKS must be >=1 for FTABLE", NULL);
  if (!ftable->initialized)
    bombElegant("Initialize FTABLE: This shouldn't happen.", NULL);
  step = ftable->length/nKicks;

  /* convert coordinate frame from local to ftable element frame. Before misalignment.*/
  if (debug) 
    dump_phase_space(&test_output, particle, np, -2, Po, 0);
  if (ftable->e1 || ftable->l1)
    ftable_frame_converter(particle, np, ftable, 0);
  if (debug) 
    dump_phase_space(&test_output, particle, np, -1, Po, 0);
  
  if (ftable->dx || ftable->dy || ftable->dz)
    offsetBeamCoordinates(particle, np, ftable->dx, ftable->dy, ftable->dz);
  if (ftable->tilt)
    rotateBeamCoordinates(particle, np, ftable->tilt);
  if (debug) 
      dump_phase_space(&test_output, particle, np, 0, Po, 0);

  s_location =step/2.;
  eomc = -particleCharge/particleMass/c_mks;
  A = (double**)czarray_2d(sizeof(double), 3, 3);
  if (debug) 
        fprintf(stdout, "ik         x       y         z   Bx  By Bz \n");
  for (ik=0; ik<nKicks; ik++) {
    for (ip=0; ip<np; ip++) {
      /* 1. get particle's coordinates */
      coord = particle[ip];
      factor = sqrt(1+sqr(coord[1])+sqr(coord[3]));
      p0 = (1.+coord[5])*Po;
      p[2] = p0/factor;
      p[0] = coord[1]*p[2];
      p[1] = coord[3]*p[2];

      /* 2. get field at the middle point */
      xyz[0] = coord[0] + coord[1]*step/2.0;
      xyz[1] = coord[2] + coord[3]*step/2.0;
      xyz[2] = s_location; 
      interpolateFTable(B, xyz, ftable);
/*      if (debug) 
        fprintf(stdout, "%5d \t %5f \t %5f \t %5f \t %10f \t %10f \t %10f \n", ik, xyz[0], xyz[1], xyz[2], B[0], B[1], B[2]);
*/      
      BA = sqrt(sqr(B[0]) + sqr(B[1]) + sqr(B[2]));
      /* 3. calculate the rotation matrix */
      A[0][0] = -(p[1]*B[2] - p[2]*B[1]);
      A[0][1] = -(p[2]*B[0] - p[0]*B[2]);
      A[0][2] = -(p[0]*B[1] - p[1]*B[0]);
      pA = sqrt(sqr(A[0][0]) + sqr(A[0][1]) + sqr(A[0][2]));
      /* When field not equal to zero or not parallel to the particles motion */
      if (BA && pA) {
        A[0][0] /= pA;
        A[0][1] /= pA;
        A[0][2] /= pA;
        A[1][0] = B[0]/BA;
        A[1][1] = B[1]/BA;
        A[1][2] = B[2]/BA;
        A[2][0] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
        A[2][1] = A[0][2]*A[1][0]-A[0][0]*A[1][2];
        A[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];
        
        /* 4. rotate coordinates from (x,y,z) to (u,v,w) with u point to BxP, v point to B */
        pz0 = p[2];
        rotate_coordinate(A, p, 0);
        if (p[2] < 0)
          bombElegant("Table function doesn't support particle going backward", NULL);
        rotate_coordinate(A, B, 0);

        /* 5. apply kick */
        rho = p[2]/(eomc*B[1]);
        theta0=theta1=theta2=0.;
        if (A[2][2]) {
          tm_a =  3.0*A[0][2]/A[2][2];
          tm_b = -6.0*A[1][2]*p[1]/p[2]/A[2][2]-6.0;
          tm_c =  6.0*step/rho/A[2][2];
#ifdef USE_GSL
          gsl_poly_solve_cubic (tm_a, tm_b, tm_c, &theta0, &theta1, &theta2);
#else
          bombElegant("gsl_poly_solve_cubic function is not available becuase this version of elegant was not built against the gsl library", NULL);
#endif
        } else if (A[0][2]) {
          tm_a = A[1][2]*p[1]/p[2]+A[2][2];
          theta0 = (tm_a-sqrt(sqr(tm_a)-2.*A[0][2]*step/rho))/A[0][2];
          theta1 = (tm_a+sqrt(sqr(tm_a)-2.*A[0][2]*step/rho))/A[0][2];
        } else {
          tm_a = A[1][2]*p[1]/p[2]+A[2][2];          
          theta0 = step/rho/tm_a;
        }
        theta=choose_theta(rho, theta0, theta1, theta2);
       if (debug) 
        fprintf(stdout, "%5d \t %5f \t %5f \t %5f \t %10f \t %10f \t %10f \t %20f \t %10f \t %10f \n", ik, xyz[0], xyz[1], xyz[2], B[0], B[1], B[2], p[2], rho, theta);

        p[0] = -p[2]*sin(theta);
        p[2] *= cos(theta);
        xyz[0] = rho*(cos(theta)-1);
        xyz[1] = (p[1]/p[2])*rho*theta;
        xyz[2] = rho*sin(theta);

        /* 6. rotate back to (x,y,z) */
        rotate_coordinate(A, xyz, 1);
        rotate_coordinate(A, p, 1);
        coord[0] += xyz[0];
        coord[2] += xyz[1];
        coord[4] += sqrt(sqr(rho*theta)+sqr(xyz[1]));
        coord[1] = p[0]/p[2];
        coord[3] = p[1]/p[2];
      } else {
        coord[0] += coord[1]*step;
        coord[2] += coord[3]*step;
        coord[4] += step*factor;         
      }
    }
    s_location += step;
  if (debug) 
    dump_phase_space(&test_output, particle, np, ik+1, Po, 0);
  }

  free_czarray_2d((void**)A,3,3);

  if (ftable->tilt)
    rotateBeamCoordinates(particle, np, -ftable->tilt);
  if (ftable->dx || ftable->dy || ftable->dz)
    offsetBeamCoordinates(particle, np, -ftable->dx, -ftable->dy, -ftable->dz);
  if (debug) 
    dump_phase_space(&test_output, particle, np, nKicks+1, Po, 0);

  /* convert coordinate frame from ftable element frame to local frame. After misalignment.*/
  if (ftable->e2 || ftable->l2)
    ftable_frame_converter(particle, np, ftable, 1);
  if (debug) 
    dump_phase_space(&test_output, particle, np, nKicks+2, Po, 0);

  return;  
}

/* 0: entrance; 1: exit */
void ftable_frame_converter(double **particle, long np, FTABLE *ftable, long entrance_exit)
{
  double *coord, x0, xp0, y0, yp0;
  double s, c, temp, dx, length;
  long ip, entrance=0, exit=1;
  
  if (entrance_exit==entrance) {
    if (!ftable->e1) {
      exactDrift(particle, np, -ftable->l1);
      return;
    }
    
    /* rotate to ftable element frame */
    s = sin(ftable->e1);
    c = cos(ftable->e1);
    for (ip=0; ip<np; ip++) {
      coord = particle[ip];
      x0  = coord[0];
      xp0 = coord[1];
      y0  = coord[2];
      yp0 = coord[3];
      length = ftable->l1-s*x0;
      temp = c-s*xp0;
      if (temp==0) 
        bombElegant("ftable_frame_converter: Particle will never get into Element.", NULL);

      coord[1] = (c*xp0+s)/temp;
      coord[3] = yp0/temp;
      coord[4] -= length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
      coord[0] = c*x0-coord[1]*length;
      coord[2] = y0-coord[3]*length;
    }
  }
  
  if (entrance_exit==exit) {
    if (!ftable->e2) {
      exactDrift(particle, np, -ftable->l2);
      return;
    }

    /* rotate to normal local frame */
    s = sin(ftable->e2);
    c = cos(ftable->e2);
    dx = ftable->l2*s;
    for (ip=0; ip<np; ip++) {
      coord = particle[ip];
      x0  = coord[0];
      xp0 = coord[1];
      y0  = coord[2];
      yp0 = coord[3];
      length = c*ftable->l2-s*x0;
      temp = c-s*xp0;
      if (temp==0) 
        bombElegant("ftable_frame_converter: Particle will never get into Element.", NULL);

      coord[1] = (c*xp0+s)/temp;
      coord[3] = yp0/temp;
      coord[4] -= length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
      coord[0] = c*x0+dx-coord[1]*length;
      coord[2] = y0-coord[3]*length;
    }
  }
  
  return;
}

void interpolateFTable(double *B, double *xyz, FTABLE *ftable) 
{
  double dummy[3];
/*
  B[0] = interpolate_bookn(ftable->Bx, dummy, xyz, 0, 0, 0, 1, ftable->verbose);
  B[1] = interpolate_bookn(ftable->By, dummy, xyz, 0, 0, 0, 1, ftable->verbose);
  B[2] = interpolate_bookn(ftable->Bz, dummy, xyz, 0, 0, 0, 1, ftable->verbose);
*/
  B[0] = interpolate_bookn(ftable->Bx, dummy, xyz, 0, 0, 0, 1, 0);
  B[1] = interpolate_bookn(ftable->By, dummy, xyz, 0, 0, 0, 1, 0);
  B[2] = interpolate_bookn(ftable->Bz, dummy, xyz, 0, 0, 0, 1, 0);
  
  return;
}

void rotate_coordinate(double **A, double *x, long inverse) {
  long i, j;
  double temp[3];

  for (i=0; i<3; i++) {
    temp[i] = 0;
    for (j=0; j<3; j++) {
      if (!inverse)
        temp[i] += A[i][j]*x[j];
      else
        temp[i] += A[j][i]*x[j];
    }
  }

  for (i=0; i<3; i++)
    x[i] =  temp[i];

  return;
}

/* choose suitable value from cubic solver */
double choose_theta(double rho, double x0, double x1, double x2)
{
  double temp;

  if (rho<0) {
    temp = -DBL_MAX;
    if (x0<0 && x0>temp) temp = x0;
    if (x1<0 && x1>temp) temp = x1;
    if (x2<0 && x2>temp) temp = x2;
  }
  if (rho>0) {
    temp = DBL_MAX;
    if (x0>0 && x0<temp) temp = x0;
    if (x1>0 && x1<temp) temp = x1;
    if (x2>0 && x2<temp) temp = x2;
  }
  return temp;
}

