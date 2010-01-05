/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: bunched_beam.c
 * purpose: Do tracking for bunched gaussian beams.
 *          See file bunched_beam.nl for input parameters.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "bunched_beam.h"

char *beam_type[N_BEAM_TYPES] = {
    "gaussian", "hard-edge", "uniform-ellipse", "shell",
    "dynamic-aperture", "line", "halo(gaussian)",
    };

static TRANSVERSE x_plane, y_plane;
static LONGITUDINAL longit;
static long x_beam_type, y_beam_type, longit_beam_type;

static SDDS_TABLE SDDS_bunch;
static long SDDS_bunch_initialized;
static long beamRepeatSeed, bunchGenerated, firstIsFiducial, beamCounter;
static long haltonID[6];

void fill_transverse_structure(TRANSVERSE *trans, double emit, double beta, 
    double alpha, double eta, double etap, long xbeam_type, double cutoff,
    double *xcentroid);
void fill_longitudinal_structure(LONGITUDINAL *xlongit, double xsigma_dp,
                                 double xsigma_s, double xcoupling_angle, 
                                 double xemit, double xbeta, double xalpha, double xchirp,
                                 long xbeam_type, double xcutoff,
                                 double *xcentroid);
void makeBucketAssignments(BEAM *beam, double Po, double frequency);
int comp_BucketNumbers(const void *coord1, const void *coord2);

void setup_bunched_beam(
    BEAM *beam,
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    OPTIM_VARIABLES *optim,
    OUTPUT_FILES *output,
    LINE_LIST *beamline,
    long n_elements,
    long save_original
    )
{
  long i, offset;
  log_entry("setup_bunched_beam");

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&bunched_beam, nltext);
  if (Po<=0)
    Po = run->p_central;
  if (use_twiss_command_values) {
    double xTwiss[5], yTwiss[5];
    long mode;
    
    if (!get_twiss_mode(&mode, xTwiss, yTwiss) && mode!=0)
      bomb("use_twiss_command_values invalid unless twiss_output command given previously with matched=0.",
           NULL);
    beta_x = xTwiss[0];
    alpha_x = xTwiss[1];
    eta_x = xTwiss[2];
    etap_x = xTwiss[3];
    beta_y = yTwiss[0];
    alpha_y = yTwiss[1];
    eta_y = yTwiss[2];
    etap_y = yTwiss[3];
  }
  
  if (echoNamelists) print_namelist(stdout, &bunched_beam);

  /* check for validity of namelist inputs */
  if (emit_nx && !emit_x)
    emit_x = emit_nx/Po;
  if (emit_ny && !emit_y)
    emit_y = emit_ny/Po;
  if (n_particles_per_bunch<=0)
    bomb("n_particles_per_bunch is invalid", NULL);
#if SDDS_MPI_IO
  beam->n_original_total = beam->n_to_track_total = n_particles_per_bunch; /* record the total number of particles being tracked */
  if (n_particles_per_bunch == 1) /* a special case for single particle tracking */
    notSinglePart = 0;
  else {
    if (isSlave) {
      long work_processors = n_processors-1;
      long my_nToTrack = n_particles_per_bunch/work_processors; 
      if (myid<=(n_particles_per_bunch%work_processors)) 
	my_nToTrack++;
      n_particles_per_bunch = my_nToTrack;
    }
    else {
      n_particles_per_bunch = 0;
    }
  }
#endif
  if (emit_x<0 || beta_x<=0)
    bomb("emit_x<=0 or beta_x<=0", NULL);
  if (emit_y<0 || beta_y<=0)
    bomb("emit_y<=0 or beta_y<=0", NULL);
  if (sigma_dp<0 || sigma_s<0) 
    bomb("sigma_dp<0 or sigma_s<0", NULL);
  if (fabs(dp_s_coupling)>1)
    bomb("|dp_s_coupling| > 1", NULL);
  if (emit_z<0)
    bomb("emit_z < 0", NULL);
  if (beta_z<0)
    bomb("beta_z < 0", NULL);
  if (emit_z && (sigma_dp || sigma_s)) 
    bomb("give emit_z or both sigma_dp and sigma_s", NULL);
  if (emit_z && beta_z<=0)
    bomb("give beta_z with emit_z", NULL);
  if ( ((alpha_z?1:0)+(dp_s_coupling?1:0)+(momentum_chirp?1:0))>1)
    bomb("give only one of alpha_z, dp_s_coupling, or momentum_chirp", NULL);
  if (emit_z && dp_s_coupling)
    bomb("give alpha_z not dp_s_coupling with emit_z", NULL);
  
  if (!distribution_type[0] || 
      (x_beam_type=match_string(distribution_type[0], beam_type, 
                                N_BEAM_TYPES, 0))<0)
    bomb("missing or unknown x distribution type", NULL);
  if (!distribution_type[1] || 
      (y_beam_type=match_string(distribution_type[1], beam_type,
                                N_BEAM_TYPES, 0))<0)
    bomb("missing or unknown y distribution type", NULL);
  if (!distribution_type[2] ||
      (longit_beam_type=match_string(distribution_type[2], beam_type,
                                     N_BEAM_TYPES, 0))<0)
    bomb("missing or unknown longitudinal distribution type", NULL);
  if (distribution_cutoff[0]<0)
    bomb("x distribution cutoff < 0", NULL);
  if (distribution_cutoff[1]<0)
    bomb("y distribution cutoff < 0", NULL);
  if (distribution_cutoff[2]<0)
    bomb("longitudinal distribution cutoff < 0", NULL);
  if ((enforce_rms_values[0] || enforce_rms_values[1] || enforce_rms_values[2]) && n_particles_per_bunch%4)
    fputs("Note: you will get better results for enforcing RMS values if n_particles_per_bunch is divisible by 4.\n", stdout);

  if (matched_to_cell) {
    unsigned long flags;
    flags = beamline->flags;
    if (!(control->cell = get_beamline(NULL, matched_to_cell, run->p_central, 0)))
      bomb("unable to find beamline definition for cell", NULL);
    if (control->cell==beamline)
      beamline->flags = flags;
  }
  else
    control->cell = NULL;

  if (bunch) {
    /* prepare dump of input particles */
    bunch = compose_filename(bunch, run->rootname);
    SDDS_PhaseSpaceSetup(&SDDS_bunch, bunch, SDDS_BINARY, 1, "bunched-beam phase space", run->runfile,
                         run->lattice, "setup_bunched_beam");
    SDDS_bunch_initialized = 1;
  }

  /* set up structures for creating initial particle distributions */ 
  if (!control->cell) {
    fill_transverse_structure(&x_plane, emit_x, beta_x, alpha_x, eta_x, etap_x,
                              x_beam_type, distribution_cutoff[0], centroid);
    fill_transverse_structure(&y_plane, emit_y, beta_y, alpha_y, eta_y, etap_y,
                              y_beam_type, distribution_cutoff[1], centroid+2);
  }
  fill_longitudinal_structure(&longit, 
                              sigma_dp, sigma_s, dp_s_coupling,
                              emit_z, beta_z, alpha_z, momentum_chirp,
                              longit_beam_type, distribution_cutoff[2], centroid+4);

  save_initial_coordinates = save_original || save_initial_coordinates;
  if (n_particles_per_bunch==1)
    one_random_bunch = 1;
  if (!one_random_bunch)
    save_initial_coordinates = save_original;
  if (isSlave || !notSinglePart)
    beam->particle = (double**)czarray_2d(sizeof(double), n_particles_per_bunch, 7);
#if SDDS_MPI_IO
  else
    beam->particle = NULL;
#endif
  if (isSlave)
    if (!save_initial_coordinates)
      beam->original = beam->particle;
    else
      beam->original = (double**)czarray_2d(sizeof(double), n_particles_per_bunch, 7);
#if SDDS_MPI_IO
  else
    beam->original = NULL;
#endif
  if(isSlave)
    if (run->acceptance) 
      beam->accepted = (double**)czarray_2d(sizeof(double), n_particles_per_bunch, 7);
    else
      beam->accepted = NULL;
#if SDDS_MPI_IO
  else
    beam->accepted = NULL;
  if (run->acceptance)
    dumpAcceptance = 1; /* indicate if the initial coordinates of transmitted particles will be dumped */
#endif
  beam->n_original = beam->n_to_track = beam->n_particle = n_particles_per_bunch;
  beam->n_accepted = beam->n_saved = 0;
  firstIsFiducial = first_is_fiducial;
  beamCounter = 0;
  
  if (one_random_bunch) {
    /* make a seed for reinitializing the beam RN generator */
    beamRepeatSeed = 1e8*random_4(1);
#if SDDS_MPI_IO
    /* All processors will have same beamRepeatSeed after here. This will make
       it easy for the serial version to generate the same sequence */
    MPI_Bcast(&beamRepeatSeed, 1, MPI_LONG, 1, MPI_COMM_WORLD);
#endif
    random_4(-beamRepeatSeed);
  }
  bunchGenerated = 0;
  
  for (i=0; i<3; i++) {
    for (offset=0; offset<2; offset++) {
      haltonID[2*i+offset] = 0;
      if (halton_sequence[i]) {
        if (halton_radix[2*i+offset]<0)
          bomb("Radix for Halton sequence can't be negative", NULL);
        if ((haltonID[2*i+offset] = startHaltonSequence(&halton_radix[2*i+offset], 0.5))<0) 
          bomb("problem starting Halton sequence", NULL);
        fprintf(stdout, "Using radix %" PRId32 " for Halton sequence for coordinate %ld\n",
                halton_radix[2*i+offset], 2*i+offset);
      }
    }
  }
  fflush(stdout);
  
  log_exit("setup_bunched_beam");
}

long new_bunched_beam(
    BEAM *beam,
    RUN *run,
    VARY *control,
    OUTPUT_FILES *output,
    long flags
    )
{    
    static long n_actual_particles;
    long i_particle, i_coord;
    unsigned long unstable;
    double s_offset, beta;
    double p_central, dummy, p, gamma;
#ifdef VAX_VMS
    char s[100];
#endif
#if USE_MPI
    double save_emit_x, save_emit_y,
           save_sigma_dp, save_sigma_s;

    if (firstIsFiducial || do_find_aperture || (beam->n_original_total==1)) {
      notSinglePart = 0;
      lessPartAllowed = 1;
    }
    else {
      notSinglePart = 1;
      lessPartAllowed = 0;
    }
#endif

    beamCounter++;
    if (firstIsFiducial && beamCounter==1) {
      beam->n_original = beam->n_to_track = 1;
#if USE_MPI
      lessPartAllowed = 1;  /* All the processors will do the same thing from now */
      if (isMaster) { /* For parallel I/O version, memory will be allocated for one particle on master */
	beam->particle=(double**)czarray_2d(sizeof(double),1,7);
	if (!save_initial_coordinates)
	  beam->original = beam->particle;
	else
	  beam->original=(double**)czarray_2d(sizeof(double),1,7);
	if (run->acceptance)
	  beam->accepted=(double**)czarray_2d(sizeof(double),1,7);
      }	
#endif
    }
    else {
      if (firstIsFiducial && beamCounter==2)
        bunchGenerated = 0;
      beam->n_original = beam->n_to_track = beam->n_particle = n_particles_per_bunch;
#if USE_MPI
/*    For DA search, we need set lessPartAllowed=1 (singlePart=1), even the particle number is larger
      than 1. Need to check if deleting these lines will affect other tests. 9/2/2009 Y. Wang
      if (beam->n_original_total != 1)
	lessPartAllowed = 0;
      else {
	lessPartAllowed = 1; 
      }
*/
      if (isMaster) {
	if (notSinglePart) {
	  beam->particle = beam->original = beam->accepted = NULL;
	  beam->n_original = beam->n_to_track = beam->n_particle = 0;
	} else {
	  beam->particle=(double**)czarray_2d(sizeof(double),n_particles_per_bunch,7);
	  if (!save_initial_coordinates)
	    beam->original = beam->particle;
	  else
	    beam->original=(double**)czarray_2d(sizeof(double),n_particles_per_bunch,7);
	  if (run->acceptance)
	    beam->accepted=(double**)czarray_2d(sizeof(double),n_particles_per_bunch,7);
	}
      }
#endif
    }

    beam->n_accepted = 0;

    p_central = beam->p0_original = run->p_central;

    if (flags&TRACK_PREVIOUS_BUNCH && beam->original==beam->particle)
      bomb("programming error: original beam coordinates needed but not saved", NULL);
    
    if (!bunchGenerated || !save_initial_coordinates) {
      if (one_random_bunch) {
            /* reseed the random number generator to get the same sequence again */
            /* means we don't have to store the original coordinates */
#if SDDS_MPI_IO 
	/* There will be no same sequence for different processors */
	if (notSinglePart)
	  random_4(-(beamRepeatSeed+2*(myid-1)));  
	else
	  random_4(-beamRepeatSeed);
#else
	random_4(-beamRepeatSeed);
#endif
      }
        if (control->cell) {
            VMATRIX *M;
	    unsigned long savedFlags;
	    savedFlags = control->cell->flags;
            M = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, &dummy,
                                   &beta_y, &alpha_y, &eta_y, &etap_y, &dummy, 
                                   (control->cell->elem_recirc?control->cell->elem_recirc:&(control->cell->elem)), 
                                   NULL, run, &unstable, NULL, NULL);
	    control->cell->flags = savedFlags;
            free_matrices(M);
            free(M);
	    fprintf(stdout, "matched Twiss parameters for beam generation:\nbetax = %13.6e m  alphax = %13.6e  etax = %13.6e m  etax' = %13.6e\n",
                beta_x, alpha_x, eta_x, etap_x);
     fflush(stdout);
            fprintf(stdout, "betay = %13.6e m  alphay = %13.6e  etay = %13.6e m  etay' = %13.6e\n",
                beta_y, alpha_y, eta_y, etap_y);
            fflush(stdout);
            fill_transverse_structure(&x_plane, emit_x, beta_x, alpha_x, eta_x, etap_x,
                    x_beam_type, distribution_cutoff[0], centroid);
            fill_transverse_structure(&y_plane, emit_y, beta_y, alpha_y, eta_y, etap_y,
                    y_beam_type, distribution_cutoff[1], centroid+2);
            }
        fprintf(stdout, "generating bunch %ld.%ld\n", control->i_step, control->i_vary);
        fflush(stdout);
#if USE_MPI
        if (firstIsFiducial && beamCounter==1) {
          /* set emittances to zero temporarily */
          save_emit_x = x_plane.emit;
          save_emit_y = y_plane.emit;
          save_sigma_dp = longit.sigma_dp;
          save_sigma_s = longit.sigma_s;
          x_plane.emit = y_plane.emit = longit.sigma_dp = longit.sigma_s = 0;
        }
#endif
#if !SDDS_MPI_IO
	if (run->random_sequence_No > 1 && one_random_bunch) {
	  /* generate the same random number sequence as the parallel version */
	  long work_processors = run->random_sequence_No-1, my_n_actual_particles;
	  long my_nToTrack, i;

	  orig_sequence_No = remaining_sequence_No = run->random_sequence_No-1;
	  n_actual_particles = 0;
	  for (i=0; i<work_processors; i++) {
	    /* This will make the serial version has the same start seed as each of the sequences 
	       in the parallel version */
	    random_4(-(beamRepeatSeed+2*i)); 
	    my_nToTrack = beam->n_to_track/work_processors;
	    if (i<(beam->n_to_track%work_processors)) 
	      my_nToTrack++;
	    if (x_plane.beam_type==DYNAP_BEAM) { 
	      /* This is a special case, generate_bunch will be called one time only */ 
	      remaining_sequence_No = 1;
	      my_nToTrack = beam->n_to_track;
	    }
	    my_n_actual_particles =
	      generate_bunch(&beam->original[n_actual_particles], my_nToTrack, &x_plane,
			     &y_plane, &longit, enforce_rms_values, limit_invariants, 
			     symmetrize, haltonID, randomize_order, limit_in_4d, Po);
	    n_actual_particles += my_n_actual_particles;
	    if (x_plane.beam_type==DYNAP_BEAM)
	      break;
	  }
	}
	else
#endif
	  n_actual_particles = 
	    generate_bunch(beam->original, beam->n_to_track, &x_plane,
			   &y_plane, &longit, enforce_rms_values, limit_invariants, 
			   symmetrize, haltonID, randomize_order, limit_in_4d, Po);
#if USE_MPI
       if (firstIsFiducial && beamCounter==1) {
          /* copy values back */
          x_plane.emit = save_emit_x;
          y_plane.emit = save_emit_y;
          longit.sigma_dp = save_sigma_dp;
          longit.sigma_s = save_sigma_s;
        }
#endif
        if (bunch) {
            if (!SDDS_bunch_initialized)
                bomb("'bunch' file is uninitialized (new_bunched_beam)", NULL);
            fprintf(stdout, "dumping bunch\n");
            fflush(stdout);
            dump_phase_space(&SDDS_bunch, beam->original, n_actual_particles, control->i_step, Po,
                             sqrt(-1.0));
            if (one_random_bunch) {
                if (!SDDS_Terminate(&SDDS_bunch)) {
                    SDDS_SetError("Problem terminating 'bunch' file (new_bunched_beam)");
                    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
                    }
                bunch = NULL;
                }
            }
        }

    beam->n_to_track = n_actual_particles;

    /* copy particles into tracking buffer, adding in the time offset 
     * and converting from deltap/p relative to bunch center (defined
     * by Po) to deltap/p relative to beamline energy (defined by
     * p_central).  Also add the particle ID.
     */
    beta = Po/sqrt(sqr(Po)+1);
    s_offset = control->bunch_frequency?
                (control->i_step-1)*(beta*c_mks)/control->bunch_frequency:0;
    for (i_particle=0; i_particle<beam->n_to_track; i_particle++) {
        if (beam->particle!=beam->original)
          for (i_coord=0; i_coord<7; i_coord++) 
            beam->particle[i_particle][i_coord] = beam->original[i_particle][i_coord];
        gamma = sqrt(sqr(p=Po*(1+beam->particle[i_particle][5]))+1);
        beta  = p/gamma;
        beam->particle[i_particle][5] = (p-p_central)/p_central;
        beam->particle[i_particle][4] += s_offset;
        }

    bunchGenerated = 1;
    return(beam->n_to_track);
    }

long track_beam(
                RUN *run,
                VARY *control,
                ERRORVAL *errcon,
                OPTIM_VARIABLES *optim,
                LINE_LIST *beamline,
                BEAM *beam,
                OUTPUT_FILES *output,
                unsigned long flags,
                long delayOutput,
                double *finalCharge
                )
{    
  double p_central;
  long n_left, effort;

  log_entry("track_beam");

  if (beam->lostOnPass)
    free(beam->lostOnPass);
  if (isSlave || !notSinglePart)
    beam->lostOnPass = tmalloc(sizeof(*(beam->lostOnPass))*beam->n_to_track);

  if (!run)
    bomb("RUN pointer is NULL (track_beam)", NULL);
  if (!control)
    bomb("VARY pointer is NULL (track_beam)", NULL);
  if (!errcon)
    bomb("ERROR pointer is NULL (track_beam)", NULL);
  if (!beamline)
    bomb("LINE_LIST pointer is NULL (track_beam)", NULL);
  if (!beam)
    bomb("BEAM pointer is NULL (track_beam)", NULL);
  if (!output)
    bomb("OUTPUT_FILES pointer is NULL (track_beam)", NULL);

  p_central = run->p_central;

  /* now track particles */
  if (!(flags&SILENT_RUNNING)) {
#if !SDDS_MPI_IO
    fprintf(stdout, "tracking %ld particles\n", beam->n_to_track);
#else
    fprintf(stdout, "tracking %ld particles\n", beam->n_to_track_total);
#endif
    fflush(stdout);
  }
  effort = 0;

#if USE_MPI
  if (isSlave && beam->n_to_track<1 && !lessPartAllowed) {
    printf("*************************************************************************************\n");
    printf("* Warning! The number of particles shouldn't be less than the number of processors! *\n");
    printf("* Less number of processors are recommended!                                        *\n");
    printf("*************************************************************************************\n");
    /*  MPI_Abort(workers, 2); */
  }
  else {  /* do tracking in parallel */ 
    notSinglePart = 1;
    if (lessPartAllowed) { /* If less number of particles is allowed, all processors will excute the same code */
      notSinglePart = 0;
    }
  }

  if (isMaster && notSinglePart) /* This is a makeup to avoid changing the execution order of the serial elegant */
    beam->lostOnPass = NULL;

  if (isSlave || !notSinglePart)
#endif
  if (control->bunch_frequency) {
    makeBucketAssignments(beam, p_central, control->bunch_frequency);
  }
  
  n_left = do_tracking(beam, NULL, 0, &effort, beamline, &p_central, 
                       beam->accepted, &output->sums_vs_z, &output->n_z_points,
                       NULL, run, control->i_step,
                       (!(run->centroid || run->sigma)?FINAL_SUMS_ONLY:0)|
                       ((control->fiducial_flag|flags)&
                        (LINEAR_CHROMATIC_MATRIX+LONGITUDINAL_RING_ONLY+FIRST_BEAM_IS_FIDUCIAL+SILENT_RUNNING
                         +FIDUCIAL_BEAM_SEEN+RESTRICT_FIDUCIALIZATION+PRECORRECTION_BEAM+IBS_ONLY_TRACKING)),
                       control->n_passes, 0, &(output->sasefel), &(output->sliceAnalysis),
		       finalCharge, beam->lostOnPass, NULL);
#if USE_MPI && (!SDDS_MPI_IO)
  notSinglePart = 0; /* All the processors will do the same thing from now */
  /* random_1(-FABS(987654321)); */
#endif  
  if (control->fiducial_flag&FIRST_BEAM_IS_FIDUCIAL && !(flags&PRECORRECTION_BEAM))
    control->fiducial_flag |= FIDUCIAL_BEAM_SEEN;
  
  if (!beam) {
    fprintf(stdout, "error: beam pointer is null on return from do_tracking (track_beam)\n");
    fflush(stdout);
    abort();
  }
  beam->p0 = p_central;
  beam->n_accepted = n_left;

  if (!(flags&SILENT_RUNNING)) {
    extern unsigned long multipoleKicksDone;
#if SDDS_MPI_IO
    if (SDDS_MPI_IO) {
      long total_effort, total_multipoleKicksDone, total_left;
      
      if (notSinglePart) {
	if (isMaster) 
	  multipoleKicksDone = effort = n_left = 0;
	MPI_Reduce (&n_left, &total_left, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce (&effort, &total_effort, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce (&multipoleKicksDone, &total_multipoleKicksDone, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	fprintf(stdout, "%ld particles transmitted, total effort of %ld particle-turns\n", total_left, total_effort);
	fflush(stdout);
	fprintf(stdout, "%lu multipole kicks done\n\n", total_multipoleKicksDone);
      }
      else {
	fprintf(stdout, "%ld particles transmitted, total effort of %ld particle-turns\n", n_left, effort);
	fflush(stdout);
	fprintf(stdout, "%lu multipole kicks done\n\n", multipoleKicksDone);
      }
      fflush(stdout);
    }
#else
    fprintf(stdout, "%ld particles transmitted, total effort of %ld particle-turns\n", n_left, effort);
    fflush(stdout);
    fprintf(stdout, "%lu multipole kicks done\n\n", multipoleKicksDone);
    fflush(stdout);
#endif
  }

  if (!delayOutput)
    do_track_beam_output(run, control, errcon, optim,
                         beamline, beam, output, flags, 
                         finalCharge?*finalCharge:0.0);

  if (!(flags&SILENT_RUNNING)) {
    report_stats(stdout, "Tracking step completed");
    fputs("\n\n", stdout);
  }
  
  log_exit("track_beam");
  return(1);
}

void do_track_beam_output(RUN *run, VARY *control,
                          ERRORVAL *errcon, OPTIM_VARIABLES *optim,
                          LINE_LIST *beamline, BEAM *beam,
                          OUTPUT_FILES *output, unsigned long flags,
                          double finalCharge)
{
  VMATRIX *M;
  long n_left;
  double p_central, p_central0;

  p_central = beam->p0;
  n_left = beam->n_accepted;
  p_central0 = beam->p0_original;
  
  if (run->output && !(flags&INHIBIT_FILE_OUTPUT)) {
    if (!output->output_initialized)
      bomb("'output' file is uninitialized (track_beam)", NULL);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping output beam data..."); fflush(stdout);
      fflush(stdout);
    dump_phase_space(&output->SDDS_output, beam->particle, n_left, control->i_step, p_central,
                     finalCharge);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "done.\n"); fflush(stdout);
      fflush(stdout);
  }

  if (run->acceptance && !(flags&INHIBIT_FILE_OUTPUT)) {
    if (!output->accept_initialized)
      bomb("'acceptance' file is uninitialized (track_beam)", NULL);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping acceptance output..."); fflush(stdout);
      fflush(stdout);
    dump_phase_space(&output->SDDS_accept, beam->accepted, beam->n_accepted, control->i_step, p_central0,
                     finalCharge);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "done.\n"); fflush(stdout);
      fflush(stdout);
  }

  if (run->losses && !(flags&INHIBIT_FILE_OUTPUT)) {
    if (!output->losses_initialized)
      bomb("'losses' file is uninitialized (track_beam)", NULL);
    if (!(flags&SILENT_RUNNING)) {
      fprintf(stdout, "Dumping lost-particle data...\n"); fflush(stdout);
#if SDDS_MPI_IO
      if (SDDS_MPI_IO) {
	long total_left;
	if (notSinglePart) {
	  if (isMaster) n_left = 0;
	  MPI_Reduce (&n_left, &total_left, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	} else 
	  total_left = n_left;
	fprintf(stdout, "n_left = %ld, n_total = %ld, n_lost = %ld\n", 
		total_left, beam->n_original_total, beam->n_original_total-total_left);
      }
#else
      fprintf(stdout, "n_left = %ld, n_total = %ld, n_lost = %ld\n", 
	      n_left, beam->n_to_track, beam->n_to_track-n_left);
#endif
      fflush(stdout);
    }
    dump_lost_particles(&output->SDDS_losses, beam->particle+n_left, beam->lostOnPass+n_left,
			beam->n_to_track-n_left, control->i_step);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "done.\n"); fflush(stdout);
    fflush(stdout);
  }

  if (run->centroid && !run->combine_bunch_statistics && !(flags&FINAL_SUMS_ONLY)  && !(flags&INHIBIT_FILE_OUTPUT)) {
    if (!output->centroid_initialized)
      bomb("'centroid' file is uninitialized (track_beam)", NULL);
    if (!output->sums_vs_z)
      bomb("missing beam sums pointer (track_beam)", NULL);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping centroid data..."); fflush(stdout);
      fflush(stdout);
    dump_centroid(&output->SDDS_centroid, output->sums_vs_z, beamline, output->n_z_points, control->i_step, p_central);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "done.\n"); fflush(stdout);
      fflush(stdout);
  }

  if (run->sigma && !run->combine_bunch_statistics && !(flags&FINAL_SUMS_ONLY) && !(flags&INHIBIT_FILE_OUTPUT)) {
    if (!output->sigma_initialized)
      bomb("'sigma' file is uninitialized (track_beam)", NULL);
    if (!output->sums_vs_z)
      bomb("missing beam sums pointer (track_beam)", NULL);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping sigma data..."); fflush(stdout);
      fflush(stdout);
    dump_sigma(&output->SDDS_sigma, output->sums_vs_z, beamline, output->n_z_points, control->i_step, p_central);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "done.\n"); fflush(stdout);
      fflush(stdout);
  }

  
  if (run->final && !run->combine_bunch_statistics) {
    if (!(M = full_matrix(&(beamline->elem), run, 1)))
      bomb("computation of full matrix for final output failed (track_beam)", NULL);
    if (!output)
      bomb("output pointer is NULL on attempt to write 'final' file (track_beam)", NULL);
    if (!output->final_initialized)
      bomb("'final' file is uninitialized (track_beam)", NULL);
    if (!output->sums_vs_z)
      bomb("beam sums array for final output is NULL", NULL);
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping final properties data..."); fflush(stdout);
    fflush(stdout);
    dump_final_properties
      (&output->SDDS_final, output->sums_vs_z+output->n_z_points, 
       control->varied_quan_value, control->varied_quan_name?*control->varied_quan_name:NULL, 
       control->n_elements_to_vary, control->n_steps*control->indexLimitProduct,
       errcon->error_value, errcon->quan_final_index, errcon->quan_final_duplicates, errcon->n_items,
       optim->varied_quan_value, optim->varied_quan_name?*optim->varied_quan_name:NULL,
       optim->n_variables?optim->n_variables+3:0,
       control->i_step, beam->particle, beam->n_to_track, p_central, M,
       finalCharge);
    if (!(flags&SILENT_RUNNING)) {
      fprintf(stdout, "done.\n"); 
      fflush(stdout);
      fflush(stdout);
    }
    free_matrices(M); free(M); M = NULL;
  }

  if (output->sasefel.active && output->sasefel.filename) {
    if (!(flags&SILENT_RUNNING)) 
      fprintf(stdout, "Dumping SASE FEL data...");
      fflush(stdout);
    doSASEFELAtEndOutput(&(output->sasefel), control->i_step);
    if (!(flags&SILENT_RUNNING)) {
      fprintf(stdout, "done.\n");
      fflush(stdout);
      fflush(stdout);
    }
  }
}



void setup_output(
    OUTPUT_FILES *output,
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    OPTIM_VARIABLES *optim,
    LINE_LIST *beamline
    )
{
    long n_elements;
    n_elements = beamline->n_elems;

    if (run->acceptance) {
        /* prepare dump of accepted particles */
        SDDS_PhaseSpaceSetup(&output->SDDS_accept, run->acceptance, SDDS_BINARY, 1, "accepted phase space", run->runfile,
                             run->lattice, "setup_output");
        output->accept_initialized = 1;
        }

    if (run->centroid) {
        /* prepare dump of centroid vs z */
        SDDS_CentroidOutputSetup(&output->SDDS_centroid, run->centroid, SDDS_BINARY, 1, "centroid output", run->runfile,
                                 run->lattice, "setup_output");
        output->centroid_initialized = 1;
        }

    if (run->sigma) {
        /* prepare dump of sigma vs z */
        SDDS_SigmaOutputSetup(&output->SDDS_sigma, run->sigma, SDDS_BINARY, 1, run->runfile,
                                 run->lattice, "setup_output");
        output->sigma_initialized = 1;
        }

    if (run->output) {
        /* prepare output dump */
        SDDS_PhaseSpaceSetup(&output->SDDS_output, run->output, SDDS_BINARY, 1, "output phase space", run->runfile,
                             run->lattice, "setup_output");
        output->output_initialized = 1;
        }

    if (run->final) {
        /* prepare 'final' dump */
        SDDS_FinalOutputSetup(&output->SDDS_final, run->final, SDDS_BINARY, 1, "final properties", run->runfile, run->lattice, 
                              control->varied_quan_name, control->varied_quan_unit, control->n_elements_to_vary, 
                              errcon->quan_name, errcon->quan_unit, errcon->n_items, 
                              errcon->quan_final_index, &errcon->quan_final_duplicates,
                              optim->varied_quan_name, optim->varied_quan_unit, 
                              (optim->n_variables?optim->n_variables+3:0),
                              "setup_output");
        output->final_initialized = 1;
        }

    if (run->losses) {
        /* prepare dump of lost particles */
        SDDS_BeamLossSetup(&output->SDDS_losses, run->losses, SDDS_BINARY, 1, 
			   "lost particle coordinates", run->runfile,
			   run->lattice, "setup_output");
        output->losses_initialized = 1;
        }

    if (output->sums_vs_z)
      free(output->sums_vs_z);
    output->sums_vs_z = NULL;
    output->n_z_points = 0;
    }

void finish_output(
    OUTPUT_FILES *output,
    RUN *run,
    VARY *control,
    ERRORVAL *errcon,
    OPTIM_VARIABLES *optim,
    LINE_LIST *beamline,
    long n_elements,
    BEAM *beam,
    double finalCharge                   
    )
{
  ELEMENT_LIST *eptr;

  log_entry("finish_output.1");
  
  if (!output)
    bomb("null OUTPUT_FILES pointer (finish_output)", NULL);
  if (!run)
    bomb("null RUN pointer (finish_output)", NULL);
  if (!control)
    bomb("null VARY pointer (finish_output)", NULL);
  if (!errcon)
    bomb("null ERROR pointer (finish_output)", NULL);
  if (!beam)
    bomb("null BEAM pointer (finish_output)", NULL);
  if (!beamline)
    bomb("null BEAMLINE pointer (finish_output)", NULL);

  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type==T_WATCH) {
      WATCH *watch;
      watch = (WATCH*)eptr->p_elem;
      if (watch->initialized) {
        if ((watch->useDisconnect && SDDS_IsDisconnected(&watch->SDDS_table) && 
             !SDDS_ReconnectFile(&watch->SDDS_table)) || !SDDS_Terminate(&watch->SDDS_table)) {
          char s[1024];
          sprintf(s, "Warning: error terminating watch file %s\n", 
                  watch->filename);
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        }
      }
    }
    eptr = eptr->succ;
  }
  
  if (run->centroid && run->combine_bunch_statistics) {
    if (!output->centroid_initialized)
      bomb("'centroid' file is uninitialized (finish_output)", NULL);
    if (!output->sums_vs_z)
      bomb("missing beam sums pointer (finish_output)", NULL);
    dump_centroid(&output->SDDS_centroid, output->sums_vs_z, beamline, output->n_z_points, 0, beam->p0);
  }
  if (run->sigma && run->combine_bunch_statistics) {
    if (!output->sigma_initialized)
      bomb("'sigma' file is uninitialized (finish_output)", NULL);
    if (!output->sums_vs_z)
      bomb("missing beam sums pointer (finish_output)", NULL);
    dump_sigma(&output->SDDS_sigma, output->sums_vs_z, beamline, output->n_z_points, 0, beam->p0);
  }
  if (run->final && run->combine_bunch_statistics) {
    VMATRIX *M;
    if (!output->final_initialized)
      bomb("'final' file is uninitialized (track_beam)", NULL);
    dump_final_properties
      (&output->SDDS_final, output->sums_vs_z+output->n_z_points, 
       control->varied_quan_value, control->varied_quan_name?*control->varied_quan_name:NULL, 
       control->n_elements_to_vary, control->indexLimitProduct*control->n_steps,
       errcon->error_value, errcon->quan_final_index, errcon->quan_final_duplicates,
       errcon->n_items,
       optim->varied_quan_value, optim->varied_quan_name?*optim->varied_quan_name:NULL,
       optim->n_variables?optim->n_variables+2:0,
       0, beam->particle, beam->n_to_track, beam->p0, M=full_matrix(&(beamline->elem), run, 1),
       finalCharge);
    free_matrices(M);
    free(M);
  }
  log_exit("finish_output.1");

  log_entry("finish_output.2");
  if (run->output) {
    if (!output->output_initialized)
      bomb("'output' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_output)) {
      SDDS_SetError("Problem terminating 'output' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    }
    output->output_initialized = 0;
  }
  if (run->acceptance) {
    if (!output->accept_initialized)
      bomb("'acceptance' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_accept)) {
      SDDS_SetError("Problem terminating 'acceptance' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    output->accept_initialized = 0;
  }
  if (run->centroid) {
    if (!output->centroid_initialized)
      bomb("'centroid' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_centroid)) {
      SDDS_SetError("Problem terminating 'centroid' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    output->centroid_initialized = 0;
  }
  if (run->sigma) {
    if (!output->sigma_initialized)
      bomb("'sigma' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_sigma)) {
      SDDS_SetError("Problem terminating 'sigma' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    output->sigma_initialized = 0;
  }
  if (run->final) {
    if (!output->final_initialized)
      bomb("'final' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_final)) {
      SDDS_SetError("Problem terminating 'final' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    output->final_initialized = 0;
  }
  if (bunch) {
    if (!SDDS_bunch_initialized)
      bomb("'bunch' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&SDDS_bunch)) {
      SDDS_SetError("Problem terminating 'bunch' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_bunch_initialized = 0;
  }
  if (run->losses) {
    if (!output->losses_initialized)
      bomb("'losses' file is uninitialized (finish_output)", NULL);
    if (!SDDS_Terminate(&output->SDDS_losses)) {
      SDDS_SetError("Problem terminating 'losses' file (finish_output)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    output->losses_initialized = 0;
  }
  log_exit("finish_output.2");

  log_entry("finish_output.3");
  if (output->sums_vs_z) {
    tfree(output->sums_vs_z);
    output->sums_vs_z = NULL;
  }
  log_exit("finish_output.3");
}

void fill_transverse_structure(
    TRANSVERSE *trans, 
    double emit, 
    double beta, 
    double alpha,
    double eta,
    double etap,
    long xbeam_type,
    double cutoff,
    double *xcentroid
    )
{
    log_entry("fill_transverse_structure");
    trans->emit      = emit;
    trans->beta      = beta;
    trans->alpha     = alpha;
    trans->eta       = eta;
    trans->etap      = etap;
    trans->beam_type = xbeam_type;
    trans->cutoff    = cutoff;
    trans->cent_posi    = xcentroid[0];
    trans->cent_slope   = xcentroid[1];
    log_exit("fill_transverse_structure");
    }

void fill_longitudinal_structure(
    LONGITUDINAL *xlongit,
    double xsigma_dp,
    double xsigma_s,
    double xdp_s_coupling,
    double xemit_z,
    double xbeta_z,
    double xalpha_z,
    double xchirp,
    long xbeam_type,
    double xcutoff,
    double *xcentroid
    )
{
    log_entry("fill_longitudinal_structure");
    xlongit->sigma_dp   = xsigma_dp;
    xlongit->sigma_s    = xsigma_s;
    xlongit->dp_s_coupling = xdp_s_coupling;
    xlongit->emit       = xemit_z;
    xlongit->beta       = xbeta_z;
    xlongit->alpha      = xalpha_z;
    xlongit->beam_type  = xbeam_type;
    xlongit->cutoff     = xcutoff;
    xlongit->cent_s     = xcentroid[0];
    xlongit->cent_dp    = xcentroid[1];
    xlongit->chirp      = xchirp;
    log_exit("fill_longitudinal_structure");
    }

#ifdef VAX_VMS
char *brief_number(double x, char *buffer)
{
    long i, n_p10;
#define N_PREFIXES 17
    struct { char letter; long power; } prefix[N_PREFIXES] = {
        {'A',  18}, {'F',  15}, 
        {'T',  12}, {'G',   9}, {'M',   6}, {'K',   3}, {'h',   2}, 
        {'D',   1}, {'|',   0}, {'d',  -1}, {'c',  -2}, {'m',  -3}, 
        {'u',  -6}, {'n',  -9}, {'p', -12}, {'f', -15}, {'a', -18}
        } ;
        
    if (x==0) {
        sprintf(buffer, "0|");
        return(buffer);
        }

    n_p10 = (long)log10(FABS(x));
    for (i=0; i<N_PREFIXES; i++) {
          if (n_p10>=prefix[i].power)
            break;
          }
    if (i==N_PREFIXES)
        sprintf(buffer, "%1.0lg|", x);
    else {
        if (n_p10>=0) {
            switch (n_p10-prefix[i].power) {
                case 0: 
                    sprintf(buffer, "%.2f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                case 1:
                    sprintf(buffer, "%.1f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                default:
                    sprintf(buffer, "%.0f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                }
            }
        else {
            switch (n_p10-prefix[i].power) {
                case 0:  
                    sprintf(buffer, "%.2f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                case 1:
                    sprintf(buffer, "%.1f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                default:
                    sprintf(buffer, "%.0f%c",
                            x/ipow(10.,prefix[i].power), prefix[i].letter);
                    break;
                }
            }
        }
    return(buffer);
    }
#endif

void makeBucketAssignments(BEAM *beam, double Po, double frequency)
{
  double tmin, tmax, *time, offset;
  long ip, bucket, buckets;

  /* Find time coordinate limits for the beam */
  time = tmalloc(sizeof(*time)*(beam->n_to_track));
  computeTimeCoordinates(time, Po, beam->particle, beam->n_to_track);
  find_min_max(&tmin, &tmax, time, beam->n_to_track);
  /* Determine how many buckets we will span.  Assumes that the bunch length is << bucket spacing */
  if ((buckets = (tmax-tmin)*frequency)>(MAX_BUCKETS-1)) {
    printf("Error: Can't have more than %d buckets at present.", MAX_BUCKETS-1);
    exit(1);
  }
  /* Determine offset to the center of the first bucket. Assumes that the bunch length is << bucket spacing and that
   * bunches have about the same length
   */
  offset = ((tmax-tmin)-buckets/frequency)/2;
  tmin += offset;
  for (ip=0; ip<beam->n_to_track; ip++) {
    /* bucket numbers are 1 or greater */
    bucket = (time[ip]-tmin)*frequency + 1.5;
    beam->particle[ip][6] += bucket/(1.0*MAX_BUCKETS);
  }
  beam->bunchFrequency = frequency;
  /* sort particles in bunch order */
  qsort(beam->particle[0], beam->n_to_track, COORDINATES_PER_PARTICLE*sizeof(double), comp_BucketNumbers);
}

int comp_BucketNumbers(const void *coord1, const void *coord2)
{
  long b1, b2;
  b1 = MAX_BUCKETS*(((double*) coord1)[6]-floor(((double*) coord1)[6]))+0.5;
  b2 = MAX_BUCKETS*(((double*) coord2)[6]-floor(((double*) coord2)[6]))+0.5;
  if (b1<b2)
    return -1;
  else if (b1>b2)
    return  1;
  else 
    return 0;
}

