/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: sdds_beam.c
 * purpose: Do tracking for beams from sdds tables.
 *          See file sdds_beam.nl for input parameters.
 *
 * Michael Borland, 1994
 */
#include "mdb.h"
#include "track.h"
#include "sdds_beam.h"

static long input_type_code;
#define ELEGANT_BEAM 0
#define SPIFFE_BEAM 1
#define N_SDDS_INPUT_TYPES 2
static char *input_type_name[N_SDDS_INPUT_TYPES] = {
  "elegant", "spiffe"
  } ;

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

  /* SDDS routines will be asked to deliver data in the order given in
   * this string, which matches the order of the #define's 
   */
static char *spiffe_columns = "r pr pz t pphi";
static char *spiffeColumn[5] = {
  "r", "pr", "pz", "t", "pphi"
  };
#define ISC_R 0
#define ISC_PR 1
#define ISC_PZ 2
#define ISC_T 3
#define ISC_PPHI 4

/* SDDS routines will be asked to deliver data in the order given in
 * this string, which matches the order of the #define's 
 */
static char *elegant_columns = "x xp y yp t p";
static char *elegantColumn[6] = {
  "x", "xp", "y", "yp", "t", "p"
  };
#define IEC_X 0
#define IEC_XP 1
#define IEC_Y 2
#define IEC_YP 3
#define IEC_T 4
#define IEC_P 5

long get_sdds_particles(double ***particle, long one_dump, long n_skip);

static char **inputFile = NULL;     /* input filenames */
static long inputFiles = 0;         /* number of input files */
static long inputFileIndex = 0;     /* index of in-use input file */
static long input_initialized = 0;  /* in-use input file has been opened */
static long has_been_read = 0;      /* data has been read from input file */

static SDDS_TABLE SDDS_input;

void setup_sdds_beam(
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
  static long initial_call = 1;

  log_entry("setup_sdds_beam");

  if (!beam)
    bombElegant("BEAM pointer is null in setup_sdds_beam", NULL);
  if (!nltext)
    bombElegant("NAMELIST_TEXT pointer is null in setup_sdds_beam", NULL);
  if (!run)
    bombElegant("RUN pointer is null in setup_sdds_beam", NULL);
  if (!control)
    bombElegant("VARY pointer is null in setup_sdds_beam", NULL);
  if (!errcon)
    bombElegant("ERROR pointer is null in setup_sdds_beam", NULL);
  if (!output)
    bombElegant("OUTPUT_FILES pointer is null in setup_sdds_beam", NULL);
  if (!beamline)
    bombElegant("beamline pointer is null in setup_sdds_beam", NULL);
  if (!run->runfile || !run->lattice)
    bombElegant("null runfile or lattice pointer in RUN structure in setup_sdds_beam", NULL);

  if (!initial_call)
    get_sdds_particles(NULL, 0, 0);
  else
    initial_call = 0;

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&sdds_beam, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &sdds_beam);
  fflush(stdout);

  /* check for validity of namelist inputs */
  if (input==NULL && input_list==NULL)
    bombElegant("no input file given in namelist sdds_beam", NULL);
  if ((selection_parameter && !selection_string) || (!selection_parameter && selection_string))
    bombElegant("must specify selection_parameter and selection_string together or not at all", NULL);

  if (inputFile)
    free(inputFile);
  inputFile = NULL;
  if (input) {
    inputFiles = 1;
    inputFile = tmalloc(sizeof(*inputFile)*1);
    cp_str(inputFile, input);
    inputFile[0] = compose_filename(inputFile[0], run->rootname);
  }
  else {
    char *ptr;
    inputFiles = 0;
    while ((ptr=get_token(input_list))) {
      inputFiles += 1;
      if (!(inputFile = SDDS_Realloc(inputFile, sizeof(*inputFile)*inputFiles)))
        bombElegant("memory allocation failure", NULL);
      cp_str(inputFile+inputFiles-1, ptr);
      inputFile[inputFiles-1] = compose_filename(inputFile[inputFiles-1], run->rootname);
    }
  }
  
  inputFileIndex = input_initialized = has_been_read = 0;

  if ((input_type_code=match_string(input_type, input_type_name, N_SDDS_INPUT_TYPES, 0))<0)
    bombElegant("unknown sdds input type", NULL);
  if (input_type_code==SPIFFE_BEAM && n_particles_per_ring<=0)
    bombElegant("n_particles_per_ring is invalid", NULL);
  if (input_type_code!=SPIFFE_BEAM)
    n_particles_per_ring = 1;
  if (n_particles_per_ring!=1 && one_random_bunch) 
    bombElegant("must have n_particles_per_ring==1 for one_random_bunch!=0", NULL);
  if (p_lower>p_upper)
    bombElegant("p_lower and p_upper are invalid", NULL);
  if (sample_interval<1)
    bombElegant("sample_interval < 1", NULL);
  if (sample_fraction>1)
    bombElegant("sample_fraction > 1", NULL);
  if (sample_fraction<1 && sample_interval>1)
    bombElegant("either sample_fraction or sample_interval must be 1", NULL);
  if (save_initial_coordinates && !reuse_bunch)
    save_initial_coordinates = 0;
/*
  if (reuse_bunch && inputFiles!=1) {
    fprintf(stdout, "*** Warning: reuse_bunch=1 with several input files.  Only the first file will be used.\n");
    while (inputFiles>1) 
      free(inputFile[--inputFiles]);
  }
*/
  
  beam->original = beam->particle = beam->accepted = NULL;
  beam->n_original = beam->n_to_track = beam->n_accepted = beam->n_saved = beam->n_particle = 0;
  save_initial_coordinates = save_original || save_initial_coordinates;
  
  log_exit("setup_sdds_beam");
}

long new_sdds_beam(
                   BEAM *beam,
                   RUN *run,
                   VARY *control,
                   OUTPUT_FILES *output,
                   long flags
                   )
{    
  double p, t_offset, gamma, beta, p_central;
#ifdef VAX_VMS
  char s[100];
#endif
  long i_store, i, j;
  static long new_particle_data, generate_new_bunch;
  double r, theta, pr, pz, pphi, delta;
  double sin_theta, cos_theta, path, *pti;

  log_entry("new_sdds_beam");

#if USE_MPI
  if (control->fiducial_flag&FIRST_BEAM_IS_FIDUCIAL) {
    lessPartAllowed = 1; /* Special for the fiducial beam */
  }
  else {
    lessPartAllowed = 0; /* This flag will control if the simulation runs in single particle mode in track_beam function */
  }
  notSinglePart = 1;  /* All CPUs will track the fiducial beam in parallel */
  partOnMaster = 0;
#endif

  if (flags&TRACK_PREVIOUS_BUNCH) {
    /* retracking bunch that has already been set up */
    if (!save_initial_coordinates)
      bombElegant("logic error---initial beam coordinates not saved", NULL);
#if SDDS_MPI_IO
    if (isSlave || !notSinglePart)
#endif
    if (beam->original==NULL)
      bombElegant("can't retrack with previous bunch--there isn't one!", NULL);
    if (n_particles_per_ring!=1) {
      fputs("Warning: can't do retracking with previous bunch when n_particles_per_ring!=1\n", stdout);
      fputs("Will use a new bunch generated from previously read data.\n", stdout);
      generate_new_bunch = 1;
    }
    else
      generate_new_bunch = 0;
  }
  else {
    if (!prebunched) {
      /* The beam in each input file is to be treated as a single bunch,
       * even though it may be spread over several pages.
       * Read in all the particles from the input file and allocate arrays
       * for storing initial coordinates and coordinates of accepted 
       * particles. */
      if ((beam->original==NULL) && !has_been_read) {
        /* no beam has been read before, or else it was purged from memory to save RAM */
        /* free any arrays we may have from previous pass */
        if (beam->particle)
          free_czarray_2d((void**)beam->particle, beam->n_particle, 7);
        if (beam->accepted)
          free_czarray_2d((void**)beam->accepted, beam->n_particle, 7);
        beam->particle = beam->accepted = beam->original = NULL;
        /* read the particle data */
        if ((beam->n_original=get_sdds_particles(&beam->original, prebunched, 0))<0) {
          bombElegant("no particles in input file", NULL);
        }
        if (reuse_bunch) {
          if (save_initial_coordinates)
            has_been_read = 1;
          else if (inputFileIndex>=inputFiles)
            inputFileIndex = 0; 
        }
        if (save_initial_coordinates || n_particles_per_ring!=1)
          beam->particle = (double**)czarray_2d
            (sizeof(double), beam->n_particle=(long)(n_particles_per_ring*beam->n_original*factor), 7);
        else {
          beam->particle = beam->original;
          beam->n_particle = beam->n_original;
        }
        if (run->acceptance)
#if USE_MPI
        if (isSlave || !notSinglePart)
#endif
          beam->accepted = (double**)czarray_2d
            (sizeof(double), (long)(beam->n_particle*factor), 7);
	if (run->losses)  
#if USE_MPI
        if (isSlave || !notSinglePart)
#endif
	  beam->lost = (double**)czarray_2d(sizeof(double),(long)(beam->n_particle*factor), 8);
        new_particle_data = 1;	
      }
      else
        /* since we've read the whole file already and saved it, there is no new data */
        new_particle_data = 0;
    }
    else {
      /* Each page in the input file is to be treated as a separate
       * bunch.  Read in the next page for tracking.
       */
      /* Free arrays from previous pass */
      if (beam->particle)
        free_czarray_2d((void**)beam->particle, beam->n_particle, 7);
      if (beam->accepted)
        free_czarray_2d((void**)beam->accepted, beam->n_particle, 7);
      if (beam->original && beam->original!=beam->particle)
        free_czarray_2d((void**)beam->original, beam->n_original, 7);
      beam->particle = beam->accepted = beam->original = NULL;

      /* read the new page */
      if ((beam->n_original=get_sdds_particles(&beam->original, prebunched, n_tables_to_skip))>=0) { 
#if SDDS_MPI_IO
        if (isSlave || !notSinglePart)  
#endif
	{
	  n_tables_to_skip = 0;    /* use the user's parameter only the first time */
	  if (save_initial_coordinates || n_particles_per_ring!=1)
	    beam->particle = (double**)czarray_2d
	      (sizeof(double), beam->n_particle=(long)(n_particles_per_ring*beam->n_original*factor), 7);
	  else {
	    beam->particle = beam->original;
	    beam->n_particle = beam->n_original;
	  }
	  if (run->acceptance) 
	    beam->accepted = (double**)czarray_2d
	      (sizeof(double), (long)(beam->n_particle*factor), 7);
	  if (run->losses)  
	    beam->lost = (double**)czarray_2d(sizeof(double), (long)(beam->n_particle*factor), 8);
	}
      }
      else { 
        log_exit("new_sdds_beam");
        return(-1);
      }
      new_particle_data = 1;
    }
  }

  t_offset = (control->bunch_frequency?(control->i_step-1)/control->bunch_frequency:0);

  p_central = beam->p0_original = run->p_central;
#if SDDS_MPI_IO
  if (isSlave && notSinglePart && new_particle_data)  { /* Compute the offset of particle ID for different processors */
      long sum=0, tmp, my_offset, *offset = tmalloc(n_processors*sizeof(*offset)),
	n_particle=(beam->n_original)/sample_interval*sample_fraction;
      MPI_Allgather(&n_particle, 1, MPI_LONG, offset, 1, MPI_LONG, workers);
      tmp = offset[0];
      for (i=1; i<n_processors; i++) {
	sum += tmp;
	tmp = offset[i];
	offset[i] = sum; 
      }
      offset[0] = 0; 
      my_offset = offset[myid-1];
      if (input_type_code==SPIFFE_BEAM) 
	particleID = my_offset+1;
      else
	particleID += my_offset;
  }
#endif
  if (new_particle_data || generate_new_bunch || 
      (input_type_code==SPIFFE_BEAM && !one_random_bunch && !(flags&TRACK_PREVIOUS_BUNCH))) {
    /* Create the initial distribution from the beam->original particle data 
     * or generate a new distribution from those data 
     */
#if SDDS_MPI_IO
    if (isSlave || (!notSinglePart)) {
#endif
    if (input_type_code==SPIFFE_BEAM) {
      if (!beam->original)
        bombElegant("beam->original array is NULL (new_sdds_beam-2)", NULL);
      if (!beam->particle)
        bombElegant("beam->particle array is NULL (new_sdds_beam-2)", NULL);
      for (i=i_store=0; i<beam->n_original; i+=sample_interval) {
        if (!beam->original[i]) {
          fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam-2)\n", i);;
          fflush(stdout);
          exitElegant(1);
        }
        if (sample_fraction!=1 && random_4(1)>sample_fraction)
          continue;
        pti = beam->original[i];
        pz = pti[ISC_PZ];
        pr = pti[ISC_PR];
        pphi = pti[ISC_PPHI];
        p = sqrt(sqr(pz) + sqr(pr) + sqr(pphi));
        gamma = sqrt(sqr(p)+1);
        if (p_lower && (p_lower>p || p_upper<p))
          continue;
        r      = pti[ISC_R];
        path   = (t_offset + pti[ISC_T])*c_mks*(beta = p/gamma) ;
        delta  = (p-p_central)/p_central;
        theta  = PIx2*random_4(1);
        for (j=0; j<n_particles_per_ring; j++, i_store++) {
          sin_theta = sin(theta);
          cos_theta = cos(theta);
          theta += PIx2/n_particles_per_ring;
          if (!beam->particle[i_store]) {
            fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam-2)\n", i_store);
            fflush(stdout);
            exitElegant(1);
          }
          beam->particle[i_store][0] = r*cos_theta;
          beam->particle[i_store][1] = (pr*cos_theta - pphi*sin_theta)/pz;
          beam->particle[i_store][2] = r*sin_theta;
          beam->particle[i_store][3] = (pr*sin_theta + pphi*cos_theta)/pz;
          beam->particle[i_store][4] = path;
          beam->particle[i_store][5] = delta;
          beam->particle[i_store][6] = particleID++;
        }
      }
      beam->n_to_track = i_store;
#ifndef VAX_VMS
      for (i_store=0; i_store<beam->n_to_track; i_store++) {
        for (i=0; i<6; i++) {
          if (!beam->particle[i_store]) {
            fprintf(stdout, "error: beam->particle[%ld] is NULL\n", i_store);
            fflush(stdout);
            exitElegant(1);
          }
          if (isnan(beam->particle[i_store][i]) || isinf(beam->particle[i_store][i])) {
            fprintf(stdout, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
            fflush(stdout);
            exitElegant(1);
          }
        }
      }
#endif
      if (center_transversely) {
        zero_centroid(beam->particle, beam->n_to_track, 0);
        zero_centroid(beam->particle, beam->n_to_track, 1);
        zero_centroid(beam->particle, beam->n_to_track, 2);
        zero_centroid(beam->particle, beam->n_to_track, 3);
      }
      if (center_arrival_time || reverse_t_sign)
        adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central, 
                                 center_arrival_time, reverse_t_sign);
    }
    else {
      /* In this case, the data is for point particles already,
       * so I just copy the data for the most part, except for sampling.
       */
      if (!beam->original && beam->n_original)
        bombElegant("beam->original array is NULL (new_sdds_beam)", NULL);
      if (!beam->particle)
        bombElegant("beam->particle array is NULL (new_sdds_beam)", NULL);
      for (i=i_store=0; i<beam->n_original; i_store++,i+=sample_interval) {
        if (sample_fraction!=1 && random_4(1)>sample_fraction) {
          i_store--;
          continue;
        }
        if (!beam->original[i]) {
          fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.2)\n", i);
          fflush(stdout);
          exitElegant(1);
        }
        gamma = sqrt(sqr(p=beam->original[i][IEC_P])+1);
        if (p_lower && (p_lower>p || p_upper<p)) {
          i_store--;
          continue;
        }
        if (!beam->particle[i_store]) {
          fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.2)\n", i_store);
          fflush(stdout);
          exitElegant(1);
        }
        for (j=0; j<4; j++)
          beam->particle[i_store][j] = beam->original[i][j];
        /* convert time to path-length */
        beam->particle[i_store][4] = (t_offset+beam->original[i][IEC_T])*c_mks*p/gamma;
        /* convert energy to dp/p */
        beam->particle[i_store][5] = (p-p_central)/p_central;
        beam->particle[i_store][6] = beam->original[i][6];
      }
      beam->n_to_track = i_store;
#ifndef VAX_VMS
      for (i_store=0; i_store<beam->n_to_track; i_store++) {
        for (i=0; i<6; i++) {
          if (isnan(beam->particle[i_store][i]) || isinf(beam->particle[i_store][i])) {
            fprintf(stdout, "error: NaN or Infinity detected in initial particle data, coordinate %ld\n", i);
            fflush(stdout);
            exitElegant(1);
          }
        }
      }
#endif
      if (center_transversely) {
        zero_centroid(beam->particle, beam->n_to_track, 0);
        zero_centroid(beam->particle, beam->n_to_track, 1);
        zero_centroid(beam->particle, beam->n_to_track, 2);
        zero_centroid(beam->particle, beam->n_to_track, 3);
      }
      if (center_arrival_time || reverse_t_sign)
        adjust_arrival_time_data(beam->particle, beam->n_to_track, p_central, 
                                 center_arrival_time, reverse_t_sign);
    }
#if SDDS_MPI_IO
    }
#endif
  }
  else {
    /* use (x, x', y, x', s, dp/p) saved in beam->original[] */
    if (!(beam->n_to_track = beam->n_saved)) {
#if SDDS_MPI_IO
      if (! beam->n_to_track_total) { /* Problem!!! In the parallel version, we need check the total number of particles */
#endif
      log_exit("new_sdds_beam");
      return -1;
    }
#if SDDS_MPI_IO
    }
    /* We need change to the singlePart mode for certain special situations */
    if (do_find_aperture || runInSinglePartMode) {
      notSinglePart = 0;
      lessPartAllowed = 1;
      partOnMaster = 1;	
    } 
    if (isSlave || (!notSinglePart) || (partOnMaster)) {
#endif
    if (!beam->original)
      bombElegant("beam->original is NULL (new_sdds_beam.3)", NULL);
    if (!beam->particle)
      bombElegant("beam->particle is NULL (new_sdds_beam.3)", NULL);
    for (i=0; i<beam->n_saved; i++) {
      if (!beam->original[i]) {
        fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.3)\n", i);
        fflush(stdout);
        exitElegant(1);
      }
      if (!beam->particle[i]) {
        fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.3)\n", i);
        fflush(stdout);
        exitElegant(1);
      }
      for (j=0; j<7; j++) 
        beam->particle[i][j] = beam->original[i][j];
      p = p_central*(1+beam->particle[i][5]);
      beta = p/sqrt(p*p+1);
      beam->particle[i][4] += t_offset*beta*c_mks;
    }
    new_particle_data = 0;
#if SDDS_MPI_IO
    }
#endif
  }
#if SDDS_MPI_IO
  if (notSinglePart) {
    if (isMaster)
      beam->n_to_track = 0;
    MPI_Allreduce (&(beam->n_to_track), &(beam->n_to_track_total), 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    beam->n_original_total = beam->n_to_track_total;
  }

  /* We need change to the singlePart mode for certain special situations */
  if (do_find_aperture || runInSinglePartMode) {
    notSinglePart = 0;
    lessPartAllowed = 1;
    partOnMaster = 1;	
  } 
  if (isSlave || (!notSinglePart) || partOnMaster) {
#endif
  if (new_particle_data && save_initial_coordinates && 
      (one_random_bunch || (reuse_bunch && input_type_code!=SPIFFE_BEAM))) {
    /* Copy the new "initial" data into original[] in case it needs to be reused,  
       but only if the stuff already in original[] is not going to be needed again
       to generate new beams.
       */
    if (beam->original==beam->particle) 
      bombElegant("logic error in new_sdds_beam: array for original coordinates is missing", NULL);
#if SDDS_MPI_IO
  if (!notSinglePart && new_particle_data) {/* Each processor will hold a copy of the whole beam */
    long np_total, np=beam->n_to_track;
    double **data_all=NULL, **data=beam->particle;
    int *offset_array= (int*)tmalloc(n_processors*sizeof(*offset_array));
    long *n_vector_long = (long*)tmalloc(n_processors*sizeof(*n_vector_long));
    int *n_vector = (int*)tmalloc(n_processors*sizeof(*n_vector));
    
    MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    data_all = (double**)resize_czarray_2d((void**)data_all, sizeof(double), (long)(np_total), 7);
    
    MPI_Allgather(&np, 1, MPI_LONG, n_vector_long, 1, MPI_LONG, MPI_COMM_WORLD);
    offset_array[0] = 0;
    for (i=0; i<n_processors-1; i++) {
      n_vector[i] = 7*n_vector_long[i]; 
      offset_array[i+1] = offset_array[i] + n_vector[i];
    }
    n_vector[n_processors-1] = 7*n_vector_long[n_processors-1]; 

    if (!beam->n_original) /* A dummy memory is allocated on Master */ 
      data = (double**)czarray_2d(sizeof(**data), 1, 7);
    MPI_Allgatherv(&data[0][0], 7*np, MPI_DOUBLE, &data_all[0][0], n_vector, offset_array, MPI_DOUBLE, MPI_COMM_WORLD );
  
    if (data)
      free_czarray_2d((void**)data, np, 7);
    if (offset_array)
      tfree(offset_array);
    if (n_vector)
      tfree(n_vector);
    if (n_vector_long)
      tfree(n_vector_long);
    if (beam->particle)
      beam->particle = data_all;
    beam->n_to_track = np_total;

    /* The memory for original will be reallocated to hold the whole beam */
    beam->original = (double**)resize_czarray_2d((void**)beam->original, sizeof(double), (long)(np_total), 7);
  }
#endif
    if (!beam->original)
      bombElegant("beam->original is NULL (new_sdds_beam.4)", NULL);
    if (!beam->particle)
      bombElegant("beam->particle is NULL (new_sdds_beam.4)", NULL);
    for (i=0; i<beam->n_to_track; i++) {
      if (!beam->original[i]) {
        fprintf(stdout, "error: beam->original[%ld] is NULL (new_sdds_beam.4)\n", i);
        fflush(stdout);
        exitElegant(1);
      }
      if (!beam->particle[i]) {
        fprintf(stdout, "error: beam->particle[%ld] is NULL (new_sdds_beam.4)\n", i);
        fflush(stdout);
        exitElegant(1);
      }
      for (j=0; j<7; j++)
        beam->original[i][j] = beam->particle[i][j];
    }
    beam->n_saved = beam->n_to_track;
    new_particle_data = 0;
  }
#if SDDS_MPI_IO
  }
#endif
  if (!save_initial_coordinates) {
    if (reuse_bunch) {
      /* close the SDDS file to free memory and read again from scratch next time */
      if (!SDDS_Terminate(&SDDS_input)) {
        SDDS_SetError("Problem terminate sdds beam input");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      input_initialized = 0;
    }
    /* free the 'original' particle data array */
    if (beam->original && beam->original!=beam->particle) {
      free_czarray_2d((void**)beam->original, beam->n_original, 7);
      beam->original = NULL;
      beam->n_original = 0;
    }
  if (beam->original==beam->particle)
      beam->original = NULL;
  }

  log_exit("new_sdds_beam");
  return(beam->n_to_track);
}

/* get_sdds_particles reads data from all of the input files and puts it into the particle array as
 * follows:
 *     for spiffe input:   (*particle)[i] = (r, pr, pz, pphi, t) for ith particle
 *     for elegant input:  (*particle)[i] = (x, xp, y, yp, t, p) for ith particle
 */
long get_sdds_particles(double ***particle, 
                        long one_dump,   /* read only one page */
                        long n_skip      /* number of pages to skip */
                        )
{
  int32_t i;
  long np_max, np, np_new, rows, dump_rejected;
  long retval, data_seen;
  double **data=NULL;
  static char s[200];
  long indexID = -1;
  double *columnData;
  long ic;
#if SDDS_MPI_IO
  long total_points = 0, total_rows; 
  int32_t active = 0; 
#endif
  log_entry("get_sdds_particles");

  if (particle==NULL) {
    /* reset for reading again */
    if (input_initialized && !SDDS_Terminate(&SDDS_input)) {
      sprintf(s, "Problem terminating SDDS beam input from file %s", inputFile[inputFileIndex]);
      SDDS_SetError(s);
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    input_initialized = 0;
    return 0;
  }
  if (inputFileIndex>=inputFiles) {
    fprintf(stdout, "No particles left in input file(s)\n");
    fflush(stdout);
    return -1;
  }
  
  retval = data_seen = np = np_new = 0;
#if SDDS_MPI_IO
  if (SDDS_input.MPI_dataset == NULL) /* For the case of one_dump, this setup should be only called once */
    SDDS_MPI_Setup(&SDDS_input, 1, n_processors, myid, MPI_COMM_WORLD, 0); /* Only slaves will read the input data */  
#endif
  while (inputFileIndex<inputFiles) {
    if (!input_initialized) {
#if SDDS_MPI_IO
      if (!SDDS_MPI_InitializeInputFromSearchPath(&SDDS_input, inputFile[inputFileIndex])) 
#else
      if (!SDDS_InitializeInputFromSearchPath(&SDDS_input, inputFile[inputFileIndex])) 
#endif
      {
        sprintf(s, "Problem opening beam input file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }

      input_initialized = 1;

      if (selection_parameter) {
        if ((i=SDDS_GetParameterIndex(&SDDS_input, selection_parameter))<0)
          fprintf(stdout, "warning: SDDS beam file %s does not contain the selection parameter %s\n",
                  inputFile[inputFileIndex], selection_parameter);
        fflush(stdout);
        if (SDDS_GetParameterType(&SDDS_input, i)!=SDDS_STRING) {
          sprintf(s, "SDDS beam file %s contains parameter %s, but parameter is not a string", 
                  inputFile[inputFileIndex], selection_parameter);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
      }
      if (input_type_code==SPIFFE_BEAM) {
        if (!check_sdds_column(&SDDS_input, "r", "m") || 
            !check_sdds_column(&SDDS_input, "pr", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "pz", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "pphi", "m$be$nc") ||
            !check_sdds_column(&SDDS_input, "t", "s")) {
          fprintf(stdout, 
                  "necessary data quantities (r, pr, pz, t) have the wrong units or are not present in %s", 
                  inputFile[inputFileIndex]);
          fflush(stdout);
          exitElegant(1);
        }
      }
      else {
        if (!check_sdds_column(&SDDS_input, "x", "m") ||
            !check_sdds_column(&SDDS_input, "y", "m") ||
            !check_sdds_column(&SDDS_input, "xp", NULL) ||
            !check_sdds_column(&SDDS_input, "yp", NULL) ||
            !check_sdds_column(&SDDS_input, "t", "s")) {
          fprintf(stdout, 
                  "necessary data quantities (x, x', y, y', t, p) have the wrong units or are not present in %s\n", 
                  inputFile[inputFileIndex]);
          fflush(stdout);
          exitElegant(1);
        }
        if (!check_sdds_column(&SDDS_input, "p", "m$be$nc")) {
          if (check_sdds_column(&SDDS_input, "p", NULL)) {
            fprintf(stdout, "Warning: p has no units.  Expected m$be$nc\n");
            fflush(stdout);
          }
        }
      }
      fprintf(stdout, "File %s opened and checked.\n", inputFile[inputFileIndex]);
      fflush(stdout);
    }
    
    np_max = np = 0;
    data = NULL;
    data_seen = 1;
    while (data_seen) {
      data_seen = 0;
#if SDDS_MPI_IO
      if ((retval=SDDS_MPI_ReadPage(&SDDS_input))==0)
#else
      if ((retval=SDDS_ReadTable(&SDDS_input))==0)
#endif
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      else if (retval!=-1) 
        data_seen = 1;
      else 
        /* end of file */
        break;
#if SDDS_MPI_IO
      particleID = (total_points + SDDS_input.MPI_dataset->start_row + 1); /* Add offset for particleID in the version with parallel I/O */
      total_points += (long) SDDS_MPI_TotalRowCount(&SDDS_input);
#endif
      printf("Read page %ld from file\n", retval);
      fflush(stdout);
      if (one_dump && n_skip>0) {
        n_skip--;
        printf("Page is being skipped.\n");
        fflush(stdout);
        continue;
      }
      dump_rejected = 0;
      if (selection_parameter) {
        char *value;
        if (!SDDS_GetParameter(&SDDS_input, selection_parameter, &value)) {
          sprintf(s, "Problem getting value of parameter %s from file %s", selection_parameter, inputFile[inputFileIndex]);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!wild_match(value, selection_string)) {
          dump_rejected = 1;
          printf("Page rejected by selection parameter value\n");
          fflush(stdout);
          break;
        }
      }
#if !SDDS_MPI_IO
      if ((rows = SDDS_RowCount(&SDDS_input))<=0) {
#else
      rows = SDDS_input.n_rows;
      total_rows = SDDS_MPI_TotalRowCount(&SDDS_input);
      if (total_rows<=0) {	
#endif
        if (rows==-1) {
          sprintf(s, "Problem counting rows of interest for file %s", inputFile[inputFileIndex]);
          SDDS_SetError(s);
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!one_dump) {
#if !SDDS_MPI_IO
          printf("Continuing to try to read file.\n");
          fflush(stdout);
#else
	  fprintf(stdout, "%ld rows in page %ld\n", total_rows, retval);
#endif
          continue;
        }
      }
#if !SDDS_MPI_IO
      fprintf(stdout, "%ld rows in page %ld\n", rows, retval);
#else
      fprintf(stdout, "%ld rows in page %ld\n", total_rows , retval);
#endif
      fflush(stdout);
      SDDS_SetColumnFlags(&SDDS_input, 0);
      if (!SDDS_SetColumnsOfInterest(&SDDS_input, SDDS_NAMES_STRING, 
                                     input_type_code==SPIFFE_BEAM?spiffe_columns:elegant_columns)) {
        sprintf(s, "Problem setting columns of interest for file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
#if SDDS_MPI_IO  
      if (isSlave) {
#endif
	if ((np_new=np+rows)>np_max) {
	  /* must reallocate to get more space */
	  np_max = np + rows;
	  data = (double**)resize_czarray_2d((void**)data, sizeof(double), (long)(np_max), 7);
	}

#if SDDS_MPI_IO
	if (total_rows<(n_processors-1)) {
	  /* if the particles number is less than n_processors, the first n_particle slaves will extract information */ 
	  if ((myid>=1) && (myid<=total_rows)) 
	    active = 1;
	}
	else if (isSlave)
	  active = 1;     	
	if (!active) /* Avoid to crash caused by returning a NULL pointer for fiducialization */ 
	  data = (double**)resize_czarray_2d((void**)data, sizeof(double), (long)(np_max), 7);
	else
	{
#endif  
	    for (ic=0; ic<(input_type_code==SPIFFE_BEAM?5:6); ic++) {
	      if (!(columnData=SDDS_GetColumnInDoubles(&SDDS_input, input_type_code==SPIFFE_BEAM?spiffeColumn[ic]:elegantColumn[ic]))) {
		sprintf(s, "Problem getting column %s for file %s", input_type_code==SPIFFE_BEAM?spiffeColumn[ic]:elegantColumn[ic], inputFile[inputFileIndex]);
		SDDS_SetError(s);
		SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	      }

	      for (i=np; i<np_new; i++) 
		data[i][ic] = columnData[i-np];
	      free(columnData);
	    }
	    if ((indexID=SDDS_GetColumnIndex(&SDDS_input, "particleID"))>=0) {
	      double *index;
	      if (!(index=SDDS_GetColumnInDoubles(&SDDS_input, "particleID"))) {
		sprintf(s, "Problem reading particleID column for file %s", inputFile[inputFileIndex]);
		SDDS_SetError(s);
		SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	      }
	      for (i=np; i<np_new; i++)
		data[i][6] = index[i-np];
	      free(index);
	    }
	    else if (input_type_code!=SPIFFE_BEAM)	
	      for (i=np; i<np_new; i++)
		data[i][6] = particleID++;
#if SDDS_MPI_IO
	  }  /* End of active */
      } /* End of isSlave */
#endif
    np = np_new;
    if (one_dump && !dump_rejected)
      break;
      }
    if (retval==-1) {
      /* go to next file */
#if SDDS_MPI_IO
      if (!SDDS_MPI_Terminate(&SDDS_input)) {
#else
      if (!SDDS_Terminate(&SDDS_input)) {
#endif
        sprintf(s, "Problem terminating SDDS beam input from file %s", inputFile[inputFileIndex]);
        SDDS_SetError(s);
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      fprintf(stdout, "File %s was used up and closed.\n", inputFile[inputFileIndex]);
      fflush(stdout);
      inputFileIndex ++;
      input_initialized = 0;
    }
#if !SDDS_MPI_IO      
    if (np)
      break;
#else
    if (SDDS_MPI_IO) {
      long np_total;
      MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (np_total)
	break;
    }
#endif
    fprintf(stdout, "Checking next file\n");
    fflush(stdout);
  }    
  
  if (input_initialized && !SDDS_ShortenTable(&SDDS_input, 1)) {
    SDDS_SetError("Problem releasing table memory when reading SDDS beam file.");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  
  if (!data_seen && one_dump)
    return -1;
#if SDDS_MPI_IO  
  fprintf(stdout, "a total of %ld data points were read\n\n", total_points);
#else
  fprintf(stdout, "a total of %ld data points were read\n\n", np);
#endif
  fflush(stdout);
  if (particle)
    *particle = data;

  log_exit("get_sdds_particles");

  return(np);
}            

void adjust_arrival_time_data(double **coord, long np, double Po, long center_t, long flip_t)
{
  long ip;
  double P, beta;
  double tave;

  if (np<1)
    return;

  if (flip_t) 
    for (ip=0; ip<np; ip++)
      coord[ip][4] *= -1;
  
  if (center_t) {
    for (ip=tave=0; ip<np; ip++) {
      P = Po*(1+coord[ip][5]);
      beta = P/sqrt(sqr(P)+1);
      tave += coord[ip][4]/beta;
    }
    tave /= np;
    for (ip=0; ip<np; ip++) {
      P = Po*(1+coord[ip][5]);
      beta = P/sqrt(sqr(P)+1);
      coord[ip][4] -= beta*tave;
    }
  }
}
