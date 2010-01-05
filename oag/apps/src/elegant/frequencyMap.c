/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: frequencyMap.c
 * purpose: Do frequency map tracking and analysis.
 *          See file frequencyMap.nl for input parameters.
 *
 * Michael Borland, 2004
 */
#include "mdb.h"
#include "track.h"
#include "frequencyMap.h"

#define IC_X 0
#define IC_Y 1
#define IC_NUX 2
#define IC_NUY 3
#define N_NOCHANGE_COLUMNS 4
#define IC_DNUX 4
#define IC_DNUY 5
#define IC_DNU 6
#define IC_DX 7
#define IC_DY 8
#define IC_DIFFUSION 9
#define N_COLUMNS 10
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"x", "&column name=x, symbol=x, units=m, type=double &end"},
    {"y", "&column name=y, symbol=y, units=m, type=double &end"},
    {"nux", "&column name=nux, symbol=$gn$r$bx$n, type=double &end"},
    {"nuy", "&column name=nuy, symbol=$gn$r$by$n, type=double &end"},
    {"dnux", "&column name=dnux, symbol=$gDn$r$bx$n, type=double &end"},
    {"dnuy", "&column name=dnuy, symbol=$gDn$r$by$n, type=double &end"},
    {"dnu", "&column name=dnu, symbol=$gDn$r, type=double &end"},
    {"dx", "&column name=dx, symbol=$gD$rx, units=m, type=double &end"},
    {"dy", "&column name=dy, symbol=$gD$ry, units=m, type=double &end"},
    {"diffusion", "&column name=diffusion, symbol=\"log$b10$n($gDn$r$bx$n$a2$n+$gDn$r$bx$n$a2$n)\", type=double &end"}
    } ;

#define IP_STEP 0
#define N_PARAMETERS 1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;

static SDDS_DATASET SDDS_fmap;

void setupFrequencyMap(
    NAMELIST_TEXT *nltext,
    RUN *run,
    VARY *control
    )
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&frequency_map, nltext);
  if (echoNamelists) print_namelist(stdout, &frequency_map);
  
  /* check for data errors */
  if (!output)
    bomb("no output filename specified", NULL);
  if (xmin>xmax)
    bomb("xmin > xmax", NULL);
  if (ymin>ymax)
    bomb("ymin > ymax", NULL);
  if (nx<1)
    nx = 1;
  if (ny<1)
    ny = 1;

  output = compose_filename(output, run->rootname);
 #if SDDS_MPI_IO
     SDDS_fmap.parallel_io = 1;
     SDDS_MPI_Setup(&SDDS_fmap, 1, n_processors, myid, MPI_COMM_WORLD, 1);
 #endif
  SDDS_ElegantOutputSetup(&SDDS_fmap, output, SDDS_BINARY, 1, "frequency map analysis",
                          run->runfile, run->lattice, parameter_definition, N_PARAMETERS,
                          column_definition, 
			  include_changes?N_COLUMNS:N_NOCHANGE_COLUMNS, 
			  "setup_frequencyMap", SDDS_EOS_NEWFILE);
  
  if (control->n_elements_to_vary) 
    if (!SDDS_DefineSimpleParameters(&SDDS_fmap, control->n_elements_to_vary,
                                     control->varied_quan_name, control->varied_quan_unit, SDDS_DOUBLE)) {
      SDDS_SetError("Unable to define additional SDDS parameters (setup_aperture_search)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
#if !SDDS_MPI_IO  
  if (!SDDS_WriteLayout(&SDDS_fmap)) {
    SDDS_SetError("Unable to write SDDS layout for aperture search");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
#endif  
}


long doFrequencyMap(
                    RUN *run,
                    VARY *control,
                    double *referenceCoord,
                    ERRORVAL *errcon,
                    LINE_LIST *beamline
                    )
{
  double firstTune[2], secondTune[2], startingCoord[6], endingCoord[6];
  double firstAmplitude[2], secondAmplitude[2];
  double dx, dy, x, y;
  long ix, iy, ip, turns;
  static double **one_part;
  double p;
  long n_part;

#if SDDS_MPI_IO
  long particles;
  /* Open file here for parallel IO */
  if (!SDDS_MPI_File_Open(SDDS_fmap.MPI_dataset, SDDS_fmap.layout.filename, SDDS_MPI_WRITE_ONLY)) 
    SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDS_fmap.MPI_dataset->MPI_file);
  if (!SDDS_MPI_WriteLayout(&SDDS_fmap))  
    SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDS_fmap.MPI_dataset->MPI_file);

  particles = nx*ny/n_processors;
  if (myid < nx*ny%n_processors) 
    particles++;
  if (!SDDS_StartPage(&SDDS_fmap, particles) || 
      !SDDS_SetParameters(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1)) {
    SDDS_SetError("Unable to start SDDS page (do_frequencyMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (verbosity) {
    verbosity = 0;
  if (myid == 0)
    printf ("Warning: In parallel version, no intermediate particle tracking information will be provided\n");
  }
#else
  if (!SDDS_StartPage(&SDDS_fmap, nx*ny) || 
      !SDDS_SetParameters(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, control->i_step, -1)) {
    SDDS_SetError("Unable to start SDDS page (do_frequencyMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
#endif

  if (control->n_elements_to_vary) {
    for (ip=0; ip<control->n_elements_to_vary; ip++)
      if (!SDDS_SetParameters(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip+1,
                              control->varied_quan_value[ip], -1)) {
        SDDS_SetError("Unable to start SDDS page (do_frequencyMap)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
  }

  /* Perform fiducialization by tracking one turn */
  if (!one_part)
    one_part = (double**)czarray_2d(sizeof(**one_part), 1, 7);
  n_part = 1;
  if (referenceCoord) {
    long i;
    for (i=0; i<6; i++)
      one_part[0][i] = referenceCoord[i];
  }
  p = run->p_central;
  if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                   NULL, run, 0, TEST_PARTICLES, 1, 0,
                   NULL, NULL, NULL, NULL, NULL)) {
    printf("Error: lost particle when fiducializing\n");
    exit(1);
  }
  
  if (nx>1)
    dx  = (xmax-xmin)/(nx-1);
  else
    dx = 0;
  if (ny>1)
    dy = (ymax-ymin)/(ny-1);
  else
    dy = 0;
  ip = 0;
  if (include_changes==0)
    turns = control->n_passes;
  else
    turns = control->n_passes/2;
  for (ix=0; ix<nx; ix++) {
    x = xmin + ix*dx;
    for (iy=0; iy<ny; iy++) {
      y = ymin + iy*dy;
      memcpy(startingCoord, referenceCoord, sizeof(*startingCoord)*6);
#if USE_MPI
      if (myid == (ix*ny+iy)%n_processors) /* Partition the job according to particle ID */
#endif
      {
      if (!computeTunesFromTracking(firstTune, firstAmplitude,
				    beamline->matrix, beamline, run,
                                    startingCoord, x, y, turns,
                                    0, endingCoord, NULL, NULL, 1) ||
	  firstTune[0]>1.0 || firstTune[0]<0 || firstTune[1]>1.0 || firstTune[1]<0) {
	if (verbosity) 
	  fprintf(stdout, "Problem with particle %ld tune determination\n", ip);
        continue;
      }
      if (!SDDS_SetRowValues(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip,
                             IC_X, x, IC_Y, y, 
                             IC_NUX, firstTune[0], 
                             IC_NUY, firstTune[1], 
			     -1)) {
        SDDS_SetError("Problem setting SDDS row values (doFrequencyMap)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (include_changes) {
	memcpy(startingCoord, endingCoord, sizeof(*startingCoord)*6);
	if (!computeTunesFromTracking(secondTune, secondAmplitude,
				      beamline->matrix, beamline, run,
				      startingCoord, 0.0, 0.0, turns,
				      0, endingCoord, NULL, NULL, 1) || 
	    secondTune[0]>1.0 || secondTune[0]<0 || secondTune[1]>1.0 || secondTune[1]<0) {
	  if (verbosity)
	    fprintf(stdout, "Problem with particle %ld tune determination\n", ip);
	  if (SDDS_fmap.n_rows) /* If the particle is lost, it will not show in the frequency map */ 
	    SDDS_fmap.n_rows--;
	  continue;
	}
	if (!SDDS_SetRowValues(&SDDS_fmap, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip,
			       IC_DNUX, fabs(secondTune[0]-firstTune[0]), 
			       IC_DNUY, fabs(secondTune[1]-firstTune[1]), 
			       IC_DNU, 
			       sqrt(sqr(secondTune[0]-firstTune[0])+sqr(secondTune[1]-firstTune[1])), 
			       IC_DX, fabs(firstAmplitude[0]-secondAmplitude[0]),
			       IC_DY, fabs(firstAmplitude[1]-secondAmplitude[1]),
                               IC_DIFFUSION, 
                               log10(sqr(secondTune[0]-firstTune[0])+SQR(secondTune[1]-firstTune[1])),
			       -1)) {
	  SDDS_SetError("Problem setting SDDS row values (doFrequencyMap)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      }
      ip++;
      if (verbosity) {
        fprintf(stdout, "Done with particle %ld of %ld\n",
                ix*ny+iy+1, nx*ny);
        fflush(stdout);
      }
      }
    }
  }

  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_fmap);
#if SDDS_MPI_IO
  if (!SDDS_MPI_WriteTable(&SDDS_fmap)) {
#else
  if (!SDDS_WriteTable(&SDDS_fmap)) {
#endif
    SDDS_SetError("Problem writing SDDS table (doFrequencyMap)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  return(1);
}

void finishFrequencyMap()
{
  if (SDDS_IsActive(&SDDS_fmap) && !SDDS_Terminate(&SDDS_fmap)) {
    SDDS_SetError("Problem terminating SDDS output (finish_aperture_search)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

