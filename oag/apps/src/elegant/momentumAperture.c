/*************************************************************************\
* Copyright (c) 2006 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2006 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: momentumAperture.c
 * purpose: Do tracking to find momentum aperture starting from the end of each element.
 * Ref: M. Belgrounne et al. PAC203, 896.
 *
 * Michael Borland, 2006
 */

#include "mdb.h"
#include "track.h"
#include "momentumAperture.h"

static SDDS_DATASET SDDSma;
static double momentumOffsetValue = 0;
static long fireOnPass = 1;

static void momentumOffsetFunction(double **coord, long np, long pass, double *pCentral)
{
  MALIGN mal;

  if (pass==fireOnPass) {
    mal.dxp = mal.dyp = mal.dz = mal.dt = mal.de = 0;
    mal.dx = x_initial;
    mal.dy = y_initial;
    mal.dp = momentumOffsetValue;
    offset_beam(coord, np, &mal, *pCentral);
  }
}

void setupMomentumApertureSearch(
                                 NAMELIST_TEXT *nltext,
                                 RUN *run,
                                 VARY *control
                                 )
{
  char description[200];
  
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&momentum_aperture, nltext);
  if (echoNamelists) print_namelist(stdout, &momentum_aperture);

  if (run->concat_order!=0)
    bomb("at present, momentum_aperture is incompatible with concatenation", NULL);
  
  /* check for data errors */
  if (!output)
    bomb("no output filename specified", NULL);
  if (delta_negative_limit>=0)
    bomb("delta_negative_limit >= 0", NULL);
  if (delta_positive_limit<=0) 
    bomb("delta_positive_limit <= 0", NULL);
  if (delta_step_size<=0)
    bomb("delta_step_size <= 0", NULL);
  if (fabs(delta_negative_limit)<=delta_step_size/2)
    bomb("|delta_negative_limit| <= delta_step_size/2", NULL);
  if (delta_positive_limit<=delta_step_size/2)
    bomb("delta_positive_limit <= delta_step_size/2", NULL);
  if (splits<0)
    bomb("splits < 0", NULL);
  if (s_start>=s_end)
    bomb("s_start >= s_end", NULL);
  if (include_name_pattern && has_wildcards(include_name_pattern) && strchr(include_name_pattern, '-'))
    include_name_pattern = expand_ranges(include_name_pattern);
  if (include_type_pattern && has_wildcards(include_type_pattern) && strchr(include_type_pattern, '-'))
    include_type_pattern = expand_ranges(include_type_pattern);
  if (skip_elements<0)
    bomb("skip_elements < 0", NULL);
  if (process_elements<=0)
    bomb("process_elements <= 0", NULL);
  
  output = compose_filename(output, run->rootname);
  sprintf(description, "Momentum aperture search");
#if SDDS_MPI_IO
  SDDS_MPI_Setup(&SDDSma, 1, n_processors, myid, MPI_COMM_WORLD, 1);
  if (!SDDS_Parallel_InitializeOutput(&SDDSma, description, "momentum aperture",  output)) {
#else
  if (!SDDS_InitializeOutput(&SDDSma, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
#endif


    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaPositiveFound", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaPositive", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "lostOnPassPositive", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "sLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "xLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "yLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaNegativeFound", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaNegative", "$gd$R$bneg$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "lostOnPassNegative", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "sLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "xLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "yLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSma, "deltaLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0){
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
#if !SDDS_MPI_IO
  /* In the version with parallel IO, the layout will be written later */
  if(!SDDS_SaveLayout(&SDDSma) || !SDDS_WriteLayout(&SDDSma)){
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
#endif
  
}

void finishMomentumApertureSearch()
{
  if (SDDS_IsActive(&SDDSma) && !SDDS_Terminate(&SDDSma)) {
    SDDS_SetError("Problem terminating SDDS output (finishMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}


long doMomentumApertureSearch(
                              RUN *run,
                              VARY *control,
                              ERRORVAL *errcon,
                              LINE_LIST *beamline,
                              double *startingCoord
                              )
{    
  double **coord;
  long nElem, iElem;
  double deltaInterval, pCentral, deltaStart;
  ELEMENT_LIST *elem, *elem0;
  long **lostOnPass, **loserFound, **survivorFound, *lostOnPass0, side;
  double deltaLimit[2], deltaLimit1, **deltaWhenLost, delta;
  double **xLost, **yLost, **deltaSurvived, **sLost, deltaLost;
  double *sStart;
  char **ElementName;
  long code;
  long processElements, skipElements, deltaSign, split;
#if USE_MPI
  long total_iElem = 0;  /* A global counter for job distribution */
  notSinglePart = 0;
#endif
#if defined(DEBUG)
  FILE *fpdeb = NULL;
  char s[1000];

  sprintf(s, "%s-Debug", output);
  fprintf(stderr, "Debug file %s\n", s);
  
  if (!(fpdeb = fopen(s, "w")))
    bomb("unable to open debug file for momentum aperture scan", NULL);
  fprintf(fpdeb, "SDDS1\n");
  fprintf(fpdeb, "&parameter name=ElementName type=string &end\n");
  fprintf(fpdeb, "&parameter name=Side type=string &end\n");
  fprintf(fpdeb, "&column name=delta type=double &end\n");
  fprintf(fpdeb, "&column name=Lost type=short &end\n");
  fprintf(fpdeb, "&column name=Loss Pass type=short &end\n");
  fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
  fflush(fpdeb);
#endif

  if (control->n_passes==1)
    fireOnPass = 0;
  else
    fireOnPass = 1;
  
  /* determine how many elements will be tracked */
  elem = &(beamline->elem);
  elem0 = NULL;
  nElem = 0;
  while (elem) {
    if (elem->end_pos>=s_start && elem->end_pos<=s_end &&
        (!include_name_pattern || wild_match(elem->name, include_name_pattern)) &&
        (!include_type_pattern || wild_match(entity_name[elem->type], include_type_pattern)) ) {
      if (!elem0)
	elem0 = elem;
      nElem++;
    } else if (elem->end_pos>s_end)
      break;
    elem = elem->succ;
  }
#if !USE_MPI
  if (nElem==0) 
    SDDS_Bomb("no elements found between s_start and s_end for momentum aperture computation");
#else
  if ((nElem<n_processors) && (myid==0)) {
    printf("Warning: The number of elements should be larger than the number of processors to avoid wasting resource.\nThe number of elements is %ld. The number of processors is %ld.\n", nElem, n_processors);
  }    
  if (verbosity) {
    verbosity = 0;
  if (myid == 0)
    printf ("Warning: In parallel version, no intermediate information will be provided\n");
  }
#endif
  /* allocate arrays for tracking */
  coord = (double**)czarray_2d(sizeof(**coord), 1, 7);

  /* allocate arrays for storing data for negative and positive momentum limits for each element */
  lostOnPass = (long**)czarray_2d(sizeof(**lostOnPass), 2, nElem);
  loserFound = (long**)czarray_2d(sizeof(**loserFound), 2, nElem);
  survivorFound = (long**)czarray_2d(sizeof(**survivorFound), 2, nElem);
  deltaSurvived = (double**)czarray_2d(sizeof(**deltaSurvived), 2, nElem);
  xLost = (double**)czarray_2d(sizeof(**xLost), 2, nElem);
  yLost = (double**)czarray_2d(sizeof(**yLost), 2, nElem);
  deltaWhenLost = (double**)czarray_2d(sizeof(**deltaWhenLost), 2, nElem);
  sLost = (double**)czarray_2d(sizeof(**sLost), 2, nElem);
  sStart = (double*)tmalloc(sizeof(*sStart)*nElem);
  ElementName = (char**)tmalloc(sizeof(*ElementName)*nElem);

  deltaLimit[0] = delta_negative_limit;
  deltaLimit[1] = delta_positive_limit;

  /* need to do this because do_tracking() in principle may realloc this pointer */
  lostOnPass0 = tmalloc(sizeof(*lostOnPass0)*1);
  
  elem = elem0;
  iElem = 0;
  processElements = process_elements;
  skipElements = skip_elements;

  if (fiducialize) {
    if (startingCoord)
      memcpy(coord[0], startingCoord, sizeof(double)*6);
    else
      memset(coord[0], 0, sizeof(**coord)*6);
    coord[0][6] = 1;
    pCentral = run->p_central;
    if (verbosity>1) 
      fprintf(stdout, "Tracking fiducial particle\n");
    code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                       NULL, NULL, NULL, NULL, run, control->i_step, 
                       FIRST_BEAM_IS_FIDUCIAL+(verbosity>1?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT, 1, 0, NULL, NULL, NULL, lostOnPass0, NULL);
    if (!code) {
      fprintf(stdout, "Fiducial particle lost. Don't know what to do.\n");
      exit(1);
    }
  }

  while (elem && processElements>0) {
    if ((!include_name_pattern || wild_match(elem->name, include_name_pattern)) &&
        (!include_type_pattern || wild_match(entity_name[elem->type], include_type_pattern))) {
      if (elem->end_pos>s_end) 
        break;
      if (skipElements>0) {
        skipElements --;
        elem = elem->succ;
        continue;
      }
#if USE_MPI
      if (myid == total_iElem%n_processors) {
#endif
      if (verbosity>0) {
        fprintf(stdout, "Searching for energy aperture for %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
        fflush(stdout);
      }
      ElementName[iElem] = elem->name;
      sStart[iElem] = elem->end_pos;
      for (side=0; side<2; side++) {
        deltaStart = 0;
        deltaSign = side==0 ? -1 : 1;
        lostOnPass[side][iElem] = -1;
        loserFound[side][iElem] = survivorFound[side][iElem] = 0;
        xLost[side][iElem] = yLost[side][iElem] = 
          deltaWhenLost[side][iElem] = sLost[side][iElem] = 
            deltaSurvived[side][iElem] =  0;
        deltaLost = deltaSign*DBL_MAX/2;
        deltaInterval = delta_step_size*deltaSign;
        
        if (verbosity>1) {
          fprintf(stdout, " Searching for %s side from 0 toward %e with interval %e\n", side==0?"negative":"positive",
                  deltaLimit[side], delta_step_size);
          fflush(stdout);
        }

        deltaLimit1 = deltaLimit[side];
        for (split=0; split<=splits; split++) {
          delta = deltaStart;
          
#if defined(DEBUG)
          fprintf(fpdeb, "%s\n%s\n",
                  elem->name, side==0?"negative":"positive");
          fflush(fpdeb);
#endif

          while (fabs(delta) <= fabs(deltaLimit1)) {
            setTrackingWedgeFunction(momentumOffsetFunction, 
                                     elem->succ?elem->succ:elem0); 
            momentumOffsetValue = delta;
            if (startingCoord)
              memcpy(coord[0], startingCoord, sizeof(double)*6);
            else
              memset(coord[0], 0, sizeof(**coord)*6);
            coord[0][6] = 1;
            pCentral = run->p_central;
            if (verbosity>3) {
              fprintf(stdout, "  Tracking with delta0 = %e (%e, %e, %e, %e, %e, %e), pCentral=%e\n", 
                      delta, coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5],
                      pCentral);
              fflush(stdout);
            }
            lostOnPass0[0] = -1;
            if (!fiducialize) {
              delete_phase_references();
              reset_special_elements(beamline, 1);
            }
            code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                               NULL, NULL, NULL, NULL, run, control->i_step, 
                               SILENT_RUNNING+INHIBIT_FILE_OUTPUT, control->n_passes, 0, NULL, NULL, NULL, lostOnPass0, NULL);
            if (!code) {
              /* particle lost */
              if (verbosity>3) {
                long i;
                fprintf(stdout, "  Particle lost with delta0 = %e at s = %e\n", delta, coord[0][4]);
                if (verbosity>4)
                  for (i=0; i<6; i++)
                    fprintf(stdout, "   coord[%ld] = %e\n", i, coord[0][i]);
                fflush(stdout);
              }
              lostOnPass[side][iElem] = lostOnPass0[0];
              xLost[side][iElem] = coord[0][0];
              yLost[side][iElem] = coord[0][2];
              sLost[side][iElem] = coord[0][4];
              deltaLost = delta;
              deltaWhenLost[side][iElem] = (coord[0][5]-pCentral)/pCentral;
              loserFound[side][iElem] = 1;
              break;
            } else {
              if (verbosity>2)
                fprintf(stdout, "  Particle survived with delta0 = %e\n", delta);
              deltaSurvived[side][iElem] = delta;
              survivorFound[side][iElem] = 1;
            }
            delta += deltaInterval;
          }
          if (split==0) {
            if (!survivorFound[side][iElem]) {
              fprintf(stdout, "Error: No survivor found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
              exit(1);
            }
            if (!loserFound[side][iElem]) {
              if (!soft_failure) {
                fprintf(stdout, "Error: No loss found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
                exit(1);
              } else {
                loserFound[side][iElem] = 1;
                split = splits;
              }
            }
          }
          deltaStart = deltaSurvived[side][iElem] - steps_back*deltaInterval;
          deltaInterval /= split_step_divisor;
          deltaStart += deltaInterval;
          deltaLimit1 = deltaLost;
          if ((deltaStart<0 && deltaSign==1) || (deltaStart>0 && deltaSign==-1))
            deltaStart = 0;
        }
#if defined(DEBUG)
        fprintf(fpdeb, "\n");
        fflush(fpdeb);
#endif
        if (verbosity>0) {
          fprintf(stdout, "Energy aperture for %s #%ld at s=%em is %e\n", elem->name, elem->occurence, elem->end_pos,
                  deltaSurvived[side][iElem]);
          fflush(stdout);
        }
      }
      
      iElem++;
#if USE_MPI
      }
      total_iElem++;
#endif
      processElements --;
    }
    elem = elem->succ;
  } 

#if SDDS_MPI_IO
  /* Open file here for parallel IO */
  if (!SDDS_LayoutWritten(&SDDSma)) {
    if (!SDDS_MPI_File_Open(SDDSma.MPI_dataset, SDDSma.layout.filename, SDDS_MPI_WRITE_ONLY)) 
      SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDSma.MPI_dataset->MPI_file);
    if (!SDDS_MPI_WriteLayout(&SDDSma))  
      SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDSma.MPI_dataset->MPI_file);
  }
#endif                                                                                       
  if (!SDDS_StartPage(&SDDSma, iElem) ||
      !SDDS_SetParameters(&SDDSma, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step",
                          control->i_step, NULL) ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, iElem, "ElementName") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, iElem, "s") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[1], iElem, "deltaPositiveFound") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[1], iElem, "deltaPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[1], iElem, "lostOnPassPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[1], iElem, "sLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[1], iElem, "xLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[1], iElem, "yLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[1], iElem, "deltaLostPositive") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], iElem, "deltaNegativeFound") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], iElem, "deltaNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], iElem, "lostOnPassNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], iElem, "sLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], iElem, "xLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], iElem, "yLostNegative") ||
      !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[0], iElem, "deltaLostNegative") ||
#if !SDDS_MPI_IO
      !SDDS_WritePage(&SDDSma)) {
#else
      !SDDS_MPI_WritePage(&SDDSma)) {
#endif
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDSma);
      

  free_czarray_2d((void**)coord, 1, 7);
  free_czarray_2d((void**)lostOnPass, 2, nElem);
  free_czarray_2d((void**)loserFound, 2, nElem);
  free_czarray_2d((void**)survivorFound, 2, nElem);
  free_czarray_2d((void**)deltaSurvived, 2, nElem);
  free_czarray_2d((void**)xLost, 2, nElem);
  free_czarray_2d((void**)yLost, 2, nElem);
  free_czarray_2d((void**)deltaWhenLost, 2, nElem);
  free_czarray_2d((void**)sLost, 2, nElem);
  free(sStart);
  free(ElementName);
  free(lostOnPass0);
  
  return 1;
}


