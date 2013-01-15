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
static double **turnByTurnCoord = NULL;
static long turnsStored = 0;

#include "fftpackC.h"
long determineTunesFromTrackingData(double *tune, double **turnByTurnCoord, long turns, double delta);

static void momentumOffsetFunction(double **coord, long np, long pass, double *pCentral)
{
  MALIGN mal;
  
  if (pass==fireOnPass) {
    mal.dxp = mal.dyp = mal.dz = mal.dt = mal.de = 0;
    mal.dx = x_initial;
    mal.dy = y_initial;
    mal.dp = momentumOffsetValue;
    offset_beam(coord, np, &mal, *pCentral);
    turnsStored = 0;
  }
  if (turnByTurnCoord) {
    /* N.B.: the arrays are ordered differently than normal, e.g.,
     * turnByTurnCoord[0] is the data for x, not particle/turn 0
     */
    turnByTurnCoord[0][turnsStored] = coord[0][0]; 
    turnByTurnCoord[1][turnsStored] = coord[0][1]; 
    turnByTurnCoord[2][turnsStored] = coord[0][2]; 
    turnByTurnCoord[3][turnsStored] = coord[0][3]; 
    turnByTurnCoord[4][turnsStored] = coord[0][5]; 
    turnsStored++;
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
  if (processNamelist(&momentum_aperture, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &momentum_aperture);

  if (run->concat_order!=0)
    bombElegant("at present, momentum_aperture is incompatible with concatenation", NULL);
  
  /* check for data errors */
  if (!output)
    bombElegant("no output filename specified", NULL);
  if (delta_negative_limit>=0)
    bombElegant("delta_negative_limit >= 0", NULL);
  if (delta_positive_limit<=0) 
    bombElegant("delta_positive_limit <= 0", NULL);
  if (delta_negative_start>0)
    bombElegant("delta_negative_start > 0", NULL);
  if (delta_positive_start<0) 
    bombElegant("delta_positive_start < 0", NULL);
  if (delta_step_size<=0)
    bombElegant("delta_step_size <= 0", NULL);
  if (fabs(delta_negative_limit-delta_negative_start)<=delta_step_size/2)
    bombElegant("|delta_negative_limit-delta_negative_start| <= delta_step_size/2", NULL);
  if (delta_positive_limit-delta_positive_start<=delta_step_size/2)
    bombElegant("delta_positive_limit-delta_positive_start <= delta_step_size/2", NULL);
  if (splits<0)
    bombElegant("splits < 0", NULL);
  if (s_start>=s_end)
    bombElegant("s_start >= s_end", NULL);
  if (include_name_pattern && has_wildcards(include_name_pattern) && strchr(include_name_pattern, '-'))
    include_name_pattern = expand_ranges(include_name_pattern);
  if (include_type_pattern && has_wildcards(include_type_pattern) && strchr(include_type_pattern, '-'))
    include_type_pattern = expand_ranges(include_type_pattern);
  if (skip_elements<0)
    bombElegant("skip_elements < 0", NULL);
  if (process_elements<=0)
    bombElegant("process_elements <= 0", NULL);
  
  output = compose_filename(output, run->rootname);
  sprintf(description, "Momentum aperture search");
#if SDDS_MPI_IO
  SDDS_MPI_Setup(&SDDSma, 1, n_processors, myid, MPI_COMM_WORLD, 1);
  if (!SDDS_Parallel_InitializeOutput(&SDDSma, description, "momentum aperture",  output)) {
#else
  if (!SDDS_InitializeOutput(&SDDSma, SDDS_BINARY, 1, description, "momentum aperture",  output)) {
#endif


    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }

  if (output_mode==0) {
    if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaPositiveFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaPositive", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPassPositive", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLostPositive", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLostPositive", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaNegativeFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaNegative", "$gd$R$bneg$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPassNegative", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLostNegative", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLostNegative", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0){
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  } else {
    if (SDDS_DefineColumn(&SDDSma, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "s", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "direction", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaFound", NULL, NULL, NULL, NULL, SDDS_SHORT, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "delta", "$gd$R$bpos$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "lostOnPass", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "sLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "xLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "yLost", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "deltaLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuxLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDSma, "nuyLost", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(&SDDSma, "Step", NULL, NULL, NULL, NULL, SDDS_LONG, NULL)<0){
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
  } 
#if !SDDS_MPI_IO
  /* In the version with parallel IO, the layout will be written later */
  if(!SDDS_SaveLayout(&SDDSma) || !SDDS_WriteLayout(&SDDSma)){
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
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
  long side;
  double **lostParticles;	
  short **loserFound, *direction=NULL, **survivorFound;
  int32_t **lostOnPass, *ElementOccurence;
  double deltaLimit1[2], deltaLimit, **deltaWhenLost, delta;
  double deltaStart1[2];
  double **xLost, **yLost, **deltaSurvived, **sLost, deltaLost;
  double **xTuneSurvived, **yTuneSurvived;
  double *sStart;
  double nominalTune[2], tune[2];
  char **ElementName, **ElementType;
  long code, outputRow, jobCounter;
  long processElements, skipElements, deltaSign, split, slot;
  char s[1000];
#if USE_MPI
  FILE *fpdMpi = NULL;
  notSinglePart = 0;

#if defined(DEBUG)
  sprintf(s, "%s-debug-%02d", output, myid);
  if (!(fpdMpi = fopen(s, "w")))
    fpdMpi = NULL;
  fprintf(fpdMpi, "ID = %d\n", myid);
  fflush(fpdMpi);
#endif
#else
#if defined(DEBUG)
  FILE *fpdeb = NULL;

  sprintf(s, "%s-Debug", output);
  fprintf(stderr, "Debug file %s\n", s);
  
  if (!(fpdeb = fopen(s, "w")))
    bombElegant("unable to open debug file for momentum aperture scan", NULL);
  fprintf(fpdeb, "SDDS1\n");
  fprintf(fpdeb, "&parameter name=ElementName type=string &end\n");
  fprintf(fpdeb, "&parameter name=Side type=string &end\n");
  fprintf(fpdeb, "&column name=delta type=double &end\n");
  fprintf(fpdeb, "&column name=Lost type=short &end\n");
  fprintf(fpdeb, "&column name=Loss Pass type=short &end\n");
  fprintf(fpdeb, "&data mode=ascii no_row_counts=1 &end\n");
  fflush(fpdeb);
#endif
#endif

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
  if ((((output_mode==0?1:2)*nElem)<n_processors) && (myid==0)) {
    printf("Warning: The number of elements should be larger than the number of processors to avoid wasting resource.\nThe number of elements is %ld. The number of processors is %d.\n", (output_mode==0?1:2)*nElem, n_processors);
    if (!output_mode) 
      printf("Tip: You can utilize more processors efficiently by setting output_mode=1\n");
  }    
  if (verbosity) {
    verbosity = 0;
    if (myid == 0)
      printf ("Warning: In parallel version, limited intermediate information will be provided\n");
  }
#endif
  /* allocate arrays for tracking */
  coord = (double**)czarray_2d(sizeof(**coord), 1, 7);

  /* allocate arrays for storing data for negative and positive momentum limits for each element */
  lostOnPass = (int32_t**)czarray_2d(sizeof(**lostOnPass), (output_mode?1:2), (output_mode?2:1)*nElem);
  loserFound = (short**)czarray_2d(sizeof(**loserFound), (output_mode?1:2), (output_mode?2:1)*nElem);
  survivorFound = (short**)czarray_2d(sizeof(**survivorFound), (output_mode?1:2), (output_mode?2:1)*nElem);
  deltaSurvived = (double**)czarray_2d(sizeof(**deltaSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  xTuneSurvived = (double**)czarray_2d(sizeof(**xTuneSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  yTuneSurvived = (double**)czarray_2d(sizeof(**yTuneSurvived), (output_mode?1:2), (output_mode?2:1)*nElem);
  xLost = (double**)czarray_2d(sizeof(**xLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  yLost = (double**)czarray_2d(sizeof(**yLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  deltaWhenLost = (double**)czarray_2d(sizeof(**deltaWhenLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  sLost = (double**)czarray_2d(sizeof(**sLost), (output_mode?1:2), (output_mode?2:1)*nElem);
  sStart = (double*)tmalloc(sizeof(*sStart)*(output_mode?2:1)*nElem);
  ElementName = (char**)tmalloc(sizeof(*ElementName)*(output_mode?2:1)*nElem);
  ElementType = (char**)tmalloc(sizeof(*ElementType)*(output_mode?2:1)*nElem);
  ElementOccurence = (int32_t*)tmalloc(sizeof(*ElementOccurence)*(output_mode?2:1)*nElem);
  if (output_mode)
    direction = (short*)tmalloc(sizeof(*direction)*2*nElem);
  deltaLimit1[0] = delta_negative_limit;
  deltaLimit1[1] = delta_positive_limit;
  deltaStart1[0] = delta_negative_start;
  deltaStart1[1] = delta_positive_start;

  if (control->n_passes==1)
    fireOnPass = 0;
  else
    fireOnPass = 1;
  turnByTurnCoord = (double**)czarray_2d(sizeof(double), 5, control->n_passes);
  
  /* need to do this because do_tracking() in principle may realloc this pointer */
  lostParticles = (double**)czarray_2d(sizeof(double),1, 8);	 
 
  elem = elem0;
  iElem = 0;
  processElements = process_elements;
  skipElements = skip_elements;

  if (fiducialize || forbid_resonance_crossing) {
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
                       FIRST_BEAM_IS_FIDUCIAL+(verbosity>1?0:SILENT_RUNNING)+INHIBIT_FILE_OUTPUT, 1, 0, NULL, NULL, NULL, lostParticles, NULL);
    if (!code) {
      fprintf(stdout, "Fiducial particle lost. Don't know what to do.\n");
      exitElegant(1);
    }
  }

  outputRow = -1;
  jobCounter = -1;

#if USE_MPI
  verbosity = 0;  
#endif

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
      if (output_mode==0) {
#if USE_MPI
        jobCounter++;
        if (myid != jobCounter%n_processors) {
          elem = elem->succ;
          continue;
        }
#endif
        outputRow++;
      }
#if !USE_MPI
      if (verbosity>0) {
        fprintf(stdout, "Searching for energy aperture for %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
        fflush(stdout);
      }
#endif
      for (side=0; side<2; side++) {
        if (output_mode==1) {
#if USE_MPI
          jobCounter++;
          if (myid!=jobCounter%n_processors)
            continue;
#endif
          outputRow++;
          slot = 0;
          direction[outputRow] = (side==0?-1:1);
        } else
          slot = side;

#if USE_MPI
        if (myid==0) {
          if (output_mode==1) {
            sprintf(s, "About %.3g%% done: ", (jobCounter*50.0)/nElem);
            report_stats(stdout, s);
            fflush(stdout);
          } else if (side==0) {
            sprintf(s, "About %.3g%% done: ", (jobCounter*100.0)/nElem);
            report_stats(stdout, s);
            fflush(stdout);
          }
        }
#else
	if (output_mode==1) {
	  sprintf(s, "About %.3g%% done: ", (outputRow*50.0)/nElem);
          report_stats(stdout, s);
          fflush(stdout);
	} else if (side==0) {
	  sprintf(s, "About %.3g%% done: ", (outputRow*100.0)/nElem);
          report_stats(stdout, s);
          fflush(stdout);
        }
        
#endif

        ElementName[outputRow] = elem->name;
        ElementType[outputRow] = entity_name[elem->type];
        ElementOccurence[outputRow] = elem->occurence;
        sStart[outputRow] = elem->end_pos;
        deltaStart = deltaStart1[side];
        deltaSign = side==0 ? -1 : 1;
        lostOnPass[slot][outputRow] = -1;
        loserFound[slot][outputRow] = survivorFound[slot][outputRow] = 0;
        xLost[slot][outputRow] = yLost[slot][outputRow] = 
          deltaWhenLost[slot][outputRow] = sLost[slot][outputRow] = 
            deltaSurvived[slot][outputRow] =  xTuneSurvived[slot][outputRow] =
              yTuneSurvived[slot][outputRow] = 0;
        deltaLost = deltaSign*DBL_MAX/2;
        deltaInterval = delta_step_size*deltaSign;
        
        if (verbosity>1) {
          fprintf(stdout, " Searching for %s side from %e toward %e with interval %e\n", side==0?"negative":"positive",
                  deltaStart, deltaLimit1[slot], delta_step_size);
          fflush(stdout);
        }

        if (forbid_resonance_crossing) {
          momentumOffsetValue = 0;
          setTrackingWedgeFunction(momentumOffsetFunction, 
                                   elem->succ?elem->succ:elem0); 
          if (startingCoord)
            memcpy(coord[0], startingCoord, sizeof(double)*6);
          else
            memset(coord[0], 0, sizeof(**coord)*6);
          coord[0][6] = 1;
          pCentral = run->p_central;
          code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                             NULL, NULL, NULL, NULL, run, control->i_step, 
                             (fiducialize?FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL:0)+SILENT_RUNNING+INHIBIT_FILE_OUTPUT, control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
          if (!code || !determineTunesFromTrackingData(nominalTune, turnByTurnCoord, turnsStored, 0.0)) {
            fprintf(stdout, "Fiducial particle tune is undefined.\n");
            exitElegant(1);
          }
          if (verbosity>3)
            fprintf(stdout, "  Nominal tunes: %e, %e\n", nominalTune[0], nominalTune[1]);
        }
        
        deltaLimit = deltaLimit1[slot];
        for (split=0; split<=splits; split++) {
          delta = deltaStart;
          
#if defined(DEBUG)
#if !USE_MPI
          fprintf(fpdeb, "%s\n%s\n",
                  elem->name, side==0?"negative":"positive");
          fflush(fpdeb);
#else
          fprintf(fpdMpi, "%s\n%s\n",
                  elem->name, side==0?"negative":"positive");
          fflush(fpdMpi);
#endif
#endif

          while (fabs(delta) <= fabs(deltaLimit)) {
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
            lostParticles[0][7] = -1;
            if (!fiducialize) {
              delete_phase_references();
              reset_special_elements(beamline, 1);
            }
            code = do_tracking(NULL, coord, 1, NULL, beamline, &pCentral, 
                               NULL, NULL, NULL, NULL, run, control->i_step, 
                               (fiducialize?FIDUCIAL_BEAM_SEEN+FIRST_BEAM_IS_FIDUCIAL:0)+SILENT_RUNNING+INHIBIT_FILE_OUTPUT, control->n_passes, 0, NULL, NULL, NULL, lostParticles, NULL);
            if (code && turnsStored>2) {
              if (!determineTunesFromTrackingData(tune, turnByTurnCoord, turnsStored, delta)) {
                if (forbid_resonance_crossing) {
                  if (verbosity>3)
                    fprintf(stdout, "   Resonance crossing detected (no tunes).  Particle lost\n");
                  code = 0; /* lost */
                }
              } else {
                if (verbosity>3) 
                  fprintf(stdout, "   Tunes: %e, %e\n", tune[0], tune[1]);
                if (forbid_resonance_crossing &&
                    ( (((long)(2*tune[0])) - ((long)(2*nominalTune[0])))!=0 ||
                     (((long)(2*tune[1])) - ((long)(2*nominalTune[1])))!=0) ) {
                  /* crossed integer or half integer */
                  if (verbosity>3)
                    fprintf(stdout, "   Resonance crossing detected (%e, %e -> %e, %e).  Particle lost\n",
                            nominalTune[0], nominalTune[1], tune[0], tune[1]);
                  code = 0;
                }
              }
            } else 
              tune[0] = tune[1] = -1;
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
              lostOnPass[slot][outputRow] = lostParticles[0][7];
              xLost[slot][outputRow] = coord[0][0];
              yLost[slot][outputRow] = coord[0][2];
              sLost[slot][outputRow] = coord[0][4];
              deltaLost = delta;
              deltaWhenLost[slot][outputRow] = (coord[0][5]-pCentral)/pCentral;
              loserFound[slot][outputRow] = 1;
              break;
            } else {
              if (verbosity>2)
                fprintf(stdout, "  Particle survived with delta0 = %e\n", delta);
              if (verbosity>3)
                fprintf(stdout, "     Final coordinates: %le, %le, %le, %le, %le, %le\n",
                        coord[0][0], coord[0][1], coord[0][2], coord[0][3], coord[0][4], coord[0][5]);
              deltaSurvived[slot][outputRow] = delta;
              survivorFound[slot][outputRow] = 1;
              xTuneSurvived[slot][outputRow] = tune[0];
              yTuneSurvived[slot][outputRow] = tune[1];
            }
            delta += deltaInterval;
          } /* delta search */
          if (split==0) {
            if (!survivorFound[slot][outputRow]) {
	      if (!soft_failure) {
		fprintf(stdout, "Error: No survivor found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
		exit(1);
	      }
              fprintf(stdout, "Warning: No survivor found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
	      deltaSurvived[slot][outputRow] = 0;
	      survivorFound[slot][outputRow] = 1;
	      split = splits;
            }
            if (!loserFound[slot][outputRow]) {
              if (!soft_failure) {
                fprintf(stdout, "Error: No loss found for initial scan for  %s #%ld at s=%em\n", elem->name, elem->occurence, elem->end_pos);
                exitElegant(1);
              }
	      loserFound[slot][outputRow] = 1;
	      split = splits;
            }
          }
          deltaStart = deltaSurvived[slot][outputRow] - steps_back*deltaInterval;
          deltaInterval /= split_step_divisor;
          deltaStart += deltaInterval;
          deltaLimit = deltaLost;
          if ((deltaStart<0 && deltaSign==1) || (deltaStart>0 && deltaSign==-1))
            deltaStart = 0;
        } /* split loop */
#if defined(DEBUG)
#if !USE_MPI
        fprintf(fpdeb, "\n");
        fflush(fpdeb);
#else
        fprintf(fpdMpi, "\n");
        fflush(fpdMpi);
#endif
#endif
        if (verbosity>0) {
          fprintf(stdout, "Energy aperture for %s #%ld at s=%em is %e\n", elem->name, elem->occurence, elem->end_pos,
                  deltaSurvived[slot][outputRow]);
          fflush(stdout);
        }
      } /* side loop */
      
      processElements --;
    } /* element loop */
    elem = elem->succ;
  } 

  outputRow++;

#if SDDS_MPI_IO
  /* Open file here for parallel IO */
  if (!SDDS_LayoutWritten(&SDDSma)) {
    if (!SDDS_MPI_File_Open(SDDSma.MPI_dataset, SDDSma.layout.filename, SDDS_MPI_WRITE_ONLY)) 
      SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDSma.MPI_dataset->MPI_file);
    if (!SDDS_MPI_WriteLayout(&SDDSma))  
      SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDSma.MPI_dataset->MPI_file);
  }
#endif                                                                                       
  if (!SDDS_StartPage(&SDDSma, outputRow) ||
      !SDDS_SetParameters(&SDDSma, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step",
                          control->i_step, NULL)) {
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if ((output_mode==0 && 
       (!SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, outputRow, "ElementName") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, outputRow, "s") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementType, outputRow, "ElementType") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementOccurence, outputRow, "ElementOccurence") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[1], outputRow, "deltaPositiveFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[1], outputRow, "deltaPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[1], outputRow, "lostOnPassPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[1], outputRow, "sLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[1], outputRow, "xLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[1], outputRow, "yLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[1], outputRow, "deltaLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], outputRow, "deltaNegativeFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], outputRow, "deltaNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], outputRow, "lostOnPassNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], outputRow, "sLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], outputRow, "xLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], outputRow, "yLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[0], outputRow, "deltaLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[0], outputRow, "nuxLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[1], outputRow, "nuxLostPositive") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[0], outputRow, "nuyLostNegative") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[1], outputRow, "nuyLostPositive"))) ||
      (output_mode==1 && 
       (!SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementName, outputRow, "ElementName") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sStart, outputRow, "s") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementType, outputRow, "ElementType") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, ElementOccurence, outputRow, "ElementOccurence") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, direction, outputRow, "direction") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, loserFound[0], outputRow, "deltaFound") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaSurvived[0], outputRow, "delta") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, lostOnPass[0], outputRow, "lostOnPass") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, sLost[0], outputRow, "sLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xLost[0], outputRow, "xLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yLost[0], outputRow, "yLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, deltaWhenLost[0], outputRow, "deltaLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, xTuneSurvived[0], outputRow, "nuxLost") ||
        !SDDS_SetColumn(&SDDSma, SDDS_SET_BY_NAME, yTuneSurvived[0], outputRow, "nuyLost")))) {
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
        
#if !SDDS_MPI_IO
  if (!SDDS_WritePage(&SDDSma)) {
#else
  if (!SDDS_MPI_WritePage(&SDDSma)) {
#endif
    SDDS_SetError("Problem writing SDDS table (doMomentumApertureSearch)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDSma);
      

  free_czarray_2d((void**)coord, 1, 7);
  free_czarray_2d((void**)lostOnPass, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)loserFound, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)survivorFound, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)deltaSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)xLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)yLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)deltaWhenLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)sLost, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)xTuneSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)yTuneSurvived, (output_mode?1:2), (output_mode?2:1)*nElem);
  free_czarray_2d((void**)turnByTurnCoord, 5, control->n_passes);
  turnByTurnCoord = NULL;
  
  free(sStart);
  free(ElementName);
  free(ElementType);
  free(ElementOccurence);
  return 1;
}


long determineTunesFromTrackingData(double *tune, double **turnByTurnCoord, long turns, double delta)
{
  double amplitude[4], frequency[4], phase[4], dummy;
  long i;

#if defined(DEBUG)
  static FILE *fpd = NULL;
  if (!fpd) {
    fpd = fopen("tbt.sdds", "w");
    fprintf(fpd, "SDDS1\n&column name=Pass type=long &end\n");
    fprintf(fpd, "&column name=x type=double units=m &end\n");
    fprintf(fpd, "&column name=xp type=double &end\n");
    fprintf(fpd, "&column name=y type=double units=m &end\n");
    fprintf(fpd, "&column name=yp type=double &end\n");
    fprintf(fpd, "&column name=delta type=double &end\n");
    fprintf(fpd, "&parameter name=delta0 type=double &end\n");
    fprintf(fpd, "&data mode=ascii &end\n");
  }
  fprintf(fpd, "%ld\n", turns);
  for (i=0; i<turns; i++) {
    fprintf(fpd, "%ld %e %e %e %e %e\n", i,
            turnByTurnCoord[0][i],
            turnByTurnCoord[1][i],
            turnByTurnCoord[2][i],
            turnByTurnCoord[3][i],
            turnByTurnCoord[4][i]);
  }
#endif

  if (PerformNAFF(&frequency[0], &amplitude[0], &phase[0], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[0], turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (x).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[1], &amplitude[1], &phase[1], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[1], turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,0, 0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (xp).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[2], &amplitude[2], &phase[2], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[2], turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (y).\n");
    return 0;
  }

  if (PerformNAFF(&frequency[3], &amplitude[3], &phase[3], 
		  &dummy, 0.0, 1.0, turnByTurnCoord[3], turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12, 0, 0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (yp).\n");
    return 0;
  }

  tune[0] = adjustTuneHalfPlane(frequency[0], phase[0], phase[1]);
  tune[1] = adjustTuneHalfPlane(frequency[2], phase[2], phase[3]);
  return 1;
}


