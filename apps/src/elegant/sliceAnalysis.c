/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/*
 * $Log: not supported by cvs2svn $
 * Revision 1.14  2010/02/05 17:48:50  soliday
 * Fixed minor issues related to compiler warnings.
 *
 * Revision 1.13  2010/02/04 15:17:27  borland
 * No longer use the bomb() routine.  Instead, use bombElegant(), which allows
 * better control of what happens when exiting.  Added "failed" semaphore option.
 * Switched from process_namelist() to processNamelist() for better response
 * to errors.
 * Includes Y. Wang's changes to parallelize shell and line beam types.
 *
 * Revision 1.12  2009/09/10 21:57:50  xiaoam
 * Add slice analysis to IBS simulation.
 *
 * Revision 1.11  2009/02/12 22:55:41  borland
 * Added ability to turn off echoing of namelists.
 *
 * Revision 1.10  2006/05/31 16:02:54  ywang25
 * The first release of Pelegant. It has passed a regression test of 100 cases.
 *
 * Revision 1.9  2005/01/27 17:39:39  borland
 * Updated calls to rpn routines.
 *
 * Revision 1.8  2002/08/14 20:23:49  soliday
 * Added Open License
 *
 * Revision 1.7  2002/07/01 16:16:02  borland
 * Slice beam properties at end of system now go into rpn for use
 * in optimization.
 *
 * Revision 1.6  2002/01/22 04:33:59  borland
 * Added alpha and beta computation to slice analysis.
 *
 * Revision 1.5  2002/01/07 20:38:25  borland
 * Trapazoid rule integratino is now optional for CSRCSBEND and CSRDRIFT.
 *
 * Revision 1.4  2002/01/07 03:42:14  borland
 * Removed purify and no-optimization from Makefile.
 * Slice analysis now uses the eta values computed from the whole beam rather
 * than the slice.  Slice analysis includes output of the eta's and the number
 * of particles in each slice.
 *
 * Revision 1.3  2002/01/04 21:38:39  borland
 * Added computation and output of the number of waists in x and y twiss parameters.
 * This is available for optimization as well.
 * CSR wake is now computed more accurately using trapazoid rule with an analytical
 * resulted use for the "divergent" term.
 *
 * Revision 1.2  2002/01/02 18:28:58  borland
 * Added slice analysis inside CSRCSBEND and CSRDRIFT elements.
 * Emittance units for tracked beam are now "m" instead of some
 * complicated sequence.
 *
 * Revision 1.1  2002/01/02 14:17:07  borland
 * First version of slice analysis code.
 *
 *
 */
#include "mdb.h"
#include "track.h"
#include "sliceAnalysis.h"
#if defined(__BORLANDC__)
#include <fdlibm.h>
#endif

static double tmp_safe_sqrt;
#define SAFE_SQRT(x) ((tmp_safe_sqrt=(x))<0?0.0:sqrt(tmp_safe_sqrt))

static SLICE_OUTPUT *sliceOutput = NULL;

long defineSliceParameters(SLICE_OUTPUT *sliceOutput, long slice);

void clearSliceAnalysis() 
{
  if (sliceOutput && sliceOutput->active && sliceOutput->filename) {
    SDDS_Terminate(&(sliceOutput->SDDSout));
    SDDS_ClearErrors();
    if (sliceOutput->enx) free(sliceOutput->enx);
    if (sliceOutput->eny) free(sliceOutput->eny);
    if (sliceOutput->ecnx) free(sliceOutput->ecnx);
    if (sliceOutput->ecny) free(sliceOutput->ecny);
    if (sliceOutput->betacx) free(sliceOutput->betacx);
    if (sliceOutput->betacy) free(sliceOutput->betacy);
    if (sliceOutput->alphacx) free(sliceOutput->alphacx);
    if (sliceOutput->alphacy) free(sliceOutput->alphacy);
    if (sliceOutput->charge) free(sliceOutput->charge);
    if (sliceOutput->particles) free(sliceOutput->particles);
    if (sliceOutput->duration) free(sliceOutput->duration);
    if (sliceOutput->Sdelta) free(sliceOutput->Sdelta);
    if (sliceOutput->Cx) free(sliceOutput->Cx);
    if (sliceOutput->Cy) free(sliceOutput->Cy);
    if (sliceOutput->Cxp) free(sliceOutput->Cxp);
    if (sliceOutput->Cyp) free(sliceOutput->Cyp);
    if (sliceOutput->Cdelta) free(sliceOutput->Cdelta);
    if (sliceOutput->Ct) free(sliceOutput->Ct);

    if (sliceOutput->enxIndex) free(sliceOutput->enxIndex);
    if (sliceOutput->enyIndex) free(sliceOutput->enyIndex);
    if (sliceOutput->ecnxIndex) free(sliceOutput->ecnxIndex);
    if (sliceOutput->ecnyIndex) free(sliceOutput->ecnyIndex);
    if (sliceOutput->betacxIndex) free(sliceOutput->betacxIndex);
    if (sliceOutput->betacyIndex) free(sliceOutput->betacyIndex);
    if (sliceOutput->alphacxIndex) free(sliceOutput->alphacxIndex);
    if (sliceOutput->alphacyIndex) free(sliceOutput->alphacyIndex);
    if (sliceOutput->chargeIndex) free(sliceOutput->chargeIndex);
    if (sliceOutput->particlesIndex) free(sliceOutput->particlesIndex);
    if (sliceOutput->durationIndex) free(sliceOutput->durationIndex);
    if (sliceOutput->SdeltaIndex) free(sliceOutput->SdeltaIndex);
    if (sliceOutput->CxIndex) free(sliceOutput->CxIndex);
    if (sliceOutput->CyIndex) free(sliceOutput->CyIndex);
    if (sliceOutput->CxpIndex) free(sliceOutput->CxpIndex);
    if (sliceOutput->CypIndex) free(sliceOutput->CypIndex);
    if (sliceOutput->CdeltaIndex) free(sliceOutput->CdeltaIndex);
    if (sliceOutput->CtIndex) free(sliceOutput->CtIndex);

    if (sliceOutput->enxMemNum) free(sliceOutput->enxMemNum);
    if (sliceOutput->enyMemNum) free(sliceOutput->enyMemNum);
    if (sliceOutput->ecnxMemNum) free(sliceOutput->ecnxMemNum);
    if (sliceOutput->ecnyMemNum) free(sliceOutput->ecnyMemNum);
    if (sliceOutput->betacxMemNum) free(sliceOutput->betacxMemNum);
    if (sliceOutput->betacyMemNum) free(sliceOutput->betacyMemNum);
    if (sliceOutput->alphacxMemNum) free(sliceOutput->alphacxMemNum);
    if (sliceOutput->alphacyMemNum) free(sliceOutput->alphacyMemNum);
    if (sliceOutput->chargeMemNum) free(sliceOutput->chargeMemNum);
    if (sliceOutput->particlesMemNum) free(sliceOutput->particlesMemNum);
    if (sliceOutput->durationMemNum) free(sliceOutput->durationMemNum);
    if (sliceOutput->SdeltaMemNum) free(sliceOutput->SdeltaMemNum);
    if (sliceOutput->CxMemNum) free(sliceOutput->CxMemNum);
    if (sliceOutput->CyMemNum) free(sliceOutput->CyMemNum);
    if (sliceOutput->CxpMemNum) free(sliceOutput->CxpMemNum);
    if (sliceOutput->CypMemNum) free(sliceOutput->CypMemNum);
    if (sliceOutput->CdeltaMemNum) free(sliceOutput->CdeltaMemNum);
    if (sliceOutput->CtMemNum) free(sliceOutput->CtMemNum);

    if (sliceOutput->sliceFound) free(sliceOutput->sliceFound);
  }
}

void setupSliceAnalysis(NAMELIST_TEXT *nltext, RUN *run, 
			OUTPUT_FILES *output_data)
{
  SDDS_DATASET *SDDSout;

#if USE_MPI
  /* This function will be parallelized in the future */
  fprintf(stdout, "Slice analysis is not supported in this version of parallel elegant.\n");
  MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
  MPI_Abort(MPI_COMM_WORLD, 1);	
#endif

  /* process namelist text */
  if (processNamelist(&slice_analysis, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &slice_analysis);

  sliceOutput = &(output_data->sliceAnalysis);

  if (!output || !strlen(output))
    bombElegant("provide a filename for slice_analysis output", NULL);
  if (n_slices<0) 
    bombElegant("invalid n_slice value", NULL);
  if (s_start<0 || s_start>s_end) 
    bombElegant("invalid s_start and/or s_end values", NULL);
  if ((s_start>0 || s_end!=DBL_MAX) && final_values_only)
    bombElegant("giving s_start or s_end with final_values_only=1 is not meaningful.", NULL);
  clearSliceAnalysis();
  
  sliceOutput->active = 1;
  sliceOutput->nSlices = n_slices;
  sliceOutput->finalValuesOnly = final_values_only;
  sliceOutput->sStart = s_start;
  sliceOutput->sEnd = s_end;

  if (!(sliceOutput->enx = malloc(sizeof(*(sliceOutput->enx))*(n_slices+2))) ||
      !(sliceOutput->eny = malloc(sizeof(*(sliceOutput->eny))*(n_slices+2))) ||
      !(sliceOutput->ecnx = malloc(sizeof(*(sliceOutput->ecnx))*(n_slices+2))) ||
      !(sliceOutput->ecny = malloc(sizeof(*(sliceOutput->ecny))*(n_slices+2))) ||
      !(sliceOutput->betacx = malloc(sizeof(*(sliceOutput->betacx))*(n_slices+2))) ||
      !(sliceOutput->betacy = malloc(sizeof(*(sliceOutput->betacy))*(n_slices+2))) ||
      !(sliceOutput->alphacx = malloc(sizeof(*(sliceOutput->alphacx))*(n_slices+2))) ||
      !(sliceOutput->alphacy = malloc(sizeof(*(sliceOutput->alphacy))*(n_slices+2))) ||
      !(sliceOutput->charge = malloc(sizeof(*(sliceOutput->charge))*(n_slices+2))) ||
      !(sliceOutput->particles = malloc(sizeof(*(sliceOutput->particles))*(n_slices+2))) ||
      !(sliceOutput->duration = malloc(sizeof(*(sliceOutput->duration))*(n_slices+2))) ||
      !(sliceOutput->Sdelta = malloc(sizeof(*(sliceOutput->Sdelta))*(n_slices+2))) ||
      !(sliceOutput->Cx = malloc(sizeof(*(sliceOutput->Cx))*(n_slices+2))) ||
      !(sliceOutput->Cy = malloc(sizeof(*(sliceOutput->Cy))*(n_slices+2))) ||
      !(sliceOutput->Cxp = malloc(sizeof(*(sliceOutput->Cxp))*(n_slices+2))) ||
      !(sliceOutput->Cyp = malloc(sizeof(*(sliceOutput->Cyp))*(n_slices+2))) ||
      !(sliceOutput->Cdelta = malloc(sizeof(*(sliceOutput->Cdelta))*(n_slices+2))) ||
      !(sliceOutput->Ct = malloc(sizeof(*(sliceOutput->Ct))*(n_slices+2))))
    bombElegant("memory allocation failure (setupSliceAnalysis)", NULL);

  if (!(sliceOutput->enxIndex = malloc(sizeof(*(sliceOutput->enxIndex))*(n_slices+2))) ||
      !(sliceOutput->enyIndex = malloc(sizeof(*(sliceOutput->enyIndex))*(n_slices+2))) ||
      !(sliceOutput->ecnxIndex = malloc(sizeof(*(sliceOutput->ecnxIndex))*(n_slices+2))) ||
      !(sliceOutput->ecnyIndex = malloc(sizeof(*(sliceOutput->ecnyIndex))*(n_slices+2))) ||
      !(sliceOutput->betacxIndex = malloc(sizeof(*(sliceOutput->betacxIndex))*(n_slices+2))) ||
      !(sliceOutput->betacyIndex = malloc(sizeof(*(sliceOutput->betacyIndex))*(n_slices+2))) ||
      !(sliceOutput->alphacxIndex = malloc(sizeof(*(sliceOutput->alphacxIndex))*(n_slices+2))) ||
      !(sliceOutput->alphacyIndex = malloc(sizeof(*(sliceOutput->alphacyIndex))*(n_slices+2))) ||
      !(sliceOutput->chargeIndex = malloc(sizeof(*(sliceOutput->chargeIndex))*(n_slices+2))) ||
      !(sliceOutput->particlesIndex = malloc(sizeof(*(sliceOutput->particlesIndex))*(n_slices+2))) ||
      !(sliceOutput->durationIndex = malloc(sizeof(*(sliceOutput->durationIndex))*(n_slices+2))) ||
      !(sliceOutput->SdeltaIndex = malloc(sizeof(*(sliceOutput->SdeltaIndex))*(n_slices+2))) ||
      !(sliceOutput->CxIndex = malloc(sizeof(*(sliceOutput->CxIndex))*(n_slices+2))) ||
      !(sliceOutput->CyIndex = malloc(sizeof(*(sliceOutput->CyIndex))*(n_slices+2))) ||
      !(sliceOutput->CxpIndex = malloc(sizeof(*(sliceOutput->CxpIndex))*(n_slices+2))) ||
      !(sliceOutput->CypIndex = malloc(sizeof(*(sliceOutput->CypIndex))*(n_slices+2))) ||
      !(sliceOutput->CdeltaIndex = malloc(sizeof(*(sliceOutput->CdeltaIndex))*(n_slices+2))) ||
      !(sliceOutput->CtIndex = malloc(sizeof(*(sliceOutput->CtIndex))*(n_slices+2))) ||
      !(sliceOutput->sliceFound = malloc(sizeof(*(sliceOutput->sliceFound))*(n_slices+2))))
    bombElegant("memory allocation failure (setupSliceAnalysis)", NULL);
  
  if (!(sliceOutput->enxMemNum = malloc(sizeof(*(sliceOutput->enxMemNum))*(n_slices+2))) ||
      !(sliceOutput->enyMemNum = malloc(sizeof(*(sliceOutput->enyMemNum))*(n_slices+2))) ||
      !(sliceOutput->ecnxMemNum = malloc(sizeof(*(sliceOutput->ecnxMemNum))*(n_slices+2))) ||
      !(sliceOutput->ecnyMemNum = malloc(sizeof(*(sliceOutput->ecnyMemNum))*(n_slices+2))) ||
      !(sliceOutput->betacxMemNum = malloc(sizeof(*(sliceOutput->betacxMemNum))*(n_slices+2))) ||
      !(sliceOutput->betacyMemNum = malloc(sizeof(*(sliceOutput->betacyMemNum))*(n_slices+2))) ||
      !(sliceOutput->alphacxMemNum = malloc(sizeof(*(sliceOutput->alphacxMemNum))*(n_slices+2))) ||
      !(sliceOutput->alphacyMemNum = malloc(sizeof(*(sliceOutput->alphacyMemNum))*(n_slices+2))) ||
      !(sliceOutput->chargeMemNum = malloc(sizeof(*(sliceOutput->chargeMemNum))*(n_slices+2))) ||
      !(sliceOutput->particlesMemNum = malloc(sizeof(*(sliceOutput->particlesMemNum))*(n_slices+2))) ||
      !(sliceOutput->durationMemNum = malloc(sizeof(*(sliceOutput->durationMemNum))*(n_slices+2))) ||
      !(sliceOutput->SdeltaMemNum = malloc(sizeof(*(sliceOutput->SdeltaMemNum))*(n_slices+2))) ||
      !(sliceOutput->CxMemNum = malloc(sizeof(*(sliceOutput->CxMemNum))*(n_slices+2))) ||
      !(sliceOutput->CyMemNum = malloc(sizeof(*(sliceOutput->CyMemNum))*(n_slices+2))) ||
      !(sliceOutput->CxpMemNum = malloc(sizeof(*(sliceOutput->CxpMemNum))*(n_slices+2))) ||
      !(sliceOutput->CypMemNum = malloc(sizeof(*(sliceOutput->CypMemNum))*(n_slices+2))) ||
      !(sliceOutput->CdeltaMemNum = malloc(sizeof(*(sliceOutput->CdeltaMemNum))*(n_slices+2))) ||
      !(sliceOutput->CtMemNum = malloc(sizeof(*(sliceOutput->CtMemNum))*(n_slices+2))))
    bombElegant("memory allocation failure (setupSliceAnalysis)", NULL);    

  sliceOutput->filename = compose_filename(output, run->rootname);
  SDDSout = &(sliceOutput->SDDSout);
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, sliceOutput->filename) ||
      !SDDS_DefineSimpleParameter(SDDSout, "Step", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ElementName", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleColumn(SDDSout, "s", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etax", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etaxp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etay", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "etayp", NULL, SDDS_DOUBLE) ) {
    fprintf(stderr, "Unable define SDDS parameter for file %s\n", sliceOutput->filename);
    fflush(stderr);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (n_slices<=1) {
    if (!defineSliceParameters(sliceOutput, 0)) {
      fprintf(stderr, "Unable define SDDS parameters for file %s\n", sliceOutput->filename);
      fflush(stderr);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else {
    long slice;
    /* slice 0 is the nominal (no slicing)
     * slice N+1 is the average over the slices 
     */
    for (slice=0; slice<=n_slices+1; slice++) {
      if (!defineSliceParameters(sliceOutput, slice)) {
	fprintf(stderr, "Unable define SDDS parameters for file %s\n", sliceOutput->filename);
	fflush(stderr);
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}


long defineSliceParameters(SLICE_OUTPUT *sliceOutput, long slice)
{
  SDDS_DATASET *SDDSout;
  char buffer[100], sliceNumString[20];
    
  SDDSout = &(sliceOutput->SDDSout);
  if (slice && sliceOutput->nSlices>1) {
    if (slice<=sliceOutput->nSlices) {
      sprintf(sliceNumString, "Slice%02ld", slice);
    } else {
      /* "slice" N+1 is the average over all slices */
      sprintf(sliceNumString, "Ave");
    }
  }
  else
    sliceNumString[0] = 0;

  sprintf(buffer, "enx%s", sliceNumString);
  if ((sliceOutput->enxIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->enxMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "eny%s", sliceNumString);
  if ((sliceOutput->enyIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->enyMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "ecnx%s", sliceNumString);
  if ((sliceOutput->ecnxIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->ecnxMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "ecny%s", sliceNumString);
  if ((sliceOutput->ecnyIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->ecnyMemNum[slice] = rpn_create_mem(buffer, 0);

  sprintf(buffer, "betacx%s", sliceNumString);
  if ((sliceOutput->betacxIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->betacxMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "betacy%s", sliceNumString);
  if ((sliceOutput->betacyIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->betacyMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "alphacx%s", sliceNumString);
  if ((sliceOutput->alphacxIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->alphacxMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "alphacy%s", sliceNumString);
  if ((sliceOutput->alphacyIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->alphacyMemNum[slice] = rpn_create_mem(buffer, 0);

  sprintf(buffer, "charge%s", sliceNumString);
  if ((sliceOutput->chargeIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "C", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->chargeMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "particles%s", sliceNumString);
  if ((sliceOutput->particlesIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->particlesMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "duration%s", sliceNumString);
  if ((sliceOutput->durationIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->durationMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Sdelta%s", sliceNumString);
  if ((sliceOutput->SdeltaIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->SdeltaMemNum[slice] = rpn_create_mem(buffer, 0);

  sprintf(buffer, "Cx%s", sliceNumString);
  if ((sliceOutput->CxIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CxMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Cy%s", sliceNumString);
  if ((sliceOutput->CyIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CyMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Cxp%s", sliceNumString);
  if ((sliceOutput->CxpIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CxpMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Cyp%s", sliceNumString);
  if ((sliceOutput->CypIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CypMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Cdelta%s", sliceNumString);
  if ((sliceOutput->CdeltaIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CdeltaMemNum[slice] = rpn_create_mem(buffer, 0);
  sprintf(buffer, "Ct%s", sliceNumString);
  if ((sliceOutput->CtIndex[slice] = 
       SDDS_DefineColumn(SDDSout, buffer, NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  sliceOutput->CtMemNum[slice] = rpn_create_mem(buffer, 0);

  return 1;
}


void performSliceAnalysisOutput(SLICE_OUTPUT *sliceOutput, double **particle, long particles, 
				long newPage, long step, double Po, double charge, 
				char *elementName, double elementPosition,
				long timeGiven)
{
  SDDS_DATASET *SDDSout;
  long slice;

  if (!sliceOutput || !sliceOutput->active) 
    return;
  if (sliceOutput->sStart<sliceOutput->sEnd &&
      (sliceOutput->sStart>elementPosition ||
       sliceOutput->sEnd<elementPosition))
    return;
  SDDSout = &(sliceOutput->SDDSout);
  if (newPage) {
    sliceOutput->rows = 0;
    if (!SDDS_StartPage(SDDSout, 1) ||
	!SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "Step", step,
			    NULL)) {
      fprintf(stderr, "Unable write data to file %s\n", sliceOutput->filename);
      fflush(stderr);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }

  performSliceAnalysis(sliceOutput, particle, particles, Po, charge, timeGiven);
  if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sliceOutput->rows,
			 "ElementName", elementName,
			 "s", elementPosition, 
                         "etax", sliceOutput->eta[0], 
                         "etaxp", sliceOutput->eta[1], 
                         "etay", sliceOutput->eta[2], 
                         "etayp", sliceOutput->eta[3], 
                         NULL)) {
      fprintf(stderr, "Unable write data to file %s\n", sliceOutput->filename);
      fflush(stderr);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (slice=0; slice<=sliceOutput->nSlices+1; slice++) {
    if (slice==1 && sliceOutput->nSlices==0)
      break;
    if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
			   sliceOutput->rows,
			   sliceOutput->enxIndex[slice], sliceOutput->enx[slice], 
			   sliceOutput->enyIndex[slice], sliceOutput->eny[slice], 
			   sliceOutput->ecnxIndex[slice], sliceOutput->ecnx[slice], 
			   sliceOutput->ecnyIndex[slice], sliceOutput->ecny[slice], 
			   sliceOutput->betacxIndex[slice], sliceOutput->betacx[slice], 
			   sliceOutput->betacyIndex[slice], sliceOutput->betacy[slice], 
			   sliceOutput->alphacxIndex[slice], sliceOutput->alphacx[slice], 
			   sliceOutput->alphacyIndex[slice], sliceOutput->alphacy[slice], 
			   sliceOutput->chargeIndex[slice], sliceOutput->charge[slice], 
			   sliceOutput->particlesIndex[slice], sliceOutput->particles[slice], 
			   sliceOutput->durationIndex[slice], sliceOutput->duration[slice], 
			   sliceOutput->SdeltaIndex[slice], sliceOutput->Sdelta[slice], 
			   sliceOutput->CxIndex[slice], sliceOutput->Cx[slice],
			   sliceOutput->CyIndex[slice], sliceOutput->Cy[slice],
			   sliceOutput->CxpIndex[slice], sliceOutput->Cxp[slice],
			   sliceOutput->CypIndex[slice], sliceOutput->Cyp[slice],
			   sliceOutput->CdeltaIndex[slice], sliceOutput->Cdelta[slice], 
			   sliceOutput->CtIndex[slice], sliceOutput->Ct[slice], 
			   -1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    rpn_store(sliceOutput->enx[slice], NULL, sliceOutput->enxMemNum[slice]);
    rpn_store(sliceOutput->eny[slice], NULL, sliceOutput->enyMemNum[slice]);
    rpn_store(sliceOutput->ecnx[slice], NULL, sliceOutput->ecnxMemNum[slice]);
    rpn_store(sliceOutput->ecny[slice], NULL, sliceOutput->ecnyMemNum[slice]);
    rpn_store(sliceOutput->betacx[slice], NULL, sliceOutput->betacxMemNum[slice]);
    rpn_store(sliceOutput->betacy[slice], NULL, sliceOutput->betacyMemNum[slice]);
    rpn_store(sliceOutput->alphacx[slice], NULL, sliceOutput->alphacxMemNum[slice]);
    rpn_store(sliceOutput->alphacy[slice], NULL, sliceOutput->alphacyMemNum[slice]);
    rpn_store(sliceOutput->charge[slice], NULL, sliceOutput->chargeMemNum[slice]);
    rpn_store(sliceOutput->particles[slice], NULL, sliceOutput->particlesMemNum[slice]);
    rpn_store(sliceOutput->duration[slice], NULL, sliceOutput->durationMemNum[slice]);
    rpn_store(sliceOutput->Sdelta[slice], NULL, sliceOutput->SdeltaMemNum[slice]);
    rpn_store(sliceOutput->Cx[slice], NULL, sliceOutput->CxMemNum[slice]);
    rpn_store(sliceOutput->Cxp[slice], NULL, sliceOutput->CxpMemNum[slice]);
    rpn_store(sliceOutput->Cy[slice], NULL, sliceOutput->CyMemNum[slice]);
    rpn_store(sliceOutput->Cyp[slice], NULL, sliceOutput->CypMemNum[slice]);
    rpn_store(sliceOutput->Cdelta[slice], NULL, sliceOutput->CdeltaMemNum[slice]);
    rpn_store(sliceOutput->Ct[slice], NULL, sliceOutput->CtMemNum[slice]);
  }
  sliceOutput->rows += 1;
  if (!SDDS_UpdatePage(SDDSout, FLUSH_TABLE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void performSliceAnalysis(SLICE_OUTPUT *sliceOutput, double **particle, long particles, 
                         double Po, double charge, long timeGiven)
{
  double *sSave=NULL;
  long i, slice, nSlices;
  double tMinAll, tMaxAll, tMin, tMax;
  long count, slicesFound=0;
  double aveCoord[6], S[6][6];
  
  if (!particles) {
    fprintf(stdout, "no particles left---can't compute slice analysis");
    fflush(stdout);
    /* fill in some dummy values */
    for (slice=0; slice<sliceOutput->nSlices+1; slice++) {
      sliceOutput->enx[slice] = sliceOutput->eny[slice] = 
        sliceOutput->ecnx[slice] = sliceOutput->ecny[slice] = 
	sliceOutput->betacx[slice] = sliceOutput->betacy[slice] = 
	sliceOutput->alphacx[slice] = sliceOutput->alphacy[slice] = 
	sliceOutput->Cx[slice] = sliceOutput->Cy[slice] =
	sliceOutput->Cxp[slice] = sliceOutput->Cyp[slice] = 
	sliceOutput->Cdelta[slice] = sliceOutput->duration[slice] = 
	sliceOutput->Sdelta[slice] = sliceOutput->charge[slice] = 
	sliceOutput->particles[slice] = DBL_MAX;
    }
    return;
  }

  if (!timeGiven) {
    /* compute the time coordinate and replace 's' with this value.
     * save s values in a new array
     */
    if (!(sSave=malloc(sizeof(*sSave)*particles)))
      SDDS_Bomb("memory allocation failure (performSliceAnalysis)");
    computeTimeCoordinates(sSave, Po, particle, particles);
    for (i=0; i<particles; i++)
      SWAP_DOUBLE(sSave[i], particle[i][4]);
  }

  /* compute normal values (over entire beam) */
  computeSliceMoments(aveCoord, S, particle, particles, 0, 0);
  sliceOutput->charge[0] = charge;  
  sliceOutput->particles[0] = particles;
  sliceOutput->Cx[0] = aveCoord[0];
  sliceOutput->Cxp[0] = aveCoord[1];
  sliceOutput->Cy[0] = aveCoord[2];
  sliceOutput->Cyp[0] = aveCoord[3];
  sliceOutput->Ct[0] = aveCoord[4];
  sliceOutput->Cdelta[0] = aveCoord[5];
  sliceOutput->Sdelta[0] = sqrt(S[5][5]);
  for (i=0; i<4; i++)
    sliceOutput->eta[i] = 0;
  if (S[5][5]) 
    for (i=0; i<4; i++)
      sliceOutput->eta[i] = S[i][5]/S[5][5];
  sliceOutput->enx[0] = Po*SAFE_SQRT(S[0][0]*S[1][1]-sqr(S[0][1]));
  sliceOutput->eny[0] = Po*SAFE_SQRT(S[2][2]*S[3][3]-sqr(S[2][3]));
  sliceOutput->ecnx[0] = Po*correctedEmittance(S, sliceOutput->eta, 0, 1,
					       &sliceOutput->betacx[0],
					       &sliceOutput->alphacx[0]);
  sliceOutput->ecny[0] = Po*correctedEmittance(S, sliceOutput->eta, 2, 3,
					       &sliceOutput->betacy[0],
					       &sliceOutput->alphacy[0]);
  /* find total bunch duration */
  tMaxAll = -(tMinAll = DBL_MAX);
  for (i=0; i<particles; i++) {
    if (tMinAll>particle[i][4])
      tMinAll = particle[i][4];
    if (tMaxAll<particle[i][4])
      tMaxAll = particle[i][4];
  }
  sliceOutput->duration[0] = tMaxAll - tMinAll;
  
  if (sliceOutput->nSlices>1) {
    /* compute values for each slice, plus average */
    nSlices = sliceOutput->nSlices;

    sliceOutput->enx[nSlices+1] = sliceOutput->eny[nSlices+1] =  0;
    sliceOutput->ecnx[nSlices+1] = sliceOutput->ecny[nSlices+1] =  0;
    sliceOutput->betacx[nSlices+1] = sliceOutput->betacy[nSlices+1] =  0;
    sliceOutput->alphacx[nSlices+1] = sliceOutput->alphacy[nSlices+1] =  0;
    sliceOutput->charge[nSlices+1] = sliceOutput->duration[nSlices+1] = 0;
    sliceOutput->particles[nSlices+1] = sliceOutput->Sdelta[nSlices+1] = 0;
    sliceOutput->Cx[nSlices+1] = sliceOutput->Cxp[nSlices+1] = 0;
    sliceOutput->Cy[nSlices+1] = sliceOutput->Cyp[nSlices+1] = 0;
    sliceOutput->Cdelta[nSlices+1] = sliceOutput->Ct[nSlices+1] = 0;

    for (slice=1; slice<=sliceOutput->nSlices; slice++) {
      /* find boundaries of slice in time */
      tMin = (slice-1)*(tMaxAll-tMinAll)/sliceOutput->nSlices + tMinAll;
      tMax = tMin + (tMaxAll-tMinAll)/sliceOutput->nSlices;

      sliceOutput->duration[slice] = tMax-tMin;
      
      count = computeSliceMoments(aveCoord, S, particle, particles, tMin, tMax);
      sliceOutput->particles[slice] = count;
      if (count<2)
        sliceOutput->sliceFound[slice] = 0;
      else {
	slicesFound++;
	sliceOutput->sliceFound[slice] = 1;

	sliceOutput->charge[slice] = (charge*count)/particles;
	sliceOutput->Cx[slice] = aveCoord[0];
	sliceOutput->Cxp[slice] = aveCoord[1];
	sliceOutput->Cy[slice] = aveCoord[2];
	sliceOutput->Cyp[slice] = aveCoord[3];
	sliceOutput->Ct[slice] = aveCoord[4];
	sliceOutput->Cdelta[slice] = aveCoord[5];
	sliceOutput->Sdelta[slice] = sqrt(S[5][5]);
        sliceOutput->enx[slice] = Po*SAFE_SQRT(S[0][0]*S[1][1]-sqr(S[0][1]));
        sliceOutput->eny[slice] = Po*SAFE_SQRT(S[2][2]*S[3][3]-sqr(S[2][3]));
        sliceOutput->ecnx[slice] = Po*correctedEmittance(S, sliceOutput->eta, 0, 1,
							 &sliceOutput->betacx[slice],
							 &sliceOutput->alphacx[slice]);
        sliceOutput->ecny[slice] = Po*correctedEmittance(S, sliceOutput->eta, 2, 3,
							 &sliceOutput->betacy[slice],
							 &sliceOutput->alphacy[slice]);
	sliceOutput->enx[nSlices+1] += sliceOutput->enx[slice];
	sliceOutput->eny[nSlices+1] += sliceOutput->eny[slice];
	sliceOutput->ecnx[nSlices+1] += sliceOutput->ecnx[slice];
	sliceOutput->ecny[nSlices+1] += sliceOutput->ecny[slice];
	sliceOutput->betacx[nSlices+1] += sliceOutput->betacx[slice];
	sliceOutput->betacy[nSlices+1] += sliceOutput->betacy[slice];
	sliceOutput->alphacx[nSlices+1] += sliceOutput->alphacx[slice];
	sliceOutput->alphacy[nSlices+1] += sliceOutput->alphacy[slice];
	sliceOutput->charge[nSlices+1] += sliceOutput->charge[slice];
	sliceOutput->particles[nSlices+1] += sliceOutput->particles[slice];
	sliceOutput->duration[nSlices+1] += sliceOutput->duration[slice];
	sliceOutput->Sdelta[nSlices+1] += sliceOutput->Sdelta[slice];
	sliceOutput->Cx[nSlices+1] += sliceOutput->Cx[slice];
	sliceOutput->Cxp[nSlices+1] += sliceOutput->Cxp[slice];
	sliceOutput->Cy[nSlices+1] += sliceOutput->Cy[slice];
	sliceOutput->Cyp[nSlices+1] += sliceOutput->Cyp[slice];
	sliceOutput->Cdelta[nSlices+1] += sliceOutput->Cdelta[slice];
	sliceOutput->Ct[nSlices+1] += sliceOutput->Ct[slice];
      }
    }
    
    if (!slicesFound)
      bombElegant("No valid slices found for slice property computation.", NULL);

    sliceOutput->enx[nSlices+1] /= slicesFound;
    sliceOutput->eny[nSlices+1] /= slicesFound;
    sliceOutput->ecnx[nSlices+1] /= slicesFound;
    sliceOutput->ecny[nSlices+1] /= slicesFound;
    sliceOutput->betacx[nSlices+1] /= slicesFound;
    sliceOutput->betacy[nSlices+1] /= slicesFound;
    sliceOutput->alphacx[nSlices+1] /= slicesFound;
    sliceOutput->alphacy[nSlices+1] /= slicesFound;
    sliceOutput->Sdelta[nSlices+1] /= slicesFound;
    sliceOutput->charge[nSlices+1] /= slicesFound;
    sliceOutput->particles[nSlices+1] /= slicesFound;
    sliceOutput->duration[nSlices+1] /= slicesFound;
    sliceOutput->Cx[nSlices+1] /= slicesFound;
    sliceOutput->Cxp[nSlices+1] /= slicesFound;
    sliceOutput->Cy[nSlices+1] /= slicesFound;
    sliceOutput->Cyp[nSlices+1] /= slicesFound;
    sliceOutput->Cdelta[nSlices+1] /= slicesFound;
    sliceOutput->Ct[nSlices+1] /= slicesFound;
  }

  if (!timeGiven) {
    for (i=0; i<particles; i++)
      SWAP_DOUBLE(sSave[i], particle[i][4]);
    free(sSave);
  }
}

double correctedEmittance(double S[6][6], double eta[4], long i1, long i2,
			  double *beta, double *alpha)
{
  double T1, T2, T3, ec2;
  T1 = S[i1][i1] + sqr(eta[i1])*S[5][5] - 2*eta[i1]*S[i1][5];
  T2 = S[i2][i2] + sqr(eta[i2])*S[5][5] - 2*eta[i2]*S[i2][5];
  T3 = S[i1][i2] - S[i2][5]*eta[i1] - S[i1][5]*eta[i2] + S[5][5]*eta[i1]*eta[i2];
  if ((ec2 = T1*T2-sqr(T3))<0)
    ec2 = 0;
  if (beta) {
    if (ec2)
      *beta = (S[i1][i1] + sqr(eta[i1])*S[5][5] - 2*eta[i1]*S[i1][5])/sqrt(ec2);
    else 
      *beta = DBL_MAX;
  }
  if (alpha) {
    if (ec2)
      *alpha = -(S[i1][i2] + eta[i1]*eta[i2]*S[5][5] -
		eta[i1]*S[i2][5] - eta[i2]*S[i1][5])/sqrt(ec2);
    else
      *alpha = DBL_MAX;
  }
  return sqrt(ec2);
}

