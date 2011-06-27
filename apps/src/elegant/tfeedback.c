/*************************************************************************\
* Copyright (c) 2003 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2003 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: tfeedback.c
 *
 * Michael Borland, 2003
 */
#include "mdb.h"
#include "track.h"

void transverseFeedbackPickup(TFBPICKUP *tfbp, double **part, long np, long pass)
{
  double position, output;
  long i, j;
  
  if (tfbp->initialized==0)
    initializeTransverseFeedbackPickup(tfbp);

  position = 0;
  if (np) {
    for (i=0; i<np; i++)
      position += part[i][tfbp->yPlane?2:0];
    position /= np;
  }
  if (tfbp->rmsNoise)
    position += gauss_rn_lim(0.0, tfbp->rmsNoise, 2, random_3);
#ifdef DEBUG
  fprintf(stdout, "TFBPICKUP: Putting value %e in slot %ld\n",
          position, pass%tfbp->filterLength);
#endif
  tfbp->data[pass%tfbp->filterLength] = position;
  if (pass<tfbp->filterLength) {
    tfbp->filterOutput = 0;
    return;
  }
  
  for (i=output=0, j=pass; i<tfbp->filterLength; i++, j--) {
#ifdef DEBUG
    fprintf(stdout, "TFBPICKUP: adding %le * %le  (data from slot %ld)\n",
            tfbp->a[i], tfbp->data[j%tfbp->filterLength], j%tfbp->filterLength);
#endif
    output += tfbp->a[i]*tfbp->data[j%tfbp->filterLength];
  }
#ifdef DEBUG
  fprintf(stdout, "TFBPICKUP: filter output is %e\n", output);
#endif
  tfbp->filterOutput = output;
}

void initializeTransverseFeedbackPickup(TFBPICKUP *tfbp)
{
  long i;
  double sum;

  if (tfbp->ID==NULL || !strlen(tfbp->ID))
    bombElegant("you must give an ID string for TFBPICKUP", NULL);

  for (i=sum=0; i<TFB_FILTER_LENGTH; i++)
    sum += tfbp->a[i];
  if (fabs(sum)>1e-6)
    fprintf(stdout, "Warning: sum of a[i] is nonzero for TFBPICKUP\n");

  for (i=TFB_FILTER_LENGTH-1; i>=0; i--) {
    if (tfbp->a[i]!=0)
      break;
  }
  if (i<0)
    bombElegant("All filter coefficients are zero for TFBPICKUP", NULL);
  tfbp->filterLength = i+1;

  if (strcmp(tfbp->plane, "x")==0 || strcmp(tfbp->plane, "X")==0) 
    tfbp->yPlane = 0;
  else {
    if (!(strcmp(tfbp->plane, "y")==0 || strcmp(tfbp->plane, "Y")==0))
      bombElegant("PLANE must be x or y for TFBPICKUP", NULL);
    tfbp->yPlane = 1;
  }

  tfbp->initialized = 1;
}

void transverseFeedbackDriver(TFBDRIVER *tfbd, double **part, long np, LINE_LIST *beamline, long pass, long nPasses, char *rootname)
{
  double kick;
  long i, j;
  
  if (tfbd->initialized==0)
      initializeTransverseFeedbackDriver(tfbd, beamline, nPasses, rootname);
  if (pass==0)
    tfbd->dataWritten = 0;
  
  if (tfbd->delay>tfbd->maxDelay) {
    if (!(tfbd->driverSignal = 
          realloc(tfbd->driverSignal, sizeof(*tfbd->driverSignal)*(tfbd->delay+1+TFB_FILTER_LENGTH))))
      bombElegant("memory allocation failure (transverseFeedbackDriver)", NULL);
    tfbd->maxDelay = tfbd->delay;
  }
    
  kick = tfbd->pickup->filterOutput*tfbd->strength;
  if (tfbd->kickLimit>0 && fabs(kick)>tfbd->kickLimit)
    kick = SIGN(kick)*tfbd->kickLimit;

#ifdef DEBUG
  fprintf(stdout, "TFBDRIVER: pass %ld\nstoring kick %e in slot %ld based on filter output of %e\n",
          pass, kick, pass%(tfbd->delay+tfbd->filterLength), tfbd->pickup->filterOutput);
#endif
  
  tfbd->driverSignal[pass%(tfbd->delay+tfbd->filterLength)] = kick;
  
  if (pass<tfbd->delay+tfbd->filterLength) {
#ifdef DEBUG
    fprintf(stdout, "TFBDRIVER: no kick applied for pass %ld due to delay of %ld\n",
            pass, tfbd->delay);
#endif
    kick = 0;
  }
  else {
    kick = 0;
    for (i=0; i<tfbd->filterLength; i++) {
#ifdef DEBUG
      fprintf(stdout, "TFBDRIVER: adding term a[%ld]=%e  *   %e\n",
              i, tfbd->a[i], tfbd->driverSignal[(pass - tfbd->delay - i)%(tfbd->delay+tfbd->filterLength)]);
#endif
      kick += tfbd->a[i]*tfbd->driverSignal[(pass - tfbd->delay - i)%(tfbd->delay+tfbd->filterLength)];
    }
    j = tfbd->pickup->yPlane?3:1;
    for (i=0; i<np; i++)
      part[i][j] += kick;
  }

  if (tfbd->outputFile) {
    if (!SDDS_SetRowValues(&tfbd->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                           pass, 
                           "Pass", pass,
                           "PickupOutput", tfbd->pickup->filterOutput, 
                           "DriverOutput", kick, NULL)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("problem writing data for TFBDRIVER output file");
    }
    if ((pass+1)==nPasses) {
      if (!SDDS_WritePage(&tfbd->SDDSout)) {
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
        SDDS_Bomb("problem writing data for TFBDRIVER output file");
      }
      tfbd->dataWritten = 1;
    }
  }
}

void flushTransverseFeedbackDriverFiles(TFBDRIVER *tfbd)
{
  if (tfbd->initialized && !(tfbd->dataWritten)) {
    if (!SDDS_WritePage(&tfbd->SDDSout)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("problem writing data for TFBDRIVER output file (flushTransverseFeedbackDriverFiles)");
    }
    tfbd->dataWritten = 1;
  }
}

void initializeTransverseFeedbackDriver(TFBDRIVER *tfbd, LINE_LIST *beamline, long nPasses, char *rootname)
{
  ELEMENT_LIST *eptr;
  long pickupFound = 0, i;

  if (tfbd->ID==NULL || !strlen(tfbd->ID))
    bombElegant("you must give an ID string for TFBDRIVER", NULL);
  
  for (i=TFB_FILTER_LENGTH-1; i>=0; i--) {
    if (tfbd->a[i]!=0)
      break;
  }
  if (i<0)
    bombElegant("All filter coefficients are zero for TFBDRIVER", NULL);
  tfbd->filterLength = i+1;

  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type==T_TFBPICKUP && strcmp(tfbd->ID, ((TFBPICKUP*)eptr->p_elem)->ID)==0) {
      pickupFound = 1;
      tfbd->pickup = ((TFBPICKUP*)eptr->p_elem);
      break;
    }
    eptr = eptr->succ;
  }
  if (!pickupFound) 
    bombElegant("pickup not found for TFBDRIVER", NULL);
  
  if (tfbd->delay<0)
    bombElegant("TFBDRIVER delay is negative", NULL);
  if (!(tfbd->driverSignal = malloc(sizeof(*(tfbd->driverSignal))*(tfbd->delay+1+TFB_FILTER_LENGTH))))
    bombElegant("memory allocation failure (TFBDRIVER)", NULL);
  tfbd->maxDelay = tfbd->delay;
  
  if (tfbd->outputFile) {
    tfbd->outputFile = compose_filename(tfbd->outputFile, rootname);
    if (!SDDS_InitializeOutput(&tfbd->SDDSout, SDDS_BINARY, 1, NULL, NULL, tfbd->outputFile) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "PickupOutput", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&tfbd->SDDSout, "DriverOutput", "rad", SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&tfbd->SDDSout) || !SDDS_StartPage(&tfbd->SDDSout, nPasses)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      SDDS_Bomb("Problem setting up TFBDRIVER output file");
    }
  }
  tfbd->dataWritten = 0;
  tfbd->initialized = 1;
}


