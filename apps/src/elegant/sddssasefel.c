/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsexpand
 * purpose: take a SDDS file and create a new file with the
 *          columns turned into parameters. Each page of the
 *          new file corresponds to a row of the old file.
 *          This is basically an inverse of sddscollapse.
 *
 * Michael Borland, 1999
 *
 $Log: not supported by cvs2svn $
 Revision 1.12  2010/08/12 15:32:21  borland
 Added exitElegant() routine so that any exit will result in creation of
 semaphore files if requested.

 Revision 1.11  2010/02/05 22:05:38  soliday
 Made some changes to reduce the number or compiler warnings on Linux.

 Revision 1.10  2003/01/21 18:49:04  borland
 Added simplexMin() arguments.

 Revision 1.9  2003/01/16 02:20:36  borland
 Added simplex_divisor parameter to optimization_setup namelist.
 Updated calls to simplexMin().

 Revision 1.8  2002/08/14 20:23:49  soliday
 Added Open License

 Revision 1.7  2002/06/15 02:21:12  borland
 Added upper and lower limits on all optimizable parameters.

 Revision 1.6  2001/10/15 20:37:07  soliday
 Cleaned up for Linux.

 Revision 1.5  2001/08/02 14:45:31  borland
 Added search_path feature to run_setup namelist and to many (all?) elements and
 commands that take input files.
 FITPOINT beam statistics data now includes x, y, and z emittances.

 Revision 1.4  2001/05/16 19:02:59  borland
 Modified calls to simplexMin() to accomodate new argument.
 Added simplified Rosenzweig/Serafini model for rf focusing in the
 body of RFCA and RFCW elements.

 Revision 1.3  2000/01/25 21:17:24  borland
 Now supports commandline setting of undulator K, period, and beta.
 Also accepts elegant .fin files as input.

 Revision 1.2  1999/10/12 21:50:02  borland
 All printouts now go to the stdout rather than stderr.  fflush statements,
 some unnecessary, were added in a mostly automated fashion.

 Revision 1.1  1999/07/01 19:19:57  borland
 First versions in repository.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

void exitElegant(long status);

#define DEBUG 0

#define SET_PIPE 0
#define SET_UNDULATORPERIOD 1
#define SET_UNDULATORK 2
#define SET_BETAX 3
#define SET_ELEGANT_INPUT 4
#define N_OPTIONS 5

char *option[N_OPTIONS] = {
  "pipe", "undulatorperiod", "undulatork", "betax", "elegantInput",
} ;

char *USAGE="sddssasefel [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
[-undulatorPeriod=<inMeters>] [-undulatorK=<value>] [-betax=<value>]\n\
[-elegantInput] \n\
Uses formulae of Ming Xie to compute performance of a SASE FEL.\n\
By default, the computations are done on the following columns of the input file:\n\
   charge (C)\n\
   rmsBunchLength (s)\n\
   pCentral (momentum in mc units)\n\
   ex0 (emittance, m-radians)\n\
   Sdelta0 (rms fractional energy spread)\n\
   undulatorPeriod (m) *\n\
   undulatorK *\n\
   betax (m) *\n\
Quantities marked with * may be omitted.  If omitted, the quantity may\n\
be specified on the commandline; otherwise, it will be optimized to \n\
minimize the gain length.\n\
If -elegantInput is given, then it is assumed that the input file is\n\
a \"final\" output file from elegant, with the data stored in the\n\
following parameters: Charge, Dt80 (for bunch length), pAverage,\n\
ex, and Sdelta.\n\
Program by Michael Borland.  (This is version 2, January 2000.)\n";

double FELScalingFunction
(double *etaDiffraction, double *etaEmittance, double *etaEnergySpread,
 double L1D, double beta, double emittance, double lightWavelength, double undulatorPeriod,
 double sigmaDelta);
void ComputeSASEFELParameters
  (double *lightWavelength, double *saturationLength, double *gainLength,  
   double *noisePower, double *saturationPower, double *PierceParameter, double *etaDiffraction, 
   double *etaEmittance, double *etaEnergySpread, double charge, double rmsBunchLength, double undulatorPeriod, 
   double undulatorK, double beta, double emittance, double sigmaDelta, double pCentral, short planar); 


/* This structure holds an array of values and a flag telling whether to
 * optimize the quantity.  If the flag is on, then the values are *set*
 * by the optimization.  Otherwise, the values are inputs.
 */
typedef struct {
  char *name, *units;
  char *elegantName;
  double *value, guess, fixedValue;
  short inFile, requiredInFile, onCommandline;
} VALUE_NAME;

void OptimizeSASEFELParameters
  (VALUE_NAME charge, VALUE_NAME rmsBunchLength, VALUE_NAME undulatorPeriod,
   VALUE_NAME undulatorK, VALUE_NAME betax, VALUE_NAME ex0, VALUE_NAME Sdelta0, VALUE_NAME pCentral,
   long row);
double SASEFELOptimFn(double *x, long *invalid);
void setOptimizationParameters(double *x0, double *dx, short *disable, long index, VALUE_NAME data, long row);
void transferOptimizationParameters(double *x0, long index, VALUE_NAME data, long row);
void ComputeSASEFELParameters
  (double *lightWavelength, double *saturationLength, double *gainLength,  double *noisePower,
   double *saturationPower, double *PierceParameter, double *etaDiffraction, double *etaEmittance,
   double *etaEnergySpread, double charge, double rmsBunchLength, double undulatorPeriod, double undulatorK, 
   double beta,  double emittance, double sigmaDelta, double pCentral, short planar);
double FELScalingFunction
  (double *etaDiffraction, double *etaEmittance,
   double *etaEnergySpread, double L1D, double beta, double emittance,
   double lightWavelength, double undulatorPeriod, double sigmaDelta);


int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long i, i_arg, noWarnings, row, rows;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  VALUE_NAME charge = { "charge", "C", "Charge", NULL, 1e-9, 0, 1, 0, 0};
  VALUE_NAME rmsBunchLength = { "rmsBunchLength", "s", "Dt80", NULL, 1e-12, 0, 1, 0, 0};
  VALUE_NAME pCentral = { "pCentral", NULL, "pAverage", NULL, 1e9, 0, 1, 0, 0};
  VALUE_NAME ex0 = { "ex0", "$gp$rm", "ex", NULL, 1e-12, 0, 1, 0, 0};
  VALUE_NAME Sdelta0 = { "Sdelta0", "", "Sdelta", NULL, 1e-4, 0, 1, 0, 0};
  VALUE_NAME undulatorPeriod = { "undulatorPeriod", "m", "undulatorPeriod", NULL, 1e-2, 0, 0, 0, 0};
  VALUE_NAME undulatorK = { "undulatorK", "", "undulatorK", NULL, 1, 0, 0, 0, 0};
  VALUE_NAME betax = {"betax", "m", "betax", NULL, 1, 0, 0, 0, 0};
  VALUE_NAME *ValueName[9] ;
  double lightWavelength, saturationLength, gainLength, noisePower, saturationPower;
  double PierceParameter, etaDiffraction, etaEmittance, etaEnergySpread, eyValue;
  long elegantInput, eyPresent;
  
  ValueName[0] = &charge;
  ValueName[1] = &rmsBunchLength;
  ValueName[2] = &pCentral;
  ValueName[3] = &ex0;
  ValueName[4] = &Sdelta0;
  ValueName[5] = &undulatorPeriod;
  ValueName[6] = &undulatorK;
  ValueName[7] = &betax;
  ValueName[8] = NULL;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = NULL;
  pipeFlags = noWarnings = elegantInput = 0;

  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_UNDULATORPERIOD:
        if (s_arg[i_arg].n_items!=2 || 
            sscanf(s_arg[i_arg].list[1], "%lf", &undulatorPeriod.fixedValue)!=1 ||
            undulatorPeriod.fixedValue<=0)
          SDDS_Bomb("invalid -undulatorPeriod syntax/values");
        undulatorPeriod.onCommandline = 1;
        break;
      case SET_UNDULATORK:
        if (s_arg[i_arg].n_items!=2 || 
            sscanf(s_arg[i_arg].list[1], "%lf", &undulatorK.fixedValue)!=1 ||
            undulatorK.fixedValue<=0)
          SDDS_Bomb("invalid -undulatorK syntax/values");
        undulatorK.onCommandline = 1;
        break;
      case SET_BETAX:
        if (s_arg[i_arg].n_items!=2 || 
            sscanf(s_arg[i_arg].list[1], "%lf", &betax.fixedValue)!=1 ||
            betax.fixedValue<=0)
          SDDS_Bomb("invalid -betax syntax/values");
        betax.onCommandline = 1;
        break;
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_ELEGANT_INPUT:
        elegantInput = 1;
        break;
      default:
        fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
        fflush(stdout);
        exitElegant(1);
        break;
      }
    }
    else {
      if (inputfile==NULL)
        inputfile = s_arg[i_arg].list[0];
      else if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  
  processFilenames("sddssasefel", &inputfile, &outputfile, pipeFlags, noWarnings, NULL);

  if (!SDDS_InitializeInput(&SDDSin, inputfile) ||
      !SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile)) 
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Check the input file for required data.  Simultaneously set up the output file
   * to reflect the input plus any optimized input parameters 
   */
  i = -1;
  while (ValueName[++i]) {
    if (!elegantInput) {
      if (SDDS_GetColumnIndex(&SDDSin, ValueName[i]->name)<0 ||
          ValueName[i]->onCommandline) {
        if (!ValueName[i]->onCommandline && ValueName[i]->requiredInFile) {
          fprintf(stdout, "Error (sddssasefel): %s is required to be in the input file\n",
                  ValueName[i]->name);
          fflush(stdout);
          exitElegant(1);
        }
        ValueName[i]->inFile = 0;
        if (!SDDS_DefineSimpleColumn(&SDDSout, ValueName[i]->name, 
                                     ValueName[i]->units, SDDS_DOUBLE))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      else {
        ValueName[i]->inFile = 1;
        if (SDDS_CheckColumn(&SDDSin, ValueName[i]->name, ValueName[i]->units, 
                             SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
          exitElegant(1);
        if (!SDDS_TransferColumnDefinition(&SDDSout, &SDDSin, ValueName[i]->name, NULL))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    } else {
      if (SDDS_GetParameterIndex(&SDDSin, ValueName[i]->elegantName)<0 ||
          ValueName[i]->onCommandline) {
        if (!ValueName[i]->onCommandline && ValueName[i]->requiredInFile) {
          fprintf(stdout, "Error (sddssasefel): %s is required to be in the input file\n",
                  ValueName[i]->name);
          fflush(stdout);
          exitElegant(1);
        }
        ValueName[i]->inFile = 0;
        if (!SDDS_DefineSimpleParameter(&SDDSout, ValueName[i]->elegantName, 
                                        ValueName[i]->units, SDDS_DOUBLE))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        ValueName[i]->inFile = 1;
        if (SDDS_CheckParameter(&SDDSin, ValueName[i]->elegantName, ValueName[i]->units, 
                                SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
          exitElegant(1);
        if (!SDDS_TransferParameterDefinition(&SDDSout, &SDDSin, ValueName[i]->elegantName, NULL))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
  
  /* Set up the output parameters in the output file. */

  eyPresent = 0;
  if (elegantInput) {
    if (SDDS_GetParameterIndex(&SDDSin, "ey")>=0) {
      if (SDDS_CheckParameter(&SDDSin, "ey", "$gp$rm", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
        exitElegant(1);
      if (!SDDS_TransferParameterDefinition(&SDDSout, &SDDSin, "ey", NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      eyPresent = 1;
    }
    if (!SDDS_DefineSimpleParameter(&SDDSout, "lightWavelength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "saturationLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "gainLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "noisePower", "W", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "saturationPower", "W", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "PierceParameter", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "etaDiffraction", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "etaEmittance", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "etaEnergySpread", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "rmsBunchLength", "s", SDDS_DOUBLE) ||
        !SDDS_SaveLayout(&SDDSout) || !SDDS_WriteLayout(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    i = -1;
    while (ValueName[++i])
      ValueName[i]->value = tmalloc(sizeof(*ValueName[i]->value)*1);
  } else if (!SDDS_DefineSimpleColumn(&SDDSout, "lightWavelength", "m", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "saturationLength", "m", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "gainLength", "m", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "noisePower", "W", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "saturationPower", "W", SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "PierceParameter", NULL, SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "etaDiffraction", NULL, SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "etaEmittance", NULL, SDDS_DOUBLE) ||
             !SDDS_DefineSimpleColumn(&SDDSout, "etaEnergySpread", NULL, SDDS_DOUBLE) ||
             !SDDS_SaveLayout(&SDDSout) || !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  while (SDDS_ReadPage(&SDDSin)>0) {
    if (!SDDS_CopyPage(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (!elegantInput) {
      if ((rows=SDDS_RowCount(&SDDSin))>0) {
        i = -1;
        while (ValueName[++i]) {
          if (ValueName[i]->inFile) {
            if (!(ValueName[i]->value 
                  = SDDS_GetColumnInDoubles(&SDDSin, ValueName[i]->name)))
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          } else
            ValueName[i]->value = tmalloc(sizeof(*ValueName[i]->value)*rows);
          if (ValueName[i]->onCommandline) 
            fill_double_array(ValueName[i]->value, rows, ValueName[i]->fixedValue);
        }
        for (row=0; row<rows; row++) {
          OptimizeSASEFELParameters(charge, rmsBunchLength, undulatorPeriod, undulatorK,
                                    betax, ex0, Sdelta0, pCentral, row);
          ComputeSASEFELParameters(&lightWavelength, &saturationLength,
                                   &gainLength, &noisePower, &saturationPower,
                                   &PierceParameter, &etaDiffraction, &etaEmittance,
                                   &etaEnergySpread,
                                   charge.value[row], rmsBunchLength.value[row],
                                   undulatorPeriod.value[row], undulatorK.value[row],
                                   betax.value[row], ex0.value[row],
                                   Sdelta0.value[row], pCentral.value[row], 1);
          if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row,
                                 "charge", charge.value[row],
                                 "rmsBunchLength", rmsBunchLength.value[row],
                                 "pCentral", pCentral.value[row],
                                 "ex0", ex0.value[row],
                                 "Sdelta0", Sdelta0.value[row],
                                 "undulatorPeriod", undulatorPeriod.value[row],
                                 "undulatorK", undulatorK.value[row],
                                 "betax", betax.value[row],
                                 "lightWavelength", lightWavelength,
                                 "saturationLength", saturationLength,
                                 "gainLength", gainLength, 
                                 "noisePower", noisePower,
                                 "saturationPower", saturationPower,
                                 "PierceParameter", PierceParameter,
                                 "etaDiffraction", etaDiffraction,
                                 "etaEmittance", etaEmittance,
                                 "etaEnergySpread", etaEnergySpread,
                                 NULL)) 
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    else {
      /* elegant input */
      i = -1;
      while (ValueName[++i]) {
        if (ValueName[i]->inFile) {
          if (!SDDS_GetParameterAsDouble(&SDDSin, ValueName[i]->elegantName,
                                         &ValueName[i]->value[0]))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        else if (ValueName[i]->onCommandline)
          ValueName[i]->value[0] = ValueName[i]->fixedValue;
      }
      row = 0;
      if (eyPresent) {
        if (!SDDS_GetParameterAsDouble(&SDDSin, "ey", &eyValue))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        ex0.value[row] = sqrt(ex0.value[row]*eyValue);
      }
      rmsBunchLength.value[0] = rmsBunchLength.value[0]/(0.8*sqrt(2*PI));
      OptimizeSASEFELParameters(charge, rmsBunchLength, undulatorPeriod, undulatorK,
                                betax, ex0, Sdelta0, pCentral, row);
      ComputeSASEFELParameters(&lightWavelength, &saturationLength,
                               &gainLength, &noisePower, &saturationPower,
                               &PierceParameter, &etaDiffraction, &etaEmittance,
                               &etaEnergySpread,
                               charge.value[row], rmsBunchLength.value[row],
                               undulatorPeriod.value[row], undulatorK.value[row],
                               betax.value[row], ex0.value[row],
                               Sdelta0.value[row], pCentral.value[row], 1);
      if (!SDDS_SetParametersFromDoubles(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                         "rmsBunchLength", rmsBunchLength.value[row],
                                         "undulatorPeriod", undulatorPeriod.value[row],
                                         "undulatorK", undulatorK.value[row],
                                         "betax", betax.value[row],
                                         "lightWavelength", lightWavelength,
                                         "saturationLength", saturationLength,
                                         "gainLength", gainLength, 
                                         "noisePower", noisePower,
                                         "saturationPower", saturationPower,
                                         "PierceParameter", PierceParameter,
                                         "etaDiffraction", etaDiffraction,
                                         "etaEmittance", etaEmittance,
                                         "etaEnergySpread", etaEnergySpread,
                                         NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return(0);
}

void OptimizeSASEFELParameters
  (VALUE_NAME charge, VALUE_NAME rmsBunchLength, VALUE_NAME undulatorPeriod,
   VALUE_NAME undulatorK, VALUE_NAME betax, VALUE_NAME ex0, VALUE_NAME Sdelta0, VALUE_NAME pCentral,
   long row)
{
  short disable[8] = { 0,0,0,0,0,0,0,0 };
#if DEBUG
  double lower[8] = { 0,0,0,0,0,0,0,0 };
  double upper[8] = { DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,
		      DBL_MAX,DBL_MAX,DBL_MAX };
#endif
  double x0[8], dx[8], result;
  long i;

  setOptimizationParameters(x0, dx, disable, 0, charge, row);
  setOptimizationParameters(x0, dx, disable, 1, rmsBunchLength, row);
  setOptimizationParameters(x0, dx, disable, 2, undulatorPeriod, row);
  setOptimizationParameters(x0, dx, disable, 3, undulatorK, row);
  setOptimizationParameters(x0, dx, disable, 4, betax, row);
  setOptimizationParameters(x0, dx, disable, 5, ex0, row);
  setOptimizationParameters(x0, dx, disable, 6, Sdelta0, row);
  setOptimizationParameters(x0, dx, disable, 7, pCentral, row);
  for (i=0; i<8; i++)
    if (!disable[i])
      break;
  if (i==8)
    return;
#if DEBUG
  fprintf(stdout, "Before optimization:\n");
  fflush(stdout);
  fprintf(stdout, " charge=%le, rmsBL=%le, uP=%le, uK=%le, betax=%le, ex0=%le, Sdelta=%le, pC=%le\n",
          x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7]);
  fflush(stdout);
  if (simplexMin(&result, x0, dx, lower, upper, disable, 8, 0.0, 1e-10,
                 SASEFELOptimFn, NULL, 500, 10, 12, 3.0, 1.0, 0)<0) {
    fprintf(stdout, "Optimization unsuccessful\n");
    fflush(stdout);
  }
#endif
  for (i=0; i<8; i++)
    dx[i] = 0.1*x0[i];
  if (simplexMin(&result, x0, dx, NULL, NULL, disable, 8, 0.0, 1e-9,
                 SASEFELOptimFn, NULL, 500, 10, 12, 3.0, 1.0, 0)<0) {
    fprintf(stdout, "Optimization unsuccessful\n");
    fflush(stdout);
  }
#if DEBUG
  fprintf(stdout, "After optimization:\n");
  fflush(stdout);
  fprintf(stdout, " charge=%le, rmsBL=%le, uP=%le, uK=%le, betax=%le, ex0=%le, Sdelta=%le, pC=%le\n",
          x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7]);
  fflush(stdout);
#endif

  transferOptimizationParameters(x0, 0, charge, row);
  transferOptimizationParameters(x0, 1, rmsBunchLength, row);
  transferOptimizationParameters(x0, 2, undulatorPeriod, row);
  transferOptimizationParameters(x0, 3, undulatorK, row);
  transferOptimizationParameters(x0, 4, betax, row);
  transferOptimizationParameters(x0, 5, ex0, row);
  transferOptimizationParameters(x0, 6, Sdelta0, row);
  transferOptimizationParameters(x0, 7, pCentral, row);
}

double SASEFELOptimFn(double *x, long *invalid)
{
  double lightWavelength, saturationLength, gainLength, noisePower, saturationPower;
  double PierceParameter, etaDiffraction, etaEmittance, etaEnergySpread;
  long i;
  *invalid = 0;
  for (i=0; i<7; i++) {
    if (x[i]<0) {
      *invalid = 1;
      return 0;
    }
  }
  ComputeSASEFELParameters(&lightWavelength, &saturationLength,
                           &gainLength, &noisePower, &saturationPower,
                           &PierceParameter, &etaDiffraction, &etaEmittance,
                           &etaEnergySpread,
                           x[0], /* charge */
                           x[1], /* bunch length */
                           x[2], /* undulator period */
                           x[3], /* undulator K */
                           x[4], /* betax */
                           x[5], /* ex0 */
                           x[6], /* Sdelta */
                           x[7], /* pCentral */
                           1);
  if (saturationLength<=0 || PierceParameter<=0 ||
      gainLength<=0 || noisePower<=0 || saturationPower<=0 ||
      etaDiffraction<=0 || etaEmittance<=0 || etaEnergySpread<=0) {
    *invalid = 1;
    return 0;
  }
#if DEBUG 
  fprintf(stdout, "opt: Q=%g, BL=%g, UP=%g, UK=%g, BX=%g, EX=%g, SD=%g, PC=%g -> %g\n",
          x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], saturationLength);
  fflush(stdout);
#endif
  return saturationLength;
}

void setOptimizationParameters(double *x0, double *dx, short *disable, long index, VALUE_NAME data, long row)
{
  if (data.inFile || data.onCommandline) {
    x0[index] = data.value[row];
    dx[index] = 0.0;
    disable[index] = 1;
  }
  else {
    x0[index] = data.guess;
    dx[index] = 0.1*x0[index];
    disable[index] = 0;
#if DEBUG
    fprintf(stdout, "optimizing %s\n", data.name);
    fflush(stdout);
#endif
  }
}

void transferOptimizationParameters(double *x0, long index, VALUE_NAME data, long row)
{
  if (!data.inFile && !data.onCommandline)
    data.value[row] = x0[index];
}
                               
