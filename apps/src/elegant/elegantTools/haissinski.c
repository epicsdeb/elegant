/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* Copyright 1994 by Louis Emery and Argonne National Laboratory,
 * all rights reserved.
 * Solves the Haissinski equation for density distribution
 * of a bunch.
 * Louis Emery, 2000
 */

#include <stdio.h>
#include "mdb.h"
#include "scan.h"
#include "match_string.h"
#include "SDDS.h"
#include "constants.h"

static char *USAGE1 = "haissinski <twissFile> <resultsFile>\n\
 {-wakeFunction=<file>,tColumn=<name>,wColumn=<name> |\n\
  -model=[L=<Henry>|Zn=<Ohms>],R=<Ohm> |\n\
  -BBResonator=Rs=<Ohm>,frequency=<Hz>,Q=<value>[,wall=<Ohms>]} \n\
 {-charge=<C>|-particles=<value>|-bunchCurrent=<A>}\n\
 {-steps=<numberOfChargeSteps>} {-outputLastStepOnly}\n\
 {-RF=Voltage=<V>,harmonic=<value>|-length=<s>}\n\
 {-harmonicCavity=Voltage=<V>,factor=<harmonicFactor>}\n\
 {-superPeriods=<number>} {-energy=<GeV>} \n\
 -integrationParameters=deltaTime=<s>,points=<number>,startTime=<s>,\n\
iterations=<number>,fraction=<value>,tolerance=<value> \n\
Calculation of steady-state longitudinal bunch density distribution \n\
in an electron storage ring. \n\
<twissFile>    Elegant twiss file that contains ring parameters required\n\
               for calculation: pCentral (beam momentum), \n\
               alphac (momentum compaction), U0 (energy loss per turn),\n\
               Sdelta0 (rms of momentum distribution), etc.\n\
<resultsFile>  Output file with bunch distribution, wake potential for the\n\
               solved bunch distribution.\n\
wakeFunction   Input the wake function for an impulse response of a 1 C charge.\n\
               Time and wake column names should be specified with \n\
               units s and V/C respectively.\n";
static char *USAGE2 = "model          Instead of a wake function, a circuit model can be entered. The inductive\n\
               part may be entered with L (inductance) or with |Z/n| value.\n\
               R is the resistance.\n\
BBResonator    Instead of a wake function or circuit model, one can specify the\n\
               parameters of a broad-band resonator, optionally with a resistive\n\
               wall term (implemented as a simple resistor as in -model).\n\
RF             RF parameters that control the length of the zero-current beam.\n\
length         Alternatively, one can specify the zero-current bunch length (rms)\n\
               in seconds directly.\n\
harmonicCavity If -RF is given, one may also specify a higher-harmonic cavity\n\
               for possible bunch shortening.  If you want to lengthen the bunch,\n\
               give a negative voltage.\n\
superPeriods   The number of superiods of lattice file in the actual machine.\n\
               Default is 1.\n\
energy         The beam energy at which to perform computations, if it is desired\n\
               that this be different than the energy in the Twiss file.\n\
intermediateSolutions   write all intermediate solutions to the results file.\n\
               Used for debugging.\n\
integrationParameters   Specifies integration parameters for solving the \n\
               Haissinski integral equation. deltaTime specifies the time \n\
               interval in seconds for evaluating wake potential and charge density.\n\
               points is the number of points requested for the charge density.\n\
               startTime is the time in ps for the first point relative to the \n\
               zero-current synchronous phase. iterationLimits is the number of \n\
               iterations for solving the integral equation. fraction is\n\
               the relaxation fraction between 0 and 1 to reduce numerical\n\
               instabilities.";

#define VERBOSE 0
#define CHARGE 1
#define PARTICLES 2
#define RF 3
#define LENGTH 4
#define INTEGRATION 5
#define STEPS 6
#define WAKE 7
#define MODEL 8
#define INTERMEDIATE_SOLUTIONS 9
#define SUPERPERIODS 10
#define ENERGY 11
#define OUTPUT_LAST_STEP_ONLY 12
#define HARMONIC_CAVITY 13
#define BUNCH_CURRENT 14
#define BBRESONATOR 15
#define N_OPTIONS 16
char *option[N_OPTIONS] = {
  "verbose",
  "charge",
  "particles",
  "rf",
  "length",
  "integrationParameters",
  "steps",
  "wakeFunction",
  "model",
  "intermediateSolutions",
  "superperiods",
  "energy",
  "outputlaststeponly",
  "harmoniccavity",
  "bunchcurrent",
  "bbresonator" };

typedef struct {
  double *y, xStart, xDelta;
  long points;
  long offset;
/* Different units from SI may be used.
   However to be compliant with elegant, the units of 
   output files should be SI units.
   These factors allows conversion to SI units. */
  double xFactor, yFactor;
} FUNCTION;

#define ITERATION_LIMITS 1000

void printFunction( char *label, FUNCTION *data);
void initializeFunction( FUNCTION *data);
void readRingParameters( char *twissFile, long superPeriods, double desiredEnergyMeV,
                        double *energyMeV, double *momentumCompaction,
                        double *U0, double *sigmaDelta,
                        double *circumference);
/* returned results is the first parameter. */
void getWakeFunction( FUNCTION *wake, char *wakeFile, char *tCol, char *wCol,
                     double delta, long points);
void integrateWakeFunction( FUNCTION *stepResponse, FUNCTION *wake);
void setupResultsFile( SDDS_TABLE *resultsPage, char *resultsFile, long points);
void getInitialDensity( FUNCTION *density, long points, double deltaTime, double charge, double length);
void copyDensity( FUNCTION *target, FUNCTION *source);
void getPotentialDistortion( FUNCTION *potentialDistortion,
                  FUNCTION *Vinduced,
                  FUNCTION *density, FUNCTION *wake);
void getPotentialDistortionFromModel( FUNCTION *potentialDistortion,
                           FUNCTION *Vinduced,
                           FUNCTION *density, double L, double R);
void calculateDistribution( FUNCTION *distribution, FUNCTION *potential,
                          FUNCTION *potentialDistortion, double length, double VrfDot);
void normalizeDensityFunction( FUNCTION *density, FUNCTION *distribution, double charge);
void writeResults( SDDS_TABLE *resultsPage, FUNCTION *density, 
                  FUNCTION *potential, FUNCTION *potentialDistortion, 
                  FUNCTION *Vind, double charge, double averageCurrent,
                  long converged);
void makeBBRWakeFunction(FUNCTION *wake, double dt, long points, 
                         double Q, double R, double omega, double rw, double T0);

long verbosity;

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  long i, j, k, outputLastStepOnly;
  char *twissFile, *resultsFile;
  SDDS_DATASET resultsPage;
  double particles, charge, finalCharge, length;
  long steps, converged, superPeriods;
  int32_t points, iterationLimits;
  long useWakeFunction=0, intermediateSolutions;
  double rfVoltage, rfHarmonic;
  double energyMeV, momentumCompaction, U0, sigmaE, circumference;
  double startTime, deltaTime;
  unsigned long dummyFlags;
  FUNCTION density, densityOld, diff, wake, stepResponse;
  FUNCTION potential, potentialDistortion, distribution, Vinduced;
  char *wakeFile, *tCol, *wCol;
  double revFrequency, syncPhase, syncTune, syncAngFrequency;
  double VrfDot, ZoverN, inductance, resistance;
  double rfHigherHarmonic=0, rfHigherHarmonicVoltage=0;
  double maxDifference, rmsDifference, madDifference, maxTolerance, fraction, lastMaxDifference;
  double maxDensity;
  double averageCurrent=0.0, bunchCurrent=0.0, desiredEnergy;
  long useBBR = 0;
  double BBR_R, BBR_Q, BBR_frequency, BBR_rw;
  
  SDDS_RegisterProgramName(argv[0]);
  argc  =  scanargs(&scanned, argc, argv);
  if (argc == 1) {
    fprintf(stderr, "%s%s\n", USAGE1, USAGE2);
    exit(1);
  }

/* initialize the FUNCTION structures */
  initializeFunction( &density );
  initializeFunction( &densityOld );
  initializeFunction( &diff );
  initializeFunction( &wake );
  initializeFunction( &stepResponse );
  initializeFunction( &potential );
  initializeFunction( &potentialDistortion );
  initializeFunction( &distribution );
  initializeFunction( &Vinduced );
  
  twissFile  =  NULL;
  resultsFile  =  NULL;
  verbosity = 0;
  particles = finalCharge = bunchCurrent = 0;
  length = 0;
  rfVoltage = rfHarmonic = 0;
  steps = 1;
  points=1000;
  iterationLimits = ITERATION_LIMITS;
  startTime=0.0;
  deltaTime=0;
  maxTolerance = 0.0001; /* fraction of maximum */
  intermediateSolutions=0;
  fraction = 0.01;
  wakeFile = tCol = wCol = NULL;
  superPeriods = 1;
  desiredEnergy = 0.0;
  outputLastStepOnly = 0;
  
  for (i = 1; i<argc; i++) {
    if (scanned[i].arg_type == OPTION) {
      delete_chars(scanned[i].list[0], "_");
      switch(match_string(scanned[i].list[0], option, N_OPTIONS, UNIQUE_MATCH)) {
      case VERBOSE:
        if(scanned[i].n_items > 1 ) {
          get_long(&verbosity, scanned[i].list[1]);
        } else {
          verbosity=1;
        }
        break;
      case CHARGE:
        if (scanned[i].n_items<2)
          bomb("invalid -charge syntax", NULL);
        get_double(&finalCharge, scanned[i].list[1]);
        break;
      case BUNCH_CURRENT:
        if (scanned[i].n_items<2)
          bomb("invalid -bunchCurrent syntax", NULL);
        get_double(&bunchCurrent, scanned[i].list[1]);
        break;
      case PARTICLES:
        if (scanned[i].n_items<2)
          bomb("invalid -particles syntax", NULL);
        get_double(&particles, scanned[i].list[1]);
        break;
      case STEPS:
        if (scanned[i].n_items<2)
          bomb("invalid -steps syntax", NULL);
        get_long(&steps, scanned[i].list[1]);
        break;
      case SUPERPERIODS:
        if (scanned[i].n_items<2)
          bomb("invalid -superPeriods syntax", NULL);
        get_long(&superPeriods, scanned[i].list[1]);
        break;
      case ENERGY:
        if (scanned[i].n_items<2)
          bomb("invalid -energy syntax", NULL);
        get_double(&desiredEnergy, scanned[i].list[1]);
        break;
      case INTERMEDIATE_SOLUTIONS:
        intermediateSolutions = 1;
        break;
      case OUTPUT_LAST_STEP_ONLY:
        outputLastStepOnly = 1;
        break;
      case LENGTH:
        if (scanned[i].n_items<2)
          bomb("invalid -length syntax", NULL);
        get_double(&length, scanned[i].list[1]);
        break;
      case RF:
        if (scanned[i].n_items<2)
          bomb("invalid -rf syntax", NULL);
        scanned[i].n_items--;
        rfVoltage = rfHarmonic = 0;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "voltage", SDDS_DOUBLE, &rfVoltage, 1, 0,
                          "harmonic", SDDS_DOUBLE, &rfHarmonic, 1, 0,
                          NULL) ||
            rfVoltage<=0 || rfHarmonic<=0)
          bomb("invalid -rf syntax/values", "-rf=voltage=<V>,harmonic=<value>");
        break;
      case HARMONIC_CAVITY:
        if (scanned[i].n_items<2)
          bomb("invalid -harmonicCavity syntax", NULL);
        scanned[i].n_items--;
        rfHigherHarmonic = rfHigherHarmonicVoltage = 0;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "voltage", SDDS_DOUBLE, &rfHigherHarmonicVoltage, 1, 0,
                          "harmonic", SDDS_DOUBLE, &rfHigherHarmonic, 1, 0,
                          NULL) ||
            rfHigherHarmonic<=0)
          bomb("invalid -harmonicCavity syntax/values", 
               "-harmonicCavity=voltage=<V>,harmonicFactor=<value>");
        break;
      case INTEGRATION:
        scanned[i].n_items--;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "points", SDDS_LONG, &points, 1, 0,
                          "startTime", SDDS_DOUBLE, &startTime, 1, 0,
                          "deltaTime", SDDS_DOUBLE, &deltaTime, 1, 0,
                          "fraction", SDDS_DOUBLE, &fraction, 1, 0,
                          "tolerance", SDDS_DOUBLE, &maxTolerance, 1, 0,
                          "iterations", SDDS_LONG, &iterationLimits, 1, 0,
                          NULL) ||
            points<=0 || iterationLimits<=0 || maxTolerance<=0 || deltaTime<0) {
          bomb("invalid -integrationParameters syntax/values", "-integrationParameters=deltaTime=<s>,points=<number>,startTime=<s>,iterations=<number>");
        }
        
        break;
      case WAKE:
        if(!strlen(wakeFile=scanned[i].list[1]))
          bomb("bad -wake syntax", "-wakeFunction=<file>");
        scanned[i].n_items -= 2;
        if (!scanItemList(&dummyFlags, scanned[i].list+2, &scanned[i].n_items, 0,
                          "tColumn", SDDS_STRING, &tCol, 1, 0,
                          "wColumn", SDDS_STRING, &wCol, 1, 0,
                          NULL) ||
            !tCol || !wCol)
          bomb("invalid -wake syntax/values", "-wakeFunction=<file>,tColumn=<name>,wColumn=<name>");
        useWakeFunction = 1;
        useBBR = 0;
        break;
      case MODEL:
        if (scanned[i].n_items<2)
          bomb("invalid -model syntax", NULL);
        scanned[i].n_items--;
        inductance = ZoverN = resistance = 0.0;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "L", SDDS_DOUBLE, &inductance, 1, 0,
                          "Zn", SDDS_DOUBLE, &ZoverN, 1, 0,
                          "R", SDDS_DOUBLE, &resistance, 1, 0,
                          NULL) ||
            inductance<0 || resistance<0 || ZoverN<0)
          bomb("invalid -model syntax/values", "-model=[L=<Henry>|Zn=<Ohms>],R=<Ohm>");
        useWakeFunction = useBBR = 0;
        break;
      case BBRESONATOR:
        if (scanned[i].n_items<3)
          bomb("invalid -BBResonator syntax", NULL);
        scanned[i].n_items--;
        BBR_Q = BBR_frequency = BBR_R = BBR_rw = 0;
        if (!scanItemList(&dummyFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "R", SDDS_DOUBLE, &BBR_R, 1, 0,
                          "Q", SDDS_DOUBLE, &BBR_Q, 1, 0,
                          "frequency", SDDS_DOUBLE, &BBR_frequency, 1, 0,
                          "wall", SDDS_DOUBLE, &BBR_rw, 1, 0,
                          NULL) ||
            inductance<0 || resistance<0 || ZoverN<0)
          bomb("invalid -model syntax/values", "-model=[L=<Henry>|Zn=<Ohms>],R=<Ohm>");
        useWakeFunction = 0;
        useBBR = 1;
        break;
      default:
        bomb("unknown option given.", NULL);  
      }
    }
    else {
      if (!twissFile)
        twissFile  =  scanned[i].list[0];
      else if (!resultsFile) 
        resultsFile =  scanned[i].list[0];
      else
        bomb("too many filenames given", NULL);
    }
  }
  if (((finalCharge!=0)+(particles!=0)+(bunchCurrent!=0))>1)
    bomb("Specify only one of -charge, -particles, or -bunchCurrent.",NULL);

  if (length && rfVoltage)
    bomb("Options length and RF cannot be both specified.",NULL);
  if (rfHigherHarmonic) {
    if (!rfVoltage)
      bomb("You must give -rf if you give -harmonicCavity", NULL);
  }
  readRingParameters( twissFile, superPeriods, desiredEnergy*1e3,
                     &energyMeV, &momentumCompaction,
                     &U0, &sigmaE, &circumference);
  revFrequency = c_mks/ circumference;

  if (bunchCurrent)
    finalCharge = bunchCurrent/revFrequency;
  if (!finalCharge) 
    finalCharge = particles * e_mks;
  if (!particles)
    particles = finalCharge/ e_mks;

  if (!length) {
    syncPhase = asin( U0/ rfVoltage);
    syncTune = sqrt( momentumCompaction * rfHarmonic * cos(syncPhase) /
                    2 / PI * rfVoltage / (energyMeV * 1e6));
    if (rfHigherHarmonic)
      syncTune 
        *= sqrt(1 + (rfHigherHarmonicVoltage/rfVoltage)*rfHigherHarmonic/cos(syncPhase));
    syncAngFrequency = syncTune * 2 * PI * revFrequency;
    length = momentumCompaction * sigmaE/ syncAngFrequency; /* length is seconds */
    /* derivative w.r.t time */
    VrfDot = rfVoltage * 2 * PI * rfHarmonic * revFrequency * 
      (cos(syncPhase) + rfHigherHarmonicVoltage/rfVoltage*rfHigherHarmonic);
  }
  else {
    syncAngFrequency = momentumCompaction * sigmaE/ length; /* length is seconds */
    syncTune = syncAngFrequency / 2 / PI / revFrequency;
    VrfDot = sqr(2 * PI) * revFrequency * sqr(syncTune) * (energyMeV * 1e6)/
      momentumCompaction;
  }
  if (deltaTime==0)
    deltaTime = length/100;
  
  if (!useWakeFunction && !useBBR) {
    if ( inductance==0 && ZoverN != 0) {
      inductance = ZoverN * circumference / 2 / PI/ c_mks;
    }
  }
  
  if (useWakeFunction || useBBR) {
    useWakeFunction = 1;
    if (useBBR) {
      /* generate the wake function for a BBR
       * W = omega*R/Q * exp(-omega*t/(2*Q)) * (cos (omega*t) - sin(omega'*t)/(2*Q'))
       * units of W are V/C
       */
#ifdef DEBUG
      FILE *fp;
      long i;
      fprintf(stderr, "Making BBR wake\n");
#endif

      makeBBRWakeFunction(&wake, deltaTime, points, BBR_Q, BBR_R, BBR_frequency*PIx2, BBR_rw, 1/revFrequency);

#ifdef DEBUG
      fp = fopen("BBRWake.sdds", "w");
      fprintf(fp, "SDDS1\n&column name=Time, type=double, units=s &end\n");
      fprintf(fp, "&column name=Wakefield, type=double, units=V/C &end\n");
      fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
      for (i=0; i<points; i++) 
        fprintf(fp, "%e %e\n", deltaTime*i, wake.y[i]);
      fclose(fp);
#endif
    }
    else {
      /*
        determine wake function from file data or from parameters.
        May require interpolation.
        */
      getWakeFunction( &wake, wakeFile, tCol, wCol, deltaTime, points);
    }
      
    /* integrate Wake function for step response
     */
    integrateWakeFunction( &stepResponse, &wake);
  } 

  setupResultsFile( &resultsPage, resultsFile, points);
  
  /* Loop over charge steps
   */
  for (i=1; i<=steps ; i++) {
    charge = 1.0 * i/ steps * finalCharge;
    if (i==1) {
      getInitialDensity( &density, points, deltaTime, charge, length);
    }
    else {
      for (j=0;j<points;j++) {
        density.y[j] *= 1.0 * i/ (i-1);
      }
    }
    /* Loop integration until convergence
     */
    converged = 0;
    if (verbosity) {
      fprintf( stdout, "Iterations for solution for charge %g\n", charge);
      fflush(stdout);
    }
    lastMaxDifference = DBL_MAX;
    for (j=0; j<iterationLimits; j++) {
      copyDensity( &densityOld, &density );
      /* Calculate potential well distortion term. 
         Induced voltage is returned as well.
       */
      if (useWakeFunction) {
        getPotentialDistortion( &potentialDistortion, &Vinduced, &density, &wake);
      }
      else {
        getPotentialDistortionFromModel( &potentialDistortion, &Vinduced, 
                              &density, inductance, resistance);
      }
      /* Calculate distribution exponential (Use K. Bane's expression in
         SLAC-PUB-5177 p. 30)
       */
      calculateDistribution( &distribution, &potential, &potentialDistortion,
                           length, VrfDot);
      /* Normalize for new density function
       */
      normalizeDensityFunction( &density, &distribution, charge);
      /* mix old density with new density to see if numerical instability
       is reduced */
      for (k=0;k<points;k++) {
        density.y[k] = fraction * density.y[k] + (1 - fraction) * densityOld.y[k];
      }

      /* Compare new density function with
         old one and check for convergence 
         */
      maxDifference = 0.0;
      madDifference = 0.0;
      rmsDifference = 0.0;
      diff = density;
      diff.y = SDDS_Malloc( sizeof(*diff.y) * diff.points);
      maxDensity = 0.0;
      for (k=0; k<diff.points;k++) {
        if (density.y[k]>maxDensity)
          maxDensity = density.y[k];
        diff.y[k] = density.y[k] - densityOld.y[k];
        if (FABS(diff.y[k]) > maxDifference) {
          maxDifference = FABS(diff.y[k]);
        }
        rmsDifference += sqr(diff.y[k]);
        madDifference += FABS(diff.y[k]);
      }
      rmsDifference /= points;
      rmsDifference = sqrt(rmsDifference);
      madDifference /= points;
      if (verbosity) {
        fprintf(stdout, "maxDifference %10.2g madDifference %10.2g rmsDifference %10.2g.\n", 
                maxDifference/maxDensity, madDifference/maxDensity, rmsDifference/maxDensity);
        fflush(stdout);
      }
      if (maxDifference/maxDensity < maxTolerance) {
        if (verbosity) {
          fprintf(stdout, "Converged for charge %g after %ld iterations.\n", charge, j);
          fflush(stdout);
        }
        converged = 1;
        break;
      }
      free(diff.y);
      if (intermediateSolutions) {
        writeResults( &resultsPage, &density, &potential, &potentialDistortion, 
                     &Vinduced, charge, averageCurrent, converged );
      }
      if (maxDifference>lastMaxDifference)
        fraction /= 2;
      lastMaxDifference = maxDifference;
    }
    /* write results whether converged or not */
    averageCurrent = charge * revFrequency;
    if (!outputLastStepOnly || i==steps) 
      writeResults( &resultsPage, &density, &potential, &potentialDistortion, 
                   &Vinduced, charge, averageCurrent, converged );
  }
  if (!SDDS_Terminate(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);

  return 0;
}

void initializeFunction( FUNCTION *data) {
  data->y = NULL;
  data->xStart = data->xDelta = 0.0;
  data->points = data->offset = 0;
  data->xFactor = data->yFactor = 0.0;
}

void readRingParameters( char *twissFile, long superPeriods, 
                        double desiredEnergyMeV, double *energyMeV, 
                        double *momentumCompaction,
                        double *U0, double *sigmaE,
                        double *circumference) {
  SDDS_TABLE twissPage;
  double pCentral, *s;
  long elements;
  
  if (verbosity) {
    fprintf( stdout, "Opening for reading \"%s\"\n", twissFile);
    fflush(stdout);
  }
  if (!SDDS_InitializeInput(&twissPage, twissFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  /* read first page of input file to get parameters 
     */
  SDDS_ReadPage(&twissPage);
  if (!SDDS_GetParameters(&twissPage,
                          "pCentral", &pCentral,
                          "alphac", momentumCompaction,
                          "U0", U0,
                          "Sdelta0", sigmaE,
                          NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  *energyMeV = sqrt(sqr(pCentral) + 1) * me_mev;
  elements = SDDS_CountRowsOfInterest(&twissPage);
  s = SDDS_GetColumnInDoubles(&twissPage, "s");
  *circumference = s[elements-1]*superPeriods;
  *U0 *= 1e6*superPeriods; /* units in eV */

  if (desiredEnergyMeV>0) {
    if (verbosity) {
      fprintf(stdout, "Scaling from %f MeV to %f MeV\n",
              *energyMeV, desiredEnergyMeV);
      fflush(stdout);
    }
    *U0 *= ipow(desiredEnergyMeV/(*energyMeV), 4);
    *sigmaE *= desiredEnergyMeV/(*energyMeV);
    *energyMeV = desiredEnergyMeV;
  }

  if (!SDDS_Terminate(&twissPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
  free(s);
}

void getWakeFunction(FUNCTION *wake, char *wakeFile, char *tCol, char *wCol,
                     double deltaTime, long points) {
  SDDS_TABLE wakePage;
  double time, *t, *V;
  long i, rows, returnCode;
  
  if (verbosity) {
    fprintf( stdout, "Opening for reading \"%s\"\n", wakeFile);
    fflush(stdout);
  }
  
  if (!SDDS_InitializeInput(&wakePage, wakeFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  SDDS_ReadPage(&wakePage);
  rows = SDDS_CountRowsOfInterest(&wakePage);
  switch(SDDS_CheckColumn(&wakePage, tCol, "s", SDDS_DOUBLE, stderr)) {
  case SDDS_CHECK_NONEXISTENT:
  case SDDS_CHECK_WRONGTYPE:
  case SDDS_CHECK_WRONGUNITS:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  switch(SDDS_CheckColumn(&wakePage, wCol, "V/C", SDDS_DOUBLE, stderr)) {
  case SDDS_CHECK_NONEXISTENT:
  case SDDS_CHECK_WRONGTYPE:
  case SDDS_CHECK_WRONGUNITS:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  t = SDDS_GetColumnInDoubles(&wakePage, tCol);
  V = SDDS_GetColumnInDoubles(&wakePage, wCol);

  /* start as close as possible to first point in file. In 
   case of wake function, there may be some negative t values
   included for some reason (even though this might correspond
   to a non-causal wake.) */
  wake->xDelta = deltaTime;
  if (t[0]<0) {
    wake->offset = ((long) (-t[0] / deltaTime)) - 1;
    wake->xStart =  wake->offset * deltaTime;
  }
  else {
    wake->offset = 0;
    wake->xStart = 0.0;
  }
  wake->points = points - wake->offset; /* add more points if there is non-causal points */
  wake->y = SDDS_Malloc( sizeof(*wake->y) * wake->points);
  wake->xFactor = 1; 
  wake->yFactor = 1; 
  /* interpolate to produce points spaced by deltaTime. */
  for (i=0; i<points; i++) {
    time = wake->xStart + i * deltaTime;
    wake->y[i] = interp( V, t, rows, time, 1, 1, &returnCode );
  }
}

void integrateWakeFunction(FUNCTION *stepResponse, FUNCTION *wake) {
  long i;
  
  *stepResponse = *wake;
  stepResponse->y = SDDS_Malloc( sizeof(*stepResponse->y) * stepResponse->points);
  stepResponse->yFactor = 1;  /* using units (V/C)-s */
  stepResponse->y[0] = 0.0;
  for (i=1; i<wake->points; i++) {
    stepResponse->y[i] = stepResponse->y[i-1] + wake->xDelta * 
      (wake->y[i] + wake->y[i-1]) / 2.0;
  }
}

void setupResultsFile( SDDS_TABLE *resultsPage, char *resultsFile, long points) {
  /* define three columns but no parameters */
  if (!SDDS_InitializeOutput(resultsPage, SDDS_BINARY, 1, 
                             "Bunch density distribution", 
                             "Bunch density distribution", 
                             resultsFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (0>SDDS_DefineParameter(resultsPage, "Charge", "Q", "C",
                             "Total Charge", NULL,
                             SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(resultsPage, "AverageCurrent", "I$bave$n", "A",
                             "Average current", NULL,
                             SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(resultsPage, "Convergence", NULL, NULL,
                             "Convergence of iterations of integral equation", NULL,
                             SDDS_STRING, NULL) ||
      0>SDDS_DefineColumn(resultsPage, "Time", "t", "s",
                          "Time relative to synchronous phase at zero current",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(resultsPage, "Density", "$gr$r", "C/m",
                          "Charge density",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(resultsPage, "Current", "I", "A",
                          "Instantaneous current",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(resultsPage, "WakeField", NULL, "V",
                          "Wake field for particular bunch distribution",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(resultsPage, "PotentialWellDistortion", NULL, NULL,
                          "Wake potential distortion term in distribution",
                          NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(resultsPage, "PotentialWell", NULL, NULL,
                          "Total Wake potential term in distribution",
                          NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteLayout(resultsPage)|| !SDDS_StartPage(resultsPage, points))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void getInitialDensity( FUNCTION *density, long points, double deltaTime, 
                       double charge, double length) {
  long i;
  double C;
  
  density->points = points;
  density->xDelta = deltaTime;
  density->y = SDDS_Malloc( sizeof(*density->y) * density->points);
  density->offset =  - ((long) points/2);
  density->xStart = density->offset * density->xDelta;
  density->xFactor = 1;  /* using seconds */
  density->yFactor = 1;  /* using C/s, i.e. density is really a current */

  C = charge/ sqrt(2*PI)/ length;
  for (i=0; i<points;i++) {
    density->y[i] = C * exp( -sqr( (i + density->offset) * density->xDelta) /
                            2.0 / sqr(length));
  }
}

void copyDensity( FUNCTION *target, FUNCTION *source) {
  long i;
  double *ptr;
  
  ptr = target->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we had to preserve. 
     */
  *target = *source;
  target->y = ptr;
  if (!target->y) {
    target->y = SDDS_Malloc( sizeof(*target->y) * source->points);
  }
  for (i=0; i<source->points;i++) {
    target->y[i] = source->y[i];
  }
}

void getPotentialDistortion( FUNCTION *potentialDistortion,
                  FUNCTION *Vind,
                  FUNCTION *density, FUNCTION *wake) {
  long i, j, index;
  double *ptr;
  double P0;
  
  ptr = potentialDistortion->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *potentialDistortion = *density;
  potentialDistortion->y = ptr;
  /* the potential well should have the same range as the density,
     i.e. the abscissa of both functions line up. */
  if (!potentialDistortion->y) {
    potentialDistortion->y = SDDS_Malloc( sizeof(*potentialDistortion->y) * potentialDistortion->points);
  }
  ptr = Vind->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *Vind = *density;
  Vind->y = ptr;
  /*  Vind should have the same range as the density,
     i.e. the abscissa of both functions line up. */
  if (!Vind->y) {
    Vind->y = SDDS_Malloc( sizeof(*Vind->y) * Vind->points);
  }

/* calculate induced voltage. */
 for (i=0; i<density->points;i++) {
    Vind->y[i] = 0.0;
    for (j=0; j<density->points;j++) {
      /* start integrating with the first element of wake->y even though
         it may correspond to negative time, i.e. non-causal data which
         may be zero, btw. */
      index = i - wake->offset - j;
      if ( index<0 ) break;
      Vind->y[i] += wake->y[index] * density->y[j] * density->xDelta;
    }
  }

/* integrate Vind to get potentialDistortion */
  if (potentialDistortion->offset > 0) {
    bomb( "problem with offset value of FUNCTION potentialDistortion.", NULL);
  }
  potentialDistortion->y[0] = 0.0;
  for (i=1; i<potentialDistortion->points; i++) {
    potentialDistortion->y[i] = potentialDistortion->y[i-1] + potentialDistortion->xDelta * 
      (Vind->y[i-1] + Vind->y[i]) / 2.0;
  }
  
/* potential at synchronous phase and middle of original
   density distribution is defined to be zero */
/*
  P0 = potentialDistortion->y[-potentialDistortion->offset];
  for (i=0; i<potentialDistortion->points; i++) {
    potentialDistortion->y[i] -= P0;
  }
*/
}

void getPotentialDistortionFromModel( FUNCTION *potentialDistortion,
                           FUNCTION *Vind,
                           FUNCTION *density, double L, double R) {
  long i;
  double *ptr, *charge;

  ptr = potentialDistortion->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *potentialDistortion = *density;
  potentialDistortion->y = ptr;
  /* the potential well should have the same range as the density,
     i.e. the abscissa of both functions line up. */
  if (!potentialDistortion->y) {
    potentialDistortion->y = SDDS_Malloc( sizeof(*potentialDistortion->y) * potentialDistortion->points);
  }

  ptr = Vind->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *Vind = *density;
  Vind->y = ptr;
  /* Vind should have the same range as the density,
     i.e. the abscissa of both functions line up. */
  if (!Vind->y) {
    Vind->y = SDDS_Malloc( sizeof(*Vind->y) * Vind->points);
  }
  
  charge = SDDS_Malloc( sizeof(*charge) * density->points);
  charge[0] = 0.0;
  for (i=1; i<potentialDistortion->points; i++) {
    charge[i] = charge[i-1] + density->xDelta *
      (density->y[i-1] + density->y[i]) / 2.0;
  }
  
  for (i=0; i<potentialDistortion->points; i++) {
    potentialDistortion->y[i] = L * density->y[i] + R * charge[i];
  }
  for (i=1; i<potentialDistortion->points - 1; i++) {
    Vind->y[i] = L * (density->y[i+1] - density->y[i-1])/ 2.0 
      /density->xDelta + R * density->y[i];
  }
  /* end points */
  Vind->y[0] = L * (density->y[1] - density->y[0])/
    density->xDelta + R * density->y[0];
  Vind->y[potentialDistortion->points - 1] = 
    L * (density->y[potentialDistortion->points - 1] - 
         density->y[potentialDistortion->points - 2])/
      density->xDelta + R * density->y[potentialDistortion->points - 1];
  free(charge);
}

void calculateDistribution( FUNCTION *distribution, FUNCTION *potential,
                          FUNCTION *potentialDistortion,
                          double length, double VrfDot) {
  long i;
  double *ptr, time;

  ptr = potential->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *potential = *potentialDistortion;
  potential->y = ptr;
  /* the potential should have the same abscissa range as the density,
     i.e. the abscissa of both functions line up. */
  if (!potential->y) {
    potential->y = SDDS_Malloc( sizeof(*potential->y) * potential->points);
  }

  ptr = distribution->y;
  /* This assignment transfers all values, including
     the pointer to the double array, which we preserved.
     */
  *distribution = *potentialDistortion;
  distribution->y = ptr;
  /* the distribution should have the same abscissa range as the density,
     i.e. the abscissa of both functions line up. */
  if (!distribution->y) {
    distribution->y = SDDS_Malloc( sizeof(*distribution->y) * distribution->points);
  }
  for (i=0; i<potential->points; i++) {
    time = distribution->xStart + i * distribution->xDelta;
    potential->y[i] = sqr(time)/ 2.0/ sqr(length) + 
      1.0/ VrfDot/ sqr(length) * potentialDistortion->y[i];
    distribution->y[i] = exp( - potential->y[i]);
  }
}

void normalizeDensityFunction( FUNCTION *density, FUNCTION *distribution, double charge) {
  long i;
  double integral;
  
  /* integrate exp(-H) to find normalization constant */
  integral = 0.0;
  for (i=0; i<distribution->points; i++) {
    integral += distribution->y[i];
  }
  integral *= distribution->xDelta;
  for (i=0; i<distribution->points; i++) {
    density->y[i] = (charge/integral) * distribution->y[i];
  }
}

void writeResults( SDDS_TABLE *resultsPage, FUNCTION *density, 
                  FUNCTION *potential, FUNCTION *potentialDistortion, 
                  FUNCTION *Vind, double charge, double averageCurrent,
                  long converged) {
  long i;
  double *current, *actualDensity, *time;

  time = SDDS_Malloc( sizeof(*time) * density->points);
  current = SDDS_Malloc( sizeof(*current) * density->points);
  actualDensity = SDDS_Malloc( sizeof(*actualDensity) * density->points);
  for (i=0;i<density->points;i++) {
    time[i] = density->xStart + i * density->xDelta;
    current[i] = density->y[i] * density->yFactor; /* convert to A units */
    actualDensity[i] = density->y[i] * density->yFactor / c_mks; /* convert to C/m units */
  }
  
  if (!SDDS_SetParameters(resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Convergence", converged?"Solution converged":"Solution did not converge",
                          "Charge", charge, 
                          "AverageCurrent", averageCurrent, NULL)||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, time, density->points, "Time") ||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, current, density->points, "Current")  ||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, actualDensity, density->points, "Density") ||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, Vind->y, density->points, "WakeField") ||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, potential->y, density->points, "PotentialWell") ||
      !SDDS_SetColumn(resultsPage, SDDS_SET_BY_NAME, potentialDistortion->y, density->points, "PotentialWellDistortion"))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if ( !SDDS_WritePage(resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free(time);
  free(current);
  free(actualDensity);
}

void printFunction( char *label, FUNCTION *data) {
  long i;

  fprintf( stderr, "Structure %s:\n", label);
  fprintf( stderr, "\t%s.xStart: %g\n", label, data->xStart);
  fprintf( stderr, "\t%s.xDelta: %g\n", label, data->xDelta);
  fprintf( stderr, "\t%s.points: %ld\n", label, data->points);
  fprintf( stderr, "\t%s.offset: %ld\n", label, data->offset);
  for (i=0;i<3;i++) {
    fprintf( stderr, "\t%s.y[%ld]: %g\n", label, i, data->y[i]);
  }
  fprintf( stderr, "\t...\n");
  for (i=data->points/2-2;i<data->points/2+3;i++) {
    fprintf( stderr, "\t%s.y[%ld]: %g\n", label, i, data->y[i]);
  }
  fprintf( stderr, "\t...\n");
  for (i=data->points-3;i<data->points;i++) {
    fprintf( stderr, "\t%s.y[%ld]: %g\n", label, i, data->y[i]);
  }
}

void makeBBRWakeFunction(FUNCTION *wake, double dt, long points, 
                         double Q, double R, double omega, double rw, double T0)
{
  double Qp, omegap, t;
  long i;
  
  wake->xStart = 0;
  wake->xDelta = dt;
  wake->points = points;
  wake->offset = 0;
  wake->xFactor = wake->yFactor = 1;
  
  Qp = sqrt(Q*Q-0.25);
  omegap = omega*Qp/Q;
  
  wake->y = tmalloc(sizeof(*(wake->y))*points);
  for (i=0; i<points; i++) {
    t = i*dt;
    wake->y[i] = (omega*R/Q)*exp(-omega*t/(2*Q))*(cos(omegap*t) - sin(omegap*t)/(2*Qp));
  }
  if (rw)
    wake->y[0] += rw/dt;
}

