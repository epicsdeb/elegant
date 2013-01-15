/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsbrightness.c
 * purpose: take twiss output from elegant and compute undulator
 *          brightness curves
 *
 * Michael Borland, 2002
 *
 $Log: not supported by cvs2svn $
 Revision 1.19  2011/02/09 20:26:33  shang
 Per Roger Dejus's fix, added statement to ensure that the emittance  down-shifted peak is less then the zero-emittance calculated energy. The down-shift is set to be  always larger than 5e-4. Errors may otherwise occur for very small emittances and long devices (large N).

 Revision 1.18  2011/01/05 22:24:48  borland
 Added ability to take error factors from a file and multiply each harmonic by the specified factor.

 Revision 1.17  2010/12/15 15:02:31  borland
 In the event that the brightness can't be computed, the corresponding rows are now filled in with 0.

 Revision 1.16  2010/01/04 23:22:08  borland
 Resolved some issues on 64 bit operating systems.  Was getting segmentation faults.

 Revision 1.15  2009/06/02 17:56:50  borland
 Fixed units for brightness.

 Revision 1.14  2009/04/14 13:19:48  borland
 Now accepts units of "m" for ex0.

 Revision 1.13  2007/03/09 19:38:02  borland
 Fixed bug in energy range calculation for helical devices.

 Revision 1.12  2006/08/24 19:20:06  soliday
 Updated so that it would compile on WIN32 again.

 Revision 1.11  2005/11/10 16:19:03  soliday
 Updated to compile with a 64bit compiler.

 Revision 1.10  2005/11/04 16:27:07  borland
 Added Xiao's code for space charge element SCMULT.

 Revision 1.8  2005/08/23 17:40:49  shang
 the pCentral units now can be m$be$nc, or NULL, or blank.

 Revision 1.7  2005/08/08 14:48:43  borland
 Changed a loop termination condition to prevent use of negative indices.

 Revision 1.6  2005/06/10 17:46:41  borland
 Changes by Shang: now optionally uses data produced by sddsanalyzebeam.
 Also, requires -coupling or -emittanceRatio if vertical emittance data is
 not given.

 Revision 1.5  2004/09/28 13:49:26  borland
 Added debugging statements.

 Revision 1.4  2004/09/28 13:21:26  borland
 Added feature whereby the parameter ey0 is accepted for the vertical emittance.
 Fixed bug with multiplication of the current by 1000 for each pass through the main
 loop.

 Revision 1.3  2004/07/06 17:18:20  borland
 Default method is now actually Dejus, consistent with the usage message.
 Previously, default method was Borland.

 Revision 1.2  2004/07/06 15:43:47  borland
 Fixed bug in Dejus method implementation: the gamma value was not being used consistently.
 Rather, in some places, a fixed energy value of 7 GeV was still used.

 Revision 1.1  2004/04/08 16:10:28  soliday
 Moved to subdirectory.

 Revision 1.8  2003/06/25 21:40:01  shang
 using pCentral as gamma instead of calculating it from the energy in dejus' method

 Revision 1.7  2003/04/17 15:20:15  soliday
 Fixed for WIN32

 Revision 1.6  2003/02/24 20:00:11  shang
 added -method option which includes dejus and walker's method for calculating
 the on-axis brightness.

 Revision 1.5  2002/08/14 20:23:48  soliday
 Added Open License

 Revision 1.4  2002/05/07 20:02:24  shang
 modified the computation of brightness by using convolution when the width of
 sinc() function is wider than that of energ spread (gaussian func)

 Revision 1.3  2002/04/24 14:26:38  borland
 Refined calculation of radiation size and divergence.  Added
 -noSpectralBroadening switch to allow turning off this part of the
 calculation.

 Revision 1.2  2002/04/22 21:01:01  borland
 Added factors of 2 for square radiation opening angle and square of
 radiation size, per Dejus.
 Refined energy-broadening effect by fitting a gaussian to the sinc function
 to get a width of xn=0.36 for sigma.

 Revision 1.1  2002/04/18 23:41:57  borland
 First version, with assistence of L. Emery.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"
#include "sddsbrightness.h"

#define SET_PIPE 0
#define SET_HARMONICS 1
#define SET_KRANGE 2
#define SET_CURRENT 3
#define SET_TOTALLENGTH 4
#define SET_PERIODLENGTH 5
#define SET_EMITTANCERATIO 6
#define SET_COUPLING 7
#define SET_NOSPECTRALBROADENING 8
#define SET_METHOD 9
#define SET_ERRORFACTORS 10
#define N_OPTIONS 11

char *option[N_OPTIONS] = {
  "pipe", "harmonics", "krange", "current", "totallength", "periodlength", 
  "emittanceratio", "coupling", "nospectralbroadening", "method", "errorFactors"
} ;

char *USAGE1="sddsbrightness [-pipe=[input][,output]] [<twissFile>] [<SDDSoutputfile>]\n\
 -harmonics=<integer> -Krange={start=<value>,end=<value>,points=<integer> | filename=<string>,column=<string>}\n\
 -current=<Amps> -totalLength=<meters> -periodLength=<meters>\n\
 [-errorFactors=<filename>,<columnName>] \n\
 [-emittanceRatio=<value> | -coupling=<value>] [-noSpectralBroadening]\n\
 [-method=<string value>,device=<string value>,neks=<value>]] \n\n\
harmonics        number of harmonics to compute\n\
Krange           range and number of undulator K parameter to evaluate or a list of K values provided by a file.\n\
current          beam current in amperes\n\
totalLength      total length of the undulator in meters\n\
periodLength     length of the undulator period in meters\n\
errorFactors     Specify SDDS file and column for factors by which to multiply\n\
                 the brightness for each harmonic, as a way of accounting for\n\
                 the effect of field errors.  Row 0 gives h=1 data, Row 1 gives h=3\n\
                 data, etc.";
char *USAGE2="emittanceRatio   ratio of y emittance to x emittance.  x emittance is\n\
                 ex0 from input file. y emittance is ratio*ex0\n\
                 Ignored if ey0 parameter or ey column is found in file.\n\
coupling         x emittance is ex0/(1+coupling), while y emittance is\n\
                 coupling*ex0/(1+coupling).\n\
                 Ignored if ey0 parameter or ey column is found in file.\n\
noSpectralBroadening\n\
                 Turns off the default inclusion of spectral broadening in\n\
                 the calculation.  Gives an over-estimate of the brightness.\n\
method           choose method for calculating brightness \n\
                 method=borland  Michael Borland's approximation method. \n\
                 method=dejus    Non-zero emittance, \n\
                                 infinite-N +convolution (Dejus' approach) (default) \n\
                 method=walkerinfinite        Non-zero emittance, \n\
                                              infinite-N +convolution (Walker's approach) \n\
                 method=walkerfinite          Non-zero emittance; finite-N (Walker's) \n\
                 device=planar or helical      undulator type. \n\
                 neks=<value> number of points for peaking search. \n\n\
Computes on-axis brightness for an undulator centered on the end of the beamline\n\
the Twiss parameters for which are in the input file.  You should generate the\n\
input file using elegant's twiss_output command with radiation_integrals=1 .\n\n\
Program by Michael Borland.  ANL(This is version 1.20, "__DATE__")\n";

void getErrorFactorData(double **errorFactor, long *errorFactors, char *filename, char *columnName);
long SetUpOutputFile(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSin, char *outputfile, long harmonics);
double ComputeBrightness(double period, long Nu, double K, long n,
                         double gamma, double ex, double ey, double Sdelta0, 
                         double current, 
                         double betax, double alphax, double etax, double etaxp,
                         double betay, double alphay, double etay, double etayp,
                         double *lambda, double *Fn, 
                         short spectralBroadening);
long GetTwissValues(SDDS_DATASET *SDDSin, 
                    double *betax, double *alphax, double *etax, double *etaxp, 
                    double *betay, double *alphay, double *etay, double *etayp, 
                    double *ex0, double *ey0, double *Sdelta0, double *pCentral, double emitRatio,
		    double coupling);
double computeFactorOfConvolution(long periods, long harmonic, double Sdelta0);
double convolutionFunc(double x);
double delta0,sincNu; /*two constants used in convolutionFunc() */

/*following functions are needed for calculating brightness using Dejus's method */
void FindPeak(double *E,double *spec,double *ep,double *sp,int32_t n);
int Gauss_Convolve(double *E,double *spec,int32_t *ns,double sigmaE);
/* note that sigmaE=Sdelta0 */
void Dejus_CalculateBrightness(double current,long nE,
                               double period_mks, long nP, long device,
                               long ihMin,long ihMax,long ihStep,double sigmaE,
                               double gamma, double ex0, double ey0, 
                               double betax, double alphax, double etax, double etaxp,
                               double betay, double alphay, double etay, double etayp,
                               long minNEKS, long maxNEKS, long neks,
			       double *Kvalue,
                               double kMin,double kMax, long method,
                               double *sigmax,double *sigmay,double *sigmaxp,double *sigmayp,
                               double **K, double ***FnOut,
                               double ***Energy, double ***Brightness, double ***LamdarOut);
void ComputeBeamSize(double period, long Nu, double ex0, double ey0, double Sdelta0, 
                     double betax, double alphax, double etax, double etaxp,
                     double betay, double alphay, double etay, double etayp,
                     double *Sx, double *Sy, double *Sxp, double *Syp);
void getKValueData(double **Kvalue, int32_t *Kpoints, char *filename, char *columnName);
/*fortran subroutine*/
void usb_();
#if defined(_WIN32)
#define usb_() USB()
#endif

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  char *Kfilename=NULL, *Kcolumn=NULL;
  double current, totalLength, periodLength, KStart, KEnd, coupling, emittanceRatio, dK=0, *Kvalue=NULL;
  int32_t KPoints;
  long harmonics, tmpFileUsed, iK, i_arg, readCode, h, ih, ij, periods;
  unsigned long dummyFlags;
  double betax, alphax, betay, alphay, etax, etaxp, etay, etayp, lambda, energy;
  double pCentral, ex0, ey0, Bn, Fn, K, Sdelta0, conFactor=1.0; /*the factor from convolution of
								  sinc() and gaussian() */
  short spectralBroadening;
  long method, device,nE,ihMin,ihMax;
  double *KK,**FnOut,**Energy,**Brightness,**LamdarOut;
  long minNEKS,maxNEKS;
  int32_t neks;
  double sigmax,sigmay,sigmaxp,sigmayp;
  char *deviceOption, *Units=NULL;
  double *errorFactor = NULL;
  long errorFactors = 0;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s%s\n", USAGE1, USAGE2);
    exit(1);
  }
  KK=NULL;
  FnOut=Energy=Brightness=LamdarOut=NULL;
  deviceOption=NULL;

  inputfile = outputfile = NULL;
  pipeFlags = dummyFlags=0;
  spectralBroadening = 1;
  current = totalLength = periodLength = 0;
  KStart = KEnd = coupling = emittanceRatio = 0;
  harmonics = KPoints = 0;
  method=1;
  device=0;
  minNEKS=100;
  maxNEKS=500;
  neks=100;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_HARMONICS:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%ld", &harmonics) ||
            harmonics<1)
          SDDS_Bomb("invalid -harmonics value");
        break;
      case SET_KRANGE:
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "start", SDDS_DOUBLE, &KStart, 1, 0,
                          "end", SDDS_DOUBLE, &KEnd, 1, 0,
                          "points", SDDS_LONG, &KPoints, 1, 0,
			  "filename", SDDS_STRING, &Kfilename, 1, 0, 
			  "column", SDDS_STRING, &Kcolumn, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -Krange syntax");
	if (!Kfilename && (KStart>=KEnd || KStart<0 || KPoints<2))
	  SDDS_Bomb("invalid -Krange start/end/points values");
	if (Kfilename && !Kcolumn)
	  SDDS_Bomb("invalid -Krange syntax, column name not provided,");
        break;
      case SET_CURRENT:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &current) ||
            current<=0)
          SDDS_Bomb("invalid -current value");
        break;
      case SET_TOTALLENGTH:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &totalLength) ||
            totalLength<=0)
          SDDS_Bomb("invalid -totalLength value");
        break;
      case SET_PERIODLENGTH:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &periodLength) ||
            periodLength<=0)
          SDDS_Bomb("invalid -periodLength value");
        break;
      case SET_EMITTANCERATIO:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &emittanceRatio) ||
            emittanceRatio<=0)
          SDDS_Bomb("invalid -emittanceRatio value");
        break;
      case SET_COUPLING:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &coupling) ||
            coupling<=0)
          SDDS_Bomb("invalid -coupling value");
        break;
      case SET_NOSPECTRALBROADENING:
        spectralBroadening = 0;
        break;
      case SET_ERRORFACTORS:
        if (s_arg[i_arg].n_items!=3)
          SDDS_Bomb("invaliud -errorFactors syntax/values");
        getErrorFactorData(&errorFactor, &errorFactors, s_arg[i_arg].list[1], s_arg[i_arg].list[2]);
        break;
      case SET_METHOD:
        if (s_arg[i_arg].n_items<2) 
          SDDS_Bomb("invalid -dejus syntax/values");
        switch (match_string(s_arg[i_arg].list[1], method_option, METHOD_OPTIONS, 0)) {
        case BORLAND:
          method=0;
          break;
        case DEJUS:
          method=1;
          break;
        case WALKER_INF:
          method=2;
          break;
        case WALKER_FIN:
          method=3;
          break;
        default:
          SDDS_Bomb("Invalid method given!");
          break;
        }
        s_arg[i_arg].n_items -=2;
        if (s_arg[i_arg].n_items>0 &&
            !scanItemList(&dummyFlags, s_arg[i_arg].list+2, &s_arg[i_arg].n_items, 0,
                          "device", SDDS_STRING, &deviceOption, 1, 0,
                          "neks", SDDS_LONG, &neks, 1, 0, NULL))
          SDDS_Bomb("invalid -method syntax/values");
        if (deviceOption) {
          switch (match_string(deviceOption, device_option, DEVICE_OPTIONS, 0)) {
          case PLANAR:
            device=0;
            break;
          case HELICAL:
            /*note that helical device has only one harmonics */
            harmonics=1;
            device=1;
            break;
          default:
            SDDS_Bomb("unknow device given!");
            break;
          }
          break;
        default:
          fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
          fflush(stdout);
          exit(1);
          break;
        }
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
#ifdef DEBUG
  fprintf(stderr, "Argument parsing done.\n");
#endif
  if (Kfilename) {
    getKValueData(&Kvalue, &KPoints, Kfilename, Kcolumn);
    /*Kvalue is sorted in decreasing order */
    KStart = Kvalue[KPoints-1];
    KEnd = Kvalue[0];
    if (KStart>=KEnd || KStart<0 || KPoints<2)
      SDDS_Bomb("Invalid K values provided.");
  }
  
  if (errorFactor && errorFactors<harmonics)
    SDDS_Bomb("too few error factors for requested number of harmonics");
  if (coupling && emittanceRatio)
    SDDS_Bomb("give only one of -coupling or -emittanceRatio");
  if (!harmonics)
    SDDS_Bomb("you must specify the number of harmonics to compute");
  if (!KPoints)
    SDDS_Bomb("you must specify the range and number of K values");
  if (!current)
    SDDS_Bomb("you must specify the current");
  if (!totalLength)
    SDDS_Bomb("you must specify the total undulator length");
  if (!periodLength)
    SDDS_Bomb("you must specify the undulator period length");
  /*poles = totalLength/periodLength; rename poles as periods, which makes sense.*/
  periods=totalLength/periodLength;
  if (periods<1)
    SDDS_Bomb("period lengths is shorter than undulator length!");
  
  processFilenames("sddsbrightness", &inputfile, &outputfile, pipeFlags, 0, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
#ifdef DEBUG
  fprintf(stderr, "Checking input file...\n");
#endif
  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (SDDS_CheckColumn(&SDDSin, "betax", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "alphax", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "etax", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etaxp", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "betay", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckColumn(&SDDSin, "alphay", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etay", "m", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etayp", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    SDDS_Bomb("something wrong with twiss parameter columns");
  }
  if (SDDS_CheckParameter(&SDDSin, "ex0", "$gp$rm", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK &&
      SDDS_CheckParameter(&SDDSin, "ex0", "m", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK) {
    if (SDDS_CheckColumn(&SDDSin,"ex", "$gp$rm", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK &&
        SDDS_CheckColumn(&SDDSin,"ex", "m", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK) {
      SDDS_Bomb("Something wrong with both ex0 parameter and ex column, one of them has to exist");
    }
  }
  if (SDDS_CheckParameter(&SDDSin, "pCentral", NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    if (!SDDS_GetParameterInformation(&SDDSin, "units", &Units, SDDS_GET_BY_NAME, "pCentral"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (Units && !is_blank(Units) && strcmp(Units, "m$be$nc"))
      SDDS_Bomb("Invalid units of pCentral parameter");
  } else if (SDDS_CheckColumn(&SDDSin, "pCentral", NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    if (!SDDS_GetColumnInformation(&SDDSin, "units", &Units, SDDS_GET_BY_NAME, "pCentral"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (Units && !is_blank(Units) && strcmp(Units, "m$be$nc"))
      SDDS_Bomb("Invalid units of pCentral column");
  } else 
    SDDS_Bomb("Something wrong with both pCentral parameter and pCentral column, one of them has to exist");
  if (SDDS_CheckParameter(&SDDSin, "Sdelta0", "", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK &&
      SDDS_CheckColumn(&SDDSin, "Sdelta", "", SDDS_ANY_FLOATING_TYPE, NULL)!=SDDS_CHECK_OK) {
    SDDS_Bomb("Something wrong with both Sdelta0 parameter and Sdelta column, one of them has to exist");
  }
  /*
  if (SDDS_CheckParameter(&SDDSin, "ex0", "$gp$rm", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK || 
      SDDS_CheckParameter(&SDDSin, "pCentral", "m$be$nc", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK ||
      SDDS_CheckParameter(&SDDSin, "Sdelta0", "", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    SDDS_Bomb("something wrong with ex0, pCentral, or Sdelta0 parameters");
  } */

#ifdef DEBUG
  fprintf(stderr, "Setting up output file...\n");
#endif
  if (!SetUpOutputFile(&SDDSout, &SDDSin, outputfile, harmonics))
    SDDS_Bomb("problem setting up output file");
  
  dK = (KEnd-KStart)/(KPoints-1);
  if (method)
    current=current*1.0e3;    /*change unit from A to mA for dejus's method */
#ifdef DEBUG
  fprintf(stderr, "Entering main loop\n");
#endif
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (!GetTwissValues(&SDDSin, 
                        &betax, &alphax, &etax, &etaxp, 
                        &betay, &alphay, &etay, &etayp, 
                        &ex0, &ey0, &Sdelta0, &pCentral, emittanceRatio, coupling))
      SDDS_Bomb("problem getting twiss parameters and other values from input file");
#ifdef DEBUG
    fprintf(stderr, "Twiss values read in.\n");
#endif
    if (!SDDS_StartPage(&SDDSout, KPoints))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
#ifdef DEBUG
    fprintf(stderr, "Started output page.\n");
#endif
    if (method) {
      ihMin=1;
      ihMax=2*(harmonics-1)+1;
      nE=KPoints;

#ifdef DEBUG
      fprintf(stderr, "Calling Dejus_CalculateBrightness\n");
#endif
      Dejus_CalculateBrightness(current,nE,periodLength, periods,device,ihMin,ihMax,2,Sdelta0,
                                pCentral,ex0,ey0,
                                betax,alphax,etax,etaxp,betay,alphay,etay,etayp,
                                minNEKS,maxNEKS,neks, Kvalue, KStart, KEnd,method,
                                &sigmax,&sigmay,&sigmaxp,&sigmayp,
                                &KK, &FnOut,&Energy,&Brightness,&LamdarOut);
#ifdef DEBUG
      fprintf(stderr, "Returned from Dejus_CalculateBrightness\n");
#endif
      if (errorFactor) {
        for (ih=0; ih<harmonics; ih++) {
          long i;
          for (i=0; i<nE; i++)
            Brightness[ih][i] *= errorFactor[ih];
        }
      }
      for (ih=0; ih<harmonics; ih++) {
        h = ih*2+1;
        if (h==1 && !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,KK,nE,0))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (!SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,Brightness[ih],nE,ih*4+1) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,FnOut[ih],nE,ih*4+2) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,LamdarOut[ih],nE,ih*4+3) ||
            !SDDS_SetColumn(&SDDSout,SDDS_SET_BY_INDEX,Energy[ih],nE,ih*4+4))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (!SDDS_SetParameters(&SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"current",current,
			      "EnergySpread",Sdelta0,"sigmax",sigmax,"sigmay",sigmay,
			      "sigmaxprime",sigmaxp,"sigmayprime",sigmayp,
			      "period",periodLength,"numberOfPeriods",periods,
			      "emitx",ex0,"emity",ey0,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {    
      for (ih=0; ih<harmonics; ih++) {
        h = ih*2+1;
        /*convolution factor is same for all K s, it varies on periods,h and Sdelta0 only*/
        conFactor=1.0; 
        if (spectralBroadening && Sdelta0!=0 && 0.36/periods/h>Sdelta0) {
          conFactor=computeFactorOfConvolution(periods,h,Sdelta0);
        }
        for (K=KStart, iK=0; iK<KPoints; iK++, K+=dK) {
	  if (Kfilename)
	    K = Kvalue[iK];
          Bn = conFactor*ComputeBrightness(periodLength, periods, K, h,
                                           pCentral, ex0, ey0, Sdelta0, current,
                                           betax, alphax, etax, etayp, 
                                           betay, alphay, etay, etaxp, 
                                           &lambda, &Fn, spectralBroadening);
          energy = 12.39/(lambda*1e10);
          if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iK,
                                 ih*4+1, Bn,
                                 ih*4+2, Fn,
                                 ih*4+3, lambda, 
                               ih*4+4, energy, -1) ||
              (h==1 &&
               !SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, iK, 0, K, -1)))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      if (!SDDS_SetParameters(&SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"current",current,
			      "EnergySpread",Sdelta0,
			      "period",periodLength,"numberOfPeriods",periods,
			      "emitx",ex0,"emity",ey0,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    
#ifdef DEBUG
    fprintf(stderr, "Writing output page.\n");
#endif
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
#ifdef DEBUG
    fprintf(stderr, "Exiting main loop.\n");
#endif

  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (Kfilename) {
    free(Kvalue);
    free(Kfilename);
    free(Kcolumn);
  }
  free_scanargs(&s_arg,argc);
  if (method) {
    for (ih=0;ih<harmonics;ih++) {
      free(Energy[ih]);
      free(FnOut[ih]);
      free(Brightness[ih]);
      free(LamdarOut[ih]);
    }
    free(Energy);
    free(FnOut);
    free(Brightness);
    free(LamdarOut);
    free(KK);
  }
  return 0;
}

long SetUpOutputFile(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSin, char *outputfile, long harmonics)
{
  long h;
  char buffer[1024];
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      !SDDS_DefineSimpleColumn(SDDSout, "K", NULL, SDDS_DOUBLE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(SDDSout,"current",NULL, "mA", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"EnergySpread",NULL, NULL, NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmax",NULL, "mm", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmay",NULL, "mm", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmaxprime",NULL,"mrad", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"sigmayprime",NULL,"mrad", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"period",NULL, "m", NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"numberOfPeriods",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"emitx","m",NULL,NULL,NULL,SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout,"emity","m",NULL,NULL,NULL,SDDS_DOUBLE, 0)<0 )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (h=1; h<2*harmonics; h+=2) {
    sprintf(buffer, "Brightness%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "photons/s/mrad$a2$n/mm$a2$n/0.1%BW", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "F%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, NULL, SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "wavelength%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "m", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sprintf(buffer, "photonEnergy%ld", h);
    if (!SDDS_DefineSimpleColumn(SDDSout, buffer, "keV", SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 1;
}

double ComputeBrightness(double period, long Nu, double K, long n, 
                         double gamma, double ex, double ey, double Sdelta0, 
                         double current, 
                         double betax, double alphax, double etax, double etaxp,
                         double betay, double alphay, double etay, double etayp,
                         double *lambdaOut, double *FnOut,
                         short spectralBroadening)
{
  double JArg, Nn;
  double Sx, Sy, Sxp, Syp, Srp2, Sr2;
  double lambda, Fn;
  double gammax, gammay, length, broadening;
  double sincWidth;
  
  JArg = n*K*K/(4+2*K*K);
  Fn = sqr(n*K/(1+K*K/2)*(jn((n+1)/2, JArg) - jn((n-1)/2, JArg)));
  
  lambda = period/(2*n*sqr(gamma))*(1+sqr(K)/2);

  /* 0.001 is for 0.1% bandwidth */
  Nn = 0.5*PI/137.04*Nu*(current/e_mks)*0.001*Fn*(1+sqr(K)/2)/n;
  length = Nu*period;
  sincWidth=0.36/Nu/n;
  
  broadening = 1;
  if (spectralBroadening) {
    /* this factor includes the loss of peak brightness due to spectral
     * broadening from the energy spread. The factor of 2 in front of Sdelta0
     * is because wavelength ~ 1/gamma^2 .  The 0.36 is from a gaussian fit
     * to the sinc function, which gives sigma(omega)/omega = 0.36/(n*Nu)
     */
    broadening = sqrt( 1 + sqr(2*Sdelta0*n*Nu/0.36) );
  }

  /* radiation divergence and size, squared.
   * numerical factors are from "Radiation Properties of an Undulator",
   * SLS-Note 4/95, by W. Joho
   */
  Srp2 = sqr(0.58)*lambda/length;   /* Joho has *broadening here, which I think is wrong */
  Sr2  = sqr(0.12/0.58)*lambda*length;

  gammax = (1+sqr(alphax))/betax;
  gammay = (1+sqr(alphay))/betay;
  Sxp = sqrt(ex*gammax + sqr(Sdelta0*etaxp) + Srp2);
  Syp = sqrt(ey*gammay + sqr(Sdelta0*etayp) + Srp2);
  Sx = sqrt(ex*betax + sqr(Sdelta0*etax) + Sr2);
  Sy = sqrt(ey*betay + sqr(Sdelta0*etay) + Sr2);
  
  *lambdaOut = lambda;
  *FnOut = Fn;
  /* 1e-12 is to convert to 1/(mm^2*mrad^2) units */
  /* Joho does not put the broadening here, which I think is wrong */
  if (Sdelta0!=0 && sincWidth>Sdelta0) 
    broadening=1.0; /*calculate the factor by convolution instead */
  return Nn/(sqr(PIx2)*Sx*Sxp*Sy*Syp)/broadening*1e-12;
}

/*this function is to get the brightness by the convolution of sinc() and
  gaussian function (the energy spreed).
  since sinc() is symmetric, the peak position of the convolution is the same
  as sinc(), that is, at x=0, thus, the convolution is calculated only at
  one point: x=0, and returned as the brightness */
/* S(x)=(sin(N*pi*x)/(N*pi*x))^2
   g(x)=Cexp(-1/8 * (x/Sdelta0)^2) */
/* g(x)=1.0*10^-7 while x=+- 4*Sdelta0*log(1e-7/C) (the range of convoluation */
/* convoluation=integral of (g(x)*S(-x)dx, using gaussianQuadrature() to calculate it */
double computeFactorOfConvolution(long periods, long harmonic, double Sdelta0)
{
  double startX,C, conFactor=1.0;
  C=1/sqrt(PI*2)*0.5/Sdelta0; /*theorectical value for C=1/sqrt(2*pi)*1/(2*Sdelat0) */
  startX=-(4*Sdelta0*sqrt(fabs(log(1e-7/C)))); 
  sincNu=harmonic*periods*PI; /*sincNu and delta0 are two globals*/
  delta0=Sdelta0;
  gaussianQuadrature(convolutionFunc,startX,fabs(startX),1,0.01,&conFactor); 
  /*fprintf(stderr,"integrated result=%.9f\n",conFactor); */
  return conFactor;
}
/*integral function of gaussianQuadrature() */
double convolutionFunc(double x)
{
  double multResult,C,sx,gx;
  C=1/sqrt(PI*2)*0.5/delta0;
  if (x==0)
    sx=1;
  else
    sx=sqr(sin(sincNu*x)/(sincNu*x));
  gx=exp(-(x*x)/(8*delta0*delta0));
  multResult=C*sx*gx;
  return multResult;
}


long GetTwissValues(SDDS_DATASET *SDDSin, 
                    double *betax, double *alphax, double *etax, double *etaxp, 
                    double *betay, double *alphay, double *etay, double *etayp, 
                    double *ex0, double *ey0, double *Sdelta0, double *pCentral, 
		    double emitRatio, double coupling)
{
  double *data;
  long rows, ey0Exist=0;

  if (!(rows=SDDS_RowCount(SDDSin)))
    return 0;

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "betax")))
    SDDS_Bomb("unable to get betax");
  *betax = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "alphax")))
    SDDS_Bomb("unable to get alphax");
  *alphax = data[rows-1];
  free(data);
  
  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etax")))
    SDDS_Bomb("unable to get etax");
  *etax = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etaxp")))
    SDDS_Bomb("unable to get etax");
  *etaxp = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "betay")))
    SDDS_Bomb("unable to get betay");
  *betay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "alphay")))
    SDDS_Bomb("unable to get alphay");
  *alphay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etay")))
    SDDS_Bomb("unable to get etay");
  *etay = data[rows-1];
  free(data);

  if (!(data=SDDS_GetColumnInDoubles(SDDSin, "etayp")))
    SDDS_Bomb("unable to get etay");
  *etayp = data[rows-1];
  free(data);
  if (SDDS_CheckParameter(SDDSin, "ex0", "$gp$rm", SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK ||
      SDDS_CheckParameter(SDDSin, "ex0", "m", SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    if (!SDDS_GetParameterAsDouble(SDDSin, "ex0", ex0))
      SDDS_Bomb("unable to get ex0 parameter from input file");
  } else {
    if (!(data=SDDS_GetColumnInDoubles(SDDSin, "ex")))
      SDDS_Bomb("unable to get ex");
    *ex0 = data[0];
    free(data);
  }
  if (SDDS_CheckParameter(SDDSin, "pCentral", NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    if (!SDDS_GetParameterAsDouble(SDDSin, "pCentral", pCentral))
      SDDS_Bomb("unable to get pCentral parameter from input file");
  } else {
    if (!(data=SDDS_GetColumnInDoubles(SDDSin, "pCentral")))
      SDDS_Bomb("unable to get pCentral");
    *pCentral = data[0];
    free(data);
  }
  if (SDDS_CheckParameter(SDDSin, "Sdelta0", "", SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    if (!SDDS_GetParameterAsDouble(SDDSin, "Sdelta0", Sdelta0))
      SDDS_Bomb("unable to get Sdelta0 parameter from input file");
  } else {
    if (!(data=SDDS_GetColumnInDoubles(SDDSin, "Sdelta")))
      SDDS_Bomb("unable to get pCentral");
    *Sdelta0 = data[0];
    free(data);
  }
  /* get ey0 if it is there */
  if (SDDS_GetParameterAsDouble(SDDSin, "ey0", ey0)) {
    ey0Exist=1;
  } else if ((data=SDDS_GetColumnInDoubles(SDDSin, "ey"))) {
    ey0Exist=1;
    *ey0 = data[0];
    free(data);
  }
  if (*ex0<=0.0)
    SDDS_Bomb("ex0 should be greater than zero.");
  if (!ey0Exist) {
    if (emitRatio==0 && coupling==0)
      SDDS_Bomb("No vertical emittance data in file: give -emittanceRatio or -coupling");
    if (coupling) {
      *ex0 = *ex0/(1+coupling);
      *ey0 = coupling*(*ex0);
    } else {
      *ey0 = *ex0*emitRatio;
    }
  }
  return 1;
}

void FindPeak(double *E,double *spec,double *ep,double *sp,int32_t n)
{
  long i;
  
  *sp=spec[0];
  *ep=E[0];
  for (i=1; i<n; i++) {
    if (*sp<spec[i]) {
      *sp=spec[i];
      *ep=E[i];
    }
  }
}

int Gauss_Convolve(double *E,double *spec,int32_t *ns,double sigmaE) 
{
  int32_t nSigma=3,nppSigma=6,ne1,ne2,ns1,np;
  int32_t i,j;
  double ep,sp,sigp,de,sum, *gs,x, *spec2;
	
  ns1=*ns;
  gs=spec2=NULL;
  
  if (!E || !spec) {
    fprintf(stderr,"No energy or spectra points!\n");
    return 1;
  }
  FindPeak(E,spec,&ep,&sp,ns1);
  
  /*generate Gaussian with correct sigma in units of x-axis */
  de=E[1]-E[0];
  sigp=2.0*sigmaE*ep/de; /*sigma in x-axis units */
  
  if (sigp < (nppSigma-1)) {
    fprintf(stderr,"too few data points for Gaussian convolution\n");
    return 1; 
  }
  np=(2*nSigma)*sigp+1;
  if (np%2==0) np=np+1; /* make odd */
  gs=(double*)calloc(np,sizeof(*gs));
  spec2=(double*)calloc(ns1,sizeof(*spec2));
  sum=0.0;
  for (i=0;i<np;i++) {
    x=i*1.0-0.5*(np-1);
    gs[i]=exp(-x*x/2.0/(sigp*sigp));
  /*  fprintf(stderr,"x=%e, gs=%e\n",x,gs[i]); */
    sum=sum+gs[i];
  }
  /*fprintf(stderr,"sigp=%e,nSigma=%d\n",sigp,nSigma); */
  
  
  /*make convolution */
  ne1=np/2;
  ne2=ns1-ne1-1;
  if (ne2<0) {
    fprintf(stderr,"Error: Check the number of peak search points\n");
    return 1;
  }
  for (i=ne1;i<=ne2;i++) {
    spec2[i]=0.0;
    for (j=0;j<np;j++)
      spec2[i]=spec2[i]+spec[i+ne1-j]*gs[j];
  }
 /* fprintf(stderr,"np=%d, sum=%e, ne1=%d, ne2=%d\n",np,sum,ne1,ne2); */
  /*retun in original array and make adjustment of array sizes */
  *ns=ne2-ne1+1;
  for (i=ne1;i<=ne2;i++) {
    E[i-ne1]=E[i];
    spec[i-ne1]=spec2[i]/sum;
  }
  free(spec2);
  free(gs);
  return 0;
}

/*compute the beam size from twiss parameters
  output: Sx --- sigmaX
          Sy --- sigmaY
          Sxp --- sigmaX prime
          Syp --- sigmaY prime
*/
void ComputeBeamSize(double period, long Nu, double ex, double ey, double Sdelta0, 
		     double betax, double alphax, double etax, double etaxp,
		     double betay, double alphay, double etay, double etayp,
		     double *Sx, double *Sy, double *Sxp, double *Syp)
{
  double gammax, gammay, length;
  
  length = Nu*period;
  
  gammax = (1+sqr(alphax))/betax;
  gammay = (1+sqr(alphay))/betay;
  if (Sxp)
    *Sxp = sqrt(ex*gammax + sqr(Sdelta0*etaxp))*1.0e3;
  if (Syp)
    *Syp = sqrt(ey*gammay + sqr(Sdelta0*etayp))*1.0e3;
  if (Sx)
    *Sx = sqrt(ex*betax + sqr(Sdelta0*etax))*1.0e3;
  if (Sy)
    *Sy = sqrt(ey*betay + sqr(Sdelta0*etay))*1.0e3;
}

/*caculate brightness by calling fortran subroutine usb(), written by Dejus */
/* input: nE number of points for calculation 
          ENERGY         Storage ring energy              (GeV)
          current        Storage ring current             (mA)
          sigmaE	 Energy spread (sigma(E)/E)
          

  output:
          sigmax:        RMS beam size (horizontal)       (mm
          sigmay:        RMS beam size (vertical)         (mm)
          sigmaxp:       RMS beam divergence (horizontal) (mrad)
          sigmayp:       RMS beam divergence (vertical)   (mrad)
         double *K,      array of K values
         double **FnOut  array of F factor fb[H][nE], H is the number of harmonics 
         double **Emergy  array of energy, eb[H][nE] in [kev]
         double **Brightness  array of brightness, sb[H][nE]
         double **LamdarOut array of lamda, LamdarOut[H][nE]
*/  

void Dejus_CalculateBrightness(double current,long nE,
                               double period_mks, long nP, long device,
                               long ihMin,long ihMax,long ihStep,double sigmaE,
                               double gamma, double ex0, double ey0, 
                               double betax, double alphax, double etax, double etaxp,
                               double betay, double alphay, double etay, double etayp,
                               long minNEKS, long maxNEKS, long neks,
			       double *Kvalue, 
                               double kMin,double kMax, long method,
                               double *sigmax,double *sigmay,double *sigmaxp,double *sigmayp,
                               double **K, double ***FnOut,
                               double ***Energy, double ***Brightness, double ***LamdarOut)
{
  double lamdar,reducedE,kx,ky,eMin,eMax,ekMin,ekMax,ep,sp,dep1,dep2,fc,fc2,de,smax;
  int32_t ih,i,j,je,errorFlag=0;
  int32_t nSigma=3,nppSigma=6,nek,ns,exitLoop=0,badPoint=0;
  double JArg,sigmaEE,gk,dek,ek;
  double *tmpE,*tmpSpec,**ei,*ptot,*pd,*kyb,**eb,**sb, eiz;
  double e,k;
  double sigmaX,sigmaX1,sigmaY,sigmaY1,period;
  double ENERGY;
  
#ifdef DEBUG
  fprintf(stderr, "In Dejus_CalculateBrightness\n");
#endif
  
  tmpE=tmpSpec=pd=ptot=NULL;
  ei=eb=sb=NULL;
  period=period_mks*1.0e2; /*use cm as units */
  if (*K) free(*K);
  if (*FnOut) free(*FnOut);
  if (*Energy) free(*Energy);
  if (*Brightness) free(*Brightness);
  if (*LamdarOut) free(*LamdarOut);
  *K=NULL;
  *FnOut=*Energy=*Brightness=*LamdarOut=NULL;
  
  if (neks<=nE) neks=nE+50;

#ifdef DEBUG
  fprintf(stderr, "Allocating memory...\n");
#endif

  tmpE=(double*)calloc(MAXIMUM_E,sizeof(*tmpE));
  tmpSpec=(double*)calloc(MAXIMUM_E,sizeof(*tmpSpec));
  kyb=(double*)calloc(neks+100,sizeof(*kyb));
  ptot=(double*)calloc(neks+100,sizeof(*ptot));
  pd=(double*)calloc(neks+100,sizeof(*pd));
 
  ei=(double**)malloc(sizeof(*ei)*(neks+100));
  eb=(double**)malloc(sizeof(*eb)*(neks+100));
  sb=(double**)malloc(sizeof(*sb)*(neks+100));
  
  for (i=0;i<(neks+100);i++) {
    ei[i]=(double*)calloc(MAXIMUM_H,sizeof(**ei));
    eb[i]=(double*)calloc(MAXIMUM_H,sizeof(**eb));
    sb[i]=(double*)calloc(MAXIMUM_H,sizeof(**sb));
  }
  
#ifdef DEBUG
  fprintf(stderr, "Done\n");
#endif

  /*gamma1=ENERGY/me_mev*1.0E3; */
  lamdar=period*1.0E8/(2.0*gamma*gamma); /*reduced wavelength A */
  ENERGY = gamma*me_mev/1e3; /* GeV */
  reducedE=c_evang/lamdar; /*reduced energy ev */
  kx=0.0;
  if (device==HELICAL) {
    kMin=kMin/sqrt(2.0);
    kMax=kMax/sqrt(2.0);
    if (Kvalue) {
      for (j=0; j<nE; j++)
	Kvalue[j] /=sqrt(2.0);
    }
  }
  eMin=reducedE/(1+kMax*kMax/2.0);
  eMax=reducedE/(1+kMin*kMin/2.0);

  /*determin peak shifts for first and second harmonics at kMin */
  ky=kMin;
  nek=maxNEKS;
  
 
  /*compute electron beam size */
#ifdef DEBUG
  fprintf(stderr, "Computing beam sizes...\n");
#endif
  ComputeBeamSize(period_mks,nP,ex0,ey0,sigmaE,
                  betax,alphax,etax,etaxp,
                  betay, alphay,etay,etayp,&sigmaX, &sigmaY, &sigmaX1, &sigmaY1);
  *sigmax=sigmaX;
  *sigmay=sigmaY;
  *sigmaxp=sigmaX1;
  *sigmayp=sigmaY1;
#ifdef DEBUG
  fprintf(stderr, "Done.\n");
#endif
  dep1 = dep2 = 0;
  for (i=1;i<3;i++) {
    eiz = i * eMax;
    if (i==1) {
      ekMin=0.95*i*eMax;
      ekMax=1.01*i*eMax;
    } else {
      ekMin=0.820*i*eMax;
      ekMax=1.002*i*eMax;
    }
    /*call fortran subroutin usb */
    usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,&nP,&kx,&ky,
        &ekMin,&ekMax,&nek,&method,tmpE,tmpSpec,&ns,&errorFlag);
    if (errorFlag) {
      fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
      exit(1);
    }
    /*find the peak */
    FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
    if (ep> 0.9995*eiz) ep = 0.9995*eiz;
    if (i==1)
      dep1=eMax*i-ep;
    else
      dep2=eMax*i-ep;
  }

#ifdef DEBUG
  fprintf(stderr, "Entering loop for harmonics\n");
#endif
  /*Main loop over harmonics and K-values */
  ih=0;
  de=(eMax-eMin)/(nE-1);
  for (i=ihMin;i<=ihMax;i=i+ihStep) {
    /*Try Kmin initially to find fc for each harmonic  
      Omit Brilliances < SPECLIM (will be set to 0.0) and peak positions will be
      calculated from the zero emittance formula. */
    smax=0;
    je=nE;
    nek=neks;
    badPoint = exitLoop = 0;
    do {
      je--;
      if (!Kvalue)
	ek=eMin+je*de;
      else
	ek = reducedE/(1+Kvalue[je]*Kvalue[je]/2.0);
      ky=sqrt(2.0*(reducedE/ek-1.0));
      if (device==HELICAL) {
        ky=ky/sqrt(2.0);
        kx=ky;
      }
      if (ih==0) {
        kyb[je]=ky;
        ptot[je]=ptot_fac*nP*(kx*kx+ky*ky)*(ENERGY*ENERGY)*current*1.0e-3/(period*1.0e-2);
        if (device==HELICAL)
          gk=HELICAK(ky);
        else
          gk=PLANARK(ky);
        pd[je]=pd_fac*nP*ky*gk*pow(ENERGY,4)*current*1.0e-3/(period*1.0e-2);
      }
      if (i%2) {
        /*odd harmonics */
        ekMin=i*ek-i*dep1;
        ekMax=i*ek+i*dep1/2.0;
        if (i==1) ekMin=ekMin-dep1;
        if (i> (ek/dep1)) {
          fprintf(stderr,"Warning: overlapping range for initial peak search for harmonic %ld\n",i);
          exitLoop=1;
          ih++;
          break;
        }
      } else {
        /*even harmonics */
        ekMin=i*ek-4.0*dep2;
        ekMax=i*ek;
      }
      /*call usb */
      usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,&nP,&kx,&ky,&ekMin,&ekMax,&nek,
          &method,tmpE,tmpSpec,&ns,&errorFlag);
      if (errorFlag) {        
        fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
        exit(1);
      }      
      smax=0;
      for (j=0;j<ns;j++)
        if (tmpSpec[j]>smax) smax=tmpSpec[j];
      ei[je][ih]=i*ek;
      eb[je][ih]=i*ek;
      sb[je][ih]=0.0;
      
    } while (smax<minB && je>0);
    if (exitLoop) {
      badPoint = 1;
    } else if (smax < minB ) {
      badPoint = 1;
    } else {
      FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
      /*define fc */
      fc=0.985*ep/ek/i;
      if (i > 1/(1-fc) ) {
	badPoint = 1;
      }
    }
    je=nE-1;
    for (j=0;j<=je;j++) {
      if (!Kvalue)
	ek=eMin+j*de;
      else
	ek = reducedE/(1+Kvalue[j]*Kvalue[j]/2.0);
      if (!badPoint) {
	ky=sqrt(2.0*(reducedE/ek-1.0));
	if (device==HELICAL) {
	  ky=ky/sqrt(2.0);
	  kx=ky;
	}
	if (ih==0) {
	  kyb[j]=ky;
	  ptot[j]=ptot_fac*nP*(kx*kx+ky*ky)*(ENERGY*ENERGY)*current*1.0e-3/(period*1.0e-2);
	  if (device==HELICAL)
	    gk=HELICAK(ky);
	  else
	    gk=PLANARK(ky);
	  pd[j]=pd_fac*nP*ky*gk*pow(ENERGY,4)*current*1.0e-3/(period*1.0e-2);
	}
      
	if (i%2) 
	  fc2=1.002;
	else
	  fc2=1.000;
	ekMin=fc*i*ek;
	ekMax=fc2*i*ek;
	/* adjust ekmin, ekmax, and number of points if beam energy spread applied */
	if (sigmaE > 0) {
	  nek=neks;
	  dek=(ekMax-ekMin)/nek;
	  sigmaEE=2.0*sigmaE*i*ek; /*  estimated width (eV) */
	  ekMin=ekMin-sigmaEE*nSigma; /* adjust for convolution */
	  ekMax=ekMax+sigmaEE*nSigma;  /* adjust for convolution */
	  if (sigmaEE/nppSigma < dek) dek=sigmaEE/nppSigma;
	  nek=(ekMax-ekMin)/dek+1; /*number of points */
	  if (nek>MAXIMUM_E) {
	    fprintf(stderr,"Energy points out of boundary (constrainted by FORTRAN usb subroutine).\n");
	    exit(1);
	  }
	}
	/* get undulator on-axis brilliance for given K value in energy range ekmin to ekmax
	   returns energy array e (eV) and spec (ph/s/mrad^2/mm^2/0.1%bw) */
	usb_(&ENERGY,&current,&sigmaX,&sigmaY,&sigmaX1,&sigmaY1,&period,
	     &nP,&kx,&ky,&ekMin,&ekMax,&nek,&method,tmpE,tmpSpec,&ns,&errorFlag);
	if (errorFlag) {
	  fprintf(stderr,"error occurred in calling fortran subroutine usb\n");
	  exit(1);
	} 
	FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
	/* fprintf(stderr,"j=%d,points %d,maxE %e,minE %e,peakE %e,peakB %e\n",j,nek,ekMax,ekMin,ep,sp); */
	if (sigmaE>0) {
	  /* gauss convolve */
	  if (Gauss_Convolve(tmpE,tmpSpec,&ns,sigmaE))
	    exit(1);
	}
	FindPeak(tmpE,tmpSpec,&ep,&sp,ns);
	/* fprintf(stderr,"after gauss convolve, peakE %e, peakB %e\n",ep,sp); */
      
	ei[j][ih]=i*ek;
	eb[j][ih]=ep;
	sb[j][ih]=sp;
      } else {
	ei[j][ih]=i*ek;
	eb[j][ih]=0;
	sb[j][ih]=0;
      }
    }
   /* fprintf(stderr,"Harmonics %d completed.\n",i); */
    ih++;
  } /*end for harmonics loop */
#ifdef DEBUG
  fprintf(stderr, "Exited loop for harmonics\n");
#endif
  
  /*output the result */
  *K=(double*)calloc(nE,sizeof(**K));
  *FnOut=(double**)malloc(sizeof(**FnOut)*ih);
  *Energy=(double**)malloc(sizeof(**Energy)*ih);
  *Brightness=(double**)malloc(sizeof(**Brightness)*ih);
  *LamdarOut=(double**)malloc(sizeof(**LamdarOut)*ih);
  for (i=0;i<ih;i++) {
    (*FnOut)[i]=calloc(nE,sizeof(***FnOut));
    (*Energy)[i]=calloc(nE,sizeof(***Energy));
    (*Brightness)[i]=calloc(nE,sizeof(***Brightness));
    (*LamdarOut)[i]=calloc(nE,sizeof(***LamdarOut));
  }
  for (j=0;j<nE;j++) {
    if (device==HELICAL)
      (*K)[j]=kyb[nE-1-j]*sqrt(2);
    else
      (*K)[j]=kyb[nE-1-j];
  }
  
  ih=0;
  for (i=ihMin;i<=ihMax;i=i+ihStep) {
    for (j=0;j<nE;j++) {
      k=(*K)[j];
      JArg = i*k*k/(4+2*k*k);
      (*FnOut)[ih][j]=pow(i*k/(1+k*k/2)*(jn((i+1)/2, JArg) - jn((i-1)/2, JArg)),2);
      (*Energy)[ih][j]=eb[nE-1-j][ih]*1.0e-3;
      (*Brightness)[ih][j]=sb[nE-1-j][ih];
      e=(*Energy)[ih][j];
      (*LamdarOut)[ih][j]=12.39/(e*1.0e10);
    }
    ih++;
  }
#ifdef DEBUG
  fprintf(stderr, "Copied results to output arrays\n");
#endif

  free(tmpE);
  free(tmpSpec);
  free(kyb);
  free(pd);
  free(ptot);
  
  for (i=0;i<(neks+100);i++) {
    free(ei[i]);
    free(eb[i]);
    free(sb[i]);
  }
  free(ei);
  free(eb);
  free(sb);
#ifdef DEBUG
  fprintf(stderr, "Freed memory.\n");
#endif
  return;
}

void getKValueData(double **Kvalue, int32_t *Kpoints, char *filename, char *columnName)
{
  SDDS_DATASET SDDSin;
  
  if (!SDDS_InitializeInput(&SDDSin, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_CheckColumn(&SDDSin, columnName, NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    fprintf(stderr, "Something wrong with column %s in %s\n", columnName, filename);
    exit(1);
  }
  if (SDDS_ReadPage(&SDDSin)<1 || (*Kpoints = SDDS_RowCount(&SDDSin))<1 || 
      !(*Kvalue=SDDS_GetColumnInDoubles(&SDDSin, columnName)))
    SDDS_Bomb("unable to read error factor file");
  /*sort kvalue in decreasing order to make the energy in increasing order */
  qsort(*Kvalue, *Kpoints, sizeof(**Kvalue), double_cmpdes);
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void getErrorFactorData(double **errorFactor, long *errorFactors, char *filename, char *columnName)
{
  SDDS_DATASET SDDSin;
  
  if (!SDDS_InitializeInput(&SDDSin, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_CheckColumn(&SDDSin, columnName, NULL, SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    fprintf(stderr, "Something wrong with column %s in %s\n", columnName, filename);
    exit(1);
  }
  if (SDDS_ReadPage(&SDDSin)<1 || (*errorFactors = SDDS_RowCount(&SDDSin))<1 || 
      !(*errorFactor=SDDS_GetColumnInDoubles(&SDDSin, columnName)))
    SDDS_Bomb("unable to read error factor file");
  if (!SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

