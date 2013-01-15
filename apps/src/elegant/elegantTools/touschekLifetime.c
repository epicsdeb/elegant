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
 * Revision 1.12  2011/10/04 21:18:45  borland
 * Added limit qualifier to -rf option, so that rf specifications will
 * limit the momentum acceptance.  Suggested by exchange with R. Bosch.
 *
 * Revision 1.11  2011/08/17 21:27:41  borland
 * Trim spaces on element names from the aperture, since they may be padded
 * if parallel version was used.
 *
 * Revision 1.10  2010/05/06 20:12:52  xiaoam
 * Fix a round off error.
 *
 * Revision 1.9  2010/02/25 03:56:42  borland
 * Fixed problems with convergence, which would sometimes result in very
 * long run times.
 *
 * Revision 1.8  2010/02/23 15:30:15  borland
 * Emittance from commandline is now treated the same as ex0 from twiss file.
 *
 * Revision 1.7  2010/01/26 03:26:33  borland
 * Fixed bug in -deltaLimit feature (used wrong length for arrays).
 *
 * Revision 1.6  2009/12/20 19:55:50  borland
 * Added -ignoreMismatch option.
 *
 * Revision 1.5  2009/12/15 03:20:14  borland
 * Improved error message.
 *
 * Revision 1.4  2009/03/17 21:16:26  borland
 * Added -deltaLimit option, which allows artificially imposing an maximum
 * momentum aperture.
 *
 * Revision 1.3  2008/05/21 21:00:15  xiaoam
 * Add check radiation integral when read twissFile. It's necessaray when user doesn't supply emitxInput.
 *
 * Revision 1.2  2007/04/12 21:11:52  soliday
 * Updated so that it will compile on WIN32.
 *
 * Revision 1.1  2007/04/12 15:57:47  soliday
 * Added touschekLifetime to elegantTools
 *
 * Revision 1.4  2006/11/07 02:02:06  borland
 * Will now work if radiation integrals are not present in the Twiss file,
 * provided the appropriate data is given on the commandline.
 *
 * Revision 1.3  2006/10/31 19:23:07  borland
 * Latest version from A. Xiao.
 *

 * sdds program to return Touschek lifetime.
 * Using A. Piwinski's formula, DESY 98-179/ISSN 0418-9833.
 * Input files are elegant twiss file with radiation integrals
 * parameters, and a momentum aperture file.
 * Calculated parameters will go in an output sdds file.
*/

#include <stdio.h>
#include "mdb.h"
#include "scan.h"
#include "match_string.h"
#include "SDDS.h"
#include "constants.h"

static char *USAGE = "touschekLifetime <resultsFile>\n\
 -twiss=<twissFile> -aperture=<momentumApertureFile>\n\
 {-charge=<nC>|-particles=<number>} {-coupling=<value>|-emityInput=<meters>}\n\
 [-deltaLimit=<percent>]\n\
 {-RF=Voltage=<MV>,harmonic=<value>,limit | -length=<mm>}\n\
 [-emitInput=<valueInMeters>] [-deltaInput=<value>] [-verbosity=<value>]\n\
 [-ignoreMismatch]\n\
Program by A. Xiao.  (This is version 5, May 2012. M. Borland)";

#define VERBOSE 0
#define CHARGE 1
#define PARTICLES 2
#define COUPLING 3
#define RF 4
#define LENGTH 5
#define EMITINPUT 6
#define DELTAINPUT 7
#define TWISSFILE 8
#define APERFILE 9
#define DELTALIMIT 10
#define IGNORE_MISMATCH 11
#define EMITXINPUT 12
#define EMITYINPUT 13
#define N_OPTIONS 14

char *option[N_OPTIONS] = {
  "verbose",
  "charge",
  "particles",
  "coupling",
  "rf",
  "length",
  "emitinput",
  "deltainput",
  "twiss",
  "aperture",
  "deltalimit",
  "ignoreMismatch",
  "emitxinput",
  "emityinput"
};

void TouschekLifeCalc();  
void FIntegral(double *tm, double *B1, double *B2, double *F, long index, long verbosity); 
double Fvalue (double t, double tm, double b1, double b2);
double linear_interpolation(double *y, double *t, long n, double t0, long i);
void limitMomentumAperture(double *dpp, double *dpm, double limit, long n);
 
/* global varibles */
long plane_orbit = 0;
double *s, *s2, *dpp, *dpm;
double *betax, *alphax, *etax, *etaxp;
double *betay, *alphay, *etay, *etayp;
double *tm, *B1, *B2, *F, *coeff;
double tLife;
long elements, elem2;
long *eOccur1;
double pCentral, sz, sigmap; 
double NP, emitx, emity; 
char **eName1, **eName2, **eType1;
long ignoreMismatch = 0;
#define NDIV 10000;

#define RF_LIMIT_APERTURE 0x01UL

int main( int argc, char **argv)
{
  SCANNED_ARG *scanned;
  char *inputfile1, *inputfile2, *outputfile;
  SDDS_DATASET twissPage, aperPage, resultsPage;
  long verbosity;
  double etaymin, etaymax;
  long i;
  unsigned long rfFlags;
  double emitInput, sigmaDeltaInput, rfVoltage, rfHarmonic, eyInput;
  double alphac, U0, circumference, EMeV;
  double coupling, emitx0, charge;
  short has_ex0 = 0, has_Sdelta0 = 0;
  double deltaLimit = 0;
  
  /****************************************************\
   * read from command line                           *
   \****************************************************/

  SDDS_RegisterProgramName(argv[0]);
  argc  =  scanargs(&scanned, argc, argv);
  if (argc == 1)
    bomb(NULL, USAGE);

  inputfile1  =  NULL;
  inputfile2  =  NULL;  
  outputfile  =  NULL;
  verbosity = 0;
  NP = 0;
  charge = 0;
  coupling = 0;
  sz = 0;
  emitInput = 0;
  sigmaDeltaInput = 0;
  rfVoltage = rfHarmonic = 0;
  rfFlags = 0;
  eyInput = 0;
  
  for (i = 1; i<argc; i++) {
    if (scanned[i].arg_type == OPTION) {
      delete_chars(scanned[i].list[0], "_");
      switch(match_string(scanned[i].list[0], option, N_OPTIONS, UNIQUE_MATCH)) {
      case VERBOSE:
        if (scanned[i].n_items > 1 ) {
          get_long(&verbosity, scanned[i].list[1]);
        } else {
          verbosity=1;
        }
        break;
      case CHARGE:
        if (scanned[i].n_items != 2 )
          bomb("invalid -charge syntax/values", "-charge=<nC>");
        get_double(&charge, scanned[i].list[1]);
        break;
      case EMITINPUT:
      case EMITXINPUT:
        if (scanned[i].n_items != 2 ) 
          bomb("invalid -emitInput syntax/values", "-emitInput=<value>");        
        get_double(&emitInput, scanned[i].list[1]);
        break;
      case DELTAINPUT:
        if (scanned[i].n_items != 2 ) 
          bomb("invalid -deltaInput syntax/values", "-deltaInput=<value>");        
        get_double(&sigmaDeltaInput, scanned[i].list[1]);
        break;
      case LENGTH:
        if (scanned[i].n_items != 2 )
          bomb("invalid -length syntax/values", "-length=<mm>");        
        get_double(&sz, scanned[i].list[1]);
        sz /= 1000;   /* convert input length from mm to m */
        break;
      case COUPLING:
        if (scanned[i].n_items != 2 )
          bomb("invalid -coupling syntax/values", "-coupling=<value>");        
        get_double(&coupling, scanned[i].list[1]);
        break;
      case EMITYINPUT:
        if (scanned[i].n_items != 2 )
          bomb("invalid -eyInput syntax/values", "-eyInput=<meters>");
        get_double(&eyInput, scanned[i].list[1]);
        break;
      case PARTICLES:
        if (scanned[i].n_items != 2 )
          bomb("invalid -particles syntax/values", "-particles=<value>");        
        get_double(&NP, scanned[i].list[1]);
        break;
      case RF:
        if (scanned[i].n_items<2)
          bomb("invalid -rf syntax", NULL);
        scanned[i].n_items--;
        if (!scanItemList(&rfFlags, scanned[i].list+1, &scanned[i].n_items, 0,
                          "voltage", SDDS_DOUBLE, &rfVoltage, 1, 0,
                          "harmonic", SDDS_DOUBLE, &rfHarmonic, 1, 0,
                          "limit", -1, NULL, 0, RF_LIMIT_APERTURE, 
                          NULL) ||
            rfVoltage<=0 || rfHarmonic<=0)
          bomb("invalid -rf syntax/values", "-rf=voltage=MV,harmonic=<value>");
        break;
      case TWISSFILE:
        if (scanned[i].n_items<2)
          bomb("invalid -twiss syntax", NULL);
        inputfile1  =  scanned[i].list[1];
        break;
      case APERFILE:
        if (scanned[i].n_items<2)
          bomb("invalid -twiss syntax", NULL);
        inputfile2  =  scanned[i].list[1];
        break;
      case DELTALIMIT:
        if (scanned[i].n_items != 2 || !get_double(&deltaLimit, scanned[i].list[1]) || deltaLimit<=0)
          bomb("invalid -deltaLimit syntax/values", "-deltaLimit=<percent>");        
        deltaLimit /= 100.0;
        break;
      case IGNORE_MISMATCH:
	ignoreMismatch = 1;
	break;
      default:
        fprintf(stderr, "unknown option \"%s\" given\n", scanned[i].list[0]);
        exit(1);
        break;
      }
    }
    else {
      if (!outputfile) 
        outputfile =  scanned[i].list[0];
      else
        bomb("too many filenames given", NULL);
    }
  }
  if (charge && NP) {
    bomb("Options charge and particles cannot be both specified.",NULL);
  }
  if (!charge) 
    charge = NP * e_mks;
  if (!NP) {
    /* command line input value is in units of nC */
    charge /= 1e9; 
    NP = charge/ e_mks;
  }
  if ((!coupling && !eyInput) || (coupling && eyInput))
    bomb("Give one and only one of coupling or eyInput",NULL);
  if (!sz && !rfVoltage) 
    bomb("Specify either the bunch length or the rf voltage.", NULL);
  
  /****************************************************\
   * Check input twissfile                            *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", inputfile1);
  if (!SDDS_InitializeInput(&twissPage, inputfile1))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  SDDS_ReadPage(&twissPage);

  switch(SDDS_CheckParameter(&twissPage, "I1", NULL, SDDS_DOUBLE, verbosity?stdout:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stdout, "\tParameter I1 not found in input file.\n");
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stdout, "Unexpected result from SDDS_CheckParameter routine while checking parameter Type.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }

  /****************************************************\
   * Check input aperturefile                         *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for checking presence of parameters.\n", inputfile2);
  if (!SDDS_InitializeInput(&aperPage, inputfile2))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* Check presence of momentum aperture */
  SDDS_ReadPage(&aperPage);
  switch(SDDS_CheckColumn(&aperPage, "deltaPositive", NULL, SDDS_DOUBLE, verbosity?stdout:NULL)) {
  case SDDS_CHECK_NONEXISTENT:
    if (verbosity)
      fprintf( stdout, "\tColumn deltaPositive is not found in input file.\n");
    exit(1);
    break;
  case SDDS_CHECK_WRONGTYPE:
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  case SDDS_CHECK_OKAY:
    break;
  default:
    fprintf( stdout, "Unexpected result from SDDS_CheckColumn routine.\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    exit(1);
    break;
  }
  /****************************************************\
   * Check output file and write file layout          *
   \****************************************************/
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", outputfile);
  if (!SDDS_InitializeOutput(&resultsPage, SDDS_BINARY, 1, "Touschek lifetime calculation",
                             "Touschek lifetime calculation", outputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Ignore error returns from these three statements */
  SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "pCentral", NULL);
  SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "Sdelta0", NULL);
  SDDS_TransferParameterDefinition(&resultsPage, &twissPage, "ex0", NULL);
  SDDS_ClearErrors();
  
  if (0>SDDS_DefineParameter(&resultsPage, "coupling", NULL, NULL, 
                             "Coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emitx", "$ge$r$bx$n", "m", 
                             "Horizontal emittance with coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "emity", "$ge$r$by$n", "m", 
                             "Vertical emittance with coupling", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Sdelta", NULL, NULL, 
                             "Fractional momentum spread", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Particles", NULL, NULL, 
                             "Particles", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "Charge", NULL, "nC", 
                             "Charge", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "sigmaz", "$gs$r$bz$n", "m", 
                             "Bunch length", NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&resultsPage, "tLifetime", NULL, "hour", 
                             "Touschek half lifetime", NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  if (!SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "s", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "betax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "alphax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etax", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etaxp", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "betay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "alphay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etay", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "etayp", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementName", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementType", NULL) ||
      !SDDS_TransferColumnDefinition(&resultsPage, &twissPage, "ElementOccurence", NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (0>SDDS_DefineColumn(&resultsPage, "tm", NULL, NULL, 
                          "Local momentum aperture (beta*dp/p)^2", NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&resultsPage, "B1", NULL, NULL, 
                          "Piwinski's parameter B1", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "B2", NULL, NULL, 
                          "Piwinski's parameter B2", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "c0", NULL, NULL, 
                          "1/T=c0*F", NULL, SDDS_DOUBLE, 0) || 
      0>SDDS_DefineColumn(&resultsPage, "F", NULL, NULL, 
                          "Piwinski's parameter F", NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&resultsPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /****************************************************\
   * read from twiss file                             *
   \****************************************************/

  if (!SDDS_GetParameters(&twissPage,
                          "pCentral", &pCentral,
                          NULL) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  elements = SDDS_CountRowsOfInterest(&twissPage);
  s = SDDS_GetColumnInDoubles(&twissPage, "s");
  elem2 = SDDS_CountRowsOfInterest(&aperPage);
  s2 = SDDS_GetColumnInDoubles(&aperPage, "s");
  if(elements<elem2)
    fprintf(stdout, "warning: Twiss file is shorter than Aperture file\n");
  if (emitInput) {
    emitx = emitInput/(1+coupling);
  } else {
    if (!SDDS_GetParameters(&twissPage, "ex0", &emitx0, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    emitx = emitx0/ ( 1 + coupling);
    has_ex0 = 1;
  }
  if (eyInput)
    emity = eyInput;
  else
    emity = emitx * coupling;
  if (sigmaDeltaInput) {
    sigmap = sigmaDeltaInput;
  } else {
    if (!SDDS_GetParameters(&twissPage, "Sdelta0", &sigmap, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    has_Sdelta0 = 1;
  }
  circumference = s[elements-1];
  EMeV = sqrt(sqr(pCentral) + 1) * me_mev;
  if (!sz) {
    /* compute length in m from rf voltage, energy spread, etc */
    if (!SDDS_GetParameters(&twissPage, "alphac", &alphac,
                            "U0", &U0, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    sz = 
      circumference*sigmap*
        sqrt(alphac*EMeV/(PIx2*rfHarmonic*sqrt(sqr(rfVoltage)-sqr(U0))));
    if (rfFlags&RF_LIMIT_APERTURE) {
      double q, rfLimit = 0;
      q = rfVoltage/U0;
      if (q<1 || (rfLimit = sqrt(2*U0/(PI*alphac*rfHarmonic*EMeV)*(sqrt(q*q-1)-acos(1/q))))==0)
        SDDS_Bomb("rf voltage too low compared to energy loss per turn");
      if (deltaLimit==0 || rfLimit<deltaLimit)
        deltaLimit = rfLimit;
    }
  }

  betax = SDDS_GetColumnInDoubles(&twissPage, "betax");
  betay = SDDS_GetColumnInDoubles(&twissPage, "betay");
  alphax = SDDS_GetColumnInDoubles(&twissPage, "alphax");
  alphay = SDDS_GetColumnInDoubles(&twissPage, "alphay");
  etax = SDDS_GetColumnInDoubles(&twissPage, "etax");
  etaxp = SDDS_GetColumnInDoubles(&twissPage, "etaxp");
  etay = SDDS_GetColumnInDoubles(&twissPage, "etay");
  etayp = SDDS_GetColumnInDoubles(&twissPage, "etayp");
  eName1 = SDDS_GetColumn(&twissPage, "ElementName");
  eType1 = SDDS_GetColumn(&twissPage, "ElementType");
  eOccur1 = SDDS_GetColumn(&twissPage, "ElementOccurence");

  eName2 = SDDS_GetColumn(&aperPage, "ElementName");
  for (i=0; i<elem2; i++)
    trim_spaces(eName2[i]);
  
  dpp = SDDS_GetColumnInDoubles(&aperPage, "deltaPositive");
  dpm = SDDS_GetColumnInDoubles(&aperPage, "deltaNegative");
  if (deltaLimit>0)
    limitMomentumAperture(dpp, dpm, deltaLimit, elem2);
  
  /****************************************************\
   * calculate Touschek Lifetime                      *
   \****************************************************/
  find_min_max(&etaymin, &etaymax, etay, elements);
  if ((etaymax-etaymin) < 1e-6)
    plane_orbit=1; 
  
  if (!(tm = SDDS_Malloc(sizeof(*tm)*elements)) ||
      !(B1 = SDDS_Malloc(sizeof(*B1)*elements)) ||
      !(B2 = SDDS_Malloc(sizeof(*B2)*elements)) ||
      !(F = SDDS_Malloc(sizeof(*F)*elements)) ||
      !(coeff = SDDS_Malloc(sizeof(*coeff)*elements)))
    bomb("memory allocation failure (integration arrays)", NULL);
  
  TouschekLifeCalc(verbosity);
  
  /****************************************************\
   * Write output file                                *
   \****************************************************/
  if (0>SDDS_StartPage(&resultsPage, elements) ||
      !SDDS_SetParameters(&resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "pCentral", pCentral, 
                          "Sdelta", sigmap,
                          "coupling", coupling, 
                          "emitx", emitx,
                          "emity", emity,
                          "Particles", NP,
                          "Charge", (1e9 * charge),
                          "tLifetime", tLife,
                          "sigmaz", sz, NULL) ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, s, elements, "s") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, betax, elements, "betax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, alphax, elements, "alphax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etax, elements, "etax") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etaxp, elements, "etaxp") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, betay, elements, "betay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, alphay, elements, "alphay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etay, elements, "etay") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, etayp, elements, "etayp") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eName1, elements, "ElementName") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eType1, elements, "ElementType") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, eOccur1, elements, "ElementOccurence") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, tm, elements, "tm") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, B1, elements, "B1") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, B2, elements, "B2") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, coeff, elements, "c0") ||
      !SDDS_SetColumn(&resultsPage, SDDS_SET_BY_NAME, F, elements, "F"))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (has_ex0 && 
      !SDDS_SetParameters(&resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "ex0", emitx0, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (has_Sdelta0 &&
      !SDDS_SetParameters(&resultsPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Sdelta0", sigmap, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WritePage(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_ReadPage(&twissPage)>0)
    fprintf( stdout, "The code doesn't support multi twiss pages.\n");
  if (SDDS_ReadPage(&aperPage)>0)
    fprintf( stdout, "The code doesn't support multi aperture pages.\n");
  
  if (!SDDS_Terminate(&twissPage) || 
      !SDDS_Terminate(&aperPage) ||
      !SDDS_Terminate(&resultsPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  return(0);
  
}

void TouschekLifeCalc(long verbosity) 
{
  long i, j;
  double sp2, sp4, beta2, betagamma2, sh2;
  double sx2, sxb2, dx2, dx_2;
  double sy2, syb2, dy2, dy_2;
  double a0, c0, c1, c2, c3;
  double pm, pp;
  
  sp2 = sqr(sigmap);
  sp4 = sqr(sp2);
  beta2 = sqr(pCentral)/(sqr(pCentral)+1);
  betagamma2 = sqr(pCentral);
  a0 = sqr(re_mks)*c_mks*NP/(4*sqrt(PI)*(sqr(pCentral)+1)*sz);
  
  i=j=0;
  for (i = 0; i<elements; i++) {
    if (verbosity>1) 
      fprintf(stderr, "Working on i=%ld, j=%ld\n", i, j);
    
/* remove zero length elements from beamline. Except the first element. */
    if(i>0) {
      while(s[i]==s[i-1]) {
        tm[i]=tm[i-1];
        B1[i]=B1[i-1];
        B2[i]=B2[i-1];
        F[i]=F[i-1];
        coeff[i]=coeff[i-1];
        if(++i==elements) break;
      }
    }
    if(i==elements) break;
    
    if(j>0) {
      while(s2[j]==s2[j-1]) {
        dpp[j]=dpp[j-1];
        dpm[j]=dpm[j-1];
        if(++j==elem2) {
          j--;
          break;
        }        
      }
    }
    /* compare two files if it's for same beam line. */       
    if(s[i]>s2[j]) j++;
    if(j==elem2) j--;
    /* Normally the first element for twiss file is _BEG_, which is not true for aperture file.
       But we need it for starting the calculation. */
    if(s[i]!=0 && s[i]==s2[j] && strcmp(eName1[i],eName2[j])!=0) {
      printf("element1 \"%s\" and elem2 \"%s\" at s1 %21.15e s2 %21.15e don't match\n",eName1[i],eName2[j],s[i],s2[j]);
      if (!ignoreMismatch)
	bomb("Twiss and Aperture file are not for same beamline",NULL);
    }

    if (verbosity>1) 
      fprintf(stderr, "Interpolating for i=%ld, j=%ld\n", i, j);
    pp = linear_interpolation(dpp, s2, elem2, s[i], j-1); 
    pm = linear_interpolation(dpm, s2, elem2, s[i], j-1);

    tm[i] = fabs(pp)>fabs(pm)?fabs(pm):fabs(pp);
    tm[i] = beta2*tm[i]*tm[i];

    sxb2 = betax[i]*emitx;
    syb2 = betay[i]*emity;
    dx2 = ipow(etax[i], 2);
    sx2  = sxb2 + dx2*sp2;
    
    dx_2 = ipow(alphax[i]*etax[i]+betax[i]*etaxp[i], 2);
    c1 = sqr(betax[i])/(2*betagamma2*sxb2);
    c2 = sqr(betay[i])/(2*betagamma2*syb2);
    
    if (plane_orbit) {
      sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2);
      B1[i] = c1*(1-sh2*dx_2/sxb2)+c2;
      c0 = sqrt(sh2)/(sigmap*betagamma2*emitx*emity);
      c3 = sx2*syb2;
      B2[i] = sqr(B1[i])-sqr(c0)*c3;   
      if (B2[i]<0) {
        if (fabs(B2[i]/sqr(B1[i]))<1e-7) {
          fprintf(stdout, "warning: B2^2<0 at \"%s\" occurence %ld. Please seek experts help.\n", eName1[i], eOccur1[i]);
        } else {
          B2[i] = 0;
        }
      }
      B2[i]=sqrt(B2[i]);
    }
    
    if (!plane_orbit) {
      dy2 = ipow(etay[i], 2);
      sy2  = syb2 + dy2*sp2;
      dy_2 = ipow(alphay[i]*etay[i]+betay[i]*etayp[i], 2);      
      sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2+(dy2+dy_2)/syb2);
      c0 = sqrt(sh2)/(sigmap*betagamma2*emitx*emity);
      c3 = sx2*sy2-sp4*dx2*dy2;
      B1[i] = c1*(1-sh2*dx_2/sxb2)+c2*(1-sh2*dy_2/syb2);
      B2[i] = sqr(B1[i])-sqr(c0)*c3;   	  
      if (B2[i]<0) {
        if (fabs(B2[i]/sqr(B1[i]))<1e-7) {
          fprintf(stdout, "warning: B2^2<0 at \"%s\" occurence %ld. Please seek experts help.\n", eName1[i], eOccur1[i]);
        } else {
          B2[i] = 0;
        }
      }
      B2[i]=sqrt(B2[i]);   	  
    }
    
    coeff[i] = a0*c0;
    if (verbosity>1) 
      fprintf(stderr, "Computing F integral\n");
    FIntegral(tm, B1, B2, F, i, verbosity);
  }
  if (verbosity>1) 
    fprintf(stderr, "Exited loop\n");
  
  tLife = 0;  
  for (i = 1; i<elements; i++) {
    if (s[i]>s[i-1]) {
      tLife += (s[i]-s[i-1])*(coeff[i]*F[i]+coeff[i-1]*F[i-1])/2;
    }
  }
  if (verbosity>1) 
    fprintf(stderr, "Exited second loop\n");

  tLife /= s[elements-1];
  tLife = 1/tLife/3600;
  return;
} 

void FIntegral(double *tm, double *B1, double *B2, double *F, long index, long verbosity) 
{
  double f0, f1, sum;
  double HPI, step;
  long converge=0, pass, f1NonzeroSeen;
  double t1, k0, k1;
  double tstart, km, b1, b2;
  
  HPI = PI/2;
  step = HPI/NDIV;
  tstart = tm[index];
  b1 = B1[index];
  b2 = B2[index];
  k0 = km = atan(sqrt(tstart));
  f0 = Fvalue(tstart, tstart, b1, b2);
  sum = 0;
  f1 = 0;
  pass = 0;
  f1NonzeroSeen = 0;
  
  while (!converge) {
    if (verbosity>2) 
      fprintf(stderr, "Computing F integral pass %ld: k0=%21.15e, f1=%21.15e\n",
              pass, k0, f1);
    pass++;
    
    k1 = k0 + step;
    t1 = sqr(tan(k1));
    if (isnan(t1)) {
      fprintf(stdout, "Integration failed at Element# %ld", index);
      break;
    }
    
    f1 = Fvalue(t1, tstart, b1, b2);
    if (f1>0)
      f1NonzeroSeen = 1;
    
    if ((f1NonzeroSeen && f1==0) || (f1>0 && f1*(HPI-k1)<1e-3*sum))
      converge = 1;
    
    sum +=(f1+f0)/2*step;
    k0 = k1;
    f0 = f1;    
  }
  F[index] = sum;
  return;
}

double Fvalue (double t, double tm, double b1, double b2)
{
  double c0, c1, c2, result;

  c0 = (sqr(2*t+1)*(t/tm/(1+t)-1)/t+t-sqrt(t*tm*(1+t))-(2+1/(2*t))*log(t/tm/(1+t)))*sqrt(1+t);
  c1 = exp(-b1*t);
  c2 = dbesi0(b2*t);
  result = c0 * c1 * c2;
  /* If overflow/underflow use approximate equation for modified bessel function. */
  if (isnan(result) || result>FLT_MAX) {
    result=c0*exp(b2*t-b1*t)/sqrt(PIx2*b2*t);
  } 
  return result;
}

/* Only linear_interpolate from i=0 to i=n-2. For i<0 using i=0. For i>n-2, using i=n-1 */
double linear_interpolation(double *y, double *t, long n, double t0, long i)
{
  if (i<0)
    i = 0;
  if (i>=n-1)
    i = n-2;
  while (i<=n-2 && t0>t[i+1])
    i++;
  if (i==n-1)
    return (y[n-1]);
  while (i>=0 && t0<t[i])
    i--;
  if (i==-1)
    return (y[0]);
  if (!(t0>=t[i] && t0<=t[i+1])) {
    fprintf(stdout, "failure to bracket point in time array: t0=%e, t[0] = %e, t[n-1] = %e\n",
            t0, t[0], t[n-1]);
    fflush(stdout);
    abort();
  }
 /* this handles zero length elements */ 
  if (t[i+1]==t[i])
    return y[i];
  
  return( y[i] + (y[i+1]-y[i])/(t[i+1]-t[i])*(t0-t[i]) );
}

void limitMomentumAperture(double *dpp, double *dpm, double limit, long n)
{
  long i;
  double limitn;
  
  limitn = -limit;
  for (i=0; i<n; i++) {
    if (dpp[i]>limit)
      dpp[i] = limit;
    if (dpm[i]<limitn)
      dpm[i] = limitn;
  }
}

