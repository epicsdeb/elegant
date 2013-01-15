#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>


#include "SDDS.h"
#include "mdb.h"
#include "scan.h"
#include "xraylib.h"
#include "xraylib-parser.h"

/* Define global variables */

#define MAXNPTS 16384
double *energy = NULL;	  	/* Photon energy of incoming beam, units = eV */	
double total[MAXNPTS];		/* Total cross section of target material for incoming photons, units = cm^2/g */	
double photo[MAXNPTS];		/* Total photoelectric cross section of target material for incoming photons, units = cm^2/g */	
double rayleigh[MAXNPTS];	/* Coherent scattering cross section of target material for incoming photons, units = cm^2/g */
double compton[MAXNPTS];	/* Incoherent scattering cross section of target material for incoming photons, units = cm^2/g */
double RefracIndexRe[MAXNPTS];	/* Real part of refractive index */
double RefracIndexIm[MAXNPTS];	/* Imaginary part of refractive index */
double delta[MAXNPTS];		/* Real part of (1 - n) */
double beta[MAXNPTS];		/* Imaginary part of (1 - n) */
double ScattFactor1[MAXNPTS];	/* Atomic scattering factor f1 */
double ScattFactor2[MAXNPTS];	/* Atomic scattering factor f2 */
double transmission[MAXNPTS];	/* Number of transmitted photons for one incoming photons hitting target, dimensionless */
double absorption[MAXNPTS];	/* Number of "absorbed" photons for one incoming photons hitting target, dimensionless */
double electronYield[MAXNPTS];	/* Number of secondary electrons emerging from the target surface for one incoming photons, units = e/ph */
double eYieldFront[MAXNPTS];	/* Number of secondary electrons emerging from the front target surface for one incoming photons, units = e/ph */
double eYieldBack[MAXNPTS];	/* Number of secondary electrons emerging from the back target surface for one incoming photons, units = e/ph */
double reflectivity[MAXNPTS];	/* Mirror reflectivity for partially polarized x-ray beam */

char   *target=NULL;		/* Common name of target material */
char   *targetFormula=NULL;	/* Chemical formula of target material */
char   shellName[64][8];	/* Designated shell name for Kissle partial cross sections. Example: K or "K,L1,L2,L3" */
long   verbose = 3;		/* Degree of details in screen output (stdout) */
int    mode = 2;		/* default mode of calculation */
int    npts = 10;		/* Number of points of input arrays: energy */
int    nShells, shellID[64];	/* number of subatomic shells and their internal ID */
char   shells[512];		/* List of shelld */
double targetThickness = 0.025;	/* Target thickness in mm */
double targetDensity = 8.9;	/* Target density in g/cm^3 */	
double thetaIn = 0.0;		/* Angle of incidence x-ray beam in degrees. 0 = normal incidence  */	
double polarization = 0.0;	/* Polarization for mirror reflection: 0 = S-polarization, 1 = P-polarization  */	
double edgeEnergy = 9000;	/* Absorption energy, units = eV */
double fluorYield = 0.44;	/* Total fluorescence yield for one excited inner shell core hole, dimensionless */
double jumpFactor = 9.0;	/* Jump factor across the absorption edge, dimensionless */
double electronConfig = 2.0;	/* Occupancy number of the shell to be excited, units = e */
double levelWidth = 1.0;	/* Energy level width, units = eV */
double teyEfficiency = 2.0e-8;	/* Proportional constant for total electron yield in Henke model, defualt for Au is 2.0E-8 */
double pi_const=3.141592653589793238462643;	/* constant */
double degToRad = 0.0174532925199433;	/* Conversion factor from degree to radian */	

static char *matTable="/home/oxygen/OAG/generalData/elementProperties.sdds";

#define modes 10
static int availableMode[modes]={0, 1, 2, 4, 6, 10, 11, 12, 14, 20};

#define SET_ENERGY 0
#define SET_MODE 1
#define SET_TARGET 2
#define SET_SHELL 3
#define SET_VERBOSE 4
#define SET_COPY 5
#define SET_POLARIZATION 6
#define N_OPTIONS 7

static char *option[N_OPTIONS]={"energy", "mode", "target", "shell", "verbose", "copy", "polarization"};

char *USAGE="sddsxra <inputFile> <outputFile> -energy=column=<colname>|begin=<value>,end=<value>,points=<integer> \n\
                     -mode=<number> -target=material=<string>|formula=<string>,thickness=<value>,density=<value>,angle=<value> \n\
                     [-shell=list of shells] [-verbose] [-copy-column|par,<list of names>] \n\
<inputFile>        optional, if provided, the energy values will be obtained from input file. \n\
<outputFile>       required, filename for data output.\n\
energy             units in eV, energy values either provided by a column from input file or specified in the command line:\n\
                   -energy=column=<colname>   Value of energy is obtained from given column in inputfile \n\
                   ‑energy=begin=xx,end=yy,points=npts energy array begin/end values and number of points \n\
mode               mode of calculation: \n\
                   0 -  X-ray cross sections of a target cm^2/g).\n\
                   1 -  Index of refraction of a target\n\
                   2 -  X-ray attenuation by a compound target\n\
                   3 -  (future development) X-ray attenuation by a multi-layer compound target\n\
                   4 -  Total electron yield from a compound film, front, back and total (e/ph)\n\
                   5 -  (future development) Total electron yield from a multi-layer compound film\n\
                   6 -  Reflectivity of a compound mirror \n\
                   7 -  (future development) Reflectivity of a multi-layer compound mirror\n\
                   10 - X-ray cross section of an element (cm2/g) \n\
                   11 - Atomic scattering factor f1 and f2 of an element \n\
                   12 - X-ray attenuation by an elemental target \n\
                   14 - Total electron yield from an elemental film, front, back and total (e/ph) \n\
                   20 - Kissel partial photoelectric cross sections of an element, summed over entire shell (cm2/g)\n\
target             target specifications)\n\
                   material = common name of the target material; or specified by formula = chemical formual of the matrix; \n\
                   thickness = thickness of target (mm)\n\
                   density = density of target (g/cm^3) \n\
                   angle = incident angle of the incoming beam in degrees (0 = normal incidence) \n\
polarization       Polarization for mirror reflection: 0 = S-polarization, 1 = P-polarizatio. \n\
shell              ‑shell= shell1, shell2,… default = “K”\n\
                   Shell names for partial photo-electric cross section data. “K” requests only K-shell partial cross section data (cm2/g)\n\
                   “K,L1,L2,L3” requests partial cross section data of K-, L1-, L2-, and L3-shells (cm2/g).\n\
                   the available shells are: K, L1, L2, L3, M1, M2, M3, M4, N1, N2, N3, N4, N5, N6, N7, O1, O2, O3, O4, O5, O6, O7, P1, P2, P3, P4, P5\n\
                   shell option is only needed for mode=20 calculation.\n\
verbose            printout results to screen.\n\
sddsxra computes x-rays aborsoption spectra using xraylib database and functions.\n\n\
Program by Bingxin Yang and Hairong Shang.  ANL(This is version 1.20, "__DATE__")\n";

/******************
  MAIN PROGRAM

  This program use xraylib functions to calculate the following quantities:
  (0)  X-ray cross sections of a compound using its formula
  (1)  Refraction index of a compound 
  (2)  Attenuation of x-ray beam through a compound target
  (4)  Total electron yield from both surfaces of a compound target
  (6)  Reflectivity of a compound mirror
  (10) X-ray cross sections of an element using its symbol
  (11) Atomic scattering factor of an elemental
  (12) Attenuation of x-ray beam through an elemental target
  (14) Total electron yield from both surfaces of an elemental target
  (20) Kissel partial x-ray cross sections of an element for specified shell

******************/

int getShellID (char *name);
void SetupOutputFile(char *outputfile, SDDS_DATASET  *SDDSout, int mode, SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar);
int getXRayCS (int mode, double massThickness);
int getAtomicXRayCS_Kissel (int shellID);
int getRefraction (int mode);
long get_material_property(char *label, char *property, char **str_value, long *int_value, double *double_value);
int main ( int argc, char *argv[] )
{
  long    i_arg, i, Z, validMode, pages, j, parMatches=0, colMatches=0, exist=0;
  int value;
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile=NULL, *outputfile=NULL, *colName=NULL, **shellList=NULL,  **columnMatch=NULL, **parMatch=NULL;
  SCANNED_ARG *s_arg;
  double eStart, eEnd, deltaE, massThickness;
  unsigned long dummyFlags=0;
  
  verbose = 0;
  eStart = 6000;
  eEnd = 13500;
  npts = 16;
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  targetThickness = 0.01; 	/* units = mm */
  targetDensity   = 2.7;  	/* units = g/cm^3 */ 
  if (argc<2) {
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }
  targetThickness = 0.01; 	/* units = mm */
  targetDensity   = 2.7;  	/* units = g/cm^3 */ 
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_ENERGY:
	s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "column", SDDS_STRING, &colName, 1, 0,
                          "begin", SDDS_DOUBLE, &eStart, 1, 0,
                          "end", SDDS_DOUBLE, &eEnd, 1, 0,
			  "points", SDDS_LONG, &npts, 1, 0, 
                          NULL))
          SDDS_Bomb("invalid -energy syntax");
	if (!colName) {
	  /*column name not provided, should use range*/
	  if (npts<=0)
	    SDDS_Bomb("invalid -energy range syntax: npts <= 0.");
	}
	break;
      case SET_MODE:
	if (s_arg[i_arg].n_items!=2 || !sscanf(s_arg[i_arg].list[1], "%d", &mode))
          SDDS_Bomb("invalid -mode value");
	validMode = 0;
	for (i=0; i<modes; i++) {
	  if (mode==availableMode[i]) {
	    validMode = 1;
	    break;
	  }
	}
	if (!validMode) {
	  fprintf(stderr, "Invalid mode %d provided, must be one of 0, 1, 2, 4, 6, 10, 11, 12, 14, 20\n", mode);
	  exit(1);
	}
        break;
      case SET_TARGET:
	s_arg[i_arg].n_items--;
	targetDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &target, 1, 0,
			  "formula", SDDS_STRING, &targetFormula, 1, 0,
                          "thickness", SDDS_DOUBLE, &targetThickness, 1, 0,
                          "density", SDDS_DOUBLE, &targetDensity, 1, 0,
			  "angle", SDDS_DOUBLE, &thetaIn, 1, 0, 
                          NULL))
          SDDS_Bomb("invalid -target syntax");
	break;
      case SET_POLARIZATION:
	if (s_arg[i_arg].n_items!=2)
	  SDDS_Bomb("Invalid -polarization syntax.");
	if (!get_double(&polarization, s_arg[i_arg].list[1]))
	  SDDS_Bomb("Invalid -polarization value provided.");
	break;
      case SET_SHELL:
	nShells = s_arg[i_arg].n_items-1;
	shellList = malloc(sizeof(*shellList)*(nShells+2));
	shellID[0] = -999;
	SDDS_CopyString(&shellList[0], "All");
	for (i=0; i<nShells; i++) {
	  shellList[i+1] = str_toupper(s_arg[i_arg].list[i+1]);
	  value = getShellID(shellList[i+1]);
	  if (value==-1) {
	    fprintf(stderr, "Invalid shell provided: %s\n", shellList[i+1]);
	    exit(1);
	  }
	  shellID[i+1] = value;
	}
	shellID[nShells+1] = 999;
	SDDS_CopyString(&shellList[nShells+1], "All_Kissel");
	/*add all and all_kissel shells */
	nShells +=2;
	break;
      case SET_VERBOSE:
 	verbose = 1;
 	if (s_arg[i_arg].n_items > 1) {
 	  if ( !get_long(&verbose, s_arg[i_arg].list[1]) )
 	    SDDS_Bomb("invalid -verbose value");
 	}
	break;
      case SET_COPY:
	if (s_arg[i_arg].n_items<3)
	  SDDS_Bomb("Invalid copy syntax provided.");
	if (strncmp_case_insensitive(s_arg[i_arg].list[1], "column", MIN(strlen(s_arg[i_arg].list[1]), 7))==0) {
	  columnMatch = tmalloc(sizeof(*columnMatch)*(colMatches=s_arg[i_arg].n_items-2));
	  for (i=0; i<colMatches; i++)
	    columnMatch[i] = s_arg[i_arg].list[i+2];
	} else if (strncmp_case_insensitive(s_arg[i_arg].list[1], "parameter", MIN(strlen(s_arg[i_arg].list[1]), 9))==0) {
	  parMatch = tmalloc(sizeof(*parMatch)*(parMatches=s_arg[i_arg].n_items-2));
	  for (i=0; i<parMatches; i++)
	    parMatch[i] = s_arg[i_arg].list[i+2];
	} else
	    SDDS_Bomb("Invalid copy syntax.");
	break;
      default:
	fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
	fflush(stdout);
	exit(1);
	break;
      }
    } else {
      if (inputfile==NULL)
        inputfile = s_arg[i_arg].list[0];
      else if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  if (!outputfile) {
    outputfile = inputfile;
    inputfile = NULL;
  }
  if (target)
    if (targetFormula) {
      exist = 1;
    } else {
      exist = get_material_property(target, "formula", &targetFormula,  NULL, NULL);
    }
  else if (targetFormula)
    exist = get_material_property(targetFormula, "name", &target,  NULL, NULL);
  else
    SDDS_Bomb("target not provided.");
  if (!exist) {
    if (target)
      SDDS_CopyString(&targetFormula, target);
    else if (targetFormula)
      SDDS_CopyString(&target, targetFormula);
    if (!targetDensity)
      // targetDensity = 2.7;
      get_material_property(target, "density", NULL, NULL, &targetDensity);
  } else if  (!targetDensity)
    get_material_property(target, "density", NULL, NULL, &targetDensity);
  
  if (mode==20 &&nShells<=2) 
    SDDS_Bomb("No shells provided for mode 20.");
  
  if (!inputfile) {
    /*get eneryg from provided range */
    if (npts>MAXNPTS || npts<1 || eEnd==eStart)
      SDDS_Bomb("Invalid energy range provided.");
    energy = malloc(sizeof(*energy)*npts);
    deltaE = (eEnd-eStart)/(npts-1);
    for (i=0; i<npts; i++)
      energy[i] = eStart + i*deltaE;
    pages = 1; 
  } else {
    if (!colName)
      SDDS_Bomb("Error: column name not provided, can not obtain energy from input file.");
    if (!SDDS_InitializeInput(&SDDSin, inputfile))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_CheckColumn(&SDDSin, colName, "eV", SDDS_ANY_FLOATING_TYPE, stderr)!=SDDS_CHECK_OK) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      fprintf(stderr, "Something wrong with %s column in input file.\n", colName);
      exit(1);
    }
    pages = 1000; /*assume a big number */
  }
  
  /* print minimum information about the calculation */
  if (verbose) {
    fprintf(stdout, "sddsxra: verbose = %ld, mode = %d, target = %s, density = %0.3f (g/cm^3)\n", verbose, mode, target, targetDensity);
    if( mode==2 || mode==4 || mode==12 || mode==14 ) {
      fprintf(stdout, "sddsxra: Target thickness = %0.4f (mm), Incident angle = %0.1f (deg)\n", targetThickness, thetaIn);
    } else if( mode==6 ) {
      fprintf(stdout, "sddsxra: Incident angle = %0.1f (deg)\n", thetaIn);
    } 
  }
  /* Set up ouput files */
  SetupOutputFile(outputfile, &SDDSout, mode, &SDDSin, colMatches, columnMatch, parMatches, parMatch);

  /* Initialize xraylib */
  XRayInit();

  /* Make calculations page-by-page */
  massThickness = 0.1 * targetThickness * targetDensity; 	/* Target mass thickness in g/cm^2 */ 
  for (j=0; j<pages; j++) {
    if (inputfile) {
      if (SDDS_ReadPage(&SDDSin)<1)
	break;
      if ((npts=SDDS_CountRowsOfInterest(&SDDSin))<1)
	continue;
      if (npts>MAXNPTS) {
	fprintf(stderr, "Error execeed maximum allowed points.\n");
	exit(1);
      }
      if (!(energy = SDDS_GetColumnInDoubles(&SDDSin, colName)))
	SDDS_Bomb("Error in reading energy value from input file.");
    }
    if (!SDDS_StartPage(&SDDSout, npts) ||
	!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "TargetMaterial", target, "sddsxra_mode", mode,
			    "TargetFormula", targetFormula,
			    "MassThickness", massThickness, 
			    "TargetThickness",  targetThickness, 
			    "TargetDensity", targetDensity,
			    "ThetaIn", thetaIn, NULL) ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, energy, npts, "PhotonEnergy"))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if ((colMatches || parMatches) && !SDDS_CopyPage(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    /* Retrieve cross sections (mode 0, 10), calculate absorption (mode 2, 12), and TEY (mode 4, 14) */
    if( mode==0 || mode==2 || mode==4 || mode==10 || mode==12 || mode==14) {
      getXRayCS(mode, massThickness);
      if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, total, npts, "TotalCS") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, photo, npts, "PhotoCS") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleigh, npts, "CoherentCS") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, compton, npts, "IncoherentCS"))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (mode==2 || mode==4 || mode==12 ||  mode==14 ) {
	if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, transmission, npts, "Transmission") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, absorption, npts, "Absorption"))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	if (mode==4 || mode==14) {
	  if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, electronYield, npts, "TotalElectronYield") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, eYieldFront, npts, "TotalElectronYieldFront") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, eYieldBack, npts, "TotalElectronYieldBack") ||
	      !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "TargetEfficiency", teyEfficiency, NULL))
	      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      }
    }
    /* Retrieve index of refraction (mode 1) and atomic scattering factor (mode 11), and calculate mirror reflectivity (mode 6) */
    if( mode==1 || mode==6 || mode==11) {
      getRefraction ( mode );
      if (mode==1 || mode==6) {
	  if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, RefracIndexRe, npts, "RefracIndex_Re") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, RefracIndexIm, npts, "RefracIndex_Im") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, delta,         npts, "delta") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, beta,          npts, "beta"))
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  if (mode==6) {
	    if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, reflectivity, npts, "Reflectivity") ||
		!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Polarization", polarization, NULL))
	      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
      } else {
	/* Fetch anomalous scattering factor */
	if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, ScattFactor1 , npts, "F1") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, ScattFactor2, npts, "F2"))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    /* Fetch partical photoelectric cross section for specified shells (Kissel)  */
    if(mode==20) {
      /* calculate photoelectric cross section shell by shell */
      if (verbose)
	fprintf(stdout, "sddsxra: Retrieve photon cross section of %d shells plus total cross sections. \n", nShells - 2);
      Z = SymbolToAtomicNumber ( targetFormula );
      /*first compute all shells */
      for ( i = 0; i < nShells; i++ ) {
	if (verbose) { fprintf(stdout, "shellIndex = %ld, shellID = %d\n", i, shellID[i] ); }
	getAtomicXRayCS_Kissel(shellID[i]);
	  if (!SDDS_StartPage(&SDDSout, npts) ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, total, npts, "TotalCS") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, photo, npts, "PhotoCS") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleigh, npts, "CoherentCS") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, compton, npts, "IncoherentCS") ||
	      !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
				  "ShellName", shellList[i], 
				  "ShellID", shellID[i],
				  "EdgeEnergy",  edgeEnergy,
				  "FluorYield", fluorYield,
				  "JumpFactor", jumpFactor,
				  "LevelWidth", levelWidth,
				  "ElectronConfig", electronConfig, NULL) ||
	      !SDDS_WritePage(&SDDSout))
	    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (inputfile)
	free(energy); energy = NULL;
	/*for mode 20, only the first page of input file is used */
	break;
    } else {
      free(energy); energy = NULL;
      if (!SDDS_WritePage(&SDDSout))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (inputfile && !SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (shellList)
    free(shellList);
  free_scanargs(&s_arg, argc);
  return 0;
}

/************************************
  getShellID (name)
************************************/
int getShellID (char *name)
{
  int value;
  if(strcmp(name, "K") == 0)  { value = K_SHELL;
  } else   if(strcmp(name, "L1") == 0)  { value = L1_SHELL;
  } else   if(strcmp(name, "L2") == 0)  { value = L2_SHELL;
  } else   if(strcmp(name, "L3") == 0)  { value = L3_SHELL;
  } else   if(strcmp(name, "M1") == 0)  { value = M1_SHELL;
  } else   if(strcmp(name, "M2") == 0)  { value = M2_SHELL;
  } else   if(strcmp(name, "M3") == 0)  { value = M3_SHELL;
  } else   if(strcmp(name, "M4") == 0)  { value = M4_SHELL;
  } else   if(strcmp(name, "M5") == 0)  { value = M5_SHELL;
  } else   if(strcmp(name, "N1") == 0)  { value = N1_SHELL;
  } else   if(strcmp(name, "N2") == 0)  { value = N2_SHELL;
  } else   if(strcmp(name, "N3") == 0)  { value = N3_SHELL;
  } else   if(strcmp(name, "N4") == 0)  { value = N4_SHELL;
  } else   if(strcmp(name, "N5") == 0)  { value = N5_SHELL;
  } else   if(strcmp(name, "N6") == 0)  { value = N6_SHELL;
  } else   if(strcmp(name, "N7") == 0)  { value = N7_SHELL;
  } else   if(strcmp(name, "O1") == 0)  { value = O1_SHELL;
  } else   if(strcmp(name, "O2") == 0)  { value = O2_SHELL;
  } else   if(strcmp(name, "O3") == 0)  { value = O3_SHELL;
  } else   if(strcmp(name, "O4") == 0)  { value = O4_SHELL;
  } else   if(strcmp(name, "O5") == 0)  { value = O5_SHELL;
  } else   if(strcmp(name, "O6") == 0)  { value = O6_SHELL;
  } else   if(strcmp(name, "O7") == 0)  { value = O7_SHELL;
  } else   if(strcmp(name, "P1") == 0)  { value = P1_SHELL;
  } else   if(strcmp(name, "P2") == 0)  { value = P2_SHELL;
  } else   if(strcmp(name, "P3") == 0)  { value = P3_SHELL;
  } else   if(strcmp(name, "P4") == 0)  { value = P4_SHELL;
  } else   if(strcmp(name, "P5") == 0)  { value = P5_SHELL;
  } else { value = -1;
  }
  return value;
}



/**********************************************************
  getRefraction ()
  
  This function retrieve the x-ray index of refraction and
  Atomic scattering factors f1 and f2
***********************************************************/
int getRefraction (int mode)
{
  int     i, Z=0;
  float   energy_keV;
  complex cosThetaIn, sinThetaIn, cosThetaOut, n, sAmpReflect, pAmpReflect;
  
  if (verbose > 3) { fprintf(stdout, "getRefraction: targetFormula = %s\n", targetFormula); }
  if (verbose > 2) { 
    if ( mode==1 ) { 
      fprintf(stdout, "\n Index  PhotonEnergy    delta        beta"); 
    } else if ( mode==6 ) { 
      fprintf(stdout, "\n Index  PhotonEnergy    delta        beta      Refelectivity"); 
    } else { fprintf(stdout, "\n Index  PhotonEnergy       f1          f2"); 
    }
    fprintf(stdout, "\n            (eV)       ");
  }

  cosThetaIn = cos(thetaIn  * degToRad);
  sinThetaIn = sin(thetaIn  * degToRad);
  if (mode > 9 )   { Z = SymbolToAtomicNumber ( targetFormula ); }
  for ( i = 0; i < npts; i++ ) {
    energy_keV  = 0.001 * energy[i];
    if ( mode < 10 ) {  /* Compound target */
      RefracIndexRe[i] = Refractive_Index_Re ( targetFormula, energy_keV, targetDensity );
      RefracIndexIm[i] = Refractive_Index_Im ( targetFormula, energy_keV, targetDensity );
      delta[i]         = 1.0 - RefracIndexRe[i];
      beta[i]          = RefracIndexIm[i];
      if (verbose>2) { fprintf(stdout, "\n %4d %12.1f %12.4g %12.4g", i, energy[i], delta[i], beta[i]); }
      if ( mode==6 ) {  /* Mirror reflectivity */
        n = RefracIndexRe[i] + RefracIndexIm[i] * _Complex_I;
        cosThetaOut     = csqrt(n * n - sinThetaIn * sinThetaIn) / n;
        sAmpReflect     = ( cosThetaIn - n * cosThetaOut ) / (cosThetaIn + n * cosThetaOut );
        pAmpReflect     = ( cosThetaOut - n * cosThetaIn ) / (cosThetaOut + n * cosThetaIn );
        reflectivity[i] = 0.5 * (1.0 + polarization) * cabs(pAmpReflect * pAmpReflect) + 0.5 * (1.0 - polarization) * cabs(sAmpReflect * sAmpReflect);
        // fprintf(stdout, "\n n = %f,%f, determ = %f, %f", n, determ); 
        // fprintf(stdout, "\n sAmpReflect = %f, %f, pAmpReflect = %f, %f, reflectivity = %f", sAmpReflect, pAmpReflect, reflectivity[i]); 
        if (verbose>2) { fprintf(stdout, "%12.4g", reflectivity[i]); }
      }
    } else { /* Atomic scattering factor */
      ScattFactor1[i] = Z + Fi( Z, energy_keV );
      ScattFactor2[i] = - Fii(Z, energy_keV );
      if (verbose>2) { fprintf(stdout, "\n %4d %12.1f %12.4g %12.4g", i, energy[i], ScattFactor1[i], ScattFactor2[i]); }
    }
  }

  fprintf(stdout, "\n");
  return 0;
}


/**********************************************************
  getXRayCS ()
  
  This function retrieve the x-ray interaction constants 
  of an element using xraylib program/database package 
***********************************************************/

int getXRayCS (int mode, double massThickness)
{
  int     i, Z=0;
  float   energy_keV;
  double  cosThetaIn;
  
  if (verbose > 3) { fprintf(stdout, "getXRayCS: targetFormula = %s\n", targetFormula); }
  if (mode > 9 )  { Z = SymbolToAtomicNumber ( targetFormula ); }
  if (verbose > 2) { 
    fprintf(stdout, "\n Index  PhotonEnergy   TotalCS    PhotoCS    coherentCS  incohrentCS");
    if ( mode==2 || mode==12 ) { fprintf(stdout, "    Transmission  Absorption "); 
    } else if ( mode==4 || mode==14 ) { fprintf(stdout, "  Transmission  Absorption  FrontTEY(QE)  BackTEY(QE)  TEY(QE) "); 
    }
    fprintf(stdout, "\n            (eV)       (cm2/g)    (cm2/g)      (cm2/g)     (cm2/g)  ");
  }

  cosThetaIn  = cos(thetaIn  * degToRad);
  for ( i = 0; i < npts; i++ ) {
    energy_keV  = 0.001 * energy[i];
    if ( mode < 10 ) {
      total[i]    = CS_Total_CP ( targetFormula, energy_keV );
      photo[i]    = CS_Photo_CP ( targetFormula, energy_keV );
      rayleigh[i] = CS_Rayl_CP  ( targetFormula, energy_keV );
      compton[i]  = CS_Compt_CP ( targetFormula, energy_keV );
    } else {
      total[i]    = CS_Total ( Z, energy_keV );
      photo[i]    = CS_Photo ( Z, energy_keV );
      rayleigh[i] = CS_Rayl  ( Z, energy_keV );
      compton[i]  = CS_Compt ( Z, energy_keV );
    }
    if (verbose>2) { fprintf(stdout, "\n %4d %10.1f %12.4g %10.4g %12.4g %12.4g", i, energy[i], total[i], photo[i], rayleigh[i], compton[i]); }
    
    if ( mode==2 || mode==4 || mode==12 || mode==14 ) {
      transmission [i] = exp(-1.0 * massThickness * total[i] / cosThetaIn ); 
      absorption[i]    = 1.0 - transmission[i];
      if (verbose>2) { fprintf(stdout, " %12.4f %12.4g", transmission[i], absorption[i]); }
      if ( mode==4 || mode==14 ) {
        eYieldFront[i]   = teyEfficiency * energy[i] * total[i] / cosThetaIn ; 	/* Front surface total electron yield */
        eYieldBack[i]    = teyEfficiency * energy[i] * total[i] / cosThetaIn * transmission [i]; /* Back surface total electron yield */
        electronYield[i] = eYieldFront[i]  + eYieldBack[i];
        if (verbose>2) { fprintf(stdout, "%12.4g %12.4g %12.4g", eYieldFront[i], eYieldBack[i], electronYield[i] ); }
      }
    }
  }
  if (verbose)
  fprintf(stdout, "\n");
  return 0;
}


/**********************************************************
  getAtomicXRayCS_Kissel ()
  
  This function retrieve the x-ray interaction constants
  of an element using xraylib program/database package,
  including partial cross sections for specified shells 
***********************************************************/

int getAtomicXRayCS_Kissel (int shellID)
{
  int   i, Z;
  float energy_keV;
  
  Z = SymbolToAtomicNumber ( targetFormula );
  edgeEnergy = 0.0;
  fluorYield = 1.0;
  jumpFactor = 1.0;
  levelWidth = 0.0;
  electronConfig  = Z;

  if (verbose > 3) { 
    fprintf(stdout, "getAtomicXRayCS: Z = %d\n", Z);
    fprintf(stdout, "getAtomicXRayCS_Kissel: shellID = %d\n", shellID);
    if (verbose>2) { fprintf(stdout, "Index  PhotonEnergy   TotalCS      PhotoCS     coherentCS   incohrentCS \n"); }
  }

  if (shellID <= 30 && shellID >= 0) {
    edgeEnergy = 1000.0 * EdgeEnergy(Z, shellID);
    fluorYield = FluorYield(Z, shellID);
    jumpFactor = JumpFactor(Z, shellID);
    electronConfig  = ElectronConfig(Z, shellID);
    levelWidth = AtomicLevelWidth(Z, shellID);
    for ( i = 0; i < npts; i++ ) {
      energy_keV  = 0.001 * energy[i];
      if ( energy[i] > edgeEnergy ) {
        photo[i]  = CS_Photo_Partial (Z, shellID, energy_keV);
      } else {
        photo[i]  = 0.0;
      }
      total[i]    = 0.0;
      rayleigh[i] = 0.0;
      compton[i]  = 0.0;
      if (verbose>2) { fprintf(stdout, "%4d %12.1f %12.4g %12.4g %12.4g %12.4g \n", i, energy[i], total[i], photo[i], rayleigh[i], compton[i]); }
    }
  } else if (shellID > 99) {
    for ( i = 0; i < npts; i++ ) {
      energy_keV  = 0.001 * energy[i];
      total[i]    = CS_Total_Kissel ( Z, energy_keV );
      photo[i]    = CS_Photo_Total  ( Z, energy_keV );
      rayleigh[i] = CS_Rayl  ( Z, energy_keV );
      compton[i]  = CS_Compt ( Z, energy_keV );
      if (verbose>2) { fprintf(stdout, "%4d %12.1f %12.4g %12.4g %12.4g %12.4g \n", i, energy[i], total[i], photo[i], rayleigh[i], compton[i]); }
    }
  } else {
    for ( i = 0; i < npts; i++ ) {
      energy_keV  = 0.001 * energy[i];
      photo[i]    = CS_Photo ( Z, energy_keV );
      total[i]    = CS_Total ( Z, energy_keV );
      rayleigh[i] = CS_Rayl  ( Z, energy_keV );
      compton[i]  = CS_Compt ( Z, energy_keV );
      if (verbose>2) { fprintf(stdout, "%4d %12.1f %12.4g %12.4g %12.4g %12.4g \n", i, energy[i], total[i], photo[i], rayleigh[i], compton[i]); }
    }
  }

  if (verbose > 1) { 
    fprintf(stdout, "edgeEnergy = %f, fluorYield = %f,  jumpFactor = %f,", edgeEnergy, fluorYield, jumpFactor);
    fprintf(stdout, " electronConfig = %f, levelWidth = %f \n", electronConfig, levelWidth);
  }
  return 0;
}

void SetupOutputFile(char *outputfile, SDDS_DATASET *SDDSout, int mode, SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar) {
  long i, cols=0, pars=0;
  char **column=NULL, **par=NULL;
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (copyCols) {
    column = getMatchingSDDSNames(SDDSin, copyCol, copyCols, &cols, SDDS_MATCH_COLUMN);
    SDDS_SetColumnFlags(SDDSin, 0);
    for (i=0; i<cols; i++) {
      if (!SDDS_TransferColumnDefinition(SDDSout, SDDSin, column[i], NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetColumnsOfInterest(SDDSin, SDDS_MATCH_STRING, column[i], SDDS_OR))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_FreeStringArray(column, cols);
  }
  if (copyPars) {
    par = getMatchingSDDSNames(SDDSin, copyPar, copyPars, &pars, SDDS_MATCH_PARAMETER);
    for (i=0; i<pars; i++) {
      if (!SDDS_TransferParameterDefinition(SDDSout, SDDSin, par[i], NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_FreeStringArray(par, pars);
  }
  if (verbose > 1 ) { fprintf(stdout, "cols %ld pars %ld\n", cols, pars); }
  if (!SDDS_DefineSimpleParameter(SDDSout, "sddsxra_mode", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetMaterial", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetFormula", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "MassThickness", "g/cm^2", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetThickness", "mm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetDensity", "g/cm^3", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "ThetaIn", "degrees", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "PhotonEnergy", "eV", SDDS_DOUBLE))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (mode!=1 && mode!=6 && mode!=11) {
    if (SDDS_DefineColumn(SDDSout, "TotalCS", NULL, "cm^2/g", "Total x-ray cross section",NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "PhotoCS", NULL, "cm^2/g", "Photoelectric cross section",NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "CoherentCS", NULL, "cm^2/g", "Coherent scattering cross section ",NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "IncoherentCS", NULL, "cm^2/g", "Incoherent scattering cross section ",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==1 || mode==6) {
    if (SDDS_DefineColumn(SDDSout, "RefracIndex_Re", NULL, NULL, "Real part of refractive index", NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "RefracIndex_Im", NULL, NULL, "Imaginary part of refractive index", NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "delta", NULL, NULL, "Real part of 1 - n", NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "beta", NULL, NULL, "Imaginary part 1 - n", NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else if (mode==11) {
    if (SDDS_DefineColumn(SDDSout, "F1", NULL, NULL, "Atomic scattering factor f1", NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "F2", NULL, NULL, "Atomic scattering factor f2", NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==6) {
    if (SDDS_DefineColumn(SDDSout, "Reflectivity", NULL, "g/cm^3", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
	!SDDS_DefineSimpleParameter(SDDSout, "Polarization", NULL, SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==2 || mode==4 || mode==12 || mode==14) {
    if (SDDS_DefineColumn(SDDSout, "Transmission", NULL, NULL, "Ratio of x-ray beam transmitted through the film target", NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "Absorption", NULL, NULL, "Ratio of x-ray beam absorbed by the film target", NULL, SDDS_DOUBLE, 0)<0 )
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==4 ||  mode==14 ) {
    if (SDDS_DefineColumn(SDDSout, "TotalElectronYieldFront", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "TotalElectronYieldBack", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
	SDDS_DefineColumn(SDDSout, "TotalElectronYield", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
	!SDDS_DefineSimpleParameter(SDDSout, "TargetEfficiency", "g/cm^2", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==20) {
    if (!SDDS_DefineSimpleParameter(SDDSout, "ShellID", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ShellName", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "EdgeEnergy", "eV", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "FluorYield", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "JumpFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "LevelWidth", "eV", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ElectronConfig", "e", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

#define Z_COL 0
#define FORMULA_COL 1
#define NAME_COL 2
#define GROUP_COL 3
#define PERIOD_COL 4
#define MASS_COL 5
#define DENSITY_COL 6
#define MELT_COL 7
#define BOIL_COL 8
#define HEAT_COL 9
#define NEG_COL 10
#define ABUND_COL 11
#define TABLE_COLS 12
char *table_column[TABLE_COLS]={"Z", "Formula", "Name", "Group", "Period", "Mass", "Density", "MeltingPoint", "BoilingPoint", "Heat", 
				"Electronegativity", "Abundance"};

long get_material_property(char *label, char *property, char **str_value, long *int_value, double *double_value) {
  long rows=0, index, col_index=-1;
  SDDS_DATASET table;
  int32_t *intValue = NULL;
  char **formula, **name;
  double *doubleValue = NULL;
  
  formula = name = NULL;
  if ((col_index=match_string(property, table_column, TABLE_COLS, 0))<0) {
    if (verbose)
      fprintf(stderr, "Property - %s does not exist in the property table.\n", property);
    return -1;
  }
  if (!SDDS_InitializeInput(&table, matTable))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_ReadPage(&table) || !(rows=SDDS_CountRowsOfInterest(&table)) ||
      !(formula=(char**)SDDS_GetColumn(&table, "Formula")) ||
      !(name=(char**)SDDS_GetColumn(&table, "Name")))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  index = match_string(label, formula, rows, EXACT_MATCH);
  if (index<0)
    index = match_string(label, name, rows, 0);
  if (index<0) {
    fprintf(stdout, "%s not found.\n", label);
    return -1;
  }
  switch (col_index) {
  case FORMULA_COL:
    SDDS_CopyString(str_value, formula[index]);
    break;
  case NAME_COL:
    SDDS_CopyString(str_value, name[index]);
    break;
  case Z_COL:
  case GROUP_COL:
  case PERIOD_COL:
    intValue = (int32_t*)SDDS_GetColumn(&table, table_column[col_index]);
    *int_value = intValue[index];
    free(intValue);
    break;
  default:
    doubleValue = (double*)SDDS_GetColumn(&table, table_column[col_index]);
    *double_value = doubleValue[index];
    free(doubleValue);
    break;
  }
  if (!SDDS_Terminate(&table))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  SDDS_FreeStringArray(formula, rows);
  SDDS_FreeStringArray(name, rows);
  free(formula);
  free(name);
  return index;
}

