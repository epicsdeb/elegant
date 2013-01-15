#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>

#include "SDDS.h"
#include "mdb.h"
#include "scan.h"
#include "xraylib.h"
#include "xraylib-parser.h"

/* Define global variables */
#define MAXNPTS 16384
double *energy=NULL;	  		/* Photon energy of incoming beam, units = eV */	
double matrixTotalCS[MAXNPTS];		/* Total cross section of matrix material for incoming photon, units = cm^2/g */	
double traceXRFCS[MAXNPTS];		/* XRF cross section of trace element for incoming photon, sum of all specified lines, units = cm^2/g */
double totalXRFYield[MAXNPTS];		/* Total XRF photon number for one incoming photon, sum of all specified lines, dimensionless. */
double xrfFluxDensity[MAXNPTS];		/* XRF photon flux density  number for one incoming photon, sum of all specified lines, units = 1/sterad */
double xrfFluxFactor[MAXNPTS];		/* Number of XRF photons through the detector window for one incoming photons, sum of all specified lines. */
double xrfPowerFactor[MAXNPTS];		/* Ratio of XRF photon power absorbed by the detector over incoming photon power, sum of all specified lines. */

char   matrix[512]="Diamond";		/* Chemical name of target (matrix) material */
char   frontFilter[512]="Diamond";	/* Common name of the material for filter upstream of the target */
char   backFilter[512]="Diamond";	/* Common name of the material for filter downstream of the target */
char   matrixFormula[512]="C";		/* Chemical formula of target material */
char   frontFilterFormula[512]="C";	/* Chemical formula for upstream filter */
char   backFilterFormula[512]="C";	/* Chemical formula for downstream filter */
char   trace[512]="Cu";			/* Chemical formula of trace element in the target */
char   detector[512]="Silicon";		/* Chemical formula of detector material */
char   detectorFormula[512]="Si";	/* Chemical formula of detector material */
char   xrfLines[512];			/* List of xrf lines */
int    nLines, xrfLineID[64];		/* number of XRF lines and their internal ID */
int    npts = 10;			/* Number of points of input arrays: energy */
long   verbose = 0;			/* Degree of details in screen output (stdout) */
long   mode = 14;			/* default mode of calculation */
long   nxDet = 20, nyDet = 20;		/* number of divisions of the detector surface for numerical integration */
long   detectDownstream;		/* Position of the detector: 0 = in front of target (reflective geometry), 1 = behind the target (transmission geometry) */
double targetThickness = 0.025;		/* Target thickness in mm */
double frontFilterThickness = 0.0;	/* Upstream filter thickness in mm */
double backFilterThickness = 0.0;	/* Downstream filter thickness in mm */
double targetDensity = 8.9;		/* Target (matrix) density in g/cm^3 */	
double frontFilterDensity = 3.51;	/* Upstream filter density in g/cm^3 */	
double backFilterDensity = 3.51;	/* Downstream filter density in g/cm^3 */
double traceConcentration = 1.0;	/* Trace element (XRF atom) concentration by weight in target */
double xrfLineEnergy[32];		/* XRF emission line energy in eV  */
double thetaIn = 0.0;			/* Angle of incidence x-ray beam in degrees. 0 = normal incidence  */	
double thetaOut = 0.0;  		/* Angle of outgoing XRF photons in degrees. 0 = normal exit */	
double detectThickness = 0.30;		/* Detector thickness in mm */
double detectDensity = 2.7;		/* Detector density in g/cm^3 */	
double detectThetaX = 90.0, detectThetaY = 90.0, detectThetaZ = 0.0; 	/* Angle between detector surface and x-, y-, z-axis.  */
double xDet = 0.0, yDet = 0.0; 		/* Transverse position of the Detector center in mm, Oxy plane is parallel to the target front surface */	
double zDet = 100;			/* Longitudinal position of the Detector center in mm, measured perpendicular to the target front surface */
double xWidth = 10, yWidth = 10;	/* Detector sizes in mm, projected to Oxy plane */
double pi_const=3.141592653589793238462643;	/* constant */
double degToRad = 0.0174532925199433;	/* Conversion factor from degree to radian */	

/****************************************************************************************
  MAIN PROGRAM

  This program use xraylib functions to calculate the following data:
  (0)  XRF cross sections in a compound matrix
  (1)  Total XRF yield in a compound target, no self absorption
  (2)  XRF flux density from a uniform compound target
  (3)  XRF detector signal from a uniform compound target
  (4)  XRF detector signal from a compound target with front and back filters
  (10) XRF cross sections of an element matrix
  (11) Total XRF yield in an elemental target, no self absorption
  (12) XRF flux density from a uniform elemental target
  (13) XRF detector signal from an elemental compound target
  (14) XRF detector signal from an elemental target with front and back filters
  
****************************************************************************************/

#define modes 11
static int availableMode[modes]={0, 1, 2, 3, 4, 10, 11, 12, 13, 14};

#define SET_ENERGY 0
#define SET_MODE 1
#define SET_TARGET 2
#define SET_VERBOSE 3
#define SET_TRACE 4
#define SET_DETECTOR 5
#define SET_FRONTFILTER 6
#define SET_BACKFILTER 7
#define SET_COPY 8
#define N_OPTIONS 9

static char *option[N_OPTIONS]={"energy", "mode", "target", "verbose", "trace", "detector", "frontfilter","backfilter", "copy"};

char *USAGE="sddsxrf <inputFile> <outputFile> -energy=column=<colname>|begin=<value>,end=<value>,points=<integer> -mode=<number> \n\
                     -trace=element=<string>,concentration=<valu>,lines=<line1,line2,line3,..> [-verbose] \n\
                     -target=material=<string>,formula=<string>,thickness=<value>,density=<value>,angle=<value> \n\
                     -detector=material=<string>,formula=<string>,thickness=<value>,density=<value>,\n\
                             xcenter=<value>,ycenter=<value>,zcenter=<value>,xwidth=<value>,ywidth=<value>, \n\
                             xtheta=<value>,ytheta=<value>,ztheta=<value>,nx=<integer>,ny=<integer>,downstream=<integer> \n\
                     -trace=element=<string>,concentration=<valu>,lines=<line1,line2,line3,..> [-verbose] \n\
                     -frontfilter=material=<string>,formula=<string>,thickness=<value>,density=<value> \n\
                     -backfilter=material=<string>,formula=<string>,thickness=<value>,density=<value> \n\
                     [-copy-column|par,<list of names>] \n\
<inputFile>        optional, if provided, the energy values will be obtained from input file. \n\
<outputFile>       required, filename for data output.\n\
energy             unites in eV, energy values either provided by input file from given column (column) or provided by a range (specified).\n\
                   -energy=column=<colname>   Value of energy is obtained from given column in inputfile \n\
                   ‑energy=begin=xx,end=yy,points=npts energy array begin/end values and number of points \n\
mode               mode of XRF (x-ray fluorescence) calculation:\n\
                   0  - XRF fluorescence cross sections of a compound matrix (cm2/g). \n\
                   1  - Total XRF yield from a compound matrix, dimensionless, no self absorption. \n\
                   2  - XRF flux density from the front or back surface of a compound target (1/sterad). \n\
                   3  - XRF detector signal from a uniform compound target\n\
                   4  - (experimental) XRF detector signal from a compound target with front and back filters \n\
                   10 - XRF fluorescence cross sections of an elemental matrix (cm2/g). \n\
                   11 - Total XRF yield from an elemental target, no self absorption. \n\
                   12 - XRF flux density from a uniform elemental target\n\
                   13 - XRF detector signal from an elemental compound target \n\
                   14 - (experimental) XRF detector signal from an elemental target with front and back filters \n\
trace              trace element specifications: \n\
                   element = chemical symbol of trace element; \n\
                   concentration = concentration by weight. \n\
                   lines = list of XRF lines, always the last part in –trace option \n\
                   KL1 or KA requests only K-line cross section  (cm2/g) \n\
                   KL1,KL2,KL2 requests Kalpha, Kbeta1 and Kbeta2 lines data (cm2/g) \n\
target             target/matrix specifications: \n\
                   material = common name of the matrix material, such as Water; \n\
                   formula = chemical formula of the matrix, such as H2O; \n\
                   thickness = thickness of target (mm) \n\
                   density = density of target (g/cm^3)  \n\
                   angle = incident angle of the incoming beam in degrees (0 = normal incidence) \n\
detctor            define detectors: \n\
                   material = common name of the detector material, such as Water; \n\
                   formula = chemical formula of the detctor; \n\
                   thickness = thickness of detector (mm) \n\
                   density = density of detector (g/cm^3) \n\
                   xcenter, ycenter, zcenter = coordinates of the detector center. \n\
                   xwidth, ywidth = width of the detector in dimensions approximately parallel to x and y axes \n\
                   nx, ny are the number of divisions of the detector surface for numerical integration.\n\
                   xtheta, ztheta = projection angle of the detector surface normal vector to coordinate planes.\n\
                   downstream = target face seen by the detector, 1 for downstream surface and 0 for upstream surface.\n\
                   center specification of detector is required for mode 2 and 12. \n\
                   full specification of detector is required for mode 3, 4, 13 and 14. \n\
frontfilter        specifications of filters upstream and downstream of the target/matrix (mode 4 and 14): \n\
backfilter         material = common name of the filter material, such as Water; \n\
                   formula = chemical formula of the filter, such as H2O; \n\
                   thickness = thickness of filter (mm) \n\
                   density = density of filter (g/cm^3)  \n\
copy               copy parameters or columns from input file.\n\
sddsxrf computes xray fluorescense yield. The bulk material receiving the beam is called the matrix or target. \n\
The element emitting fluorescence x-ray photons is called the “trace,” its concentration may be as low as 1 ppm or \n\
as high as 100%. The concentration of the trace element is given by its weight fraction.\n\n\
Program by Bingxin Yang and Hairong Shang.  ANL(This is version 1.20, "__DATE__")\n";

int getLineID (char *name);
int getXRFCS (int mode);
void SetupOutputFile(char *outputfile, SDDS_DATASET  *SDDSout, int mode, SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar);
long get_material_property(char *label, char *property, char **str_value, long *int_value, double *double_value);

int main ( int argc, char *argv[] )
{
  char  *material, *formula, *str_value=NULL, **columnMatch=NULL, **parMatch=NULL;
  long  i_arg, i, validMode, pages, len=0, upstream=0, Z, colMatches=0, parMatches=0;
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile=NULL, *outputfile=NULL, *colName=NULL, *traceLineList[126], *valuePtr=NULL, *prevKeyword=NULL;
  SCANNED_ARG *s_arg;
  double eStart, eEnd, deltaE, defaultDensity;
  unsigned long dummyFlags=0;
  
  eStart = eEnd = 0;
  npts = 0;
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  
  if (argc<2) {
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }
  
  /* Process command line options */
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
	  if (eStart==eEnd  || npts==0)
	    SDDS_Bomb("invalid -energy range syntax");
	}
	break;
      case SET_MODE:
	if (s_arg[i_arg].n_items!=2 ||
            !get_long(&mode, s_arg[i_arg].list[1]))
          SDDS_Bomb("invalid -mode value");
	validMode = 0;
	for (i=0; i<modes; i++) {
	  if (mode==availableMode[i]) {
	    validMode = 1;
	    break;
	  }
	}
	if (!validMode) {
	  fprintf(stderr, "Invalid mode %ld provided, must be one of 0, 1, 2, 3, 4, 10, 11, 12, 13, 14\n", mode);
	  exit(1);
	}
        break;
      case SET_FRONTFILTER:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	frontFilterDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &frontFilterThickness, 1, 0,
                          "density", SDDS_DOUBLE, &frontFilterDensity, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -target syntax");
        if (material!=NULL) { 
  	  strcpy(frontFilter, material);
  	  if (formula!=NULL) {
    	    strcpy(frontFilterFormula, formula);
  	  } else {
	    if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	      strcpy(frontFilterFormula, material);
	    else 
	      strcpy(frontFilterFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
  	  strcpy(frontFilterFormula, formula);
	  if (get_material_property(formula, "name", &str_value, NULL, NULL)<0) 
	    strcpy(frontFilter, formula);
	  else
	    strcpy(frontFilter, str_value);
  	} else {
          SDDS_Bomb("Front filter material is not specified.");
  	}
	if (frontFilterDensity<=0) {
	  if (get_material_property(frontFilter, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: front filter density can not be obtained from database, has to be provided.");
	  frontFilterDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value = NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;
      case SET_BACKFILTER:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	backFilterDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &backFilterThickness, 1, 0,
                          "density", SDDS_DOUBLE, &backFilterDensity, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -target syntax");
        if (material!=NULL) { 
  	  strcpy(backFilter, material);
  	  if (formula!=NULL) {
    	    strcpy(backFilterFormula, formula);
  	  } else {
	    if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	      strcpy(backFilterFormula, material);
	    else 
	      strcpy(backFilterFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
  	  strcpy(backFilterFormula, formula);
	  if (get_material_property(formula, "hame", &str_value, NULL, NULL)<0) 
	    strcpy(backFilter, formula);
	  else 
	    strcpy(backFilter, str_value );
  	} else {
          SDDS_Bomb("Back filter material is not specified.\n");
  	}
	if (backFilterDensity<=0) {
	  if (get_material_property(backFilter, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: back filter density can not be obtained from database, has to be provided.");
	  backFilterDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;
      case SET_TARGET:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	targetDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &targetThickness, 1, 0,
                          "density", SDDS_DOUBLE, &targetDensity, 1, 0,
			  "angle", SDDS_DOUBLE, &thetaIn, 1, 0, 
			  "thetaIn", SDDS_DOUBLE, &thetaIn, 1, 0, 
                          NULL))
          SDDS_Bomb("invalid -target syntax");
        if (material!=NULL) { 
          if (verbose > 3) printf(" User provides target material name. \n");
  	  strcpy(matrix, material);
  	  if (formula!=NULL) {
    	    strcpy(matrixFormula, formula);
  	  } else {
	     if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	       strcpy(matrixFormula, material);
	    else 
	      strcpy(matrixFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
          if (verbose > 3) printf(" User provides target formula. \n");
  	  strcpy(matrixFormula, formula);
	  if (get_material_property(formula, "formula", &str_value, NULL, NULL)<0) 
	    strcpy(matrix, formula);
	  else 
	    strcpy(matrix, str_value );
  	} else {
          SDDS_Bomb("target formula not specified.\n");
  	}
	if (targetDensity<=0) {
	  if (get_material_property(matrix, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: back filter density can not be obtained from database, has to be provided.");
	  targetDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;
      case SET_TRACE:
	nLines =0;
	prevKeyword = NULL;
	for (i=1; i<s_arg[i_arg].n_items; i++) {
	  if ((valuePtr=strstr(s_arg[i_arg].list[i], "="))) {
	    len = valuePtr - s_arg[i_arg].list[i];
	    valuePtr++ ;
	    prevKeyword =  s_arg[i_arg].list[i];
	    if (strncmp_case_insensitive(prevKeyword, "element", len)==0) {
	      strcpy(trace, valuePtr);
	    } else if (strncmp_case_insensitive(prevKeyword, "concentration", len)==0) {
	      if (!get_double(&traceConcentration, valuePtr))
		SDDS_Bomb("Invalid value provided for trace concentration!");
	      /*traceConcentration = atof(valuePtr);*/
	    } else if (strncmp_case_insensitive(prevKeyword, "lines", len)==0) {
	      SDDS_CopyString(&traceLineList[nLines], valuePtr);
	      nLines ++;
	    } else {
	      if (!prevKeyword || strncmp_case_insensitive(prevKeyword, "lines", len)) {
		fprintf(stderr, "Unknown trace option keyword - %s provided.\n", s_arg[i_arg].list[i]);
		exit(1);
	      }
	    }
	  } else {
	    if (!prevKeyword || strncmp_case_insensitive(prevKeyword, "lines", len)) {
	      fprintf(stderr, "Unknown trace option keyword - %s provided.\n", s_arg[i_arg].list[i]);
	      exit(1);
	    }
	    SDDS_CopyString(&traceLineList[nLines], s_arg[i_arg].list[i]);
	    nLines ++;
	  }
	}
	strcpy(xrfLines,"");
	if (nLines > 0) {
	  for (i=0; i<nLines; i++) {
	    if ((xrfLineID[i] = getLineID(traceLineList[i]))>0) {
	      fprintf(stderr, "Error: invalid trace line provided - %s\n", traceLineList[i]);
	      exit(1);
	    }
	    strcat(xrfLines, traceLineList[i]);
	    strcat(xrfLines, " ");
	  }
	} else 
	  SDDS_Bomb("No trace line provided.");
	break;
      case SET_DETECTOR:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	detectDensity=0;
	upstream=0;
	detectThetaX = 90.0;
	detectThetaY = 90.0;	
	detectThetaZ = 00.0;	
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &detectThickness, 1, 0,
                          "density", SDDS_DOUBLE, &detectDensity, 1, 0,
			  "xcenter", SDDS_DOUBLE, &xDet, 1, 0,
			  "ycenter", SDDS_DOUBLE, &yDet, 1, 0,
			  "zcenter", SDDS_DOUBLE, &zDet, 1, 0,
			  "xwidth", SDDS_DOUBLE, &xWidth, 1, 0,
			  "ywidth", SDDS_DOUBLE, &yWidth, 1, 0,
			  "xtheta", SDDS_DOUBLE, &detectThetaX, 1, 0,
			  "ytheta", SDDS_DOUBLE, &detectThetaY, 1, 0,
			  "ztheta", SDDS_DOUBLE, &detectThetaZ, 1, 0,
			  "nx", SDDS_LONG, &nxDet, 1, 0,
			  "ny", SDDS_LONG, &nyDet, 1, 0,
			  "downstream", SDDS_LONG, &detectDownstream, 1, 0,
			  "upstream", SDDS_LONG, &upstream, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -detector syntax");
	if (upstream)  detectDownstream = 0;
	if (detectThetaZ < 0.0)  detectThetaZ = fabs(detectThetaZ);	
	if (detectThetaX < 89.999999) {
  	  detectThetaY = acos(1.0 - pow(cos(detectThetaZ*degToRad), 2) -  pow(cos(detectThetaX*degToRad), 2)) / degToRad;
	} else {
  	  detectThetaX = acos(1.0 - pow(cos(detectThetaZ*degToRad), 2) -  pow(cos(detectThetaY*degToRad), 2)) / degToRad;
	}
        if (material!=NULL) { 
  	  strcpy(detector, material);
  	  if (formula!=NULL) {
    	    strcpy(detectorFormula, formula);
  	  } else {
	     if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	       strcpy(detectorFormula, material);
	    else 
	      strcpy(detectorFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
  	  strcpy(detectorFormula, formula);
	  if (get_material_property(formula, "name", &str_value, NULL, NULL)<0) 
	    strcpy(detector, formula);
	  else 
	    strcpy(detector, str_value );
  	} else {
  	  strcpy(detector, "Silicon");
  	  strcpy(detectorFormula, "Si");
  	}
	if (detectDensity<=0) {
	  if (get_material_property(detector, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: back filter density can not be obtained from database, has to be provided.");
	  detectDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
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
  if (!outputfile)
    SDDS_Bomb("output file not provided.");
  /* Validate element symbol of trace */
  Z = SymbolToAtomicNumber ( trace ); 
  if (Z <= 0 || Z >= 120 )
    SDDS_Bomb("trace element not valid.");
  if (!nLines)
    SDDS_Bomb("No trace lines provided.");
  if (verbose > 2) fprintf(stderr, "ok here\n");
  if (mode<0 || mode>14 || (mode>=5 && mode<10)) {
    fprintf(stderr, "sddsxrf: mode is not recognized: %ld \n Exiting ....\n", mode);
    exit(1);
  }
  if (!inputfile) {
    /*get eneryg from provided range */
    if (npts>MAXNPTS || npts<1 || eEnd==eStart)
      SDDS_Bomb("Invalid or none energy range provided.");
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
  if (verbose) {
    /* Confirm minimum information about the calculation */
    printf("sddsxrf:  mode = %ld\n",  mode);
    printf("trace = %s, concentration = %f, lines = ", trace, traceConcentration);
    for ( i = 0; i < nLines; i++ ) printf( "%s,", traceLineList[i]);
    printf("\n");
    printf("matrix = %s, matrixFormula = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3), thetaIn = %0.1f (deg),\n", 
	   matrix, matrixFormula, targetThickness, targetDensity, thetaIn );
    if (mode == 2 || mode == 12) {
      printf("detectDownstream = %ld\n", detectDownstream);
    } else if (mode == 3 || mode == 4 || mode == 13 || mode == 14) {
      printf("detector = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3), xwidth = %0.1f (mm), ywidth = %0.1f (mm)\n", 
	     detector, detectThickness, detectDensity, xWidth, yWidth );
      printf("detectorFormula = %s, xcenter = %0.1f (mm), ycenter = %0.1f (mm), zcenter = %0.1f (mm), nx = %ld, ny = %ld\n", 
	     detectorFormula, xDet, yDet, zDet, nxDet, nyDet);
      printf("detectThetaX = %0.2f (deg), detectThetaY = %0.2f (deg), detectThetaZ = %0.2f (deg), detectDownstream = %ld\n", 
	     detectThetaX, detectThetaY, detectThetaZ, detectDownstream);
      if (mode == 4 || mode == 14) {
	printf("frontFilter = %s, frontFilterFormula = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3)\n", 
	       frontFilter, frontFilterFormula, frontFilterThickness, frontFilterDensity);
	printf("backFilter = %s, backFilterFormula = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3)\n", 
	       backFilter, backFilterFormula, backFilterThickness, backFilterDensity);
      }
    }
  }
  /* Initialize xraylib */
  XRayInit();
  SetupOutputFile(outputfile, &SDDSout, mode,  &SDDSin, colMatches, columnMatch, parMatches, parMatch);
  for  (i=0; i<pages; i++) {
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
    /* Calculations for matrixTotalCS cross section: 
       0  <= mode <= 4,  x-ray fluorescence cross sections from an compound matrix
       10 <= mode <= 14, x-ray fluorescence cross sections from a elemental matrix    */
    getXRFCS(mode);
    if (!SDDS_StartPage(&SDDSout, npts) ||
	!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sddsxrf_mode", mode, "Matrix", matrix, "MatrixFormula", matrixFormula, 
	                      "Trace", trace, "xrfLines", xrfLines, "xrfLineEnergy0", xrfLineEnergy[0], NULL) ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, energy, npts, "PhotonEnergy") ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, matrixTotalCS, npts, "MatrixTotalCS") ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, traceXRFCS, npts, "TraceXRFCS"))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (inputfile && (colMatches || parMatches) && !SDDS_CopyPage(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if ((mode>0 && mode<5) || mode>=11) {
      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "TargetThickness", targetThickness,
			      "TargetDensity", targetDensity,
			      "TraceConcentration", traceConcentration, "ThetaIn", thetaIn,
			      "ThetaOut", thetaOut, NULL) ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, totalXRFYield, npts, "TotalXRFYield"))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==2 || mode==12 ) {
      if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, xrfFluxDensity , npts, "xrfFluxDensity"))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==3 || mode==4 || mode==13 || mode==14) {
      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Detector", detector,
			      "xDet", xDet, "yDet", yDet, "zDet", zDet, "xWidth", xWidth, "yWidth", yWidth,
			      "nx", nxDet, "ny", nyDet, NULL) ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, xrfFluxFactor, npts, "xrfFluxFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, xrfPowerFactor, npts, "xrfPowerFactor"))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==4 || mode==14) {
      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                              "frontFilterFormula",frontFilterFormula,
			      "frontFilterThickness", frontFilterThickness,
			      "frontFilterDensity", frontFilterDensity,
			      "backFilterFormula", backFilterFormula,
			      "backFilterThickness", backFilterThickness,
			      "backFilterDensity", backFilterDensity, NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    free(energy); energy=NULL;
  }
  if (inputfile && !SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (nLines)
    for (i=0; i<nLines; i++)
      free(traceLineList[i]);
  free_scanargs(&s_arg, argc);
  return 0;
}




/************************************
  getLineID (name)
************************************/
int getLineID (char *name)
{
  int value;
  
  if(strcmp(name, "KA") == 0)  { value = KA_LINE;
  } else   if(strcmp(name, "KB") == 0)  { value = KB_LINE;
  } else   if(strcmp(name, "LA") == 0)  { value = LA_LINE;
  } else   if(strcmp(name, "LB") == 0)  { value = LB_LINE;
  /* Single lines */
  } else   if(strcmp(name, "KA1") == 0) { value = KA1_LINE;
  } else   if(strcmp(name, "KA2") == 0) { value = KA2_LINE;
  } else   if(strcmp(name, "KB1") == 0) { value = KB1_LINE;
  } else   if(strcmp(name, "KB2") == 0) { value = KB2_LINE;
  } else   if(strcmp(name, "KB3") == 0) { value = KB3_LINE;
  } else   if(strcmp(name, "KB4") == 0) { value = KB4_LINE;
  } else   if(strcmp(name, "KB5") == 0) { value = KB5_LINE;
  } else   if(strcmp(name, "LA1") == 0) { value = LA1_LINE;
  } else   if(strcmp(name, "LA2") == 0) { value = LA2_LINE;
  } else   if(strcmp(name, "LB1") == 0) { value = LB1_LINE;
  } else   if(strcmp(name, "LB2") == 0) { value = LB2_LINE;
  } else   if(strcmp(name, "LB3") == 0) { value = LB3_LINE;
  } else   if(strcmp(name, "LB4") == 0) { value = LB4_LINE;
  } else   if(strcmp(name, "LB5") == 0) { value = LB5_LINE;
  } else   if(strcmp(name, "LB6") == 0) { value = LB6_LINE;
  } else   if(strcmp(name, "LB7") == 0) { value = LB7_LINE;
  } else   if(strcmp(name, "LB9") == 0) { value = LB9_LINE;
  } else   if(strcmp(name, "LB10") == 0){ value = LB10_LINE;
  } else   if(strcmp(name, "LB15") == 0){ value = LB15_LINE;
  } else   if(strcmp(name, "LB17") == 0){ value = LB17_LINE;
  } else   if(strcmp(name, "LG1") == 0) { value = LG1_LINE;
  } else   if(strcmp(name, "LG2") == 0) { value = LG2_LINE;
  } else   if(strcmp(name, "LG3") == 0) { value = LG3_LINE;
  } else   if(strcmp(name, "LG4") == 0) { value = LG4_LINE;
  } else   if(strcmp(name, "LG5") == 0) { value = LG5_LINE;
  } else   if(strcmp(name, "LG6") == 0) { value = LG6_LINE;
  } else   if(strcmp(name, "LG8") == 0) { value = LG8_LINE;
  } else   if(strcmp(name, "LE") == 0)  { value = LE_LINE;
  } else   if(strcmp(name, "LL") == 0)  { value = LL_LINE;
  } else   if(strcmp(name, "LS") == 0)  { value = LS_LINE;
  } else   if(strcmp(name, "LT") == 0)  { value = LT_LINE;
  } else   if(strcmp(name, "LU") == 0)  { value = LU_LINE;
  } else   if(strcmp(name, "LV") == 0)  { value = LV_LINE;
  } else   if(strcmp(name, "MA1") == 0) { value = MA1_LINE;
  } else   if(strcmp(name, "MA2") == 0) { value = MA2_LINE;
  } else   if(strcmp(name, "MB") == 0)  { value = MB_LINE;
  } else   if(strcmp(name, "MG") == 0)  { value = MG_LINE;
  } else   if(strcmp(name, "KL1") == 0)  { value = KL1_LINE;
  } else   if(strcmp(name, "KL2") == 0)  { value = KL2_LINE;
  } else   if(strcmp(name, "KL3") == 0)  { value = KL3_LINE;
  } else   if(strcmp(name, "KM1") == 0)  { value = KM1_LINE;
  } else   if(strcmp(name, "KM2") == 0)  { value = KM2_LINE;
  } else   if(strcmp(name, "KM3") == 0)  { value = KM3_LINE;
  } else   if(strcmp(name, "KM4") == 0)  { value = KM4_LINE;
  } else   if(strcmp(name, "KM5") == 0)  { value = KM5_LINE;
  } else   if(strcmp(name, "KN1") == 0)  { value = KN1_LINE;
  } else   if(strcmp(name, "KN2") == 0)  { value = KN2_LINE;
  } else   if(strcmp(name, "KN3") == 0)  { value = KN3_LINE;
  } else   if(strcmp(name, "KN4") == 0)  { value = KN4_LINE;
  } else   if(strcmp(name, "KN5") == 0)  { value = KN5_LINE;
  } else   if(strcmp(name, "KN6") == 0)  { value = KN6_LINE;
  } else   if(strcmp(name, "KN7") == 0)  { value = KN7_LINE;
  } else   if(strcmp(name, "KO1") == 0)  { value = KO1_LINE;
  } else   if(strcmp(name, "KO2") == 0)  { value = KO2_LINE;
  } else   if(strcmp(name, "KO3") == 0)  { value = KO3_LINE;
  } else   if(strcmp(name, "KO4") == 0)  { value = KO4_LINE;
  } else   if(strcmp(name, "KO5") == 0)  { value = KO5_LINE;
  } else   if(strcmp(name, "KO6") == 0)  { value = KO6_LINE;
  } else   if(strcmp(name, "KO7") == 0)  { value = KO7_LINE;
  } else   if(strcmp(name, "KP1") == 0)  { value = KP1_LINE;
  } else   if(strcmp(name, "KP2") == 0)  { value = KP2_LINE;
  } else   if(strcmp(name, "KP3") == 0)  { value = KP3_LINE;
  } else   if(strcmp(name, "KP4") == 0)  { value = KP4_LINE;
  } else   if(strcmp(name, "KP5") == 0)  { value = KP5_LINE;
  } else   if(strcmp(name, "L1L2") == 0)  { value = L1L2_LINE;
  } else   if(strcmp(name, "L1L3") == 0)  { value = L1L3_LINE;
  } else   if(strcmp(name, "L1M1") == 0)  { value = L1M1_LINE;
  } else   if(strcmp(name, "L1M2") == 0)  { value = L1M2_LINE;
  } else   if(strcmp(name, "L1M3") == 0)  { value = L1M3_LINE;
  } else   if(strcmp(name, "L1M4") == 0)  { value = L1M4_LINE;
  } else   if(strcmp(name, "L1M5") == 0)  { value = L1M5_LINE;
  } else   if(strcmp(name, "L1N1") == 0)  { value = L1N1_LINE;
  } else   if(strcmp(name, "L1N2") == 0)  { value = L1N2_LINE;
  } else   if(strcmp(name, "L1N3") == 0)  { value = L1N3_LINE;
  } else   if(strcmp(name, "L1N4") == 0)  { value = L1N4_LINE;
  } else   if(strcmp(name, "L1N5") == 0)  { value = L1N5_LINE;
  } else   if(strcmp(name, "L1N6") == 0)  { value = L1N6_LINE;
  } else   if(strcmp(name, "L1N7") == 0)  { value = L1N7_LINE;
  } else   if(strcmp(name, "L1N67") == 0) { value = L1N67_LINE;
  } else   if(strcmp(name, "L1O1") == 0)  { value = L1O1_LINE;
  } else   if(strcmp(name, "L1O2") == 0)  { value = L1O2_LINE;
  } else   if(strcmp(name, "L1O3") == 0)  { value = L1O3_LINE;
  } else   if(strcmp(name, "L1O4") == 0)  { value = L1O4_LINE;
  } else   if(strcmp(name, "L1O45") == 0) { value = L1O45_LINE;
  } else   if(strcmp(name, "L1O5") == 0)  { value = L1O5_LINE;
  } else   if(strcmp(name, "L1O6") == 0)  { value = L1O6_LINE;
  } else   if(strcmp(name, "L1O7") == 0)  { value = L1O7_LINE;
  } else   if(strcmp(name, "L1P1") == 0)  { value = L1P1_LINE;
  } else   if(strcmp(name, "L1P2") == 0)  { value = L1P2_LINE;
  } else   if(strcmp(name, "L1P23") == 0) { value = L1P23_LINE;
  } else   if(strcmp(name, "L1P3") == 0)  { value = L1P3_LINE;
  } else   if(strcmp(name, "L1P4") == 0)  { value = L1P4_LINE;
  } else   if(strcmp(name, "L1P5") == 0)  { value = L1P5_LINE;
  } else   if(strcmp(name, "L2L3") == 0)  { value = L2L3_LINE;
  } else   if(strcmp(name, "L2M1") == 0)  { value = L2M1_LINE;
  } else   if(strcmp(name, "L2M2") == 0)  { value = L2M2_LINE;
  } else   if(strcmp(name, "L2M3") == 0)  { value = L2M3_LINE;
  } else   if(strcmp(name, "L2M4") == 0)  { value = L2M4_LINE;
  } else   if(strcmp(name, "L2M5") == 0)  { value = L2M5_LINE;
  } else   if(strcmp(name, "L2N1") == 0)  { value = L2N1_LINE;
  } else   if(strcmp(name, "L2N2") == 0)  { value = L2N2_LINE;
  } else   if(strcmp(name, "L2N3") == 0)  { value = L2N3_LINE;
  } else   if(strcmp(name, "L2N4") == 0)  { value = L2N4_LINE;
  } else   if(strcmp(name, "L2N5") == 0)  { value = L2N5_LINE;
  } else   if(strcmp(name, "L2N6") == 0)  { value = L2N6_LINE;
  } else   if(strcmp(name, "L2N7") == 0)  { value = L2N7_LINE;
  } else   if(strcmp(name, "L2O1") == 0)  { value = L2O1_LINE;
  } else   if(strcmp(name, "L2O2") == 0)  { value = L2O2_LINE;
  } else   if(strcmp(name, "L2O3") == 0)  { value = L2O3_LINE;
  } else   if(strcmp(name, "L2O4") == 0)  { value = L2O4_LINE;
  } else   if(strcmp(name, "L2O5") == 0)  { value = L2O5_LINE;
  } else   if(strcmp(name, "L2O6") == 0)  { value = L2O6_LINE;
  } else   if(strcmp(name, "L2O7") == 0)  { value = L2O7_LINE;
  } else   if(strcmp(name, "L2P1") == 0)  { value = L2P1_LINE;
  } else   if(strcmp(name, "L2P2") == 0)  { value = L2P2_LINE;
  } else   if(strcmp(name, "L2P23") == 0) { value = L2P23_LINE;
  } else   if(strcmp(name, "L2P3") == 0)  { value = L2P3_LINE;
  } else   if(strcmp(name, "L2P4") == 0)  { value = L2P4_LINE;
  } else   if(strcmp(name, "L2P5") == 0)  { value = L2P5_LINE;
  } else   if(strcmp(name, "L2Q1") == 0)  { value = L2Q1_LINE;
  } else   if(strcmp(name, "L3M1") == 0)  { value = L3M1_LINE;
  } else   if(strcmp(name, "L3M2") == 0)  { value = L3M2_LINE;
  } else   if(strcmp(name, "L3M3") == 0)  { value = L3M3_LINE;
  } else   if(strcmp(name, "L3M4") == 0)  { value = L3M4_LINE;
  } else   if(strcmp(name, "L3M5") == 0)  { value = L3M5_LINE;
  } else   if(strcmp(name, "L3N1") == 0)  { value = L3N1_LINE;
  } else   if(strcmp(name, "L3N2") == 0)  { value = L3N2_LINE;
  } else   if(strcmp(name, "L3N3") == 0)  { value = L3N3_LINE;
  } else   if(strcmp(name, "L3N4") == 0)  { value = L3N4_LINE;
  } else   if(strcmp(name, "L3N5") == 0)  { value = L3N5_LINE;
  } else   if(strcmp(name, "L3N6") == 0)  { value = L3N6_LINE;
  } else   if(strcmp(name, "L3N7") == 0)  { value = L3N7_LINE;
  } else   if(strcmp(name, "L3O1") == 0)  { value = L3O1_LINE;
  } else   if(strcmp(name, "L3O2") == 0)  { value = L3O2_LINE;
  } else   if(strcmp(name, "L3O3") == 0)  { value = L3O3_LINE;
  } else   if(strcmp(name, "L3O4") == 0)  { value = L3O4_LINE;
  } else   if(strcmp(name, "L3O45") == 0) { value = L3O45_LINE;
  } else   if(strcmp(name, "L3O5") == 0)  { value = L3O5_LINE;
  } else   if(strcmp(name, "L3O6") == 0)  { value = L3O6_LINE;
  } else   if(strcmp(name, "L3O7") == 0)  { value = L3O7_LINE;
  } else   if(strcmp(name, "L3P1") == 0)  { value = L3P1_LINE;
  } else   if(strcmp(name, "L3P2") == 0)  { value = L3P2_LINE;
  } else   if(strcmp(name, "L3P23") == 0) { value = L3P23_LINE;
  } else   if(strcmp(name, "L3P3") == 0)  { value = L3P3_LINE;
  } else   if(strcmp(name, "L3P4") == 0)  { value = L3P4_LINE;
  } else   if(strcmp(name, "L3P45") == 0) { value = L3P45_LINE;
  } else   if(strcmp(name, "L3P5") == 0)  { value = L3P5_LINE;
  } else   if(strcmp(name, "L3Q1") == 0)  { value = L3Q1_LINE;
  } else   if(strcmp(name, "M1M2") == 0)  { value = M1M2_LINE;
  } else   if(strcmp(name, "M1M3") == 0)  { value = M1M3_LINE;
  } else   if(strcmp(name, "M1M4") == 0)  { value = M1M4_LINE;
  } else   if(strcmp(name, "M1M5") == 0)  { value = M1M5_LINE;
  } else   if(strcmp(name, "M1N1") == 0)  { value = M1N1_LINE;
  } else   if(strcmp(name, "M1N2") == 0)  { value = M1N2_LINE;
  } else   if(strcmp(name, "M1N3") == 0)  { value = M1N3_LINE;
  } else   if(strcmp(name, "M1N4") == 0)  { value = M1N4_LINE;
  } else   if(strcmp(name, "M1N5") == 0)  { value = M1N5_LINE;
  } else   if(strcmp(name, "M1N6") == 0)  { value = M1N6_LINE;
  } else   if(strcmp(name, "M1N7") == 0)  { value = M1N7_LINE;
  } else   if(strcmp(name, "M1O1") == 0)  { value = M1O1_LINE;
  } else   if(strcmp(name, "M1O2") == 0)  { value = M1O2_LINE;
  } else   if(strcmp(name, "M1O3") == 0)  { value = M1O3_LINE;
  } else   if(strcmp(name, "M1O4") == 0)  { value = M1O4_LINE;
  } else   if(strcmp(name, "M1O5") == 0)  { value = M1O5_LINE;
  } else   if(strcmp(name, "M1O6") == 0)  { value = M1O6_LINE;
  } else   if(strcmp(name, "M1O7") == 0)  { value = M1O7_LINE;
  } else   if(strcmp(name, "M1P1") == 0)  { value = M1P1_LINE;
  } else   if(strcmp(name, "M1P2") == 0)  { value = M1P2_LINE;
  } else   if(strcmp(name, "M1P3") == 0)  { value = M1P3_LINE;
  } else   if(strcmp(name, "M1P4") == 0)  { value = M1P4_LINE;
  } else   if(strcmp(name, "M1P5") == 0)  { value = M1P5_LINE;
  } else   if(strcmp(name, "M2M3") == 0)  { value = M2M3_LINE;
  } else   if(strcmp(name, "M2M4") == 0)  { value = M2M4_LINE;
  } else   if(strcmp(name, "M2M5") == 0)  { value = M2M5_LINE;
  } else   if(strcmp(name, "M2N1") == 0)  { value = M2N1_LINE;
  } else   if(strcmp(name, "M2N2") == 0)  { value = M2N2_LINE;
  } else   if(strcmp(name, "M2N3") == 0)  { value = M2N3_LINE;
  } else   if(strcmp(name, "M2N4") == 0)  { value = M2N4_LINE;
  } else   if(strcmp(name, "M2N5") == 0)  { value = M2N5_LINE;
  } else   if(strcmp(name, "M2N6") == 0)  { value = M2N6_LINE;
  } else   if(strcmp(name, "M2N7") == 0)  { value = M2N7_LINE;
  } else   if(strcmp(name, "M2O1") == 0)  { value = M2O1_LINE;
  } else   if(strcmp(name, "M2O2") == 0)  { value = M2O2_LINE;
  } else   if(strcmp(name, "M2O3") == 0)  { value = M2O3_LINE;
  } else   if(strcmp(name, "M2O4") == 0)  { value = M2O4_LINE;
  } else   if(strcmp(name, "M2O5") == 0)  { value = M2O5_LINE;
  } else   if(strcmp(name, "M2O6") == 0)  { value = M2O6_LINE;
  } else   if(strcmp(name, "M2O7") == 0)  { value = M2O7_LINE;
  } else   if(strcmp(name, "M2P1") == 0)  { value = M2P1_LINE;
  } else   if(strcmp(name, "M2P2") == 0)  { value = M2P2_LINE;
  } else   if(strcmp(name, "M2P3") == 0)  { value = M2P3_LINE;
  } else   if(strcmp(name, "M2P4") == 0)  { value = M2P4_LINE;
  } else   if(strcmp(name, "M2P5") == 0)  { value = M2P5_LINE;
  } else   if(strcmp(name, "M3M4") == 0)  { value = M3M4_LINE;
  } else   if(strcmp(name, "M3M5") == 0)  { value = M3M5_LINE;
  } else   if(strcmp(name, "M3N1") == 0)  { value = M3N1_LINE;
  } else   if(strcmp(name, "M3N2") == 0)  { value = M3N2_LINE;
  } else   if(strcmp(name, "M3N3") == 0)  { value = M3N3_LINE;
  } else   if(strcmp(name, "M3N4") == 0)  { value = M3N4_LINE;
  } else   if(strcmp(name, "M3N5") == 0)  { value = M3N5_LINE;
  } else   if(strcmp(name, "M3N6") == 0)  { value = M3N6_LINE;
  } else   if(strcmp(name, "M3N7") == 0)  { value = M3N7_LINE;
  } else   if(strcmp(name, "M3O1") == 0)  { value = M3O1_LINE;
  } else   if(strcmp(name, "M3O2") == 0)  { value = M3O2_LINE;
  } else   if(strcmp(name, "M3O3") == 0)  { value = M3O3_LINE;
  } else   if(strcmp(name, "M3O4") == 0)  { value = M3O4_LINE;
  } else   if(strcmp(name, "M3O5") == 0)  { value = M3O5_LINE;
  } else   if(strcmp(name, "M3O6") == 0)  { value = M3O6_LINE;
  } else   if(strcmp(name, "M3O7") == 0)  { value = M3O7_LINE;
  } else   if(strcmp(name, "M3P1") == 0)  { value = M3P1_LINE;
  } else   if(strcmp(name, "M3P2") == 0)  { value = M3P2_LINE;
  } else   if(strcmp(name, "M3P3") == 0)  { value = M3P3_LINE;
  } else   if(strcmp(name, "M3P4") == 0)  { value = M3P4_LINE;
  } else   if(strcmp(name, "M3P5") == 0)  { value = M3P5_LINE;
  } else   if(strcmp(name, "M3Q1") == 0)  { value = M3Q1_LINE;
  } else   if(strcmp(name, "M4M5") == 0)  { value = M4M5_LINE;
  } else   if(strcmp(name, "M4N1") == 0)  { value = M4N1_LINE;
  } else   if(strcmp(name, "M4N2") == 0)  { value = M4N2_LINE;
  } else   if(strcmp(name, "M4N3") == 0)  { value = M4N3_LINE;
  } else   if(strcmp(name, "M4N4") == 0)  { value = M4N4_LINE;
  } else   if(strcmp(name, "M4N5") == 0)  { value = M4N5_LINE;
  } else   if(strcmp(name, "M4N6") == 0)  { value = M4N6_LINE;
  } else   if(strcmp(name, "M4N7") == 0)  { value = M4N7_LINE;
  } else   if(strcmp(name, "M4O1") == 0)  { value = M4O1_LINE;
  } else   if(strcmp(name, "M4O2") == 0)  { value = M4O2_LINE;
  } else   if(strcmp(name, "M4O3") == 0)  { value = M4O3_LINE;
  } else   if(strcmp(name, "M4O4") == 0)  { value = M4O4_LINE;
  } else   if(strcmp(name, "M4O5") == 0)  { value = M4O5_LINE;
  } else   if(strcmp(name, "M4O6") == 0)  { value = M4O6_LINE;
  } else   if(strcmp(name, "M4O7") == 0)  { value = M4O7_LINE;
  } else   if(strcmp(name, "M4P1") == 0)  { value = M4P1_LINE;
  } else   if(strcmp(name, "M4P2") == 0)  { value = M4P2_LINE;
  } else   if(strcmp(name, "M4P3") == 0)  { value = M4P3_LINE;
  } else   if(strcmp(name, "M4P4") == 0)  { value = M4P4_LINE;
  } else   if(strcmp(name, "M4P5") == 0)  { value = M4P5_LINE;
  } else   if(strcmp(name, "M5N1") == 0)  { value = M5N1_LINE;
  } else   if(strcmp(name, "M5N2") == 0)  { value = M5N2_LINE;
  } else   if(strcmp(name, "M5N3") == 0)  { value = M5N3_LINE;
  } else   if(strcmp(name, "M5N4") == 0)  { value = M5N4_LINE;
  } else   if(strcmp(name, "M5N5") == 0)  { value = M5N5_LINE;
  } else   if(strcmp(name, "M5N6") == 0)  { value = M5N6_LINE;
  } else   if(strcmp(name, "M5N7") == 0)  { value = M5N7_LINE;
  } else   if(strcmp(name, "M5O1") == 0)  { value = M5O1_LINE;
  } else   if(strcmp(name, "M5O2") == 0)  { value = M5O2_LINE;
  } else   if(strcmp(name, "M5O3") == 0)  { value = M5O3_LINE;
  } else   if(strcmp(name, "M5O4") == 0)  { value = M5O4_LINE;
  } else   if(strcmp(name, "M5O5") == 0)  { value = M5O5_LINE;
  } else   if(strcmp(name, "M5O6") == 0)  { value = M5O6_LINE;
  } else   if(strcmp(name, "M5O7") == 0)  { value = M5O7_LINE;
  } else   if(strcmp(name, "M5P1") == 0)  { value = M5P1_LINE;
  } else   if(strcmp(name, "M5P2") == 0)  { value = M5P2_LINE;
  } else   if(strcmp(name, "M5P3") == 0)  { value = M5P3_LINE;
  } else   if(strcmp(name, "M5P4") == 0)  { value = M5P4_LINE;
  } else   if(strcmp(name, "M5P5") == 0)  { value = M5P5_LINE;
  } else   if(strcmp(name, "N1N2") == 0)  { value = N1N2_LINE;
  } else   if(strcmp(name, "N1N3") == 0)  { value = N1N3_LINE;
  } else   if(strcmp(name, "N1N4") == 0)  { value = N1N4_LINE;
  } else   if(strcmp(name, "N1N5") == 0)  { value = N1N5_LINE;
  } else   if(strcmp(name, "N1N6") == 0)  { value = N1N6_LINE;
  } else   if(strcmp(name, "N1N7") == 0)  { value = N1N7_LINE;
  } else   if(strcmp(name, "N1O1") == 0)  { value = N1O1_LINE;
  } else   if(strcmp(name, "N1O2") == 0)  { value = N1O2_LINE;
  } else   if(strcmp(name, "N1O3") == 0)  { value = N1O3_LINE;
  } else   if(strcmp(name, "N1O4") == 0)  { value = N1O4_LINE;
  } else   if(strcmp(name, "N1O5") == 0)  { value = N1O5_LINE;
  } else   if(strcmp(name, "N1O6") == 0)  { value = N1O6_LINE;
  } else   if(strcmp(name, "N1O7") == 0)  { value = N1O7_LINE;
  } else   if(strcmp(name, "N1P1") == 0)  { value = N1P1_LINE;
  } else   if(strcmp(name, "N1P2") == 0)  { value = N1P2_LINE;
  } else   if(strcmp(name, "N1P3") == 0)  { value = N1P3_LINE;
  } else   if(strcmp(name, "N1P4") == 0)  { value = N1P4_LINE;
  } else   if(strcmp(name, "N1P5") == 0)  { value = N1P5_LINE;
  } else   if(strcmp(name, "N2N3") == 0)  { value = N2N3_LINE;
  } else   if(strcmp(name, "N2N4") == 0)  { value = N2N4_LINE;
  } else   if(strcmp(name, "N2N5") == 0)  { value = N2N5_LINE;
  } else   if(strcmp(name, "N2N6") == 0)  { value = N2N6_LINE;
  } else   if(strcmp(name, "N2N7") == 0)  { value = N2N7_LINE;
  } else   if(strcmp(name, "N2O1") == 0)  { value = N2O1_LINE;
  } else   if(strcmp(name, "N2O2") == 0)  { value = N2O2_LINE;
  } else   if(strcmp(name, "N2O3") == 0)  { value = N2O3_LINE;
  } else   if(strcmp(name, "N2O4") == 0)  { value = N2O4_LINE;
  } else   if(strcmp(name, "N2O5") == 0)  { value = N2O5_LINE;
  } else   if(strcmp(name, "N2O6") == 0)  { value = N2O6_LINE;
  } else   if(strcmp(name, "N2O7") == 0)  { value = N2O7_LINE;
  } else   if(strcmp(name, "N2P1") == 0)  { value = N2P1_LINE;
  } else   if(strcmp(name, "N2P2") == 0)  { value = N2P2_LINE;
  } else   if(strcmp(name, "N2P3") == 0)  { value = N2P3_LINE;
  } else   if(strcmp(name, "N2P4") == 0)  { value = N2P4_LINE;
  } else   if(strcmp(name, "N2P5") == 0)  { value = N2P5_LINE;
  } else   if(strcmp(name, "N3N4") == 0)  { value = N3N4_LINE;
  } else   if(strcmp(name, "N3N5") == 0)  { value = N3N5_LINE;
  } else   if(strcmp(name, "N3N6") == 0)  { value = N3N6_LINE;
  } else   if(strcmp(name, "N3N7") == 0)  { value = N3N7_LINE;
  } else   if(strcmp(name, "N3O1") == 0)  { value = N3O1_LINE;
  } else   if(strcmp(name, "N3O2") == 0)  { value = N3O2_LINE;
  } else   if(strcmp(name, "N3O3") == 0)  { value = N3O3_LINE;
  } else   if(strcmp(name, "N3O4") == 0)  { value = N3O4_LINE;
  } else   if(strcmp(name, "N3O5") == 0)  { value = N3O5_LINE;
  } else   if(strcmp(name, "N3O6") == 0)  { value = N3O6_LINE;
  } else   if(strcmp(name, "N3O7") == 0)  { value = N3O7_LINE;
  } else   if(strcmp(name, "N3P1") == 0)  { value = N3P1_LINE;
  } else   if(strcmp(name, "N3P2") == 0)  { value = N3P2_LINE;
  } else   if(strcmp(name, "N3P3") == 0)  { value = N3P3_LINE;
  } else   if(strcmp(name, "N3P4") == 0)  { value = N3P4_LINE;
  } else   if(strcmp(name, "N3P5") == 0)  { value = N3P5_LINE;
  } else   if(strcmp(name, "N4N5") == 0)  { value = N4N5_LINE;
  } else   if(strcmp(name, "N4N6") == 0)  { value = N4N6_LINE;
  } else   if(strcmp(name, "N4N7") == 0)  { value = N4N7_LINE;
  } else   if(strcmp(name, "N4O1") == 0)  { value = N4O1_LINE;
  } else   if(strcmp(name, "N4O2") == 0)  { value = N4O2_LINE;
  } else   if(strcmp(name, "N4O3") == 0)  { value = N4O3_LINE;
  } else   if(strcmp(name, "N4O4") == 0)  { value = N4O4_LINE;
  } else   if(strcmp(name, "N4O5") == 0)  { value = N4O5_LINE;
  } else   if(strcmp(name, "N4O6") == 0)  { value = N4O6_LINE;
  } else   if(strcmp(name, "N4O7") == 0)  { value = N4O7_LINE;
  } else   if(strcmp(name, "N4P1") == 0)  { value = N4P1_LINE;
  } else   if(strcmp(name, "N4P2") == 0)  { value = N4P2_LINE;
  } else   if(strcmp(name, "N4P3") == 0)  { value = N4P3_LINE;
  } else   if(strcmp(name, "N4P4") == 0)  { value = N4P4_LINE;
  } else   if(strcmp(name, "N4P5") == 0)  { value = N4P5_LINE;
  } else   if(strcmp(name, "N5N6") == 0)  { value = N5N6_LINE;
  } else   if(strcmp(name, "N5N7") == 0)  { value = N5N7_LINE;
  } else   if(strcmp(name, "N5O1") == 0)  { value = N5O1_LINE;
  } else   if(strcmp(name, "N5O2") == 0)  { value = N5O2_LINE;
  } else   if(strcmp(name, "N5O3") == 0)  { value = N5O3_LINE;
  } else   if(strcmp(name, "N5O4") == 0)  { value = N5O4_LINE;
  } else   if(strcmp(name, "N5O5") == 0)  { value = N5O5_LINE;
  } else   if(strcmp(name, "N5O6") == 0)  { value = N5O6_LINE;
  } else   if(strcmp(name, "N5O7") == 0)  { value = N5O7_LINE;
  } else   if(strcmp(name, "N5P1") == 0)  { value = N5P1_LINE;
  } else   if(strcmp(name, "N5P2") == 0)  { value = N5P2_LINE;
  } else   if(strcmp(name, "N5P3") == 0)  { value = N5P3_LINE;
  } else   if(strcmp(name, "N5P4") == 0)  { value = N5P4_LINE;
  } else   if(strcmp(name, "N5P5") == 0)  { value = N5P5_LINE;
  } else   if(strcmp(name, "N6N7") == 0)  { value = N6N7_LINE;
  } else   if(strcmp(name, "N6O1") == 0)  { value = N6O1_LINE;
  } else   if(strcmp(name, "N6O2") == 0)  { value = N6O2_LINE;
  } else   if(strcmp(name, "N6O3") == 0)  { value = N6O3_LINE;
  } else   if(strcmp(name, "N6O4") == 0)  { value = N6O4_LINE;
  } else   if(strcmp(name, "N6O5") == 0)  { value = N6O5_LINE;
  } else   if(strcmp(name, "N6O6") == 0)  { value = N6O6_LINE;
  } else   if(strcmp(name, "N6O7") == 0)  { value = N6O7_LINE;
  } else   if(strcmp(name, "N6P1") == 0)  { value = N6P1_LINE;
  } else   if(strcmp(name, "N6P2") == 0)  { value = N6P2_LINE;
  } else   if(strcmp(name, "N6P3") == 0)  { value = N6P3_LINE;
  } else   if(strcmp(name, "N6P4") == 0)  { value = N6P4_LINE;
  } else   if(strcmp(name, "N6P5") == 0)  { value = N6P5_LINE;
  } else   if(strcmp(name, "N7O1") == 0)  { value = N7O1_LINE;
  } else   if(strcmp(name, "N7O2") == 0)  { value = N7O2_LINE;
  } else   if(strcmp(name, "N7O3") == 0)  { value = N7O3_LINE;
  } else   if(strcmp(name, "N7O4") == 0)  { value = N7O4_LINE;
  } else   if(strcmp(name, "N7O5") == 0)  { value = N7O5_LINE;
  } else   if(strcmp(name, "N7O6") == 0)  { value = N7O6_LINE;
  } else   if(strcmp(name, "N7O7") == 0)  { value = N7O7_LINE;
  } else   if(strcmp(name, "N7P1") == 0)  { value = N7P1_LINE;
  } else   if(strcmp(name, "N7P2") == 0)  { value = N7P2_LINE;
  } else   if(strcmp(name, "N7P3") == 0)  { value = N7P3_LINE;
  } else   if(strcmp(name, "N7P4") == 0)  { value = N7P4_LINE;
  } else   if(strcmp(name, "N7P5") == 0)  { value = N7P5_LINE;
  } else { value = -999;
  }

  return value;
}


/**********************************************************
  getXRFCS ()
  
  This function retrieve the x-ray interaction constants 
  of an element using xraylib program/database package 
***********************************************************/

int getXRFCS (int mode)
{
  int    i, ju, jv, k, Z, Z0, Z1;
  float  energy_keV=0.0, xrfEnergy_keV=0.0;
  double detectMassThickness, targetMassThickness, cosThetaIn, cosThetaOut;
  double sumXRFCS, xrfCS, matrixXRFCS[32], detectXRFCS[32], backFilterCS[32];
  double sumXRFFluxDensity, effectiveCS, angleFactor, prefactor, pixelPrefactor;
  double pixelXRFFlux, sumXRFFlux, sumXRFAbsorb, du, dv, uStart, vStart, up, vp, xp, yp, zp, rSq;
  double cosThetaX, cosThetaY, cosThetaZ, sinThetaY, cosDetectThetaIn, detectorAbs;
  double frontFilterMassThickness, frontFilterCS, firstFilterTrans=1.0;
  double backFilterMassThickness,  secondFilterTrans=1.0;
  
  /* Calculate incoming and outgoing angle at the detector */
  cosThetaIn  = cos (thetaIn   * degToRad);
  cosThetaOut = zDet / sqrt(xDet * xDet + yDet * yDet + zDet * zDet);
  thetaOut    = acos( cosThetaOut ) / degToRad;
  if ( cosThetaIn < 1.0E-6 ) { printf("Incident angle too close to target surface. thetaIn = %f (deg)\n", thetaIn); }
  if ( cosThetaOut< 1.0E-6 ) { printf("Outgoing angle too close to target surface. thetaOut = %f (deg)\n", thetaOut); }
  
  /* Calculate directional cosine of incoming x-ray at detector surface */
  cosThetaX   = cos (detectThetaX * degToRad);
  cosThetaZ   = cos (detectThetaZ * degToRad);
  sinThetaY   = sqrt (cosThetaX * cosThetaX + cosThetaZ * cosThetaZ);
  cosThetaY   = sqrt (1.0 - sinThetaY * sinThetaY);
  
  Z  = SymbolToAtomicNumber ( trace );
  if (mode > 9 ) { 
    Z0 = SymbolToAtomicNumber ( frontFilterFormula ); 
    Z1 = SymbolToAtomicNumber ( matrixFormula ); 
  } else { 
    Z0 = Z1 = 0; 
  }
  targetMassThickness            = 0.1 * targetThickness * targetDensity;	// Target mass thickness in g/cm^2
  frontFilterMassThickness = 0.1 * frontFilterThickness * frontFilterDensity;	// Upstream filter mass thickness in g/cm^2
  detectMassThickness      = 0.1 * detectThickness * detectDensity;	// Detector mass thickness in g/cm^2
  angleFactor              = 0.25 / pi_const;				// Inverse of 4-PI solid angle
  if (detectDownstream < 1) {
    backFilterMassThickness = frontFilterMassThickness;			// XRF photon goes through front filter in reflective geometry
  } else {
    backFilterMassThickness = 0.1 * backFilterThickness * backFilterDensity;	// downstream filter mass thickness in g/cm^2
  }

  /* Integration parameters for detector */
  du = xWidth / nxDet;
  dv = yWidth / nyDet;
  uStart = -0.5 * (nxDet - 1.0) * du;
  vStart = - 0.5 * (nyDet - 1.0) * dv;

  if (verbose > 3) { printf( "getXRFCS: trace-Z = %d, targetMassThickness = %f\n", Z, targetMassThickness); }
  if (verbose > 2) { 
    printf( "\n Index  PhotonEnergy     TotalCS        XRFCS");
    if ( mode == 1 || mode == 11 ) { printf("     xrfYield"); }
    if ( mode == 2 || mode == 12 ) { printf("   xrfFluxDensity"); }
    if ( mode == 3 || mode == 4 || mode == 13 || mode == 14 ) { printf("     xrfFluxFactor  xrfPowerFactor"); }
    printf( "\n             (eV)        (cm^2/g)     (cm^2/g)");
    if ( mode == 2 || mode == 12 ) { printf("     (1/sterad)"); }
  }

  /* Retrieve matrix (self) absorption data for XRF */
  for ( k = 0; k < nLines; k++ ) { 
    xrfEnergy_keV = LineEnergy(Z, xrfLineID[k]);			/* trace fluorescence line energy in keV */
    xrfLineEnergy[k] = 1000.0 * xrfEnergy_keV;				/* trace fluorescence line energy in eV */
    matrixXRFCS[k] = CS_Total_CP (matrixFormula,   xrfEnergy_keV);	/* matrix absorption of trace fluorescence */
    detectXRFCS[k] = CS_Total_CP (detectorFormula, xrfEnergy_keV);	/* detector absorption of trace fluorescence */
    if (detectDownstream >= 1) {
      backFilterCS[k]    = CS_Total_CP (backFilterFormula, xrfEnergy_keV);	/* XRF cross section behind the target */
    } else {
      backFilterCS[k]    = CS_Total_CP (frontFilterFormula, xrfEnergy_keV);		/* XRF cross section in front of the target */
    }
    if (verbose > 4) { printf( "\ngetXRFCS: Index = %d, xrfEnergy_keV = %f, matrixXRFCS = %f", k, xrfEnergy_keV, matrixXRFCS[k]); }
  }

  /* MAIN LOOP */
  for ( i = 0; i < npts; i++ ) {
    energy_keV  = 0.001 * energy[i];	/* Excitation energy in keV */
    /* Retrieve total cross section of matrix and filter */
    if (mode < 10 ) {
      frontFilterCS = CS_Total_CP (frontFilterFormula, energy_keV);	/* filter absorption of incoming x-ray */
      matrixTotalCS[i]   = CS_Total_CP (matrixFormula, energy_keV);		/* matrix absorption of incoming x-ray */
    } else {
      frontFilterCS = CS_Total (Z0, energy_keV);	/* filter absorption of incoming x-ray */
      matrixTotalCS[i]   = CS_Total (Z1, energy_keV);	/* matrix absorption of incoming x-ray */
    }
    /* Prepare for loop covering each XRF line */
    sumXRFCS          = 0.0;
    sumXRFFluxDensity = 0.0;
    sumXRFFlux        = 0.0;
    sumXRFAbsorb      = 0.0;
    for ( k = 0; k < nLines; k++ ) { 
      xrfCS = CS_FluorLine(Z, xrfLineID[k], energy_keV);
      sumXRFCS = sumXRFCS + xrfCS;
      if (verbose > 5) { printf( "\n Sum lines: k = %d, xrfLineID = %d, xrfCS = %f, sumXRFCS = %f", k, xrfLineID[k], xrfCS, sumXRFCS); }
      
      /* Calculate upstream and downstream filter transmission factor */
      if ( mode == 4 || mode == 14 ) { 
        firstFilterTrans  = exp( -1.0 * frontFilterMassThickness * frontFilterCS  / cosThetaIn ); 
        secondFilterTrans = exp( -1.0 * backFilterMassThickness * backFilterCS[k] / cosThetaOut);
      } else {
        firstFilterTrans  = 1.0; 
        secondFilterTrans = 1.0;
      }

      /* Calculate xrf flux density, self absorption included */
      prefactor   = angleFactor * traceConcentration * xrfCS / cosThetaIn;
      if ( xrfCS < 1e-9 ) {
        /* Do nothing */
      } else if ( mode == 2 || mode == 12 ) {  /* Calculate XRF angular distribution */
        if (detectDownstream < 1) { /* Detector in front of target */
          effectiveCS = matrixTotalCS[i] / cosThetaIn + matrixXRFCS[k] / cosThetaOut;
        } else { /* Detector behind target */
          effectiveCS = matrixTotalCS[i] / cosThetaIn - matrixXRFCS[k] / cosThetaOut;
          prefactor   = prefactor * exp(-1.0 * matrixXRFCS[k] * targetMassThickness / cosThetaOut);
        }
        sumXRFFluxDensity = sumXRFFluxDensity + prefactor / effectiveCS * (1.0 - exp(-1.0 * targetMassThickness * effectiveCS));
        if (verbose > 5) { printf( "\ngetXRFCS: prefactor = %f, effectiveCS = %f", prefactor, effectiveCS); }
      } else if ( mode == 3 || mode == 4 || mode == 13 || mode == 14 ) { 
        /* Calculate XRF flux through a detector window, self absorption included */
        if (verbose > 9) { printf( "\n Detector Integration: k   ju  jv       xp       yp    zp    cosThetaOut cosDetectThetaIn prefactor FilterTrans1 FilterTrans1"); }
        for ( ju = 0; ju < nxDet; ju++ ) { 
          up = uStart + ju * du;
          for ( jv = 0; jv < nyDet; jv++ ) { 
            vp = vStart + jv * dv;
            if (fabs(sinThetaY) < 0.0001) {
              xp  = xDet + up ;
              yp  = yDet;
              zp  = zDet - vp;    
            } else {
              xp = xDet + (up * cosThetaZ - vp * cosThetaX * cosThetaY) / sinThetaY ;
              yp = yDet +  vp * sinThetaY;            
              zp = zDet - (up * cosThetaX + vp * cosThetaZ * cosThetaY) / sinThetaY ;    
            }
            rSq = xp * xp + yp * yp + zp * zp;
            cosThetaOut   = zp / sqrt(rSq);
            cosDetectThetaIn = (cosThetaX * xp + cosThetaY * yp + cosThetaZ * zp) / sqrt(rSq);
            pixelPrefactor   = prefactor * du * dv * cosDetectThetaIn / rSq;
            if ( detectDownstream < 1 ) { /* Calculate XRF flux through detector window in reflective geometry */
              effectiveCS  = matrixTotalCS[i] / cosThetaIn + matrixXRFCS[k] / cosThetaOut;
            } else {  /* Calculate XRF flux through detector window in transmission geometry */
              effectiveCS    = matrixTotalCS[i]  / cosThetaIn - matrixXRFCS[k] / cosThetaOut;
              pixelPrefactor = pixelPrefactor * exp(-1.0 * targetMassThickness * matrixXRFCS[k] / cosThetaOut);
            } 
            /* Calculate pixel dependent downstream filter transmission factor (angle change) */
            if ( mode == 4 || mode == 14 ) { secondFilterTrans = exp( -1.0 * backFilterMassThickness * backFilterCS[k] / cosThetaOut); }
            pixelXRFFlux = pixelPrefactor / effectiveCS * (1.0 - exp(-1.0 * targetMassThickness * effectiveCS));
            if (verbose > 9) { printf( "\n Detector Integration: %3d %3d %3d %8.2f %8.2f  %8.2f %8.3g %10.4g %10.4g %10.4g %10.4g",
               k, ju, jv, xp, yp, zp, cosThetaOut, cosDetectThetaIn, prefactor, firstFilterTrans, secondFilterTrans);
            }
            detectorAbs  = 1.0 - exp (-1.0 * detectMassThickness * detectXRFCS[k] / cosDetectThetaIn);
            sumXRFFlux   = sumXRFFlux   + firstFilterTrans * pixelXRFFlux * secondFilterTrans ;
            sumXRFAbsorb = sumXRFAbsorb + firstFilterTrans * pixelXRFFlux * secondFilterTrans * detectorAbs * 0.001 * xrfLineEnergy[k] / energy_keV;
          }
        }
      }
    }
    traceXRFCS[i]     = sumXRFCS;		/* Total XRF cross section */ 
    xrfFluxDensity[i] = sumXRFFluxDensity;	/* XRF angular distribution in reflective or transmission geometry */
    xrfFluxFactor[i]  = sumXRFFlux;		/* XRF signal, flux through a detector window */
    xrfPowerFactor[i] = sumXRFAbsorb; 		/* Ratio of power absorbed by the detector over power of incoming photon */
    
    /* Total XRF yield, if self absorption not included */
    totalXRFYield[i]  = traceConcentration * sumXRFCS / matrixTotalCS[i] * (1.0 - exp (-1.0 * targetMassThickness * matrixTotalCS[i] / cosThetaIn ));

    if (verbose > 2) { 
      printf( "\n %4d %12.1f %12.4g %12.4g", i, energy[i], matrixTotalCS[i], traceXRFCS[i]); 
      if ( mode == 1 || mode == 11 ) { printf( "%12.4f", totalXRFYield[i]); 
      } else if ( mode == 2 || mode == 12 ) { printf( "%16.4g", xrfFluxDensity[i]); 
      } else if ( mode == 3 ||  mode == 4 || mode == 13 ||  mode == 14 ) { printf( "%16.4g, %14.4g", xrfFluxFactor[i], xrfPowerFactor[i]); 
      }
    }
  }

  if (verbose > 2) { printf( "\n"); }
  return 0;
}


void SetupOutputFile(char *outputfile, SDDS_DATASET *SDDSout, int mode,  SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar) {
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
  if (!SDDS_DefineSimpleParameter(SDDSout, "sddsxrf_mode", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(SDDSout, "Matrix", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "MatrixFormula", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "Trace", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "xrfLines", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "xrfLineEnergy0", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "PhotonEnergy", "eV", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "MatrixTotalCS", "cm^2/g", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "TraceXRFCS", "cm^2/g", SDDS_DOUBLE) )
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if ((mode>0 && mode<5) || mode>=11) {
    if (!SDDS_DefineSimpleParameter(SDDSout, "TargetThickness", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "TargetDensity", "g/cm^3", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "TraceConcentration", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ThetaIn", "degree", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ThetaOut", "degree", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "TotalXRFYield", NULL, SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (mode==2 || mode==12 ) {
    if (!SDDS_DefineSimpleColumn(SDDSout, "xrfFluxDensity", "1/sterad", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (mode==3 || mode==4 || mode==13 || mode==14) {
    if (!SDDS_DefineSimpleParameter(SDDSout, "Detector", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "xDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "yDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "zDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "xWidth", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "yWidth", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "nx", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ny", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleColumn(SDDSout, "xrfFluxFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "xrfPowerFactor", NULL, SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==4 || mode==14) {
    if (!SDDS_DefineSimpleParameter(SDDSout, "frontFilterFormula", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "frontFilterThickness", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "frontFilterDensity", "g/cm^3", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "backFilterFormula", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "backFilterThickness", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "backFilterDensity", "g/cm^3", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

static char *matTable="/home/oxygen/OAG/generalData/elementProperties.sdds";
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
