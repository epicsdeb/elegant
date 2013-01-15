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
double energy[MAXNPTS];	  		/* Photon energy of incoming beam, units = eV */	
double thetaOut[MAXNPTS];  		/* Inclination angle of outgoing Compton photons in degrees. 0 = normal exit */	
double phiOut[MAXNPTS];  		/* Azmuth angle of outgoing Compton photons in degrees. 0 = along x-axis */
double rayleighDCS[MAXNPTS];		/* differential cross section of coherent scattering off target material, units = cm^2/g/sterad */
double totalRayleighYield[MAXNPTS];	/* Total coherently scattered photon number for one incoming photon, dimensionless. */
double rayleighFluxDensity[MAXNPTS];	/* Coherently scattered photon flux density for one incoming photons, units = 1/sterad */
double rayleighFluxFactor[MAXNPTS];	/* Number of coherently scattered photons through the detector window for one incoming photon. */
double rayleighPowerFactor[MAXNPTS];	/* Ratio of coherently scattered photon power absorbed by the detector over incoming power. */
double comptonDCS[MAXNPTS];		/* differential cross section of incoherent scattering off target material, units = cm^2/g/sterad */
double comptonEnergy[MAXNPTS];		/* Energy of the Compton-scattered x-ray photons, units = eV */
double totalComptonYield[MAXNPTS];	/* Total incoherently scattered photon number for one incoming photon, dimensionless. */
double comptonFluxDensity[MAXNPTS];	/* Incoherently scattered photon flux density  for one incoming photon, units = 1/sterad */
double comptonFluxFactor[MAXNPTS];	/* Number of incoherently scattered photons through the detector window for one incoming photon. */
double comptonPowerFactor[MAXNPTS];	/* Ratio of incoherently scattered photon power absorbed by the detector over incoming power. */
double scatterFluxFactor[MAXNPTS];	/* Number of scattered photons through the detector window for one incoming photon. */
double scatterPowerFactor[MAXNPTS];	/* Ratio of scattered photon power absorbed by the detector over incoming power. */

char   target[512]="Beryllium";		/* Common name of target (matrix) material */
char   filter[512]="Diamond";		/* Common name of the material for filter upstream of the target */
char   detector[512]="Silicon";		/* Common name of detector material */
char   detectFilter[512]="Aluminum";	/* Common name of detector filter material */
char   targetFormula[512]="Be";		/* Chemical formula of target material */
char   filterFormula[512]="C";		/* Chemical formula for upstream filter */
char   detectorFormula[512]="Si";	/* Chemical formula of detector material */
char   detectFilterFormula[512]="Al";	/* Chemical formula of detector filter material */
int    npts = 8;			/* Number of points of input arrays: energy */
long   verbose = 1;			/* Degree of details in screen output (stdout) */
long   mode = 5;			/* default mode of calculation */
int    nxWidth = 4, nyWidth = 4;	/* number of divisions of the detector surface for numerical integration */
int    detectBehindTarget=1;		/* Position of the detector: 0 = in front of target (reflective geometry), 1 = behind the target (transmission geometry) */
double polarization = 0.0;		/* Incoming beam polarization: 1 = fully polarized along x-axis, 0 = completely unpolarized */
double targetThickness = 0.25;		/* Target thickness in mm */
double filterThickness = 0.0;		/* Upstream filter thickness in mm */
double detectThickness = 0.3;		/* Detector thickness in mm */
double detectFilterThickness = 0.25;	/* Detector filter thickness in mm */
double targetDensity = 1.85;		/* Target (matrix) density in g/cm^3 */	
double filterDensity = 3.51;		/* Upstream filter density in g/cm^3 */	
double detectDensity = 2.33;		/* Detector density in g/cm^3 */	
double detectFilterDensity = 2.33;	/* Detector filter density in g/cm^3 */	
double targetThetaX = 90.0, targetThetaY = 90.0, targetThetaZ = 0.0;	/* Projection angle of target normal vector in the beam coordinates */	
double detectThetaX = 90.0,  detectThetaY = 90.0, detectThetaZ = 0.0;	/* Angle between detector surface normal and x-, y-, z-axis.  */
double phThetaX = 90.0, phThetaY = 90.0, phThetaZ = 0.0;	/* Projection angle of pinhole aperture in the beam coordinates */	
double xDet = 0.0, yDet = 0.0; 		/* Transverse position of the Detector center in mm, Oxy plane is parallel to the target front surface */	
double zDet = 100.0;			/* Longitudinal position of the Detector center in mm, measured perpendicular to the target front surface */
double xDetWidth = 10, yDetWidth = 10;	/* Detector sizes in mm, projected to Oxy plane */	
double xPH = 0.0, yPH = 0.0, zPH = 100;	/* Center position of the pinhole aperture for the detector */
double xPHWidth = 10, yPHWidth = 10;	/* Pinhole aperture sizes in mm, projected to Oxy plane */
double degToRad = 0.0174532925199433;		/* Conversion factor from degree to radian */
double pi_const=3.141592653589793238462643;	/* constant */

#define MODES 15
static int availableMode[MODES]={0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 20, 21, 22};

#define FRONT_FLAG 0x0001U
#define BACK_FLAG  0x0002U

#define SET_ENERGY 0
#define SET_MODE 1
#define SET_THETA 2
#define SET_PHI 3
#define SET_POLARIZATION 4
#define SET_TARGET 5
#define SET_DETECTOR 6
#define SET_COPY 7
#define SET_VERBOSE 8
#define SET_FILTER 9
#define SET_DETECTORFILTER 10
#define N_OPTIONS 11

#define UNDEFINED   0
#define BYCOLUMN    1
#define BYPARAMETER 2
#define BYVALUE     3
#define BYRANGE     4
#define INVALID_ANGLE 9999

static char *option[N_OPTIONS]={"energy", "mode", "theta", "phi", "polarization", "target", "detector", "copy", "verbose", "filter", "detectfilter"};

char *USAGE="sddscompton <inputFile> <outputFile> -mode=<number> -verbose=<number> -polarization=<value> \n\
                -energy=column=<colname>|begin=<value>,end=<value>,points=<integer>|specified=<value>|parameter=<name> \n\
                -theta=column=<colname>|begin=<value>,end=<value>,points=<integer>|specified=<value>|parameter=<name> \n\
                -phi=column=<colname>|begin=<value>,end=<value>,points=<integer>|specified=<value>|parameter=<name> \n\
                -target=material=<string>,formula=<string>,thickness=<value>,density=<value>,thetax=<value>,thetaz=<value> \n\
                -filter=material=<string>,formula=<string>,thickness=<value>,density=<value> \n\
                -detector=material=<string>,formula=<string>,thickness=<value>,density=<value>,\n\
                    xcenter=<value>,ycenter=<value>,zcenter=<value>,xwidth=<value>,ywidth=<value>, \n\
                    thetax=<value>,thetaz=<value>,nx=<integer>,ny=<interger>,[back|front] \n\
                -detectorfilter=material=<string>,formula=<string>,thickness=<value>,density=<value> \n\
                -copy=column|par,<list of names> \n\
<inputFile>     optional, if provided, the energy, theta, or phi values will be obtained from input file. \n\
<outputFile>    required, filename for data output.\n\
mode            mode of x-ray scattering calculation:\n\
                0 -  Cross sections of a target using its formula \n\
                1 -  Compton scattering yield in a compound target, no self absorption \n\
                2 -  Scattered photon flux density from a compound target in reflective or transmission geometry \n\
                3 -  Detector signals from a compound target in reflective or transmission geometry \n\
                4 -  Filtered detector signals from a double layer compound target in reflective or transmission geometry \n\
                5 -  (Future development) Compton detector signal from a double layer compound target with a downstream aperture \n\
                10 - Compton cross sections of an element using its symbol \n\
                11 - Total Compton yield in an elemental target, no self absorption \n\
                12 - Scattered photon flux density from an elemental target in reflective or transmission geometry \n\
                13 - Detector signals from an elemental target in reflective or transmission geometry \n\
                14 - Filtered detector signals from a double layer elemental target in reflective or transmission geometry \n\
                15 - (Future development) Compton detector signal from a double layer elemental target with a downstream aperture. \n\
                20 - Klein-Nishina Compton cross section with the target's total electron density. \n\
                21 - Klein-Nishina Compton scattering yield in a compound target, no self absorption \n\
                22 - Klein-Nishina Compton photon flux density from a compound target in reflective or transmission geometry \n\
energy          unites in eV, energy values either provided by input file from given column (parameter) or provided by a range (specified).\n\
                  -energy=column=<colname>   Value of energy is obtained from given column in inputfile \n\
                  ‑energy=begin=xx,end=yy,points=npts energy array begin/end values and number of points \n\
                  -energy=parameter=<parametername>   Value of energy is obtained from given parameter in inputfile \n\
                  ‑energy=specified=<value> energy value specified in command line \n\
                only one of the four options should be provided.\n\
theta           units in degree, the sub-options are similar to energy, for an array of theta values.\n\
                  -theta=column=<colname>   Value of theta is obtained from given column in inputfile \n\
                  ‑theta=begin=xx,end=yy,points=npts theta array begin/end values and number of points \n\
                  -theta=parameter=<parametername>   Value of theta is obtained from given parameter in inputfile \n\
                  ‑theta=specified=<value> theta value specified in command line \n\
                only one of the four options should be provided.\n\
phi             units in degree, the sub-options are similar to energy, for an array of phi values.\n\
                  -phi=column=<colname>   Value of phi is obtained from given column in inputfile \n\
                  ‑phi=begin=xx,end=yy,points=npts  phi array begin/end values and number of points \n\
                  -phi=parameter=<parametername>   Value of phi is obtained from given parameter in inputfile \n\
                  ‑phi=specified=<value> phi value specified in command line \n\
                only one of the four options should be provided.\n\
target          target/matrix specifications: \n\
                material = common name of the matrix material, such as Water; \n\
                formula = chemical formula of the matrix, such as H2O; \n\
                thickness = thickness of target (mm) \n\
                density = density of target (g/cm^3)  \n\
                thetax,  thetaz = angles between the target surface normal vector and coordinate axes. \n\
filter          filter specifications. \n\
                material = common name of the filter material, such as Water; \n\
                formula = chemical formula of the filter, such as H2O; \n\
                thickness = thickness of filter (mm) \n\
                density = density of filter (g/cm^3)  \n\
detctor         define detectors: \n\
                material = common name of the detector material, such as Water; \n\
                formula = chemical formula of the detctor; \n\
                thickness = thickness of detector (mm) \n\
                density = density of detector (g/cm^3) \n\
                xcenter, ycenter, zcenter = coordinates of the detector center. \n\
                xwidth, ywidth = width of the detector in dimensions approximately parallel to x and y axes \n\
                thetax,  thetaz = angles between the detector surface normal vector and coordinate axes. \n\
                \"back\" or \"front\" provides the position of the detector, behind the target -- back (transmission geometry) \n\
                or in front of target -- front (reflective geometry). \n\
detectorfilter  detector filter specifications. \n\
                material = common name of the detector filter material, such as Water; \n\
                formula = chemical formula of the detector filter, such as H2O; \n\
                thickness = thickness of detector filter (mm) \n\
                density = density of detector filter (g/cm^3)  \n\
polarization    Incoming beam polarization: 1 = fully polarized along x-axis, 0 = completely unpolarized.\n\
copy            copy parameters or columns from input file.\n\
sddscompton computes coherent and incoherent x ray scattering. However, since coherent scattering from polycrystaline material samples \n\
are often dominated by diffraction, only the incoherent (Compton) scattering results are likely to be useful in most cases. \n\
The coherent scattering results are useful for gases, liquids and amorphous solids. The program performs the following tasks \n\
1) Retrieve the coherent (Rayleigh) and incoherent (Compton) photon interaction cross sections of the target material. \n\
2) Calculate the total number of scattered photons. \n\
3) Calculate the angular distributions of the scattered photons, considering the targets’ self absorption.\n\
4) Calculate the integrate flux of scattered photons over a detector area. \n\
5) Calculate the integrate power of scattered photons absorbed by a detector. \n\n\
Program by Bingxin Yang and Hairong Shang.  ANL(This is version 1.20, "__DATE__")\n";

/**************************************************************************************
  MAIN PROGRAM

  This program use xraylib functions to calculate the following data:
  (0)  Cross sections of a target using its formula
  (1)  Compton yield in a compound target, no self absorption
  (2)  Compton flux density from a compound target in reflective or transmission geometry
  (3)  Compton detector signal from a compound target in reflective or transmission geometry
  (4)  Compton detector signal from a double layer compound target in reflective or transmission geometry
  (5)  Compton detector signal from a double layer compound target with a downstream aperture
  (10) Compton cross sections of an element using its symbol
  (11) Total Comptonyield in an elemental target, no self absorption
  (12) Compton flux density from an elemental target in reflective or transmission geometry
  (13) Compton detector signal from an elemental target in reflective or transmission geometry
  (14) Compton detector signal from a double layer elemental target in reflective or transmission geometry
  (15) Compton detector signal from a double layer elemental target with a downstream aperture
  (20) Klein-Nishina cross section of a target using its nominal electron density
  (21) Klein-Nishina Compton yield in a compound target, no self absorption
  (22) Klein-Nishina Compton flux density from a compound target in reflective or transmission geometry
  
****************************************************************************************/

void SetupOutputFile(char *outputfile, SDDS_DATASET  *SDDSout, long mode, SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar);
long get_material_property(char *label, char *property, char **str_value, long *int_value, double *double_value);
double getRayleighDCS (long mode, char *targetFormula, int Z, double polarization, float energy_keV, float thetaOut_Rad, float phiOut_Rad);
/**********************************************************
  getComptonEnergy ()
  
  This function retrieve the compton energy 
***********************************************************/
double getComptonEnergy (float energy_keV, float thetaOut_Rad);
double getComptonDCS (long mode, char *targetFormula, int Z, double polarization, float energy_keV, float thetaOut_Rad, float phiOut_Rad);
int getDetectorSignal (long mode);
int getScatterCS (long mode);
double getRayleighDCS (long mode, char *targetFormula, int Z, double polarization, float energy_keV, float thetaOut_Rad, float phiOut_Rad);
char *modeString (long mode);

int main ( int argc, char *argv[] )
{
  char  *material, *formula, *str_value, **columnMatch, **parMatch, *inputfile, *outputfile, *energyCol, *energyPar, *phiCol, *phiPar, *thetaCol, *thetaPar;
  int    i_arg, i, j,  validMode, energySpecSource, thetaSpecSource, phiSpecSource;
  long   calc_mode, energyPoints, phiPoints, thetaPoints, pages=0, rows=0, parMatches=0, colMatches=0;
  double deltaE, deltaPhi, deltaTheta, testvalue;
  double eStart, eEnd, defaultDensity, phiStart, phiEnd, thetaStart, thetaEnd, energyValue, thetaValue, phiValue, *energyInput, *phiInput, *thetaInput;
  unsigned long dummyFlags=0;
  SCANNED_ARG *s_arg;
  SDDS_DATASET SDDSin, SDDSout;
  
  material = formula = str_value = inputfile = outputfile = energyCol = energyPar = phiCol = phiPar = thetaCol = thetaPar = NULL;
  energyInput = phiInput = thetaInput = NULL;
  energySpecSource = thetaSpecSource = phiSpecSource = UNDEFINED;
  npts = energyPoints = phiPoints = thetaPoints = 0;
  eStart = eEnd = energyValue = 0;
  thetaStart = thetaEnd = thetaValue = phiStart = phiEnd = phiValue = INVALID_ANGLE + 1;
  columnMatch = parMatch = NULL;
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  
  /* If no arguments print the info message */ 
  if (argc<2) {
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
  }
  
  /* Process command line arguments: Check for syntex errors */
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {

      case SET_ENERGY:
	s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "column", SDDS_STRING, &energyCol, 1, 0,
			  "parameter", SDDS_STRING, &energyPar, 1, 0,
                          "begin", SDDS_DOUBLE, &eStart, 1, 0,
                          "end", SDDS_DOUBLE, &eEnd, 1, 0,
			  "points", SDDS_LONG, &energyPoints, 1, 0, 
			  "specified", SDDS_DOUBLE, &energyValue, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -energy syntax");
	if (energyCol) energySpecSource = 10 * energySpecSource + BYCOLUMN;
	if (energyPar) energySpecSource = 10 * energySpecSource + BYPARAMETER;
	if (energyValue>=10) energySpecSource = 10 * energySpecSource + BYVALUE;
        if (eStart>10 && eEnd>10 && energyPoints>0) energySpecSource = 10 * energySpecSource + BYRANGE;
        if (energySpecSource>=10)
	  SDDS_Bomb("Invalid range or duplicate energy specification. Use only ONE of the four specifications: column, parameter, specified or range.");
	break;

      case SET_THETA:
	s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "column", SDDS_STRING, &thetaCol, 1, 0,
			  "parameter", SDDS_STRING, &thetaPar, 1, 0,
                          "begin", SDDS_DOUBLE, &thetaStart, 1, 0,
                          "end", SDDS_DOUBLE, &thetaEnd, 1, 0,
			  "points", SDDS_LONG, &thetaPoints, 1, 0,
			  "specified", SDDS_DOUBLE, &thetaValue, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -theta syntax");
	if (thetaCol) thetaSpecSource = 10 * thetaSpecSource + BYCOLUMN;
	if (thetaPar) thetaSpecSource = 10 * thetaSpecSource + BYPARAMETER;
	if (thetaValue<INVALID_ANGLE) thetaSpecSource = 10 * thetaSpecSource + BYVALUE;
        if (thetaStart<INVALID_ANGLE && thetaEnd<INVALID_ANGLE && thetaPoints>0) thetaSpecSource = 10 * thetaSpecSource + BYRANGE;
        if (thetaSpecSource>=10)
	  SDDS_Bomb("Invalid range or duplicate theta specification. Use only ONE of the four specifications: column, parameter, specified or range.");
	break;

      case SET_PHI:
	s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "column", SDDS_STRING, &phiCol, 1, 0,
			  "parameter", SDDS_STRING, &phiPar, 1, 0,
                          "begin", SDDS_DOUBLE, &phiStart, 1, 0,
                          "end", SDDS_DOUBLE, &phiEnd, 1, 0,
			  "points", SDDS_LONG, &phiPoints, 1, 0, 
			  "specified", SDDS_DOUBLE, &phiValue, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -phi syntax");
	if (phiCol) phiSpecSource = 10 * phiSpecSource + BYCOLUMN;
	if (phiPar) phiSpecSource = 10 * phiSpecSource + BYPARAMETER;
	if (phiValue<INVALID_ANGLE) phiSpecSource = 10 * phiSpecSource + BYVALUE;
        if (phiStart<INVALID_ANGLE && phiEnd<INVALID_ANGLE && phiPoints>0) phiSpecSource = 10 * phiSpecSource + BYRANGE;
        if (phiSpecSource>=10)
	  SDDS_Bomb("Invalid range or duplicate phi specification. Use only ONE of the four specifications: column, parameter, specified or range.");
	break;

      case SET_POLARIZATION:
	if (s_arg[i_arg].n_items !=2)
	  SDDS_Bomb("Invalid -polarization syntax.");
	if (!get_double(&polarization, s_arg[i_arg].list[1]))
	  SDDS_Bomb("Invalid -polarization format.");
	if (polarization < -1.0 || polarization > 1.0)
	  SDDS_Bomb("Invalid -polarization value provided.");
	break;

      case SET_VERBOSE:
	verbose =1;
	if (s_arg[i_arg].n_items==2 && !get_long(&verbose,  s_arg[i_arg].list[1]))
	  SDDS_Bomb("Invalid verbose value provided.");
	break;

      case SET_MODE:
	if (s_arg[i_arg].n_items!=2 ||
            !get_long(&mode, s_arg[i_arg].list[1]))
          SDDS_Bomb("invalid -mode value");
	validMode = 0;
	for (i=0; i<MODES; i++) {
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

      case SET_TARGET:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	targetDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &targetThickness, 1, 0,
                          "density", SDDS_DOUBLE, &targetDensity, 1, 0,
			  "thetax", SDDS_DOUBLE, &targetThetaX, 1, 0,
			  "thetaz", SDDS_DOUBLE, &targetThetaZ, 1, 0,
			  "angle", SDDS_DOUBLE, &targetThetaZ, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -target syntax");
        if (material!=NULL) { 
          if (verbose > 3) fprintf(stdout," User provides target material name. \n");
  	  strcpy(target, material);
  	  if (formula!=NULL) {
    	    strcpy(targetFormula, formula);
  	  } else {
	    if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	       strcpy(targetFormula, material);
	    else 
	      strcpy(targetFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
          if (verbose > 3) fprintf(stdout," User provides target formula. \n");
  	  strcpy(targetFormula, formula);
	  if (get_material_property(formula, "name", &str_value, NULL, NULL)<0) 
	    strcpy(target, formula);
	  else 
	    strcpy(target, str_value );
  	} else {
          SDDS_Bomb("target formula not specified.\n");
  	}
	if (targetDensity<=0) {
	  if (get_material_property(target, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: target density can not be obtained from database, has to be provided.");
	  targetDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;

      case SET_FILTER:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	filterDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &filterThickness, 1, 0,
                          "density", SDDS_DOUBLE, &filterDensity, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -filter syntax");
        if (material!=NULL) { 
          if (verbose > 3) fprintf(stdout," User provides filter material name. \n");
  	  strcpy(filter, material);
  	  if (formula!=NULL) {
    	    strcpy(filterFormula, formula);
  	  } else {
	    if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	      strcpy(filterFormula, material);
	    else 
	      strcpy(filterFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
          if (verbose > 3) fprintf(stdout," User provides filter formula. \n");
  	  strcpy(filterFormula, formula);
	  if (get_material_property(formula, "name", &str_value, NULL, NULL)<0) 
	    strcpy(filter, formula);
	  else 
	    strcpy(filter, str_value);
	} else {
	  SDDS_Bomb("filter formula not specified.\n");
	}
	if (filterDensity<=0) {
	  if (get_material_property(filter, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: filter density can not be obtained from database, has to be provided.");
	  filterDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;

      case SET_DETECTORFILTER:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	filterDensity = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &detectFilterThickness, 1, 0,
                          "density", SDDS_DOUBLE, &detectFilterDensity, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -filter syntax");
        if (material!=NULL) { 
  	  strcpy(detectFilter, material);
          if (verbose > 3) fprintf(stdout," User provides detector filter material name. \n");
  	  if (formula!=NULL) {
    	    strcpy(detectFilterFormula, formula);
  	  } else {
	    if (get_material_property(material, "formula", &str_value, NULL, NULL)<0) 
	      strcpy(detectFilterFormula, material);
	    else 
	      strcpy(detectFilterFormula, str_value );
  	  }
	} else if (formula!=NULL) { 
  	  strcpy(detectFilterFormula, formula);
          if (verbose > 3) fprintf(stdout," User provides detector filter formula. \n");
	  if (get_material_property(formula, "name", &str_value, NULL, NULL)<0) 
	    strcpy(detectFilter, formula);
	  else 
	    strcpy(detectFilter, str_value);
	} else {
	  SDDS_Bomb("detector filter formula not specified.\n");
	}
	if (detectFilterDensity<=0) {
	  if (get_material_property(detectFilter, "density", NULL, NULL, &defaultDensity)<0) 
	    SDDS_Bomb("Error: filter density can not be obtained from database, has to be provided.");
	  detectFilterDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
	break;

      case SET_DETECTOR:
	s_arg[i_arg].n_items--;
        material = formula = str_value = NULL;
	detectDensity=0;
	dummyFlags = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "material", SDDS_STRING, &material, 1, 0,
                          "formula", SDDS_STRING, &formula, 1, 0,
                          "thickness", SDDS_DOUBLE, &detectThickness, 1, 0,
                          "density", SDDS_DOUBLE, &detectDensity, 1, 0,
			  "xcenter", SDDS_DOUBLE, &xDet, 1, 0,
			  "ycenter", SDDS_DOUBLE, &yDet, 1, 0,
			  "zcenter", SDDS_DOUBLE, &zDet, 1, 0,
			  "xwidth", SDDS_DOUBLE, &xDetWidth, 1, 0,
			  "ywidth", SDDS_DOUBLE, &yDetWidth, 1, 0,
			  "thetax", SDDS_DOUBLE, &detectThetaX, 1, 0,
			  "thetaz", SDDS_DOUBLE, &detectThetaZ, 1, 0,
			  "nx", SDDS_LONG, &nxWidth, 1, 0,
			  "ny", SDDS_LONG, &nyWidth, 1, 0,
			  "back", -1, NULL, 0, BACK_FLAG,
			  "front", -1, NULL, 0, FRONT_FLAG, 
                          NULL))
          SDDS_Bomb("invalid -detector syntax");
	if ((dummyFlags & BACK_FLAG) && (dummyFlags & FRONT_FLAG))
	  SDDS_Bomb("Invalid detector position given, has to back or front, can not be both.");
	if (dummyFlags & BACK_FLAG)
	  detectBehindTarget =1;
	if (dummyFlags & FRONT_FLAG)
	  detectBehindTarget = 0;
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
	    SDDS_Bomb("Error: detector density can not be obtained from database, has to be provided.");
	  detectDensity = defaultDensity;
	}
	if (str_value) free(str_value); str_value=NULL;
	if (material) free(material); material=NULL;
	if (formula) free(formula); formula=NULL;
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
	fprintf(stderr, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
	fflush(stderr);
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
 
  /* Confirm info about independent variables are available */
  calc_mode = mode % 10;
  if (!inputfile) {
    if (energySpecSource!=BYVALUE && energySpecSource!=BYRANGE) SDDS_Bomb("error: energy specifications are invalid or nonexisting!");
    if (calc_mode==0 || calc_mode==2) {
      if (thetaSpecSource!=BYVALUE && thetaSpecSource!=BYRANGE)
        SDDS_Bomb("error: theta specifications are invalid or nonexisting!");
      if (phiSpecSource!=BYVALUE && phiSpecSource!=BYRANGE && fabs(polarization)<=0.001 ) 
        SDDS_Bomb("error: phi specifications are invalid or nonexisting!");
    }
  } else {
    if (energySpecSource<=0) SDDS_Bomb("error: energy specifications are invalid or nonexisting!");
    if (calc_mode==0 || calc_mode==2) {
      if (thetaSpecSource<=0)    SDDS_Bomb("error: theta specifications are invalid or nonexisting!");
      if (phiSpecSource<=0 && fabs(polarization)<=0.001) SDDS_Bomb("error: phi specifications are invalid or nonexisting!");
    }
  }
    
  /* Calculate derived parameters in the geometric spec */
  testvalue = pow(cos(targetThetaX * degToRad), 2) + pow(cos(targetThetaZ * degToRad), 2);
  if (testvalue > 1.0) {
    fprintf(stdout,"\nsddscompton: Input target angles invalid, targetThetaX = %f, targetThetaZ = %f.\n", targetThetaX, targetThetaZ);
    return -1;
  }
  targetThetaY  = acos( sqrt (1.0 - testvalue) ) / degToRad;
  if (verbose) {
    /* Print minimum information about the calculation */
    fprintf(stdout,"\nsddscompton: verbose = %ld, mode = %ld: %s\n", verbose, mode, modeString(mode));
    fprintf(stdout,"polarization = %f, energySpecSource = %d, thetaSpecSource = %d, phiSpecSource = %d\n", polarization, energySpecSource, thetaSpecSource, phiSpecSource);
    fprintf(stdout,"target = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3)\n", target, targetThickness, targetDensity );
    fprintf(stdout,"targetFormula = %s, targetThetaX = %0.1f (deg), targetThetaY = %0.1f (deg), targetThetaZ = %0.1f (deg)\n", 
            targetFormula, targetThetaX, targetThetaY, targetThetaZ );
    if (calc_mode==4 || calc_mode==5) {
      fprintf(stdout,"filter = %s, formula = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3)\n", filter, filterFormula, filterThickness, filterDensity);
      fprintf(stdout,"detectFilter = %s, formula = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3)\n", 
	      detectFilter, detectFilterFormula, detectFilterThickness, detectFilterDensity);
    }
    if (calc_mode==3 || calc_mode==4 || calc_mode==5) {
      fprintf(stdout,"detector = %s, thickness = %0.3f (mm), density = %0.3f (g/cm^3), xwidth = %0.1f (mm), ywidth = %0.1f (mm)\n", 
	      detector, detectThickness, detectDensity, xDetWidth, yDetWidth );
      fprintf(stdout,"detectorFormula = %s, xcenter = %0.1f (mm), ycenter = %0.1f (mm), zcenter = %0.1f (mm), nx = %d, ny = %d\n", 
	      detectorFormula, xDet, yDet, zDet, nxWidth, nyWidth);
      fprintf(stdout,"detectThetaX = %0.2f (deg), detectThetaY = %0.2f (deg), detectThetaZ = %0.2f (deg), detectBehindTarget = %d\n", 
	      detectThetaX, detectThetaY, detectThetaZ, detectBehindTarget);
    }
    if (calc_mode==2) 
      fprintf(stdout,"detectBehindTarget = %d\n", detectBehindTarget);
    if (calc_mode==5) {
      fprintf(stdout,"pinhole: xcenter = %0.1f (mm), ycenter = %0.1f (mm), zcenter = %0.1f (mm), xPHWidth = %0.1f (mm), yPHWidth = %0.1f (mm)", 
	      xPH, yPH, zPH, xPHWidth, yPHWidth);
      fprintf(stdout,"phThetaX = %0.2f (deg), phThetaY = %0.2f (deg), phThetaZ = %0.2f (deg)\n", phThetaX, phThetaY, phThetaZ);
    }
  }
  if (energyPoints<=0) energyPoints = 1;
  if (phiPoints<=0)    phiPoints = 1;
  if (thetaPoints<=0)  thetaPoints = 1;
  if (!inputfile) {
    npts = MAX(energyPoints, MAX(phiPoints, thetaPoints));
    if (energySpecSource==BYVALUE) {
      for (i=0; i<npts; i++)
	energy[i] = energyValue;
    } else {
      deltaE = 0;
      if (npts>1) deltaE = (eEnd - eStart)/(npts - 1);
      for (i=0; i<npts; i++)
	energy[i] = eStart + deltaE * i;
      energyValue = energy[0];
    }
    if (thetaSpecSource==BYVALUE) {
      for (i=0; i<npts; i++)
	thetaOut[i] = thetaValue;
    } else {
      deltaTheta = 0.0;
      if (npts>1) deltaTheta = (thetaEnd - thetaStart)/(npts - 1);
      for (i=0; i<npts; i++) 
	thetaOut[i] = thetaStart + i * deltaTheta;
      thetaValue = thetaOut[0];
    }
    if (phiSpecSource==BYVALUE) {
      for (i=0; i<npts; i++)
	phiOut[i] = phiValue;
    } else {
      deltaPhi = 0.0;
      if (npts>1) deltaPhi = (phiEnd - phiStart)/(npts - 1);
      for (i=0; i<npts; i++) 
	phiOut[i] = phiStart + i * deltaPhi;
      phiValue = phiOut[0];
    }
    pages = 1;
  } else {
    if (!SDDS_InitializeInput(&SDDSin, inputfile))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    pages = 10000; /*assume a big number */
    /* if (!energyCol && !energyPar)
      SDDS_Bomb("Error: energy column or parameter name not provided in input file.\n"); */
  }

  /* Initialize xraylib and output file */
  XRayInit();
  SetupOutputFile(outputfile, &SDDSout, mode,  &SDDSin, colMatches, columnMatch, parMatches, parMatch);
  
  /* Calculations for total cross section: 
        0 <= mode <=  5, x-ray fluorescence cross sections from an elemental target
       10 <= mode <= 15, x-ray fluorescence cross sections from a compound target    */
  for (i=0; i<pages; i++) {
    if (inputfile) {
      if (SDDS_ReadPage(&SDDSin)<1)
	break;
      rows=SDDS_CountRowsOfInterest(&SDDSin);
      // if (!rows && (energyCol || phiCol || thetaCol)) continue;
      if (energyCol) {
	if (!(energyInput = SDDS_GetColumnInDoubles(&SDDSin, energyCol)))
	  SDDS_Bomb("Error in reading energy column value from input file.");
      } else if (energyPar) {
	if (!SDDS_GetParameterAsDouble(&SDDSin, energyPar, &energyValue))
	  SDDS_Bomb("Error in reading energy parameter value from input file.");
      } else if (energySpecSource==BYRANGE) 
        energyValue = eStart;     /* When range is specified with inputfile, only use the start value */
      if (thetaCol) {
	if (!(thetaInput = SDDS_GetColumnInDoubles(&SDDSin, thetaCol)))
	  SDDS_Bomb("Error in reading theta column value from input file.");
      } else if (thetaPar) {
	if (!SDDS_GetParameterAsDouble(&SDDSin, thetaPar, &thetaValue))
	  SDDS_Bomb("Error in reading theta parameter value from input file.");
      } else if (thetaSpecSource==BYRANGE) 
        thetaValue = thetaStart;  /* When range is specified with inputfile, only use the start value */
      if (phiCol) {
	if (!(phiInput = SDDS_GetColumnInDoubles(&SDDSin, phiCol)))
	  SDDS_Bomb("Error in reading phi column value from input file.");
      } else if (phiPar) {
	if (!SDDS_GetParameterAsDouble(&SDDSin, phiPar, &phiValue))
	  SDDS_Bomb("Error in reading phi parameter value from input file.");
      } else if (phiSpecSource==BYRANGE) 
        phiValue = phiStart;     /* When range is specified with inputfile, only use the start value */
      if (rows>0)
	npts = rows;
      if (npts>MAXNPTS) {
	fprintf(stderr, "Error execeed maximum allowed points.\n");
	exit(1);
      }
      for (j=0; j<npts; j++) {
	if (energyCol)
	  energy[j] = energyInput[j];
	else
	  energy[j] = energyValue;
	if (thetaCol)
	  thetaOut[j]=thetaInput[j];
	else
	  thetaOut[j]= thetaValue;
	if (phiCol)
	  phiOut[j] = phiInput[j];
	else 
	  phiOut[j] = phiValue;
      }
      if (energyInput) free(energyInput); energyInput = NULL;
      if (phiInput) free(phiInput); phiInput = NULL;
      if (thetaInput) free(thetaInput); thetaInput = NULL;
    }
    
    if (!SDDS_StartPage(&SDDSout, npts) ||
	!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sddscompton_mode", mode, "Polarization", polarization,
	      "Target", target, "TargetFormula", targetFormula, "TargetThickness", targetThickness, "TargetDensity", targetDensity, 
	      "TargetThetaX", targetThetaX, "TargetThetaY", targetThetaY, "TargetThetaZ", targetThetaZ, NULL) ||
	!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, energy, npts, "PhotonEnergy"))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (inputfile && (colMatches || parMatches) && !SDDS_CopyPage(&SDDSout, &SDDSin))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      
    if ( calc_mode==0 || calc_mode==1 || calc_mode==2 ) { 
      getScatterCS (mode);
      if ( calc_mode==0 || calc_mode==2 )
        if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonEnergy, npts, "ComptonEnergy") ||
  	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, thetaOut, npts, "Theta") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, phiOut, npts, "Phi") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleighDCS, npts, "RayleighDCS") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonDCS, npts, "ComptonDCS"))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (calc_mode==1 || calc_mode==2) {
	if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, totalRayleighYield, npts, "TotalRayleighYield") ||
	    !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, totalComptonYield, npts, "TotalComptonYield"))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	if (calc_mode==2) {
	  if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleighFluxDensity, npts, "RayleighFluxDensity") ||
	      !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonFluxDensity, npts, "ComptonFluxDensity") ||
	      !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "DetectBehindTarget", detectBehindTarget, NULL))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      }
    } else {
      /* mode =3, 4, 5, 13, 14, 15 */
      getDetectorSignal (mode);
      if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonEnergy, npts, "ComptonEnergy") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, totalRayleighYield, npts, "TotalRayleighYield") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, totalComptonYield, npts, "TotalComptonYield") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleighFluxFactor, npts, "RayleighFluxFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, rayleighPowerFactor, npts, "RayleighPowerFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonFluxFactor, npts, "ComptonFluxFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, comptonPowerFactor, npts, "ComptonPowerFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, scatterFluxFactor, npts, "ScatterFluxFactor") ||
	  !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, scatterPowerFactor, npts, "ScatterPowerFactor") ||
	  !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Detector", detector, 
			      "DetectorFormula", detectorFormula, "DetectorThickness", detectThickness,
			      "DetectorDensity", detectDensity, "DetectBehindTarget", detectBehindTarget,
			      "xDet", xDet, "yDet", yDet, "zDet", zDet, "xDetWidth", xDetWidth, "yDetWidth", yDetWidth, 
	                      "DetectorThetaX", detectThetaX, "DetectorThetaY", detectThetaY, "DetectorThetaZ", detectThetaZ,
			      "nx", nxWidth, "ny", nyWidth, NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if ( calc_mode==4 || calc_mode==5) { 
        if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Filter", filter, 
			      "FilterFormula", filterFormula, "FilterThickness", filterThickness,
			      "FilterDensity", filterDensity, "DetectorFilter", detectFilter, 
			      "DetectorFilterFormula", detectFilterFormula, "DetectorFilterThickness", detectFilterThickness,
			      "DetectorFilterDensity", detectFilterDensity, NULL))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if ( calc_mode==5 ) { 
        if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
			      "xPinhole", xPH, "yPinhole", yPH, "zPinhole", zPH, "xPinholeWidth", xPHWidth, "yPinholeWidth", yPHWidth,
	                      "PinholeThetaX", phThetaX, "PinholeThetaY", phThetaY, "PinholeThetaZ", phThetaZ, NULL))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (!SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (inputfile && !SDDS_Terminate(&SDDSin))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free_scanargs(&s_arg, argc);
  return 0;
}


/**********************************************************
  getRayleighDCS ()
  
  This function retrieve the compton scattering constant 
  using xraylib program/database package 
***********************************************************/
double getRayleighDCS (long mode, char *targetFormula, int Z, double polarization, float energy_keV, float thetaOut_Rad, float phiOut_Rad)
{
  double halfPi = 3.141592653589793238462643/2.0;
  double thetaOut0, value;
  
  if ( detectBehindTarget == 0 ) {
    thetaOut0 = pi_const - thetaOut_Rad;
  } else {
    thetaOut0 = thetaOut_Rad;
  }
  if ( fabs(thetaOut0) < 1.0E-12 ) { thetaOut0 = 1.0E-12; };

  if ( fabs(polarization) < 0.001 ) {
    if ( mode < 10 ) {
      value = DCS_Rayl_CP (targetFormula, energy_keV, thetaOut0);
    } else {
      value = DCS_Rayl (Z, energy_keV, thetaOut0);
    }
  } else if ( polarization >= 1.0 ) {
    if ( mode < 10 ) {
      value = DCSP_Rayl_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad);
    } else {
      value = DCSP_Rayl (Z, energy_keV, thetaOut0, phiOut_Rad);
    }
  } else if ( polarization <= -1.0 ) {
    if ( mode < 10 ) {
      value = DCSP_Rayl_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else {
      value = DCSP_Rayl (Z, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    }
  } else {
    if ( mode < 10 ) {
      value = 0.5 * (1.0 + polarization) * DCSP_Rayl_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad) +
              0.5 * (1.0 - polarization) * DCSP_Rayl_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else {
      value = 0.5 * (1.0 + polarization) * DCSP_Rayl (Z, energy_keV, thetaOut0, phiOut_Rad) +
              0.5 * (1.0 - polarization) * DCSP_Rayl (Z, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    }
  }

  if (verbose > 9) {
    fprintf(stdout, "\n getRayleighDCS: mode=%ld, targetFormula=%s, Z=%d, polarization=%f ", mode, targetFormula, Z, polarization);
    fprintf(stdout, "\n getRayleighDCS: energy_keV=%f, thetaOut_Rad=%f, phiOut=%f, rayleighDCS = %f", energy_keV, thetaOut_Rad, phiOut_Rad, value);
  }
  return value;    
}


/**********************************************************
  getComptonEnergy ()
  
  This function retrieve the compton energy 807

***********************************************************/
double getComptonEnergy (float energy_keV, float thetaOut_Rad)
{
  double value;
  /* In reflective mode, the angle theta is given from the 180 back scattering direction */
  if ( detectBehindTarget == 0 ) {
    value = ComptonEnergy(energy_keV, pi_const - thetaOut_Rad);
  } else {
    value = ComptonEnergy(energy_keV, thetaOut_Rad);		
  }
  return value;
}

/**********************************************************
  KNElectronDensity ()
  
  This function calculates the unit conversion factor for Klein-Nishina cross sections

***********************************************************/
double KNElectronDensity (char *targetFormula)
{
  int    Z;
  double value;

  /* Calculate conversion factor for Klein-Nishina. Future upgrade to use compound parser data struct compoundData */
  Z     = SymbolToAtomicNumber(targetFormula); 
  value = Z * CS_Total_CP ( targetFormula, 10.0 ) / CSb_Total_CP ( targetFormula, 10.0 );
  return value;
}

/**********************************************************
  getComptonDCS ()
  
  This function retrieve the compton scattering constant 
  using xraylib program/database package 
***********************************************************/
double getComptonDCS (long mode, char *targetFormula, int Z, double polarization, float energy_keV, float thetaOut_Rad, float phiOut_Rad)
{
  double halfPi = 3.141592653589793238462643/2.0, thetaOut0, value;
  
  /* In reflective mode, the angle theta is given from the 180 back scattering direction */
  if ( detectBehindTarget == 0 ) {
    thetaOut0 = pi_const - thetaOut_Rad;
  } else {
    thetaOut0 = thetaOut_Rad;
  }
  if (fabs(thetaOut0) < 1.0E-12 ) { thetaOut0 = 1.0E-12; };
  
  if ( fabs(polarization) < 0.001 ) {
    if ( mode < 10 ) {
      value = DCS_Compt_CP (targetFormula, energy_keV, thetaOut0);
    } else if ( mode < 20 ) {
      value = DCS_Compt (Z, energy_keV, thetaOut0);
    } else {
      value = KNElectronDensity(targetFormula) * DCS_KN (energy_keV, thetaOut0);
    }
    // fprintf(stdout, "\ngetComptonDCS: polarization = %f, thetaOut0 = %f, value = %f\n", polarization, thetaOut0, value);
  } else if ( polarization > 0.999 ) {
    // fprintf(stdout, "getComptonDCS: polarization = %f, phiOut_Rad = %f\n", polarization, phiOut_Rad);
    if ( mode < 10 ) {
      value = DCSP_Compt_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad);
    } else if ( mode < 20 ) {
      value = DCSP_Compt (Z, energy_keV, thetaOut0, phiOut_Rad);
    } else {
      value = KNElectronDensity(targetFormula) * DCSP_KN (energy_keV, thetaOut0, phiOut_Rad);
    }
  } else if ( polarization <= -0.999 ) {
    if ( mode < 10 ) {
      value = DCSP_Compt_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else if ( mode < 20 ) {
      value = DCSP_Compt (Z, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else {
      value = KNElectronDensity(targetFormula) * DCSP_KN (energy_keV, thetaOut0, phiOut_Rad + halfPi);
    }
  } else {
    if ( mode < 10 ) {
      value = 0.5 * (1.0 + polarization) * DCSP_Compt_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad) +
              0.5 * (1.0 - polarization) * DCSP_Compt_CP (targetFormula, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else if ( mode < 20 ) {
      value = 0.5 * (1.0 + polarization) * DCSP_Compt (Z, energy_keV, thetaOut0, phiOut_Rad) +
              0.5 * (1.0 - polarization) * DCSP_Compt (Z, energy_keV, thetaOut0, phiOut_Rad + halfPi);
    } else {
      value = 0.5 * (1.0 + polarization) * KNElectronDensity(targetFormula) * DCSP_KN (energy_keV, thetaOut0, phiOut_Rad) +
              0.5 * (1.0 - polarization) * KNElectronDensity(targetFormula) * DCSP_KN (energy_keV, thetaOut0, phiOut_Rad + halfPi);
    }
  }

  if (verbose > 9) {
    fprintf(stdout, "\n getComptonDCS: mode=%ld, targetFormula=%s, Z=%d, polarization=%f ", mode, targetFormula, Z, polarization);
    fprintf(stdout, "\n getComptonDCS: energy_keV=%f, thetaOut_Rad=%f, phiOut=%f, comptonDCS = %f", energy_keV, thetaOut_Rad, phiOut_Rad, value);
  }
  return value;    
}


/******************************************************************************
  getDetectorSignal ()
  
  This function calculates Compton scattering signal of a double layer target
*******************************************************************************/
int getDetectorSignal (long mode)
{
  int    i, ju, jv, Z, Z1;
  long   calc_mode;
  float  energy_keV=0.0, ComptonEnergy_keV=0.0, thetaOut_Rad, phiOut_Rad, pixelPhiOut_Rad;
  double totalCS=0, rayleighCS=0, comptonCS=0, targetCompTotalCS;
  double filterTotalCS, filterRayleighCS, filterComptonCS, filterCompTotalCS, filterRayleighDCS, filterComptonDCS;
  double detectTotalCS, detectCompTotalCS;
  double massThickness, cosThetaX, cosThetaY, cosThetaZ, cosTargetThetaIn, cosTargetThetaOut;
  double cosDetectThetaX, cosDetectThetaY, cosDetectThetaZ, sinDetectThetaY, cosDetectThetaIn, detectAbs;
  double detMassThickness, du, dv, uStart, vStart, up, vp, xp, yp, zp, R;
  double detFilterMassThickness, detectFilterRayleighCS, detectFilterCompTotalCS, detectFilterRayleighTrans, detectFilterComptonTrans;
  double pixelRayleighDCS, pixelComptonDCS, cosPixelThetaOut;
  double prefactor, prefactorRayl, prefactorComp, effRayleighCS, effComptonCS;
  double pixelRayleighFlux, pixelComptonFlux, sumRayleighFlux, sumComptonFlux, sumRayleighPower, sumComptonPower;
  double filterMassThickness, filterTrans0, targetTrans0, phTrans;
  double filterPrefactorRayl, filterPrefactorComp, filterEffRayleighCS, filterEffComptonCS;

  /* calc_mode designate the output quantities */
  calc_mode = mode % 10;

  /* Calculate beam entry angle into the target. Assume known quantities: targetThetaX, targetThetaY, targetThetaZ */
  cosThetaX  = cos (targetThetaX * degToRad);
  cosThetaY  = cos (targetThetaY * degToRad);
  cosThetaZ  = cos (targetThetaZ * degToRad);
  cosTargetThetaIn = cosThetaZ;

  if (mode > 9 ) { 
    Z  = SymbolToAtomicNumber ( targetFormula ); 
    Z1 = SymbolToAtomicNumber ( filterFormula ); 
  } else { 
    Z = Z1 = 0; 
  }
  filterMassThickness    = 0.1 * filterThickness * filterDensity;	// Upstream filter mass thickness in g/cm^2
  massThickness          = 0.1 * targetThickness * targetDensity;	// Target mass thickness in g/cm^2
  detMassThickness       = 0.1 * detectThickness * detectDensity;	// Detector mass thickness in g/cm^2
  detFilterMassThickness = 0.1 * detectFilterThickness * detectFilterDensity;	// Detector filter mass thickness in g/cm^2
  filterTotalCS          = 0;
  filterCompTotalCS      = 0;
  if ( calc_mode==3 ) { 
    filterMassThickness = 0.0;
    detFilterMassThickness = 0.0; 
  }
  /* Calculate detector surface normal vector */
  cosDetectThetaX = cos (detectThetaX * degToRad);
  cosDetectThetaZ = cos (detectThetaZ * degToRad);
  sinDetectThetaY = sqrt (cosDetectThetaX * cosDetectThetaX + cosDetectThetaZ * cosDetectThetaZ);
  cosDetectThetaY = sqrt (1.0 - sinDetectThetaY * sinDetectThetaY);

  /* Calculate scattered photon angle of exit from the target to center of detector */
  R = sqrt(xDet * xDet + yDet * yDet + zDet * zDet);
  cosTargetThetaOut = (xDet * cosThetaX  + yDet * cosThetaY + zDet * cosThetaZ ) / R;
  cosDetectThetaIn  = (xDet * cosDetectThetaX + yDet * cosDetectThetaY + zDet * cosDetectThetaZ) / R;
  thetaOut_Rad      = acos (zDet / R);	/* Outgoing photon scattering angle to detector center */
  phiOut_Rad        = atan2(yDet, xDet);		/* Outgoing photon scattering angle to detector center */
  if (verbose >3) { fprintf(stdout, "getDetectorSignal: cosTargetThetaIn = %f, cosTargetThetaOut = %f\n", cosTargetThetaIn, cosTargetThetaOut); }
  /* Integration parameters for detector */
  du = xDetWidth / nxWidth;
  dv = yDetWidth / nyWidth;
  uStart = -0.5 * (nxWidth - 1.0) * du;
  vStart = -0.5 * (nyWidth - 1.0) * dv;

  /* Print verbose mode information */
  if (verbose > 3) { fprintf(stdout, "getDetectorSignal: Z = %d, massThickness = %f\n", Z, massThickness); }
  if (verbose > 2) { 
    fprintf(stdout, "\n Index PhotonEnergy  TotalCS  rayleighCS  comptonCS RayleighYield ComptonYield RayleighFluxFactor ComptonFluxFactor RayleighPowerFactor ComptonPowerFactor");
    fprintf(stdout, "\n             (eV)   (cm^2/g)   (cm^2/g)   (cm^2/g) \n");
  } 

  /* Main loop */
  for ( i = 0; i < npts; i++ ) {
    /* Calculate quantities at detector center */
    energy_keV        = 0.001 * energy[i];	  /* Incoming photon energy */
    rayleighDCS[i]    = getRayleighDCS(mode, targetFormula, Z, polarization, energy_keV, thetaOut_Rad, phiOut_Rad);
    comptonDCS[i]     = getComptonDCS (mode, targetFormula, Z, polarization, energy_keV, thetaOut_Rad, phiOut_Rad);
    ComptonEnergy_keV = getComptonEnergy(energy_keV, thetaOut_Rad);	/* Compton-scattered x-ray energy */
    comptonEnergy[i]  = 1000.0 * ComptonEnergy_keV;
    /* mode = 0, 10, and others: Retriece absorption and patial cross section of (Rayleigh, Compton) scattering */
    if ( mode < 10 ) {
      totalCS           = CS_Total_CP ( targetFormula, energy_keV );
      rayleighCS        = CS_Rayl_CP  ( targetFormula, energy_keV );
      comptonCS         = CS_Compt_CP ( targetFormula, energy_keV );
      targetCompTotalCS = CS_Total_CP ( targetFormula, ComptonEnergy_keV );
      filterTotalCS     = CS_Total_CP ( filterFormula, energy_keV );
      filterRayleighCS  = CS_Rayl_CP  ( filterFormula, energy_keV );
      filterComptonCS   = CS_Compt_CP ( filterFormula, energy_keV );
      filterCompTotalCS = CS_Total_CP ( filterFormula, ComptonEnergy_keV );
    } else {
      totalCS           = CS_Total ( Z, energy_keV );
      rayleighCS        = CS_Rayl  ( Z, energy_keV );
      comptonCS         = CS_Compt ( Z,  energy_keV );
      targetCompTotalCS = CS_Total ( Z, ComptonEnergy_keV );
      filterTotalCS     = CS_Total ( Z1, energy_keV );
      filterRayleighCS  = CS_Rayl  ( Z1, energy_keV );
      filterComptonCS   = CS_Compt ( Z1, energy_keV );
      filterCompTotalCS = CS_Total ( Z1, ComptonEnergy_keV );
    }
    detectTotalCS           = CS_Total_CP ( detectorFormula, energy_keV );
    detectCompTotalCS       = CS_Total_CP ( detectorFormula, ComptonEnergy_keV );
    detectFilterRayleighCS  = CS_Total_CP ( detectFilterFormula, energy_keV );
    detectFilterCompTotalCS = CS_Total_CP ( detectFilterFormula, ComptonEnergy_keV );
   
    /* Total Compton yield, self absorption is not included */
    filterTrans0          = exp (-1.0 * filterMassThickness * filterTotalCS / cosTargetThetaIn );
    targetTrans0          = exp (-1.0 * massThickness * totalCS / cosTargetThetaIn );
    totalRayleighYield[i] = filterRayleighCS / filterTotalCS * (1 - filterTrans0) + filterTrans0 * rayleighCS / totalCS * (1.0 - targetTrans0);
    totalComptonYield[i]  = filterComptonCS  / filterTotalCS * (1 - filterTrans0) + filterTrans0 * comptonCS  / totalCS * (1.0 - targetTrans0);
    /* Calculate Comptonflux through a detector window, self absorption included */
    if (verbose > 9) { fprintf(stdout, "\n Detector Integration: i  ju  jv       xp       yp  cosTargetThetaOut cosDetectThetaIn prefactor"); }
    sumRayleighFlux  = 0.0;
    sumComptonFlux   = 0.0;
    sumRayleighPower = 0.0;
    sumComptonPower  = 0.0;
    for ( ju = 0; ju < nxWidth; ju++ ) { 
      up = uStart + ju * du;
      for ( jv = 0; jv < nyWidth; jv++ ) { 
        vp  = vStart + jv * dv;
        if (fabs(sinDetectThetaY) < 1e-6) {
          xp  = xDet + up ;
          yp  = yDet;
          zp  = zDet - vp;    
        } else {
          xp  = xDet + (up * cosDetectThetaZ - vp * cosDetectThetaX * cosDetectThetaY) / sinDetectThetaY ;
          yp  = yDet +  vp * sinDetectThetaY;            
          zp  = zDet - (up * cosDetectThetaX + vp * cosDetectThetaZ * cosDetectThetaY) / sinDetectThetaY ;    
        }
        R = sqrt(xp * xp + yp * yp + zp * zp);
        cosPixelThetaOut  = (xp * cosThetaX  + yp * cosThetaY + zp * cosThetaZ ) / R;
        cosDetectThetaIn  = (xp * cosDetectThetaX + yp * cosDetectThetaY + zp * cosDetectThetaZ) / R;
        /* Outgoing photon scattering angle for pixel */
 	thetaOut_Rad      = acos(zp / R);
        pixelPhiOut_Rad   = atan2(yp, xp);
        /* pixel specific differential scattering cross sections */
        pixelRayleighDCS  = getRayleighDCS (mode, targetFormula, Z,  polarization, energy_keV, thetaOut_Rad, pixelPhiOut_Rad);
        pixelComptonDCS   = getComptonDCS  (mode, targetFormula, Z,  polarization, energy_keV, thetaOut_Rad, pixelPhiOut_Rad);
        filterRayleighDCS = getRayleighDCS (mode, filterFormula, Z1, polarization, energy_keV, thetaOut_Rad, pixelPhiOut_Rad);
        filterComptonDCS  = getComptonDCS  (mode, filterFormula, Z1, polarization, energy_keV, thetaOut_Rad, pixelPhiOut_Rad);
        if (verbose > 10) { fprintf(stdout, "\n pixelRayleighDCS = %f, pixelComptonDCS = %f, cosTargetThetaIn = %f", pixelRayleighDCS, pixelComptonDCS, cosTargetThetaIn); }

        /* pixel specific scattering absoption coefficient for Compton-scattered photons, including filter attenuation of incoming beam */
        prefactor           = du * dv * cosDetectThetaIn / (R * R) / cosTargetThetaIn ;
        filterPrefactorRayl = prefactor * filterRayleighDCS;
        filterPrefactorComp = prefactor * filterComptonDCS;
        prefactorRayl       = prefactor * pixelRayleighDCS * filterTrans0;
        prefactorComp       = prefactor * pixelComptonDCS  * filterTrans0;
        if (verbose > 10) { fprintf(stdout, "\n prefactorRayl = %f, prefactorComp = %f", prefactorRayl, prefactorComp); }
        if (detectBehindTarget < 1) { /* Detector in front of target */ 
          /* Set up self absorption terms for reflective geometry */
          effRayleighCS       = totalCS       / cosTargetThetaIn + totalCS           / cosPixelThetaOut;
          effComptonCS        = totalCS       / cosTargetThetaIn + targetCompTotalCS / cosPixelThetaOut;
          filterEffRayleighCS = filterTotalCS / cosTargetThetaIn + filterTotalCS     / cosPixelThetaOut;
          filterEffComptonCS  = filterTotalCS / cosTargetThetaIn + filterCompTotalCS / cosPixelThetaOut;
          /* Filter attenuates scattered beam from the target */
          prefactorRayl       = prefactorRayl * exp(-1.0 * filterMassThickness * filterTotalCS     / cosPixelThetaOut);
          prefactorComp       = prefactorComp * exp(-1.0 * filterMassThickness * filterCompTotalCS / cosPixelThetaOut);
        } else { /* Detector behind target */
          /* Set up self absorption terms for transmission geometry */
          prefactorRayl       = prefactorRayl       * exp(-1.0 * massThickness       * totalCS           / cosPixelThetaOut); 
          prefactorComp       = prefactorComp       * exp(-1.0 * massThickness       * targetCompTotalCS / cosPixelThetaOut); 
          filterPrefactorRayl = filterPrefactorRayl * exp(-1.0 * filterMassThickness * filterTotalCS     / cosPixelThetaOut);
          filterPrefactorComp = filterPrefactorComp * exp(-1.0 * filterMassThickness * filterCompTotalCS / cosPixelThetaOut);
          effRayleighCS       = totalCS        / cosTargetThetaIn - totalCS           / cosPixelThetaOut;
          effComptonCS        = totalCS        / cosTargetThetaIn - targetCompTotalCS / cosPixelThetaOut;
          filterEffRayleighCS = filterTotalCS  / cosTargetThetaIn - filterTotalCS     / cosPixelThetaOut;
          filterEffComptonCS  = filterTotalCS  / cosTargetThetaIn - filterCompTotalCS / cosPixelThetaOut;
          /* Target attenuates scattered beam from the filter */
          filterPrefactorRayl = filterPrefactorRayl * exp(-1.0 * massThickness * totalCS           / cosPixelThetaOut);
          filterPrefactorRayl = filterPrefactorRayl * exp(-1.0 * massThickness * targetCompTotalCS / cosPixelThetaOut);
        } 
        /* Rayleigh flux for one pixel in the filter */
        if ( fabs(filterEffRayleighCS * filterMassThickness) < 0.001 ) {
          pixelRayleighFlux = filterPrefactorRayl * filterMassThickness;
        } else { 
          pixelRayleighFlux = filterPrefactorRayl / filterEffRayleighCS * (1.0 - exp(-1.0 * filterMassThickness * filterEffRayleighCS));
        }
        /* Compton flux for one pixel in the filter */
        if ( fabs(filterEffComptonCS * filterMassThickness) < 0.001 ) {
          pixelComptonFlux  = filterPrefactorComp * filterMassThickness;
        } else { 
          pixelComptonFlux  = filterPrefactorComp / filterEffComptonCS  * (1.0 - exp(-1.0 * filterMassThickness * filterEffComptonCS));
        }
        /* Add Rayleigh flux for one pixel in the target */
        if ( fabs(effRayleighCS * massThickness) < 0.001 ) {
          pixelRayleighFlux = pixelRayleighFlux + prefactorRayl * massThickness;
        } else {
          pixelRayleighFlux = pixelRayleighFlux + prefactorRayl / effRayleighCS * (1.0 - exp(-1.0 * massThickness * effRayleighCS));
        }
        /* Add Compton flux for one pixel in the target */
        if ( fabs(effComptonCS * massThickness) < 0.001 ) {
          pixelComptonFlux  = pixelComptonFlux  + prefactorComp * massThickness;
        } else {
          pixelComptonFlux  = pixelComptonFlux  + prefactorComp / effComptonCS  * (1.0 - exp(-1.0 * massThickness * effComptonCS));
        }
        if (verbose > 9) { 
          fprintf(stdout, "\n Detector Integration: %4d%4d%4d%9.2f%9.2f%9.3g%11.4g%11.4g%11.4g", 
                  i,ju,jv,xp,yp,cosPixelThetaOut,cosDetectThetaIn,prefactorRayl,prefactorComp);
        }
        
        phTrans = 1.0;
        if (calc_mode==5) {
          /* Calculate pinhole aperture transmission */
          phTrans = 1.0;
        }
        detectFilterRayleighTrans = exp (-1.0 * detFilterMassThickness * detectFilterRayleighCS  / cosDetectThetaIn);
        detectFilterComptonTrans  = exp (-1.0 * detFilterMassThickness * detectFilterCompTotalCS / cosDetectThetaIn);
        sumRayleighFlux   = sumRayleighFlux + pixelRayleighFlux * phTrans * detectFilterRayleighTrans;
        sumComptonFlux    = sumComptonFlux  + pixelComptonFlux  * phTrans * detectFilterComptonTrans;
        detectAbs         = 1.0 - exp (-1.0 * detMassThickness * detectTotalCS / cosDetectThetaIn);
        sumRayleighPower  = sumRayleighPower + pixelRayleighFlux * detectAbs * phTrans * detectFilterRayleighTrans;
        detectAbs         = 1.0 - exp (-1.0 * detMassThickness * detectCompTotalCS / cosDetectThetaIn);
        ComptonEnergy_keV = getComptonEnergy(energy_keV, thetaOut_Rad);	/* Compton-scattered x-ray energy */
        sumComptonPower  = sumComptonPower + pixelComptonFlux * phTrans  * detectFilterComptonTrans * detectAbs * ComptonEnergy_keV / energy_keV;
      }
    }
    /* Record the integration results */
    rayleighFluxFactor[i]  = sumRayleighFlux;	/* Rayleight signal flux, through a detector window */
    rayleighPowerFactor[i] = sumRayleighPower;	/* Ratio of power absorbed by the detector over power of incoming photon */
    comptonFluxFactor[i]   = sumComptonFlux;	/* Comptonsignal, flux through a detector window */
    comptonPowerFactor[i]  = sumComptonPower; 	/* Ratio of power absorbed by the detector over power of incoming photon */\
    scatterFluxFactor[i]   = sumComptonFlux  + sumRayleighFlux;		/* Scatter signal, flux through a detector window */
    scatterPowerFactor[i]  = sumComptonPower + sumRayleighPower; 	/* Ratio of power absorbed by the detector over power of incoming photon */\
    
    if (verbose > 2) { 
      fprintf(stdout, " %4d %12.1f %10.4g %9.4g %10.4g %12.4g %12.4g", i, energy[i], totalCS, rayleighCS, comptonCS, totalRayleighYield[i], totalComptonYield[i]); 
      fprintf(stdout, "%16.4g %18.4g %18.4g %18.4g\n", rayleighFluxFactor[i], comptonFluxFactor[i], rayleighPowerFactor[i], comptonPowerFactor[i]);
    }
  }
  if (verbose > 2) { fprintf(stdout, "\n"); }
  return 0;
}


/**********************************************************
  getScatterCS ()
  
  This function retrieve the x-ray interaction constants 
  of an element using xraylib program/database package 
***********************************************************/
int getScatterCS (long mode)
{
  int    i, Z=0;
  long   calc_mode;
  float  energy_keV=0.0, ComptonEnergy_keV=0.0, thetaOut_Rad, phiOut_Rad;
  double totalCS, rayleighCS, comptonCS, targetComptonCS;
  double massThickness, cosThetaX, cosThetaY, cosThetaZ, cosTargetThetaIn, cosTargetThetaOut;
  double prefactorRayl, prefactorComp, effRayleighCS, effComptonCS;
  
  /* calc_mode designate the output quantities */
  calc_mode    = mode % 10;
  if (mode > 9 ) Z = SymbolToAtomicNumber  ( targetFormula ); 

  /* Calculate beam entry angle into the target. Assume known quantities: targetThetaX, targetThetaY, targetThetaZ */
  cosThetaX  = cos (targetThetaX * degToRad);
  cosThetaY  = cos (targetThetaY * degToRad);
  cosThetaZ  = cos (targetThetaZ * degToRad);
  cosTargetThetaIn = cosThetaZ;
  massThickness = 0.1 * targetThickness * targetDensity;	// Target mass thickness in g/cm^2

  /* Calculate scattered photon angle of exit from the target */
  cosTargetThetaOut  = (xDet*cosThetaX + yDet*cosThetaY + zDet*cosThetaZ ) / sqrt(xDet*xDet + yDet*yDet + zDet*zDet);
  if (verbose > 3) { fprintf(stdout, "getScatterCS: cosTargetThetaIn = %f, cosTargetThetaOut = %f\n", cosTargetThetaIn, cosTargetThetaOut); }
  /* Print verbose mode information */
  if (verbose > 3) { fprintf(stdout, "getScatterCS: Z = %d, massThickness = %f\n", Z, massThickness); }
  if (verbose > 2) { 
    fprintf(stdout, "\n Index PhotonEnergy thetaOut  phiOut comptonEnergy TotalCS  rayleighCS  rayleighDCS   comptonCS   comptonDCS ");
    if ( calc_mode==1 ) { fprintf(stdout," RayleighYield  ComptonYield"); }
    if ( calc_mode==2 ) { fprintf(stdout,"RayleighFluxDensity ComptonFluxDensity"); }
    fprintf(stdout, "\n             (eV)     (deg)   (deg)       (eV)    (cm^2/g)   (cm^2/g) (cm^2/g/sterad)  (cm^2/g) (cm^2/g/sterad)");
    if ( calc_mode==2 ) { fprintf(stdout,"   1/sterad           1/sterad"); }
  } 

  /* Main loop */
  for ( i = 0; i < npts; i++ ) {
    energy_keV        = 0.001 * energy[i];				/* Incoming photon energy */
    thetaOut_Rad      = thetaOut[i] * degToRad;				/* Outgoing photon scattering angle */
    phiOut_Rad        = phiOut[i]   * degToRad;				/* Outgoing photon scattering angle */
    ComptonEnergy_keV = getComptonEnergy(energy_keV, thetaOut_Rad);	/* Compton-scattered x-ray energy */
    comptonEnergy[i]  = 1000.0 * ComptonEnergy_keV;
    
    /* mode = 0, 10, and others: Retrieve absorption and patial cross section of (Rayleigh, Compton) scattering */
    if ( mode < 10 ) {
      totalCS         = CS_Total_CP ( targetFormula, energy_keV );
      rayleighCS      = CS_Rayl_CP  ( targetFormula, energy_keV );
      comptonCS       = CS_Compt_CP ( targetFormula, energy_keV );
      targetComptonCS = CS_Total_CP ( targetFormula, ComptonEnergy_keV );
    } else {
      totalCS         = CS_Total ( Z, energy_keV );
      rayleighCS      = CS_Rayl  ( Z, energy_keV );
      targetComptonCS = CS_Total ( Z, ComptonEnergy_keV );
      if ( mode < 20 )
        comptonCS     = CS_Compt ( Z, energy_keV );
      else
        comptonCS     = KNElectronDensity(targetFormula) * CS_KN ( energy_keV );
    }
    if (calc_mode!=1) {
      rayleighDCS[i]  = getRayleighDCS(mode, targetFormula, Z, polarization, energy_keV, thetaOut_Rad, phiOut_Rad);
      comptonDCS[i]   = getComptonDCS (mode, targetFormula, Z, polarization, energy_keV, thetaOut_Rad, phiOut_Rad);
    }
    
    /* mode = 1 or 11: Total Compton yield, self absorption is not included */
    totalRayleighYield[i] = rayleighCS / totalCS  * (1.0 - exp (-1.0 * massThickness * totalCS / cosTargetThetaIn ));
    totalComptonYield[i]  = comptonCS  / totalCS  * (1.0 - exp (-1.0 * massThickness * totalCS / cosTargetThetaIn ));
    
    /* Calculate Compton flux density, self absorption included */
    if ( calc_mode==2 ) {        
      prefactorRayl = rayleighDCS[i] / cosTargetThetaIn;
      prefactorComp = comptonDCS[i]  / cosTargetThetaIn;
      /* Calculate Rayleigh angular distribution with self absorption: single target */
      if (detectBehindTarget < 1) { /* Detector in front of a single target */
        effRayleighCS = totalCS / cosTargetThetaIn + totalCS         / cosTargetThetaOut;
        effComptonCS  = totalCS / cosTargetThetaIn + targetComptonCS / cosTargetThetaOut;
      } else { /* Detector behind a single target */
        prefactorRayl = prefactorRayl * exp(-1.0 * totalCS         * massThickness / cosTargetThetaOut);
        prefactorComp = prefactorComp * exp(-1.0 * targetComptonCS * massThickness / cosTargetThetaOut);
        effRayleighCS = totalCS / cosTargetThetaIn - totalCS         / cosTargetThetaOut;
        effComptonCS  = totalCS / cosTargetThetaIn - targetComptonCS / cosTargetThetaOut;
      }
      if (verbose > 3) 
        fprintf(stdout, "\n effRayleighCS = %f, massThickness = %f, fabs(effRayleighCS * massThickness) = %g", 
                            effRayleighCS, massThickness, fabs(effRayleighCS * massThickness));
      if ( fabs(effRayleighCS * massThickness) < 0.001 ) {
        rayleighFluxDensity[i] = prefactorRayl * massThickness;
      } else {
        rayleighFluxDensity[i] = prefactorRayl / effRayleighCS * (1.0 - exp(-1.0 * massThickness * effRayleighCS));
      }
      if ( fabs(effComptonCS * massThickness) < 0.001 ) {
        comptonFluxDensity[i]  = prefactorComp * massThickness;
      } else {
        comptonFluxDensity[i]  = prefactorComp / effComptonCS  * (1.0 - exp(-1.0 * massThickness * effComptonCS));
      }
      if (verbose > 5) 
        fprintf(stdout, "\ngetScatterCS: prefactorRayl = %f, effRayleighCS = %f, prefactorComp = %f, effComptonCS = %f", 
                           prefactorRayl, effRayleighCS, prefactorComp, effComptonCS);
    }
    
    if (verbose > 2) { 
      fprintf(stdout, "\n %4d %12.1f %8.3f %8.3f %10.1f %10.4g", i, energy[i], thetaOut[i], phiOut[i], comptonEnergy[i], totalCS);
      fprintf(stdout, "%10.4g %12.4g %12.4g %12.4g", rayleighCS, rayleighDCS[i], comptonCS, comptonDCS[i]); 
      if ( calc_mode==1 ) { fprintf(stdout, "%14.4g %14.4g", totalRayleighYield[i], totalComptonYield[i]); }
      if ( calc_mode==2 ) { fprintf(stdout, "%16.4g %18.4g", rayleighFluxDensity[i], comptonFluxDensity[i]); }
    }
  }

  if (verbose > 2) { fprintf(stdout, "\n"); }
  return 0;
}


/**********************************************************
  modeString ()
  
  This function returns a message string to explain program action in this mode
***********************************************************/
char *modeString (long mode)
{
  char *msg;
  
  /* Use constant string to store message */
  msg = "unknown mode.";
  if (mode==0)         { msg = "Cross sections of a target using its formula";
  } else if (mode==1)  { msg = "Compton yield in a compound target, no self absorption";
  } else if (mode==2)  { msg = "Compton flux density from a compound target";
  } else if (mode==3)  { msg = "Compton detector signal from a compound target";
  } else if (mode==4)  { msg = "Compton detector signal from a two-layer compound target: filter-target";
  } else if (mode==5)  { msg = "Compton detector signal from a two-layer compound target with pinhole aperture";
  } else if (mode==10) { msg = "Compton cross sections of an element using its symbol";
  } else if (mode==11) { msg = "Total Comptonyield in an elemental target, no self absorption";
  } else if (mode==12) { msg = "Compton flux density from an elemental target";
  } else if (mode==13) { msg = "Compton detector signal from an elemental target";
  } else if (mode==14) { msg = "Compton detector signal from a two-layer elemental target: filter-target";
  } else if (mode==15) { msg = "Compton detector signal from a two-layer elemental target with pinhole aperture";
  } else if (mode==20) { msg = "Klein-Nishina cross sections of an elemental target using its nominal electron density";
  } else if (mode==21) { msg = "Klein-Nishina Compton yield in an elemental target, no self absorption";
  } else if (mode==22) { msg = "Klein-Nishina Compton flux density from an elemental target";
  }

  return msg;
}

void SetupOutputFile(char *outputfile, SDDS_DATASET *SDDSout, long mode,  SDDS_DATASET *SDDSin, long copyCols, char **copyCol, long copyPars, char **copyPar) {
  long i, cols=0, pars=0, calc_mode;
  char **column=NULL, **par=NULL;

  /* calc_mode designate the output quantities */
  calc_mode    = mode % 10;

  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
 
  if (!SDDS_DefineSimpleParameter(SDDSout, "sddscompton_mode", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleParameter(SDDSout, "Polarization", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "Target", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetFormula", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetThickness", "mm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetDensity", "g/cm^3", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetThetaX", "degree", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetThetaY", "degree", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "TargetThetaZ", "degree", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "ThetaIn", "degree", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "PhotonEnergy", "eV", SDDS_DOUBLE))
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (calc_mode==0 || calc_mode==1 || calc_mode==2 ) {
    if (calc_mode==0 || calc_mode==2 )
      if (!SDDS_DefineSimpleColumn(SDDSout, "ComptonEnergy", "eV", SDDS_DOUBLE)   ||
          !SDDS_DefineSimpleColumn(SDDSout, "Theta", "degree", SDDS_DOUBLE)    ||
	  !SDDS_DefineSimpleColumn(SDDSout, "Phi", "degree", SDDS_DOUBLE)      ||
	  !SDDS_DefineSimpleColumn(SDDSout, "RayleighDCS", "cm^2/g", SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleColumn(SDDSout, "ComptonDCS", "cm^2/g", SDDS_DOUBLE))
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (calc_mode==1 || calc_mode==2 ) {
      if (!SDDS_DefineSimpleColumn(SDDSout, "TotalRayleighYield", NULL, SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleColumn(SDDSout, "TotalComptonYield", NULL, SDDS_DOUBLE))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (calc_mode==2) {
	if (!SDDS_DefineSimpleColumn(SDDSout, "RayleighFluxDensity", "1/stera?", SDDS_DOUBLE) ||
	    !SDDS_DefineSimpleColumn(SDDSout, "ComptonFluxDensity", "1/stera", SDDS_DOUBLE) ||
	    !SDDS_DefineSimpleParameter(SDDSout, "DetectBehindTarget", NULL, SDDS_LONG))
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  } else {
    /* for mode=3, 4, 5, 13, 14 ,15 */
    if (!SDDS_DefineSimpleColumn(SDDSout, "ComptonEnergy", "eV", SDDS_DOUBLE)   ||
        !SDDS_DefineSimpleColumn(SDDSout, "TotalRayleighYield", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "TotalComptonYield", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "RayleighFluxFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "RayleighPowerFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "ComptonFluxFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "ComptonPowerFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "ScatterFluxFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(SDDSout, "ScatterPowerFactor", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "Detector", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorFormula", NULL, SDDS_STRING) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorThickness", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorDensity", "g/cm^3", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectBehindTarget", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorThetaX", "degree", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorThetaY", "degree", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "DetectorThetaZ", "degree", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "xDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "yDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "zDet", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "xDetWidth", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "yDetWidth", "mm", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDSout, "nx", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleParameter(SDDSout, "ny", NULL, SDDS_LONG))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if ( calc_mode==4 || calc_mode==5 ) { 
      if (!SDDS_DefineSimpleParameter(SDDSout, "Filter", NULL, SDDS_STRING) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "FilterFormula", NULL, SDDS_STRING) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "FilterThickness", "mm", SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "FilterDensity", "g/cm^3", SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "DetectorFilter", NULL, SDDS_STRING) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "DetectorFilterFormula", NULL, SDDS_STRING) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "DetectorFilterThickness", "mm", SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleParameter(SDDSout, "DetectorFilterDensity", "g/cm^3", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if ( calc_mode==5 ) { 
      if (!SDDS_DefineSimpleParameter(SDDSout, "xPinhole", "mm", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "yPinhole", "mm", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "zPinhole", "mm", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "xPinholeWidth", "mm", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "yPinholeWidth", "mm", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "PinholeThetaX", "degree", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "PinholeThetaY", "degree", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(SDDSout, "PinholeThetaZ", "degree", SDDS_DOUBLE))
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
   if (copyCols) {
    column = getMatchingSDDSNames(SDDSin, copyCol, copyCols, &cols, SDDS_MATCH_COLUMN);
    SDDS_SetColumnFlags(SDDSin, 0);
    for (i=0; i<cols; i++) {
      if (SDDS_GetColumnIndex(SDDSout, column[i])<0 && !SDDS_TransferColumnDefinition(SDDSout, SDDSin, column[i], NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetColumnsOfInterest(SDDSin, SDDS_MATCH_STRING, column[i], SDDS_OR))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_FreeStringArray(column, cols);
  }
  if (copyPars) {
    par = getMatchingSDDSNames(SDDSin, copyPar, copyPars, &pars, SDDS_MATCH_PARAMETER);
    for (i=0; i<pars; i++) {
      if (SDDS_GetParameterIndex(SDDSout, par[i])<0 && !SDDS_TransferParameterDefinition(SDDSout, SDDSin, par[i], NULL))
	SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    SDDS_FreeStringArray(par, pars);
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
  long rows=0,  index, col_index=-1;
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
  /* if (index<0) {
    fprintf(stdout, "%s not found.\n", label);
    return -1; 
  } */
  switch (col_index) {
  case FORMULA_COL:
    if (index>=0)
      SDDS_CopyString(str_value, formula[index]);
    else
      strcpy(*str_value, label);
    break;
  case NAME_COL:
    if (index>=0)
      SDDS_CopyString(str_value, name[index]);
    else
      strcpy(*str_value, label);
    break;
  case Z_COL:
  case GROUP_COL:
  case PERIOD_COL:
    if (index>=0) {
      intValue = (int32_t*)SDDS_GetColumn(&table, table_column[col_index]);
      *int_value = intValue[index];
      free(intValue);
    }
    break;
  default:
    if (index>=0) {
      doubleValue = (double*)SDDS_GetColumn(&table, table_column[col_index]);
      *double_value = doubleValue[index];
      free(doubleValue);
    }
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


