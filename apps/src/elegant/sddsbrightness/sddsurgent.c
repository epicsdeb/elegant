/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsurgent.c
 * transfer from urgent.f by R.P. Walker and B. Divacco, Sincrotron Trieste
 * URGENT is A program for calculating Undulator Radiation properties.
 * URGENT is designed fro the accurate and efficient calculation of the basic
 *        properties (angular, spectral, polarization, power density) of the 
 *        radiation generated in ideal plane, helical or elliptical undulators,
 *        and also the crossed-undulators scheme [Nikitin,Kim].
 * Hairong Shang, May 2005

$Log: not supported by cvs2svn $
Revision 1.19  2011/02/17 18:33:41  shang
changed the default value of nSigma to 3.

Revision 1.18  2010/04/18 02:35:25  borland
Correctly computes spatial distribution for harmonic > 1.

Revision 1.17  2010/04/17 22:08:54  borland
Now correctly computes K when harmonic number is not 1 and energy is given.

Revision 1.16  2009/05/01 21:53:39  shang
fixed a segmentation bug

Revision 1.15  2008/09/09 18:28:51  shang
fixed a type conversion problem in linux 64 bit machine

Revision 1.14  2007/03/13 17:39:23  shang
fixed the bug in computing on-axis power density for spatial distribution, modified the code
to read electron beam energy from each page of the input file and update the undulator's Kx and Ky if
they are computed from electron beam energy.

Revision 1.13  2007/02/02 19:20:38  shang
replaced the parameter name "TotalPowerDensity" by "OnAxisPowerDensity" for US output.

Revision 1.12  2007/02/01 18:31:58  shang
added checking input parameters; program prints error message and exits if checking failed.

Revision 1.11  2007/01/09 21:10:45  shang
fixed a typo (replaced spectrom by spectrum)

Revision 1.10  2007/01/09 21:03:31  shang
corrected the usage for calculation method.

Revision 1.9  2007/01/09 20:55:46  shang
modified to be able to calculation method and mode through description text.

Revision 1.8  2006/10/23 19:49:43  soliday
Updated to fix an issue with linux-x86_64

Revision 1.7  2006/08/31 23:14:07  borland
Fixed bug in last change: confused kx and ky when user specifies the
photon energy.  Clarified usage message.

Revision 1.6  2006/08/31 17:45:46  borland
Removed printout statements.

Revision 1.5  2006/08/30 22:37:35  borland
Now computes the energy range if none is given.
Now also computes the K value for the undulator if a photon energy is given for
the -undulator option.

Revision 1.4  2006/08/24 19:20:06  soliday
Updated so that it would compile on WIN32 again.

Revision 1.3  2006/07/31 23:49:26  jiaox
Removed the -1 option for mode and method. Increased the NXP and NYP limits from 50 to 500. Fixed the bug of pipe option. Clean up the code to asuure the proper the mode/method combinations.

*/


#include "mdb.h"
#include "scan.h"
#include "oagphy.h"

#define SET_CALCULATION 0
#define SET_UNDULATOR 1
#define SET_ELECTRON_BEAM 2
#define SET_PINHOLE 3
#define SET_ALPHA 4
#define SET_OMEGA 5
#define SET_US 6
#define SET_PHOTONENERGY 7
#define SET_NPHI 8
#define SET_PIPE 9
#define SET_NOWARNINGS 10
#define SET_COUPLING 11
#define SET_EMITTANCERATIO 12
#define N_OPTIONS 13
char *option[N_OPTIONS] = {
  "calculation", "undulator", "electronbeam","pinhole", "alpha", "omega",
  "us", "photonenergy", "nphi", "pipe", "nowarnings", "coupling", "emittanceratio",
};

typedef struct {
  long nPeriod, itype;
  double period, kx, ky, phase, energy;
} UNDULATOR_PARAM;

typedef struct {
  long nsigma;
  double energy, current, sigmax, sigmay, sigmaxp, sigmayp, energySpread;
} ELECTRON_BEAM_PARAM;

typedef struct {
  double xPC, yPC, xPS, yPS, distance;
  long nXP, nYP;
} PINHOLE_PARAM;


#define UNDULATOR_KX_GIVEN 0x0001U
#define UNDULATOR_KY_GIVEN 0x0002U
#define UNDULATOR_PERIOD_GIVEN 0x0004U
#define UNDULATOR_NUMBER_OF_PERIODS_GIVEN 0x0008U
#define ELECTRON_SIGMAX_GIVEN 0x0010U
#define ELECTRON_SIGMAY_GIVEN 0x0020U
#define ELECTRON_SIGMAXP_GIVEN 0x0040U
#define ELECTRON_SIGMAYP_GIVEN 0x0080U
#define PINHOLE_XPC_GIVEN 0x0100U
#define PINHOLE_YPC_GIVEN 0x0200U
#define PINHOLE_XPS_GIVEN 0x0400U
#define PINHOLE_YPS_GIVEN 0x0800U
#define PINHOLE_NXP_GIVEN 0x1000U
#define PINHOLE_NYP_GIVEN 0x2000U
#define UNDULATOR_ENERGY_GIVEN 0x4000U

#define MAXIMUM_E 50000  /*maximum energy points */
#define MAXIMUM_H 750000  /*maximum harmonics */
/*fortran subroutine*/
void urgent_();
void us_();
#if defined(_WIN32)
#define urgent_() URGENT()
#define us_() US()
#endif

#define CLO_DEJUS_METHOD 0
#define CLO_WALKER_INFINITE_METHOD 1
#define CLO_WALKER_FINITE_METHOD 2
#define CLO_METHODS 3
static char *method_options[CLO_METHODS]={
  "dejus", "walkerinfinite", "walkerfinite"
  };

#define CLO_FLUX_DISTRIBUTION_MODE 0
#define CLO_FLUX_SPECTRUM_MODE 1
#define CLO_BRIGHTNESS_MODE 2
#define CLO_BRILLIANCE_MODE 3
#define CLO_PINHOLE_SPECTRUM_MODE 4
#define CLO_INTEGRATED_SPECTRUM_MODE 5
#define CLO_POWER_DENSITY_MODE 6
#define CLO_MODES 7
static char *mode_options[CLO_MODES]={
  "fluxDistribution", "fluxSpectrum", "brightness", "brilliance", "pinholeSpectrum", 
  "integratedSpectrum", "powerDensity"
  };

void SetupUSOutput(SDDS_DATASET *SDDSout, char *outputfile, long mode, long iang, long isub);
void DefineParameters(SDDS_DATASET *SDDSout, long mode, long iang, long us);
void WriteUSResultsToOutput(SDDS_DATASET *SDDSout, UNDULATOR_PARAM  undulator_param, 
                            ELECTRON_BEAM_PARAM electron_param, PINHOLE_PARAM pinhole_param,
                            long nPhi, long nAlpha, double dAlpha, long nOmega, double dOmega,
                            long mode, long icalc, int32_t iharm, long nEE, long isub, long iang,
                            double ptot, double pd, double ptot1, double ftot, long imax, long imin,
                            double emin, double emax, double E1, double lamda1, double *xPMM, double *yPMM,
                            double *L1, double *L2, double *L3, double *L4, double *EE, double *irradiance,
                            double *spec1);
void SpecifyDescription(char **description, long isub, long mode, long iang, long i_max, int32_t iharm);
SDDS_DATASET *SetupOutputFile(char *outputfile, SDDS_DATASET *SDDSout, long mode, long itype, long icalc, 
                              long isub, long iang, long idebug, int32_t iharm, long special, SDDS_DATASET *SDDSout1);
void WriteToOutput(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSout2, char *description, char *undulatorType,
                   long idebug,
                   long itype, double period, double kx, double ky, double phase,
                   long nPeriod, double emin, double emax, long nE,
                   double energy, double energySpread, double current, double sigx, double sigy, double sigxp,
                   double sigyp, double distance, double xPC, double yPC, 
                   double xPS, double yPS,long nXP, long nYP,
                   long mode, long icalc, int32_t iharm, long nPhi, long nSig, long nAlpha,
                   double dAlpha, long nOmega, double dOmega,
                   double E1, double lamda1, double ptot, double pd, long isub, long iang,
                   double *EE, double *lamda, long min_harmonic, long max_harmonic, long nEE,
                   double *xPMM, double *yPMM, double *irradiance, 
                   double *L1, double *L2, double *L3, double *L4, double max_irradiance, double *power,
                   long *I1, long *I2, long i_max, double *EI, double *spec1, double *spec2,
                   double *spec3, double pdtot, double ptot1, double ftot, long *harmonics,
                   SDDS_DATASET *SDDSout1, long special, int32_t iharm5, double *x5, double *y5, double *E5,
                   double *power5, double *flux5);
long GetISub(long mode, long icalc);
void check_input_parameters(UNDULATOR_PARAM *undulator_param, ELECTRON_BEAM_PARAM *electron_param, 
                            PINHOLE_PARAM *pinhole_param, long nE, long nPhi, long nAlpha, long nOmega, double dOmega,
                            long mode, long icalc, int32_t iharm, long inputSupplied, long us);

char *USAGE1="sddsurgent <inputFile> <outputFile>\n\
    [-calculation=mode={1|2|3|4|5|6|-6|fluxDistribution|fluxSpectrum|brightness|brilliance|pinholeSpectrum|integratedSpectrum|powerDensity|,method={1|2|3|4|14|dejus|walkerinfinite|walkerfinite},harmonics=<integer>] \n\
    [-undulator=period=<value>,numberOfPeriods=<integer>,kx=<value>,ky=<value>,phase=<value>,energy=<eV>] \n\
    [-electronBeam=current=<value>,energy=<value>,spread=<value>,xsigma=<value>,ysigma=<value>,xprime=<value>,yprime=<value>,nsigma=<number>] \n\
    [-pinhole=distance=<value>,xposition=<value>,yposition=<value>,xsize=<value>,ysize=<value>,xnumber=<integer>,ynumber=<integer>]\n\
    [-alpha=steps=<integer>,delta=<value>] \n\
    [-omega=steps=<integer>,delta=<value>] [-nphi] \n\
    [-us] [-photonEnergy=maximum=<value>,minimum=<value>,points=<number>]\n\
    [-nowarnings] [-coupling=<value>] [-emittanceRatio=<value>] \n\
<inputfile>   if both <inputfile> and <outputfile> are provided, \n\
              the inputfile is the same input file as sddsbrightness, it can be a \n\
              twiss file or output from sddsanalyzebeam. <inputfile> provides \n\
              electron beam size (sigmax, sigmay, sigmaxp and sigmayp). \n\
              If <outputfile> is not provided, the <inputfile> will be taken as \n\
              <outputfile>, the electron beam size should be provided through \n\
              -electronBeam option. \n\
nowarnings    If provided, no warning messages will be printed.\n\
coupling       same as the coupling in sddsbrightness;\n\
              x emittance is ex0/(1+coupling), while y emittance is\n\
              coupling*ex0/(1+coupling).\n\
              Ignored if ey0 parameter or ey column is found in file.\n";
char *USAGE2="emittanceRatio   ratio of y emittance to x emittance.  x emittance is\n\
                 ex0 from input file. y emittance is ratio*ex0\n\
                 Ignored if ey0 parameter or ey column is found in file.\n\
us            If provided, use Roger's us program for spectral calculation.\n\
              otherwise, use Walker's urgent program for spectral calculation.\n\
calculation   specifies calculation method and mode. \n\
              method: both urgent and us have following methods\n\
                1:                    Non-zero emittance; finite-N. \n\
                2:                    Non-zero emittance; infinite-N. \n\
                3 | WalkerFinite:   Zero emittance;     finite-N.   \n\
                us has additional methods as following. note that urgent method\n\
                already includes convolution.\n\
                4  | Dejus:           Non-zero emittance; infinite-N + convolution (Dejus) \n\
                14 | WalkerInfinite:    Non-zero emittance; infinite-N + convolution (Walker) \n\
                \n\
              mode: both urgent and us has following mode: \n\
                1 | fluxDistribution:        Angular/spatial flux density distribution \n\
                                             Flux distribution at the energy chosen as minimum energy.\n\
                2 | fluxSpectrum:            Angular/spatial flux density spectrum \n\
                                             Spectrum at any given point in space as selected by the X and Y\n\
                                             coordinate for the center of the pinhole. X is horizontal and Y is vertical \n\
                3 | brightness | brilliance: On-axis brilliance spectrum \n\
                4 | pinholeSpectrum:         Flux spectrum through a pinhole \n\
                                             Spectrum through a pinhole centered at X-center and Y-center with \n\
                                             size X-size and Y-size.  The energy range is from the minimum to the \n\
                                             maximum energy. \n\
                5 | integratedSpectrum:      Flux spectrum integrated over all angles \n\
                                             The pinhole parameters have no significance here. \n\
                6 | powerDensity:            Power density and integrated power \n\
                                             Integrated over all energies, thus the energy parameters have no significance here.\n\
                \n";
char *USAGE3="  urgent can have -6 mode (only valid for harmonics<=0), \n\
                which does everything mode=6 does, plus \n\
                the angular/spatial distribution of power density \n\
                for each harmonic. For mode=6(urgent), two output files, \n\
                <outputfile> and <outputfile>.total will be created;\n\
                while for mode=-6(urgent), three output files will be created: \n\
                <outputfile>, <outputfile>.total, and <outputfile>.harmonic which contains\n\
                the flux distribution for each harmonics.\n\
                              also note that urgent mode =6 the method 2 is actully \n\
                Zero emittance. So mode=6,method=2 for urgent corresponds to \n\
                mode=6,method=3 -us.\n\
                harmonics  specifies harmonic number. \n\
                 =0, include all harmonics for calculation.\n\
                 =-1, Lowest order harmonic \n\
                 =I (I>0), Ith harmonic; \n\
                 =-I (I>0), include harmonics from 1 to I.\n\
                 For urgent mode=1,2,3,4 with method=3, harmonics is not relevent,while the above \n\
                 numbers are valid in us calculation.\n\
                \n";
char *USAGE4="undulator     Specifies the undulator parameters: period in m units, number of periods, \n\
               photon energy (in eV, results in calculation of ky and forces kx=0),\n\
               horizontal-plane (vertical field) deflection parameter (ky)\n\
               vertical-plane (horizontal field) deflection parameter (kx)\n\
               and phase difference (degree) of canted undulator.\n\
               The default value of phase is 0, i.e., single undulator.\n\
electronBeam  Specifies the electron beam parameters (which can also be provided by input file): \n\
               current  electron beam current in A. (default is 0.1A). \n\
               energy   electron energy in Gev. (default is 7.0Gev).\n\
               spread   electron energy spread. \n\
               xsigma   horizontal RMS beam size (mm) \n\
               ysigma   vertical RMS beam size (mm) \n\
               xprime   horizontal RMS divergence (mrad) \n\
               yprime   vertical RMS divergence (mrad) \n\
               nsigma   no. of standard deviations of electron beam dimensions \n\
                        (size and divergence) to be included. \n\
pinhole    Specifies pinhole parameters: \n\
            distance   distances from the source (m) \n\
                       (distance=0.0  angular units) \n\
            xposition  X-coordinate for center of pinhole (mm) or (mrad) \n\
            yposition  Y-coordinate for center of pinhole (mm) or (mrad) \n\
            xsize      X-size of pinhole (full width)	(mm) or (mrad) \n\
            ysize      y-size of pinhole (full width)	(mm) or (mrad) \n\
            (for angular units (d=0.0) values are entered in mrad) \n\
            xnumber    Number of subdivisions of pinhole in X (max 50) \n\
            ynumber    Number of subdivisions of pinhole in Y (max 50) \n\
            pinhole parameters are not needed for computing on-axis brilliance (i.e., mode=3). \n";
char *USAGE5="nphi       specifies number of steps in phi between o and pi/2.0. nphi<100 \n\
            used in (calculation mode=1,2,3,4,5 calculation method=1,2) \n\
alpha      steps specifies no. of steps in angle alpha (gamma*theta) (<100). \n\
            used in mode=1, method=1.\n\
            delta specifies range of angles in alpha**2 to be used, in units of \n\
            the angular equivalent to 1/N [2.0].\n\
            used in  (mode=1, method=1) and method=3.\n\
omega      steps specifies no. of steps in photon energy for the natural lineshape. \n\
            steps < 5000. used in mode=2,3,4, method=1. \n\
            delta specifies range of photon energies to be included in the natural lineshape \n\
            in units (energy of fundamental/N) [2.0] i.e. the default value covers the range +/- \n\
            2/N of the natural lineshape. used in mode=2,3,4,5, method=1. \n\
photonEnergy  specifies the maximum and minimum photon energy in eV, \n\
               and the number of energy points to be computed.\n\n\
URGENT (converted to sddsurgent) is program for calculating undulator radiation spectra. It is designed for \n\
the accurate and efficient calculation of the basic properties (angular, spectral, polarization, power density) \n\
of the radiation generated in ideal plane, helical or elliptical undulators, and also the crossed-undulators scheme.\n\n\
sddsurgent also incorporated US program by Roger Dejusto to calculate undulator spectra within the Bessel function \n\
approximation for an ideal planar undulator or an ideal elliptical undulator (including polarization in both cases).\n\n";

int main(int argc, char **argv) {
  char *inputfile, *outputfile, *undulatorType, *description, *output, *method_str, *mode_str;
  SDDS_DATASET sddsin, sddsout, *sddsout2=NULL, sddsout1;
  unsigned long pipeFlags=0,dummyFlags=0;
  long i_arg, tmpFileUsed, i, j, special=0, mode_index=-1, nXP0, nYP0;
  long hUndulator = 1;

  SCANNED_ARG *s_arg;

  /*input parameters */
  UNDULATOR_PARAM undulator_param;
  ELECTRON_BEAM_PARAM electron_param;
  PINHOLE_PARAM pinhole_param;
  double emax,  emin, dAlpha, dOmega, coupling, emittanceRatio, period, current;
  long nE, mode, icalc,  nAlpha, nOmega, idebug, us, nPhi, points, nowarnings;
  int32_t iharm;
  /* input parameters from input file */
  TWISS_PARAMETER twiss;
  long page, inputPages=1, inputSupplied=0;
  
  /*output parameters*/
  double lamda1, E1, ptot, pd, *EE, *lamda, *xPMM, *yPMM, *irradiance;
  double *L1,*L2,*L3,*L4,*power,*EI, *spec1, *spec2;
  double *spec3, max_irradiance, totalPD, totalPower, totalFlux, pdtot, ptot1, ftot;
  long isub,iang,min_harmonic,max_harmonic,nEE,*I1,*I2,i_max, *harmonics,iharm5;
  double *E5,*flux5,*power5,*x5,*y5;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s%s%s%s%s", USAGE1, USAGE2, USAGE3, USAGE4, USAGE5);
    exit(1);
  }
  method_str = mode_str = NULL;
  output= NULL;
  lamda1=E1=ptot=pd=max_irradiance=totalPD=totalPower=totalFlux=pdtot=ptot1=ftot=0;
  EE=lamda=xPMM=yPMM=irradiance=L1=L2=L3=L4=power=EI=spec1=spec2=spec3=NULL;
  isub=iang=min_harmonic=max_harmonic=nEE=i_max=idebug=iharm5=us=0;
  I1=I2=harmonics=NULL;
  E5=flux5=power5=x5=y5=NULL;
  undulatorType=description=NULL;
  inputfile=outputfile=NULL;
  pipeFlags=0;
  tmpFileUsed=0;
  us = 0;
  nowarnings = 1;
  electron_param.current=0.1;
  electron_param.energy = 7.0;
  electron_param.energySpread = 0.0;
  electron_param.nsigma = 3;
  electron_param.sigmaxp = electron_param.sigmayp = 0.0;
  electron_param.sigmax = electron_param.sigmay =0.0;
  undulator_param.itype=1;
  undulator_param.phase=0.0;
  pinhole_param.distance=0;
  nPhi=20;
  pinhole_param.xPC=pinhole_param.yPC=pinhole_param.xPS=pinhole_param.yPS=0;
  pinhole_param.nXP=pinhole_param.nYP=0;
  pinhole_param.distance = 1;
  pinhole_param.nXP = pinhole_param.nYP = 0;

  nAlpha=15;
  dAlpha=2.0;
  nOmega=16;
  dOmega=2.0;
  nE=100;
  mode=-1;
  icalc=-1;
  iharm=0;
  emin=emax=0;
  
  coupling = emittanceRatio = 0;
  undulator_param.energy = undulator_param.kx = undulator_param.ky = 0;
  undulator_param.period = 0; 
  undulator_param.nPeriod = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_NOWARNINGS:
        nowarnings = 1;
        break;
      case SET_COUPLING:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &coupling) ||
            coupling<=0)
          SDDS_Bomb("invalid -coupling value");
        break;
      case SET_EMITTANCERATIO:
        if (s_arg[i_arg].n_items!=2 ||
            !sscanf(s_arg[i_arg].list[1], "%le", &emittanceRatio) ||
            emittanceRatio<=0)
          SDDS_Bomb("invalid -emittanceRatio value");
        break;
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_CALCULATION:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -calculation syntax");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "mode", SDDS_STRING, &mode_str, 1, 0,
                          "method", SDDS_STRING, &method_str, 1, 0,
                          "harmonics", SDDS_LONG, &iharm, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -calculation syntax");
        if (method_str) {
          if (!get_long(&icalc, method_str)) {
            switch(match_string(method_str, method_options, CLO_METHODS, 0)) {
            case CLO_DEJUS_METHOD:
              icalc = 4;
              break;
            case CLO_WALKER_INFINITE_METHOD:
              icalc = 14;
              break;
            case CLO_WALKER_FINITE_METHOD:
              icalc = 3;
              break;
            default:
              fprintf(stderr,"Unknown calculation method - %s provided.\n", method_str);
              exit(1);
              break;
            }
          }
        }
        if (mode_str) {
          if (!get_long(&mode, mode_str)) {
            if ((mode_index=match_string(mode_str, mode_options, CLO_MODES, 0))==-1) {
              fprintf(stderr, "Unknow mode - %s provided.\n", mode_str);
              exit(1);
            }
            if (mode_index<3)
              mode = mode_index+1;
            else  
              mode = mode_index;
          }
        }
        if ( (mode<1 && mode != -6 )  || mode>6)
          SDDS_Bomb("invalid mode given for calculation, it should be 1,2,3,4,5,6,-6");
        if (icalc<1 || (icalc>4 && icalc!=14))
          SDDS_Bomb("invalid method given for calculation, it should be -1, or 1,2,3,4,14");
        if (iharm<-1000 || iharm>1000)
          SDDS_Bomb("harmonics exceeds the range of which from -1000 to 1000.");
	if (iharm>0)
	  hUndulator = iharm;
        s_arg[i_arg].n_items++; 
        if (method_str) free(method_str);
        if (mode_str) free(mode_str);
        method_str = mode_str = NULL;
        break;
      case SET_UNDULATOR:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -undulator syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "period", SDDS_DOUBLE, &undulator_param.period, 1, UNDULATOR_PERIOD_GIVEN,
                          "numberOfPeriods", SDDS_LONG, &undulator_param.nPeriod, 1, UNDULATOR_NUMBER_OF_PERIODS_GIVEN,
                          "kx", SDDS_DOUBLE, &undulator_param.kx, 1, UNDULATOR_KX_GIVEN,
                          "ky", SDDS_DOUBLE, &undulator_param.ky, 1, UNDULATOR_KY_GIVEN,
                          "energy", SDDS_DOUBLE, &undulator_param.energy, 1, UNDULATOR_ENERGY_GIVEN,
                          "phase", SDDS_DOUBLE, &undulator_param.phase, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -undulator syntax");
        if (!(dummyFlags&UNDULATOR_ENERGY_GIVEN) && !(dummyFlags&UNDULATOR_KX_GIVEN) &&
            !(dummyFlags&UNDULATOR_KY_GIVEN)) {
          SDDS_Bomb("invalid -undulator syntax, give kx, ky or energy!");
        }
        if (!(dummyFlags&UNDULATOR_NUMBER_OF_PERIODS_GIVEN))
          SDDS_Bomb("invalid -undulator syntax, numberOfPeriods is not given!");
        if (undulator_param.phase>0) undulator_param.itype=2;
        s_arg[i_arg].n_items++;
        break;
      case SET_ELECTRON_BEAM:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -electronBeam syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "current", SDDS_DOUBLE, &electron_param.current, 1, 0,
                          "energy", SDDS_DOUBLE, &electron_param.energy, 1, 0,
                          "xsigma", SDDS_DOUBLE, &electron_param.sigmax, 1, ELECTRON_SIGMAX_GIVEN,
                          "xprime", SDDS_DOUBLE, &electron_param.sigmaxp, 1, ELECTRON_SIGMAXP_GIVEN,
                          "ysigma", SDDS_DOUBLE, &electron_param.sigmay, 1, ELECTRON_SIGMAY_GIVEN,
                          "yprime", SDDS_DOUBLE, &electron_param.sigmayp, 1, ELECTRON_SIGMAYP_GIVEN,
                          "spread", SDDS_DOUBLE, &electron_param.energySpread, 1, 0,
                          "nsigma", SDDS_LONG, &electron_param.nsigma, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -electronbeam syntax");
        s_arg[i_arg].n_items++;
        break; 
      case SET_PINHOLE:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -pinhole syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "distance", SDDS_DOUBLE, &pinhole_param.distance, 1, 0,
                          "xposition", SDDS_DOUBLE, &pinhole_param.xPC, 1, 0,
                          "yposition", SDDS_DOUBLE, &pinhole_param.yPC, 1, 0,
                          "xsize", SDDS_DOUBLE, &pinhole_param.xPS, 1, 0,
                          "ysize", SDDS_DOUBLE, &pinhole_param.yPS, 1, 0,
                          "xnumber", SDDS_LONG, &pinhole_param.nXP, 1, 0,
                          "ynumber", SDDS_LONG, &pinhole_param.nYP, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -pinhole syntax.");
        if (pinhole_param.nXP>500 || pinhole_param.nYP>500) 
          SDDS_Bomb("Number of intervals of x/y acceptance (nXP/nYP) should be between 0 and 500");
        s_arg[i_arg].n_items++;
        break;
      case SET_ALPHA:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -alpha syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "steps", SDDS_LONG, &nAlpha, 1, 0,
                          "delta", SDDS_DOUBLE, &dAlpha, 1, 0, NULL))
          SDDS_Bomb("invalid -alpha syntax.");
        s_arg[i_arg].n_items++;
        break;
      case SET_OMEGA:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -omega syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "steps", SDDS_LONG, &nOmega, 1, 0,
                          "delta", SDDS_DOUBLE, &dOmega, 1, 0, NULL))
          SDDS_Bomb("invalid -omega syntax.");
        s_arg[i_arg].n_items++;
        break;
      case SET_NPHI:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -nphi syntax.");
        if (!get_long(&nPhi, s_arg[i_arg].list[1]))
          SDDS_Bomb("invalid -nphi syntax.");
        break;
      case SET_US:
        us=1;
        break;
      case SET_PHOTONENERGY:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -photonenergy syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "maximum", SDDS_DOUBLE, &emax, 1, 0,
                          "minimum", SDDS_DOUBLE, &emin, 1, 0,
                          "points", SDDS_LONG, &nE, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -photonEnergy syntax");
        s_arg[i_arg].n_items++;
        break;
      default:
        fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
        fflush(stdout);
        exit(1);
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
  if (coupling && emittanceRatio)
    SDDS_Bomb("give only one of -coupling or -emittanceRatio");   
  
  if( !pipeFlags) {
      processFilenames("sddsurgent", &inputfile, &outputfile, pipeFlags, nowarnings, &tmpFileUsed); 
      if (tmpFileUsed)
          outputfile = inputfile;
  } else {
       if ( pipeFlags&USE_STDIN) {
          outputfile = inputfile;
          inputfile = NULL;
       }
       if (pipeFlags&USE_STDOUT) {
          if(inputfile)
             processFilenames("sddsurgent",&inputfile,&outputfile,pipeFlags,nowarnings, &tmpFileUsed);
       }
  }   
  
  if (emin<0 || emax<0) 
    SDDS_Bomb("The minimum or maximum photon energy is less than 0.");
  if (nE<2)
    nE = 2;
  
  SDDS_ZeroMemory(&twiss, sizeof(TWISS_PARAMETER));
  twiss.beams = 0;
  if ((inputfile && !tmpFileUsed) || (pipeFlags&USE_STDIN) ) {
    if (!(ReadTwissInput(inputfile, &twiss, coupling, emittanceRatio,
                         undulator_param.period, undulator_param.nPeriod)))
      SDDS_Bomb("Unable to read twiss parameter from input.");                     
    inputPages = twiss.beams;
    inputSupplied = 1;
    if (twiss.pCentral)
      electron_param.energy = twiss.pCentral[0]*me_mev*0.001;
    
  } else {
    if (undulator_param.energy) {
      /* calculate the Kx value for the undulator */
      double gamma, lambda, x1;
      gamma = electron_param.energy*1e3/me_mev;
      lambda = h_mks*c_mks/(undulator_param.energy*e_mks);
      undulator_param.kx = 0;
      x1 = 4*gamma*gamma*lambda/undulator_param.period;
      hUndulator = (long)(2/x1 + 0.5);
      if (hUndulator%2==0)
	hUndulator += 1;
      undulator_param.ky = sqrt(hUndulator*x1-2);
    }
    if (emin==emax && emin==0) {
      double gamma, lambda;
      gamma = electron_param.energy*1e3/me_mev;
      lambda = undulator_param.period/(2*gamma*gamma)*(1+0.5*(sqr(undulator_param.kx)+sqr(undulator_param.ky)));
      emin = emax = hUndulator*h_mks*c_mks/lambda/e_mks;
      nE = 2;
    }
  }
  check_input_parameters(&undulator_param, &electron_param, &pinhole_param, nE, nPhi, nAlpha, nOmega, dOmega,
                         mode, icalc, iharm, inputSupplied, us);
  
  
  if (undulator_param.itype==1)
    SDDS_CopyString(&undulatorType,"Undulator");
  else
    SDDS_CopyString(&undulatorType,"Canted undulator");
  
  special = 0;
  /* nXP and nYP are changed by US and urgent program
     the return value of nXP = nXP + 1, nYP = nYP +1
     so nXP, nYp is changing when page number increasing
     nXP0, nYP0 are to remember the input value of nXP and nYP
   */
  nXP0 = pinhole_param.nXP;
  nYP0 = pinhole_param.nYP;
  
  if (nE>=(pinhole_param.nXP+1)*(pinhole_param.nYP+1))
    points=nE+100;
  else
    points=(pinhole_param.nXP+1)*(pinhole_param.nYP+1)+100;
  
  EE=(double*)malloc(sizeof(*EE)*points);
  if (!us) 
    lamda=(double*)malloc(sizeof(*lamda)*points);
  xPMM=(double*)malloc(sizeof(*xPMM)*points); 
  yPMM=(double*)malloc(sizeof(*yPMM)*points); 
  irradiance=(double*)malloc(sizeof(*irradiance)*points);
  L1=(double*)malloc(sizeof(*L1)*points);
  L2=(double*)malloc(sizeof(*L2)*points);
  L3=(double*)malloc(sizeof(*L3)*points);
  L4=(double*)malloc(sizeof(*L4)*points);
  if (!us) { 
    power=(double*)malloc(sizeof(*power)*points);
    EI=(double*)malloc(sizeof(*EI)*points);
  }
  spec1=(double*)malloc(sizeof(*spec1)*points);
  if (!us) {
    spec2=(double*)malloc(sizeof(*spec2)*points);
    spec3=(double*)malloc(sizeof(*spec3)*points);
    I1=(long*)malloc(sizeof(*I1)*points);
    I2=(long*)malloc(sizeof(*I2)*points);
    harmonics=(long*)malloc(sizeof(*harmonics)*points);
  }
  if (mode==-6) {
    special=1;
    E5=(double*)malloc(sizeof(*E5)*MAXIMUM_H);
    x5=(double*)malloc(sizeof(*x5)*MAXIMUM_H);
    y5=(double*)malloc(sizeof(*y5)*MAXIMUM_H);
    power5=(double*)malloc(sizeof(*power5)*MAXIMUM_H);
    flux5=(double*)malloc(sizeof(*flux5)*MAXIMUM_H);
  }
  for (page=0; page<inputPages; page++) {
    /* change the pinhole nXP and nYP back to the input value */
    pinhole_param.nXP = nXP0;
    pinhole_param.nYP = nYP0; 
    if (inputSupplied) {
      electron_param.sigmax = twiss.sigmax[page];
      electron_param.sigmaxp = twiss.sigmaxp[page];
      electron_param.sigmay = twiss.sigmay[page];
      electron_param.sigmayp = twiss.sigmayp[page];
      electron_param.energySpread = twiss.Sdelta0[page];
      electron_param.energy = twiss.pCentral[page]*me_mev*0.001;
      if (undulator_param.energy) {
        /* calculate the Kx value for the undulator */
        double gamma, lambda;
        gamma = electron_param.energy*1e3/me_mev;
        lambda = h_mks*c_mks/(undulator_param.energy*e_mks);
        undulator_param.kx = 0;
        undulator_param.ky = sqrt(4*gamma*gamma*lambda/undulator_param.period-2);
      }
      if (emin==emax && emin==0) {
        double gamma, lambda;
        gamma = electron_param.energy*1e3/me_mev;
        lambda = undulator_param.period/(2*gamma*gamma)*(1+0.5*(sqr(undulator_param.kx)+sqr(undulator_param.ky)));
        emin = emax = h_mks*c_mks/lambda/e_mks;
        nE = 2;
      }
      check_input_parameters(&undulator_param, &electron_param, &pinhole_param, nE, nPhi, nAlpha, nOmega, dOmega,
                             mode, icalc, iharm, inputSupplied, us);
    }
    
    if (!us) {
      urgent_(&undulator_param.itype,&undulator_param.period,&undulator_param.kx,
              &undulator_param.ky,&undulator_param.phase,&undulator_param.nPeriod,
              &emin,&emax,&nE,
              &electron_param.energy,&electron_param.current,
              &electron_param.sigmax,&electron_param.sigmay,
              &electron_param.sigmaxp,&electron_param.sigmayp,
              &pinhole_param.distance,&pinhole_param.xPC,&pinhole_param.yPC,
              &pinhole_param.xPS,&pinhole_param.yPS,&pinhole_param.nXP,&pinhole_param.nYP,
              &mode, &icalc, &iharm, &nPhi, &electron_param.nsigma, &nAlpha,
              &dAlpha,&nOmega, &dOmega, /*end of input parameters */
              &E1,&lamda1,&ptot,&pd,&isub,&iang,
              EE,lamda,&min_harmonic,&max_harmonic,&nEE,
              xPMM,yPMM,irradiance,L1,L2,L3,L4,&max_irradiance,power,
              I1,I2,&i_max,EI,spec1,spec2,
              spec3,&pdtot,&ptot1,&ftot,harmonics,
              &iharm5,x5,y5,E5,power5,flux5);
    } else {
      /*in us, undulator period uses cm units, while urgent uses m units */
      /*in us, current is in mA units, while urgent is in A units */
      period=undulator_param.period*100;
      current=electron_param.current*1000;
      us_(&electron_param.energy, &current, &electron_param.sigmax,
          &electron_param.sigmay, &electron_param.sigmaxp, &electron_param.sigmayp,
          &period, &undulator_param.nPeriod, &undulator_param.kx, &undulator_param.ky,
          &emin, &emax, &nE, 
          &pinhole_param.distance, &pinhole_param.xPC, &pinhole_param.yPC, 
          &pinhole_param.xPS, &pinhole_param.yPS, &pinhole_param.nXP, &pinhole_param.nYP,
          &mode, &icalc, &iharm, &nPhi, &nAlpha, &dAlpha, &nOmega, &dOmega, 
          &electron_param.nsigma, /*end of input parameters */
          &E1, &lamda1, &ptot, &pd, &ptot1, &isub, &iang, 
          &ftot, &max_harmonic, &min_harmonic, L1, L2, L3, L4, &nEE,
          xPMM, yPMM, irradiance, spec1, EE);
     }
    
    
    output = outputfile;
    if (!us) {
      if (!page) {
        sddsout2=SetupOutputFile(output, &sddsout, mode, undulator_param.itype, icalc, isub, iang, idebug, iharm, special, &sddsout1);
      }
      
      SpecifyDescription(&description, isub, mode, iang, i_max, iharm);
      if (!iang)
        pd = pd / (pinhole_param.distance * pinhole_param.distance);
      
      WriteToOutput(&sddsout, sddsout2, description, undulatorType, idebug,
                    undulator_param.itype, undulator_param.period, undulator_param.kx,
                    undulator_param.ky, undulator_param.phase, undulator_param.nPeriod,
                    emin, emax, nE,
                    electron_param.energy, electron_param.energySpread, electron_param.current,
                    electron_param.sigmax, electron_param.sigmay, electron_param.sigmaxp,
                    electron_param.sigmayp, pinhole_param.distance, pinhole_param.xPC,
                    pinhole_param.yPC, pinhole_param.xPS, pinhole_param.yPS, pinhole_param.nXP,
                    pinhole_param.nYP, mode, icalc, iharm, nPhi, electron_param.nsigma, nAlpha,
                    dAlpha,nOmega, dOmega,
                    E1,lamda1,ptot,pd,isub,iang,
                    EE,lamda,min_harmonic,max_harmonic,nEE,
                    xPMM,yPMM,irradiance,L1,L2,L3,L4,max_irradiance,power,
                    I1,I2,i_max,EI,spec1,spec2,
                    spec3,pdtot,ptot1,ftot,harmonics, &sddsout1,
                    special, iharm5,x5,y5,E5,power5,flux5);
      free(description);
    } else {
      if (!page) 
        SetupUSOutput(&sddsout, outputfile, mode, iang, isub); 
      if (!iang)
        pd = pd /(pinhole_param.distance * pinhole_param.distance); 
       WriteUSResultsToOutput(&sddsout, undulator_param, electron_param, pinhole_param,
                              nPhi, nAlpha, dAlpha, nOmega, dOmega,
                              mode, icalc, iharm, nEE, isub, iang,
                              ptot, pd, ptot1, ftot, max_harmonic, min_harmonic,
                              emin, emax, E1, lamda1, xPMM, yPMM,
                              L1, L2, L3, L4, EE, irradiance, spec1);
    }
  } /*end of for page */
  /* terminate sdds file */
  if (!us) {
    if (special && !SDDS_Terminate(&sddsout1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (sddsout2) {
      if (!SDDS_Terminate(sddsout2))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      free(sddsout2);
      sddsout2 = NULL;
    }
    special=0;
  }
  if (!SDDS_Terminate(&sddsout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      
  
  if (special) {
    free(x5);
    free(y5);
    free(E5);
    free(flux5);
    free(power5);
  }
  free(undulatorType);
  free(EE);
  if (lamda) free(lamda);
  free(xPMM);
  free(yPMM);
  free(irradiance);
  free(L1);
  free(L2);
  free(L3);
  free(L4);
  if (I1) free(I1);
  if (I2) free(I2);
  if (harmonics) free(harmonics);
  if (EI) free(EI);
  if (power) free(power);
  free(spec1);
  if (spec2) free(spec2);
  if (spec3) free(spec3);
  if (inputSupplied) {
    for (i=0; i<TWISS_DATA_TYPES; i++)
      free(twiss.data[i]);
    free(twiss.data);
  }
  
  free_scanargs(&s_arg, argc); 
  return 0;
}

SDDS_DATASET *SetupOutputFile(char *outputfile, SDDS_DATASET *SDDSout, long mode, long itype, long icalc, 
                              long isub, long iang, long idebug, int32_t iharm, long special, SDDS_DATASET *SDDSout1)
{
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, outputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  DefineParameters(SDDSout, mode, iang, 0);
  if (isub==5) {
    SDDS_DATASET *SDDSout2=NULL;
    if (iang==0) {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (iharm>0) {
      if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (SDDS_DefineColumn(SDDSout, "PowerDensity", NULL, "Watts/mm^2 or mrad^2", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "Flux", NULL, "photons/s/0.1%bandwidth", NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (outputfile) {
      if (iharm<=0) {
        char out2[1024];
        SDDSout2=(SDDS_DATASET *)malloc(sizeof(*SDDSout2));
        sprintf(out2, "%s.total", outputfile);
        if (!SDDS_InitializeOutput(SDDSout2, SDDS_BINARY, 1, NULL, NULL, out2))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (SDDS_DefineParameter(SDDSout2, "Description",  NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
            SDDS_DefineParameter(SDDSout2, "TotalPowerDensity",  NULL, 
                                 "Watts/mm^2 or mrad^2", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineParameter(SDDSout2, "TotalPower",  NULL, "Watts", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineParameter(SDDSout2, "TotalFlux",  NULL, 
                                 "photons/s/0.1%bandwidth", NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (SDDS_DefineColumn(SDDSout2, "Harmonic", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
            SDDS_DefineColumn(SDDSout2, "Energy", NULL, "eV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout2, "PowerDensity", NULL, "Watts/mm^2 or mrad^2", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout2, "Power", NULL, "Watts", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout2, "Flux", NULL, "photons/s/0.1%bandwidth", NULL, NULL, SDDS_DOUBLE, 0)<0 )
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (!SDDS_WriteLayout(SDDSout2))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (special) {
        char out1[1024];
        sprintf(out1,"%s.harmonic",outputfile);
        if (!SDDS_InitializeOutput(SDDSout1,SDDS_BINARY, 0, NULL, NULL, out1))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (SDDS_DefineParameter(SDDSout1, "Description",  NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (SDDS_DefineColumn(SDDSout1, "X", NULL, NULL, "mm", NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout1, "Y", NULL, NULL, "mm", NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout1, "Energy", NULL, "eV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout1, "PowerDensity", NULL, "Watts/mm^2 or mrad^2", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout1, "Flux", NULL, "photons/s/0.1%bandwidth", NULL, NULL, SDDS_DOUBLE, 0)<0 )
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        if (!SDDS_WriteLayout(SDDSout1))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (!SDDS_WriteLayout(SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    return SDDSout2;
  } else {
    if (mode==1) {
      if (iang==0) {
        if (SDDS_DefineColumn(SDDSout, "X", NULL, "mm", "X position",NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout, "Y", NULL, "mm", "Y position",NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (SDDS_DefineColumn(SDDSout, "X", NULL, "mrad", "X position",NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineColumn(SDDSout, "Y", NULL, "mrad", "Y position",NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    
    if (isub==2 || isub==3 || (isub==4 && mode!=1)) {
      /*isub=2*/
      if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", "photon energy",NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Lambda", NULL, "Angstrom", "photon wavelength",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    /*isub=1,2*/
    if (mode==1 || (mode==2 && (isub==2 || isub==4))) {
      if (iang==0) {
        if (SDDS_DefineColumn(SDDSout, "Irradiance", NULL, "photons/s/mm^2/0.1%bandwidth", 
                              "spatial flux density",NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (SDDS_DefineColumn(SDDSout, "AngularFluxDensity", NULL, "photons/s/mrad^2/0.1%bandwidth", 
                              "angular flux density",NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (mode==1 && isub==4) {
      if (SDDS_DefineColumn(SDDSout, "I", NULL, NULL, NULL,NULL, SDDS_LONG, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (isub==2 || isub==4) {
      switch (mode) {
      case 2:
        if (SDDS_DefineColumn(SDDSout, "SpectralPowerDensity", NULL, "Watts/mm^2 or mrad^2/eV bandwidth", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      case 3:
        if (SDDS_DefineColumn(SDDSout, "Brightness", NULL, "photons/s/mm^2/mrad^2/0.1%bandwidth", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      case 4:
        if (SDDS_DefineColumn(SDDSout, "Flux", NULL, "photons/s/0.1%bandwidth", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0 || 
            SDDS_DefineColumn(SDDSout, "SpectralPower", NULL, "Watts/eV bandwidth", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      }
    }
    if (isub==3) {
      if (SDDS_DefineColumn(SDDSout, "Flux", NULL, "photons/s/0.1%bandwidth", 
                            NULL, NULL, SDDS_DOUBLE, 0)<0 || 
          SDDS_DefineColumn(SDDSout, "SpectralPower", NULL, "Watts/eV bandwidth", 
                            NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode!=1) {
      if (SDDS_DefineColumn(SDDSout, "IMIN", NULL, NULL, NULL,NULL, SDDS_LONG, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "IMAX", NULL, NULL, NULL,NULL, SDDS_LONG, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode!=3) {
      if (SDDS_DefineColumn(SDDSout, "L1", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "L2", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "L3", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "L4", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return NULL;
}

void SpecifyDescription(char **description, long isub, long mode, long iang, long i_max, int32_t iharm)
{             
  char desc[2046];
  
  switch (isub) {
  case 1:
    if (iang==0)
      SDDS_CopyString(description,"Spatial flux density distribution at minimum photon energy.");
    else
      SDDS_CopyString(description,"Angular flux density distribution at minimum photon energy.");
    break;
  case 2:
    switch (mode)  {
    case 2:
      if (iang==0) 
        SDDS_CopyString(description,"Spatial flux density at position xPC, yPC.");
      else
        SDDS_CopyString(description,"Angular flux density at position xPC, yPC.");
      break;
    case 3:
      SDDS_CopyString(description, "Spectrum of on-axis brightness. Distance, XPC, YPC, XPS, YPS, NXP and NYP are irrelevant.");
      break;
    case 4:
      SDDS_CopyString(description, "Spectrum of flux integrated over all angles.");
      break;
    }
    break;
  case 3:
    SDDS_CopyString(description, "Spectrum of flux integrated over all angles. Distance, XPC, YPC, XPS, YPS, NXP and NYP are irrelevant.");
    break;
  case 4:
    switch (mode) {
    case 1:
      if (iang==0)
        SDDS_CopyString(description,"Spatial flux density distribution at minimum photon energy with zero emittance.");
      else
        SDDS_CopyString(description,"Angular flux density distribution at minimum photon energy with zero emittance.");
      break;
    case 2:
      if (iang==0) 
        SDDS_CopyString(description,"Spatial flux density at position xPC, yPC with zero emittance.");
      else
        SDDS_CopyString(description,"Angular flux density at position xPC, yPC with zero emittance.");
      break;
    case 3:
      SDDS_CopyString(description, "Spectrum of on-axis brightness with zero emittance.");
      break;
    case 4:
      SDDS_CopyString(description, "Spectrum of flux integrated over all angles with zero emittance.");
      break;
    }
    break;
  case 5:
    if (iharm>0) 
      sprintf(desc,"Angular distribution for harmonic %d.",iharm);
    else
      sprintf(desc,"Angular distribution for harmonic from 1 to %ld.",i_max);
    SDDS_CopyString(description, desc);
    break;
  default:
    fprintf(stderr,"Invalid isub obtained from urgent()\n");
    exit(1);
    break;
  }
}

long GetISub(long mode, long icalc)
{
  long isub=0;
  if (mode==1) isub=1;
  if (mode>=2 && mode<=4) isub=2;
  if (mode==5) isub=3;
  if (icalc==3) isub=4;
  if (mode==6 || mode==-6) isub=5;
  if (isub==0) {
    fprintf(stderr,"invalid mode and calculation method given.\n");
    exit(1);
  }
  return isub;
}

void WriteToOutput(SDDS_DATASET *SDDSout, SDDS_DATASET *SDDSout2, char *description, char *undulatorType, long idebug,
                   long itype, double period, double kx, double ky, double phase,
                   long nPeriod, double emin, double emax, long nE,
                   double energy, double energySpread, double current, double sigx, double sigy, double sigxp,
                   double sigyp, double distance, double xPC, double yPC, 
                   double xPS, double yPS,long nXP, long nYP,
                   long mode, long icalc, int32_t iharm, long nPhi, long nSig, long nAlpha,
                   double dAlpha, long nOmega, double dOmega,
                   double E1, double lamda1, double ptot, double pd, long isub, long iang,
                   double *EE, double *lamda, long min_harmonic, long max_harmonic, long nEE,
                   double *xPMM, double *yPMM, double *irradiance, 
                   double *L1, double *L2, double *L3, double *L4, double max_irradiance, double *power,
                   long *I1, long *I2, long i_max, double *EI, double *spec1, double *spec2,
                   double *spec3, double pdtot, double ptot1, double ftot, long *harmonics,
                   SDDS_DATASET *SDDSout1, long special, int32_t iharm5, double *x5, double *y5, double *E5,
                   double *power5, double *flux5) 
{
  long i, n1=0, n2=0;
  char desc[2046];
  if (special) {
    for (i=0;i<iharm5;i++) {
      if (!SDDS_StartPage(SDDSout1, nEE))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      sprintf(desc,"Angular Distribution for harmonic %ld",i+1);
      if (!SDDS_SetParameters(SDDSout1,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"Description", desc,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetColumn(SDDSout1, SDDS_SET_BY_NAME, x5+i*nEE, nEE, "X") ||
          !SDDS_SetColumn(SDDSout1, SDDS_SET_BY_NAME, y5+i*nEE, nEE, "Y") ||
          !SDDS_SetColumn(SDDSout1, SDDS_SET_BY_NAME, E5+i*nEE, nEE, "Energy") ||
          !SDDS_SetColumn(SDDSout1, SDDS_SET_BY_NAME, power5+i*nEE, nEE, "PowerDensity") ||
          !SDDS_SetColumn(SDDSout1, SDDS_SET_BY_NAME, flux5+i*nEE, nEE, "Flux"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_WritePage(SDDSout1))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (isub!=5 && energySpread>0) {
    /* does guass-convolve for energy-spectrum*/
    if (isub==2 || isub==3 || (isub==4 && mode!=1)) {
      Gauss_Convolve(EE, irradiance, &nEE, energySpread, &n1, &n2);
      for (i=n1;i<=n2;i++) 
        lamda[i-n1]=lamda[i];
      if (isub==2 || isub==4) {
        switch (mode) {
        case 2:
        case 4:
          for (i=n1;i<=n2;i++)
            power[i-n1] = power[i];
          break;
        }
      } else if (isub==3) {
        for (i=n1;i<=n2;i++)
          power[i-n1] = power[i];
      }
      if (mode!=1) {
        for (i=n1;i<=n2;i++) {
          I1[i-n1] = I1[i];
          I2[i-n1] = I2[i];
        }
      }
      if (mode!=3) {
        for (i=n1;i<=n2;i++) {
          L1[i-n1] = L1[i];
          L2[i-n1] = L2[i];
          L3[i-n1] = L3[i];
          L4[i-n1] = L4[i];
        }
      }
    }
  }
  if (!SDDS_StartPage(SDDSout, nEE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_SetParameters(SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"Description", description,
                          "UndulatorType", undulatorType, "Period", period, "Phase", phase,
                          "Kx", kx, "Ky", ky, "NPeriod", nPeriod, "MinEnergy", emin, "MaxEnergy", emax,
                          "NE", nEE, "ElectronEnergy", energy, "ElectronCurrent", current,
                          "SigmaX", sigx, "SigmaY", sigy, "SigmaXPrime", sigxp, "SigmaYPrime", sigyp,
                          "AcceptanceDistance", distance, "XCenter", xPC, "YCenter", yPC,
                          "XSize", xPS, "YSize", yPS, "NXP", nXP, "NYP", nYP,
                          "Mode", mode, "ICalc", icalc, "Harmonics", iharm, "NPhi", nPhi,
                          "NSig", nSig, "NAlpha", nAlpha, "DAlpha", dAlpha, "NOmega", nOmega,
                          "DOmega", dOmega, "E1", E1, "Lambda1", lamda1, "TotalPower", ptot, "OnAxisPowerDensity", pd,NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (isub==5) {
    if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, xPMM, nEE, "X") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, yPMM, nEE, "Y") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, power, nEE, "PowerDensity") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Flux"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (iharm>0) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, EE, nEE, "Energy"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (SDDSout2) {
      if (!SDDS_StartPage(SDDSout2, i_max))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetParameters(SDDSout2, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                              "Description", "Central power density and integrated power over the range of acceptance.",
                              "TotalPowerDensity", pdtot, "TotalPower", ptot1, "TotalFlux", ftot, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetColumn(SDDSout2, SDDS_SET_BY_NAME, harmonics, i_max, "Harmonic") ||
          !SDDS_SetColumn(SDDSout2, SDDS_SET_BY_NAME, EI, i_max, "Energy") ||
          !SDDS_SetColumn(SDDSout2, SDDS_SET_BY_NAME, spec1, i_max, "PowerDensity") ||
          !SDDS_SetColumn(SDDSout2, SDDS_SET_BY_NAME, spec2, i_max, "Power") ||
          !SDDS_SetColumn(SDDSout2, SDDS_SET_BY_NAME, spec3, i_max, "Flux"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_WritePage(SDDSout2))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else {
    if (mode==1) {
      if (!SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Energy", EE[0], 
                              "Lambda", lamda[0],
                              "MaxIrradiance", max_irradiance, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "MinHarmonics", min_harmonic, 
                              "MaxHarmonics", max_harmonic, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, xPMM, nEE, "X") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, yPMM, nEE, "Y"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (isub==4) {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, I1, nEE, "I"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (isub==2 || isub==3 || (isub==4 && mode!=1)) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, EE, nEE, "Energy") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, lamda, nEE, "Lambda"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==2) {
      if (!SDDS_SetParameters(SDDSout,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                              "TotalPowerDensity1", pdtot,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==4 || mode==5) {
      if (!SDDS_SetParameters(SDDSout,SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                              "TotalPower1", ptot1,NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    if (mode==1 || (mode==2 && (isub==2 || isub==4))) {
      if (iang==0) {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Irradiance"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "AngularFluxDensity"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (isub==2 || isub==4) {
      switch (mode) {
      case 2:
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, power, nEE, "SpectralPowerDensity"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      case 3:
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Brightness"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      case 4:
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Flux") ||
            SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, power, nEE, "SpectralPower"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        break;
      }
    }
    if (isub==3) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Flux") || 
          SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, power, nEE, "SpectralPower"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode!=1) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, I1, nEE, "IMIN") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, I2, nEE, "IMAX"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode!=3) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L1, nEE, "L1") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L2, nEE, "L2") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L3, nEE, "L3") ||
          !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L4, nEE, "L4"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WritePage(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void DefineParameters(SDDS_DATASET *SDDSout, long mode, long iang, long us) {
  if (SDDS_DefineParameter(SDDSout, "Description", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "UndulatorType", NULL, NULL,
                           "corresponding to itype=1 (Undulator) or itype=2(canted undulator)", 
                           NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Kx", NULL, NULL, "manget field", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Ky", NULL, NULL, "manget field", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Period", NULL, "m", "length of one period of undulator", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Phase", NULL, "degree",
                           "phase difference between the two undulators in the canted undulator",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NPeriod", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MinEnergy", NULL, "eV", "Minimum photon energy", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MaxEnergy", NULL, "eV", "Maximum photon energy",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NE", NULL, NULL, "Number of Energy points",
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "ElectronEnergy", NULL, "GeV", "Electron beam energy", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "ElectronCurrent", NULL, "A", "Electron beam current",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "SigmaX", NULL, "mm", "horizontal rms electron beam size",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "SigmaY", NULL, "mm", "vertical rms electron beam size",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "SigmaXPrime", NULL, "mrad", 
                           "horizontal rms electron beam divergence", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "SigmaYPrime", NULL, "mrad", 
                           "vertical rms electron beam divergence", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "AcceptanceDistance", NULL, "mm", "acceptance distance",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "XCenter", NULL, "mm", "acceptance center x position",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "YCenter", NULL, "mm", "acceptance center y position", 
                           NULL, SDDS_DOUBLE, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (iang) {
    if (SDDS_DefineParameter(SDDSout, "XSize", NULL, "mrad", "total horizontal acceptance",
                             NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "YSize", NULL, "mrad", "total vertical acceptance", 
                             NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else { 
    if (SDDS_DefineParameter(SDDSout, "XSize", NULL, "mm", "total horizontal acceptance",
                             NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "YSize", NULL, "mm", "total vertical acceptance", 
                             NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (SDDS_DefineParameter(SDDSout, "NXP", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NYP", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Mode",  NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "ICalc",  NULL, NULL, "calculation method", 
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Harmonics",  NULL, NULL,
                           "number of harmonics to be included", NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NPhi",  NULL, NULL, 
                           "no. of steps in phi between 0 and pi/2.0 used in (MODE=1,2,3,4,5 ICALC=1,2)",
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NSig",  NULL, NULL, 
                           "no. of standard deviations of electron beam dimensions (size and divergence) to be included used in (MODE=1,2,3,4 ICALC=1,2) and (MODE=6 ICALC=1)", 
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NAlpha",  NULL, NULL, 
                           "no. of steps in angle alpha (gamma*theta) used in (MODE=1 ICALC=1)", 
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "DAlpha",  NULL, NULL, 
                           "range of angles in alpha**2 to be used, in units of the angular equivalent of 1/N, used in (MODE=1 ICALC=1) and ICALC=3 ", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NOmega",  NULL, NULL, 
                           "no. of steps in photon energy for the natural lineshape",
                           NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "DOmega",  NULL, NULL, 
                           "range of photon energies to be included in the natural lineshape", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "E1",  NULL, "eV", "energy of first harmonic", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Lambda1",  NULL, "Angstrom", "wavelenth of first harmonic", 
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "TotalPower",  NULL, "Watts", NULL, 
                           NULL, SDDS_DOUBLE, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  if (iang) {
    if (SDDS_DefineParameter(SDDSout, "OnAxisPowerDensity", NULL, "Watts/mrad^2", NULL,
                             NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else {
    if (SDDS_DefineParameter(SDDSout, "OnAxisPowerDensity", NULL, "Watts/mm^2", NULL,
                             NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!us) {
    if (mode==1) {
      if (SDDS_DefineParameter(SDDSout, "Energy",  NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineParameter(SDDSout, "Lambda",  NULL, "Angstrom", NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      
      if (SDDS_DefineParameter(SDDSout, "MinHarmonics",  NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
          SDDS_DefineParameter(SDDSout, "MaxHarmonics",  NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      
      if (SDDS_DefineParameter(SDDSout, "MaxIrradiance",  NULL, "photons/s/mm^2/0.1%bandwidth", NULL, NULL,
                               SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);    
    }
    
    if (mode==2) {
      if (SDDS_DefineParameter(SDDSout, "TotalPowerDensity1", NULL, "Watts/mrad^2", NULL,
                               NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }  
    if (mode==4 || mode==5) {
      if (SDDS_DefineParameter(SDDSout, "TotalPower1",  NULL, "Watts", NULL, 
                               NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else {
    if (SDDS_DefineParameter(SDDSout, "MaximumHarmonics", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "MinimumHarmonics", NULL, NULL, NULL, NULL, SDDS_LONG, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (mode==1) {
      if (SDDS_DefineParameter(SDDSout, "TotalFluxDensity", NULL, "photons/s/0.1%bw", 
                               NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (mode==2) {
      if (iang) {
        if (SDDS_DefineParameter(SDDSout, "IntegFluxDensity", NULL, "photons/s/mrad^2", NULL,
                                 NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineParameter(SDDSout, "IntegPowerDensity", NULL, "watts/s/mrad^2", NULL,
                                 NULL, SDDS_DOUBLE, 0)<0 )
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (SDDS_DefineParameter(SDDSout, "IntegFluxDensity", NULL, "photons/s/mm^2", NULL,
                                 NULL, SDDS_DOUBLE, 0)<0 ||
            SDDS_DefineParameter(SDDSout, "IntegPowerDensity", NULL, "watts/s/mm^2", NULL,
                                 NULL, SDDS_DOUBLE, 0)<0 )
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if (mode==4 || mode==5) {
      if (SDDS_DefineParameter(SDDSout, "IntegFlux", NULL, "photons/s", NULL,
                               NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineParameter(SDDSout, "IntegPower", NULL, "watts/s", NULL,
                               NULL, SDDS_DOUBLE, 0)<0 )
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  
}

void SetupUSOutput(SDDS_DATASET *SDDSout, char *outputfile, long mode, long iang, long isub)
{
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, outputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  DefineParameters(SDDSout, mode, iang, 1);
  if (isub==1 || (isub==5 && mode==1)) {
    if (iang) {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0  ||
          SDDS_DefineColumn(SDDSout, "AngularFluxDensity", NULL, "photons/s/mrad^2/0.1%bandwidth", 
                            "angular flux density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Irradiance", NULL, "photons/s/mm^2/0.1%bandwidth", 
                            "spatial flux density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else if (isub==2 || isub==3 || isub==5) { 
    if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (mode==2) {
      if (iang) {
        if (SDDS_DefineColumn(SDDSout, "AngularFluxDensity", NULL, "photons/s/mrad^2/0.1%bw", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (SDDS_DefineColumn(SDDSout, "Irradiance", NULL, "photons/s/mm^2/0.1%bw", 
                              NULL, NULL, SDDS_DOUBLE, 0)<0)
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      
    } else if (mode==3) {
      if (SDDS_DefineColumn(SDDSout, "Brilliance", NULL, "photons/s/mrad^2/mm^2/0.1%bw",
                            NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else if (mode==4 || mode==5) {
      if (SDDS_DefineColumn(SDDSout, "Flux", NULL, "photons/s/0.1%bw",
                            NULL, NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else if (isub==4) {
    if (iang) {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0  ||
          SDDS_DefineColumn(SDDSout, "PowerDensity", NULL, "watts/mrad^2", 
                            "angular power density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (SDDS_DefineColumn(SDDSout, "X", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
          SDDS_DefineColumn(SDDSout, "Y", NULL, "mm", NULL, NULL, SDDS_DOUBLE, 0)<0  ||
          SDDS_DefineColumn(SDDSout, "PowerDensity", NULL, "watts/mm^2", 
                            "power density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (mode!=3)
    if (SDDS_DefineColumn(SDDSout, "L1", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "L2", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "L3", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "L4", NULL, NULL, NULL,NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void WriteUSResultsToOutput(SDDS_DATASET *SDDSout, UNDULATOR_PARAM  undulator_param, 
                            ELECTRON_BEAM_PARAM electron_param, PINHOLE_PARAM pinhole_param,
                            long nPhi, long nAlpha, double dAlpha, long nOmega, double dOmega,
                            long mode, long icalc, int32_t iharm, long nEE, long isub, long iang,
                            double ptot, double pd, double ptot1, double ftot, long imax, long imin,
                            double emin, double emax, double E1, double lamda1, double *xPMM, double *yPMM,
                            double *L1, double *L2, double *L3, double *L4, double *EE, double *irradiance,
                            double *spec1)
{
  char desc[2054];
  long i, n1=0, n2=0;
  
  if (isub==1 || (isub==5 && mode==1)) {    
    if (iang)
      sprintf(desc,"Angular flux density distribution for %lf ev (ph/s/mr^2/0.1%%bw).", emin);
    else
      sprintf(desc, "Irradiance for %lf ev (ph/s/mm^2/0.1%%bm).", emin);
  } else if (isub==2 || isub==3 || isub==5) { 
    if (mode==2) 
      sprintf(desc, "Angular flux density spectrum for x=%f mrad and y=%lf mrad.", xPMM[0], yPMM[0]);
    else if (mode==3)
      sprintf(desc,"%s", "On-axis brilliance (ph/s/mrad^2/mm^2/0.1%%bw).");
    else if (mode==4) {
      if (iang)
        sprintf(desc, "Flux through %lf mrad, %lf mrad pinhole at %lf mrad, %lf mrad.", 
                pinhole_param.xPS, pinhole_param.yPS, pinhole_param.xPC, pinhole_param.yPC);
      else
        sprintf(desc, "Flux through %f mrad, %f mm pinhole at %f mm, %f mm.", 
                pinhole_param.xPS, pinhole_param.yPS, pinhole_param.xPC, pinhole_param.yPC);
    } else if (mode==5) 
      sprintf(desc, "%s", "Angle-integrated spectrum.");
  } else if (isub==4) {
    if (iang) 
      sprintf(desc, "%s", "Power density distribution (W/mrad^2)");
    else
      sprintf(desc, "%s", "Power density distribution (W/mm^2)");
  }
  /* write values */
  if (icalc==3)
    strcat(desc, "  Zero emittance. US program");
  else
    strcat(desc, "  US program");
  if (electron_param.energySpread>0) {
    if (mode!=1 && mode!=6) {
      Gauss_Convolve(EE, spec1, &nEE, electron_param.energySpread, &n1, &n2);
      if (mode!=3) {
        for (i=n1; i<=n2; i++) {
          L1[i-n1] = L1[i];
          L2[i-n1] = L2[i];
          L3[i-n1] = L3[i];
          L4[i-n1] = L4[i];
        }
      }
    }
  }
  if (!SDDS_StartPage(SDDSout, nEE))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_SetParameters(SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"Description", desc,
                          "UndulatorType", "Undulator", 
                          "Period", undulator_param.period, "Phase", 0.0,
                          "Kx", undulator_param.kx, "Ky", undulator_param.ky, "NPeriod", undulator_param.nPeriod, 
                          "MinEnergy", emin, "MaxEnergy", emax,
                          "NE", nEE, "ElectronEnergy", electron_param.energy, "ElectronCurrent", electron_param.current,
                          "SigmaX", electron_param.sigmax, "SigmaY", electron_param.sigmay, 
                          "SigmaXPrime", electron_param.sigmaxp, "SigmaYPrime", electron_param.sigmayp,
                          "AcceptanceDistance", pinhole_param.distance, "XCenter", pinhole_param.xPC, 
                          "YCenter", pinhole_param.yPC, "XSize", pinhole_param.xPS, "YSize", pinhole_param.yPS, 
                          "NXP", pinhole_param.nXP, "NYP", pinhole_param.nYP,
                          "Mode", mode, "ICalc", icalc, "Harmonics", iharm, "NPhi", nPhi,
                          "NSig", electron_param.nsigma, "NAlpha", nAlpha, "DAlpha", dAlpha, "NOmega", nOmega,
                          "DOmega", dOmega, "E1", E1, "Lambda1", lamda1, "TotalPower", ptot, "OnAxisPowerDensity", pd,
                          "MaximumHarmonics", imax, "MinimumHarmonics", imin, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (mode==1) {
    if (!SDDS_SetParameters(SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "TotalFluxDensity", ftot, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else if (mode==2) {
    if (!SDDS_SetParameters(SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "IntegFluxDensity", ftot, "IntegPowerDensity", ptot1, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } else if (mode==4 || mode==5) {
    if (!SDDS_SetParameters(SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "IntegFlux", ftot, "IntegPower", ptot1, NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode!=3) {
    if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L1, nEE, "L1") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L2, nEE, "L2") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L3, nEE, "L3") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, L4, nEE, "L4"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (mode==1 || mode==6) {
    if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, xPMM, nEE, "X") ||
        !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, yPMM, nEE, "Y"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (mode==1) {
      if (iang) {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "AngularFluxDensity"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "Irradiance"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    } else {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, irradiance, nEE, "PowerDensity"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  } else {
    if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, EE, nEE, "Energy"))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (mode==2) {
      if (iang) {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, spec1, nEE, "AngularFluxDensity"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } else {
        if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, spec1, nEE, "Irradiance"))
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    } else if (mode==3) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, spec1, nEE, "Brilliance"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else if (mode==4 || mode==5) {
      if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, spec1, nEE, "Flux"))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_WritePage(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors); 
}

void check_input_parameters(UNDULATOR_PARAM *undulator_param, ELECTRON_BEAM_PARAM *electron_param, 
                            PINHOLE_PARAM *pinhole_param, long nE, long nPhi, long nAlpha, long nOmega, double dOmega,
                            long mode, long icalc, int32_t iharm, long inputSupplied, long us)
{
  long error=0;
  if (nPhi>100 || nAlpha>100) {
    fprintf(stderr, "Number of steps of phi or alpha (nPhi, nAlpha) should not exceed 100.\n");
    error ++;
  }
  
  if (us) {
    if (mode<1 || mode>6) {
      fprintf(stderr, "The calculation mode for US has to be an integer from 1 to 6.\n");
      error++;
    }
    if (icalc==2) {
      fprintf(stderr, "us calculation method 2 is for test purposes only, use 1,3,4,14 instead.\n"); 
      error++;
    }
    
    if (icalc<1 || (icalc>3 && icalc!=4 && icalc!=14)) {
      fprintf(stderr, "The calculation method for US has to be 1, 3, 4, or 14.\n");
      error++;
    }
    if (mode==1 && (icalc==2 || icalc==4 || icalc==14)) {
      fprintf(stderr, "Method=%ld, not valid for the flux density distribution.\n", icalc);
      error++;
    }
    if (mode==5 && (icalc==1 || icalc==3)) {
      fprintf(stderr, "Method=%ld, not valid for angle-integrated spectrum.\n", icalc);
      error++;
    }
    if (undulator_param->itype!=1) {
      fprintf(stderr, "US computation only works for regular undulators.\n");
      error++;
    }
  } else {
    if (undulator_param->itype<1 || undulator_param->itype>2) {
      fprintf(stderr, "Invalid undulator type given, it can only be 1 (regular undulator) or 2 (canted undulator whose phase is greater than 0).\n");
      error++;
    }
    if (mode>6 || (mode<1 && mode!=-6)) {
      fprintf(stderr, "Invalid calculation mode given for urgent, it has to be an integer from 1 to 6 or -6.\n");
      error++;
    }
    if (icalc<1 || icalc>3) {
      fprintf(stderr, "Invalid calculation method give for urgent, it has to be an integer from 1 to 3.\n");
      error++;
    }
    if (mode==5 && icalc==3) {
      fprintf(stderr, "In URGENT, integrated flux spectrum (mode=5) is not available for WalkerFinite method (method=3).\n");
      error++;
    }
    if (mode==-6 && iharm>0) {
      fprintf(stderr, "In URGENT, mode=-6 is only valid for harmonics<=0.\n");
      error++;
    }
    
    if (abs(mode)==6 && icalc==3) {
      fprintf(stderr, "In URGENT, mode=6 or mode=-6 (power density and integrated power calculation) is not available for method=3 (WalkerFinite method).");
      error++;
    }
    if (undulator_param->itype==2 && (mode==3 || mode==5 || abs(mode)==6)) {
      fprintf(stderr, "In URGENT, mode=3, 5, 6 or -6 is not available for canted undulator (itype=2, phase>0).\n");
      error++;
    }
    if (iharm<-1000 || iharm>1000) {
      fprintf(stderr, "The harmonics number %ld is out of URGENT harmonics range (-1000, 1000).\n", iharm);
      error++;
    }
    if (nE<1 || nE>5001) {
      fprintf(stderr, "The number of photon energy points (%ld) is out of URGENT range (1, 5001).\n", nE);
      error++;
    }
    
    if (undulator_param->itype==2 && nPhi>25) {
      fprintf(stderr, "The value of nPhi (%ld) can not be greater than 25 for canted undulator (phase>0) calculation.\n", nPhi);
      error++;
    }
    
    if (icalc==1 && (mode==2 || mode==3 || mode==5) && undulator_param->itype==2) {
      if (nOmega>5000) {
        fprintf(stderr, "Too many omega points %ld (exceeds 5000).\n", nOmega);
        error++;
      }
      if (nOmega/dOmega<4.0) {
        fprintf(stderr, "Rule nOmega/dOmega >= 4 was expected for itype=1 (non-canted undulator).\n");
        error++;
      }
    }
  }
  
  if (pinhole_param->nXP<0 || pinhole_param->nXP>500) {
    fprintf(stderr, "The xnumber of pinhole parameter - %ld is out of range (0, 500).\n", pinhole_param->nXP);
    error++;
  }
  if (pinhole_param->nYP<0 || pinhole_param->nYP>500) {
    fprintf(stderr, "The ynumber of pinhole parameter - %ld is out of range (0, 500).\n", pinhole_param->nYP);
    error++;
  }
  if (undulator_param->period==0) {
    fprintf(stderr,"The undulator period is zero, you need provide the period of undulator!\n");
    error++;
  }
  if (undulator_param->nPeriod==0) {
    fprintf(stderr,"The number of undulator periods is zero, you have to provide a positive integer!\n");
    error++;
  }
  if (electron_param->energy==0 || electron_param->current==0) {
    fprintf(stderr, "The energy and current of electron beam can not be zero!");
    error++;
  }
  
  if (pinhole_param->xPC<0 || pinhole_param->yPC<0) {
    fprintf(stderr, "The center (x/y position) of pinhole can not be less than zero, it must lie in the first quadrant.\n");
    error++;
  }
  if (mode==1 && us) {
    if (pinhole_param->nXP<0 || pinhole_param->nYP<0) {
      fprintf(stderr, "The subdivision number of pinhole in X and Y (xnumber, ynumber) can not be less than zero for US flux distribution calculation!\n");
      error++;
    }
    if (pinhole_param->distance==0.0 && pinhole_param->xPC==0.0 && pinhole_param->yPC==0.0) { 
      fprintf(stderr, "Invalid distance, xPC and yPC values given for mode=1, they should not be zero for mode=1 calculation.\n");
      error++;
    }
  }
  if (!inputSupplied && (electron_param->sigmaxp==0.0 || electron_param->sigmayp==0.0)) {
    if (us) {
      if (mode==6 && icalc!=3) {
        fprintf(stderr, "The rms x and y divergence should not be zero for US mode=6, and method=1,4,14 for US calculation)\n");
        error++;
      }
    } else {
      if ((abs(mode)<5 && icalc!=3) || (abs(mode)!=6 && icalc==1)) {
        fprintf(stderr, "The rms x and y divergence  should not be zero for URGENT mode 1,2,3,4 icalc=1,2 or mode=+-6, and method=1 calculation.\n");
        error++;
      }
    }
  }  
  if (error) {
    fprintf(stderr, "Input parameters checking failed.\n");
    exit(1);
  }
}
