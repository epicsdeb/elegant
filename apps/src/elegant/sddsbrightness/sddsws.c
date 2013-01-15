/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/
/* program: sddsus.c
   transfer by Shang from ws.f by Roger Dejus 
   for calculating wiggler and bending magnet spectra using the bessel function approximation.

$Log: not supported by cvs2svn $
Revision 1.14  2012/01/08 21:26:02  borland
Fixed units for Kx and Ky in output file.

Revision 1.13  2011/03/14 15:06:15  shang
fixed the current parameter units which should be mA instead of A

Revision 1.12  2011/03/09 20:51:54  shang
added the description for computing bending magnet spectra.

Revision 1.11  2011/03/09 20:44:36  shang
added description

Revision 1.10  2009/06/05 15:03:21  shang
added warning message for possible inaccurate results.

Revision 1.9  2009/05/01 20:53:50  shang
fixed a bug in compute gk

Revision 1.8  2009/04/30 22:30:34  shang
optimized the code to reduce the pow() calls and improved the computation speed by about one times.

Revision 1.7  2009/04/29 20:16:11  shang
removed printing statements

Revision 1.6  2009/04/29 20:06:32  shang
to be consistent with ws, assume it is bending magnet when the number of periods is set to 0.5

Revision 1.5  2009/04/29 19:20:40  shang
fixed another bug in computing the pinhole ydelta when the start position is not at zero  that xsize instead of ysize was used.

Revision 1.4  2009/04/29 19:02:03  shang
fixed a typo in syntax text and fixed a bug in allocating memory for cy.

Revision 1.3  2009/04/21 20:46:16  shang
removed qromb8 statement

Revision 1.2  2009/04/21 20:35:52  shang
removed all numerical recipe routines and now uses qromb() for integration, which was converted from netlib.org; and removed the inputfile option.

Revision 1.1  2009/04/21 15:23:02  shang
first version, converted from ws.f by Roger Dejus

*/

#include "scan.h"
#include "SDDS.h"
#include "mdb.h"
#include "ws_constants.h"

#define CLO_ELECTRON_BEAM 0
#define CLO_PHOTON_ENERGY 1
#define CLO_MAGNET 2
#define CLO_PHIN_HOLE 3
#define CLO_CALCULATION_MODE 4
#define CLO_PIPE 5
#define CLO_NO_WARNINGS 6
#define N_OPTIONS 7

static char *option[N_OPTIONS] = {
  "electronBeam", "photonEnergy", "magnet", "pinhole", "calculationMode", "pipe", "nowarnings" };

#define BENDING_MAGNET_FLAG 0x0001

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

char *USAGE="sddsws <outputFile> [-pipe[=out]] [-nowarnings] \n\
     [-electronBeam=current=<value>(mA),energy=<value>(GeV)] \n\
     [-photoEnergy=maximum=<value>(eV),minimum=<value>(eV),points=<number>] \n\
     [-magnet=period=<value>(cm),<numberOfPeriods>=<number>,kx=<value>,ky=<valu>[,bendingMagnet]] \n\
     [-pinhole=distance=<value>,xposition=<value>,yposition=<value>,xsize=<value>,ysize=<value>,xnumber=<integer>,ynumber=<integer>]\n\
     [-calculationMode=fluxDistribution|fluxSpectrum|brightness|pinholeSpectrum|integratedSpectrum|powerDensity|1-6] \n\
electronBeam     Specifies the electron beam (storage ring) parameters: \n\
                 current  electron beam current in mA. (default is 100mA). \n\
                 energy   electron energy in Gev. (default is 7.0Gev).\n\
photonEnergy     specifies the maximum and minimum photon energy in eV, \n\
                 and the number of energy points to be computed.\n\
magnet           specifies the wigger or bending magnets parameters, \n\
                 period, in cm units \n\
                 number of periods, ignored if bendingMagnet specified.\n\
                 kx, horizontal-plane (vertical field) deflection parameter (ky)\n\
                 ky, vertical-plane (horizontal field) deflection parameter (kx)\n\
                 currently kx!=0 not available (force kx=0) \n\
                 bendingMagnet (optional), specifies that it is bendingMagnet, otherwise, wiggler.\n\
pinhole          Specifies pinhole parameters: \n\
                 distance   distances from the source (m) \n\
                            (if distance=0.0 use angular units) \n\
                 xposition  X-coordinate for center of pinhole (mm) or (mrad) \n\
                 yposition  Y-coordinate for center of pinhole (mm) or (mrad) \n\
                 xsize      X-size of pinhole (full width)	(mm) or (mrad) \n\
                 ysize      y-size of pinhole (full width)	(mm) or (mrad) \n\
                 (for angular units (d=0.0) values are entered in mrad) \n\
                 xnumber    Number of subdivisions of pinhole in X (max 50) \n\
                 ynumber    Number of subdivisions of pinhole in Y (max 50) \n\
                 pinhole parameters are not needed for computing on-axis brilliance (i.e., mode=3). \n\
calculationMode specifies calculation mode: \n\
                1 | fluxDistribution:        Angular/spatial flux density distribution \n\
                                             Flux distribution at the energy chosen as minimum energy. \n\
                2 | fluxSpectrum:            Angular/spatial flux density spectrum \n\
                                             Spectrum at any given point in space as selected by the X and Y \n\
                                             coordinate for the center of the pinhole. X is horizontal and Y is vertical.\n\
                3 | brightness | brilliance: On-axis brilliance spectrum (not implemented) \n\
                4 | pinholeSpectrum:         Flux spectrum through a pinhole \n\
                                             Spectrum through a pinhole centered at X-center and Y-center with \n\
                                             size X-size and Y-size.  The energy range is from the minimum to the \n\
                                             maximum energy. \n\
                5 | integratedSpectrum:      Flux spectrum integrated over all angles \n\
                                             The pinhole parameters have no significance here. \n\
                6 | powerDensity:            Power density and integrated power \n\
                                             Integrated over all energies, thus the energy parameters have no significance here.\n\
transformed from ws.f by Roger Dejus, Hairong Shang (April, 2009).\n\n\
NOTE: THE POLARIZATION PARAMETERS ARE PROVIDED (P1, P2, P3 and P4 in the output file.) \n\
ALTHOUGH NOT THOROUGHLY TESTED - USE WITH CAUTION.\n\
sddsws calculates wiggler spectra using the Bessel function approximation. \n\
The input parameters are divided into sections related to the storage ring, the wiggler device, and the quantity to be calculated.\n\n\
Note: For a bending magnet source: set N=0.5, and make Ky large and adjust the period length accordingly. \n\
For example, put Ky=9.34 and calculate the period length from, Period (cm) = 10.0/B0(T), where B0 is the known\n\
strength of the magnetic field (in Tesla) for the bending magnet.  The calculated power density (pd) is \n\
correct, but the total power (ptot) is irrelevant. Typically make the extend of the pinhole small in the\n\
horizontal direction (theta << Ky/gamma) as the intensity should not depend on the horizontal offset.\n\n";
         
void checkWSInput(long mode, double *xpc, double *ypc, double xsize, double ysize, long nE, double kx, double ky, 
                  long bendingMagnet, long *isAngular, double *pdistance,
                  double *xps, double *yps, long *nxp, long *nyp, double emin, double emax);

void compute_constants(long nE, long nxp, long nyp, double nPeriod, 
                       double energy, double current, double kx, double ky, double period, double pdistance,
                       double emax, double emin, double xpc, double ypc, double xsize, double ysize);

void space_distribution(long mode, long bendingMagnet, long nxp, long nyp, long nE, double nPeriod, double kx, double ky, long isAngular, double emin,
                        double **xpp, double **ypp, double **irradiance, 
                        double **p1, double **p2, double **p3, double **p4, double *flux, double *power);
void SetupOutputFile(char *filename, SDDS_DATASET *SDDSout, long mode, long isAngular);
double trapz2(double *ra, long nxp, long nyp);
void compute_irradiance(long nxp, long nyp, long bendingMagnet,  double kx, double ky, double photonE, long mode, double *ra0, double *ra1,  double *ra3);
void spectral_distribution(long mode, long nE, long nxp, long nyp, long bendingMagnet, double kx, double ky, 
                           double **irradiance, double **p1, double **p2, double **p3, double **p4, double *flux, double *power);
void angle_integration(long lopt, long nE, 
                       double **irradiance, double **p1, double **p2, double **p3, double **p4, double *flux, double *power) ;
void wiggler_power_distribution(long nxp, long nyp, long nE, double **xpp, double **ypp,
                                double **irradiance, double **p1, double **p2, double **p3, double **p4, double *power);
void bendingMagnet_power_distribution(long nxp, long nyp, long nE, double ky, double **xpp, double **ypp,
                                      double **irradiance, double **p1, double **p2, double **p3, double **p4, double *power);

double fkh_integrand(double alpha);
double fkv_integrand(double alpha);

/*global variables */
double gamma_, g2, lamdar, er, k2, k3, lamda1, e1z, d2, len, gk, k_magnet, totalPower, totalPowerDensity, b0, ec0, psi0, facs, faca, facp, facp1, fac, dE, dxp, dyp;
double *xp, *yp, *photonE, *cx, *cy;
void fk(double xg, double yg, double k_magnet, double *s0, double *s1, double *s2, double *s3);
void fkl(double xg, double yg, double k_magnet, double *s0, double *s1, double *s2, double *s3);

typedef struct {
  double xc, yc, kc;
} INTEG_PARAMETER; /*the parameters used in the integration function*/

INTEG_PARAMETER integ_par;

int main(int argc, char **argv)
{
  char  *outputfile=NULL;
  double energy=7.0, current=100.0, emax=-1, emin=-1, xpc=0, ypc=0, xsize=-1, ysize=-1, kx=0, ky=-1, period=-1, pdistance=0;
  long mode=-1, nxp=20, nyp=20, nE=500, bendingMagnet=0, i_arg, nowarnings=0, mode_index=-1, isAngular=0, total_rows=0, index=0;
  double nPeriod;
  SDDS_DATASET SDDSout;
  SCANNED_ARG *s_arg;
  unsigned long dummyFlags=0, pipeFlags=0;
  
  double *xpp, *ypp, *irradiance, *p1, *p2, *p3, *p4, flux, power;
  char desc[2048];
  
  photonE = xp = yp = cx = cy = NULL;
  xpp = ypp = irradiance = p1 = p2 = p3 = p4 = NULL;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s", USAGE);
    exit(1);
  }
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case CLO_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        if (!(pipeFlags&USE_STDOUT))
          SDDS_Bomb("invalid -pipe flag provided.");
        break;
      case CLO_ELECTRON_BEAM:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -electronBeam syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "current", SDDS_DOUBLE, &current, 1, 0,
                          "energy", SDDS_DOUBLE, &energy, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -electronbeam syntax");
        s_arg[i_arg].n_items++;
        break;
      case CLO_PHOTON_ENERGY:
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
      case CLO_MAGNET:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -undulator syntax.");
        s_arg[i_arg].n_items--;
        dummyFlags = 0;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "period", SDDS_DOUBLE, &period, 1, 0,
                          "numberOfPeriods", SDDS_DOUBLE, &nPeriod, 1, 0,
                          "kx", SDDS_DOUBLE, &kx, 1, 0,
                          "ky", SDDS_DOUBLE, &ky, 1, 0,
                          "bendingMagnet", -1, NULL, 0, BENDING_MAGNET_FLAG, 
                          NULL))
          SDDS_Bomb("invalid -undulator syntax");
        s_arg[i_arg].n_items++;
        if (dummyFlags & BENDING_MAGNET_FLAG)
          bendingMagnet = 1;
        if (nPeriod==0.5)
          bendingMagnet = 1;
        if (bendingMagnet)
          nPeriod = 0.5;
        break;
      case CLO_PHIN_HOLE:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -pinhole syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "distance", SDDS_DOUBLE, &pdistance, 1, 0,
                          "xposition", SDDS_DOUBLE, &xpc, 1, 0,
                          "yposition", SDDS_DOUBLE, &ypc, 1, 0,
                          "xsize", SDDS_DOUBLE, &xsize, 1, 0,
                          "ysize", SDDS_DOUBLE, &ysize, 1, 0,
                          "xnumber", SDDS_LONG, &nxp, 1, 0,
                          "ynumber", SDDS_LONG, &nyp, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -pinhole syntax.");
        s_arg[i_arg].n_items++;
        break;
      case CLO_CALCULATION_MODE:
        if (!get_long(&mode, s_arg[i_arg].list[1])) {
          if ((mode_index=match_string(s_arg[i_arg].list[1], mode_options, CLO_MODES, 0))==-1) {
            fprintf(stderr, "Unknow mode - %s provided.\n", s_arg[i_arg].list[1]);
            exit(1);
          }
          if (mode_index<3)
            mode = mode_index+1;
          else  
            mode = mode_index;
        } else {
          if (mode<3)
            mode_index = mode-1;
          else
            mode_index = mode;
        }
        break;
      case CLO_NO_WARNINGS:
        nowarnings = 1;
        break;
      default:
        fprintf(stderr, "Invalid option - %s provided.\n",  s_arg[i_arg].list[1]);
        exit(1);
      }
    } else {
      if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  if (!outputfile && !pipeFlags)
    SDDS_Bomb("output file not provided.");
  if (outputfile && pipeFlags) 
    SDDS_Bomb("Too many files provided.");
  if (!nowarnings) {
    fprintf(stdout, "Warning: the polarization parameters (P1, P2, P3 and P4 in the output file) are provided although not thoroughly tested, use with caution.\n");
    if (bendingMagnet) 
      fprintf(stdout, "Warning: the result may be inaccurate for bending magnet, use with caution.\n");
  }
  checkWSInput(mode, &xpc, &ypc, xsize, ysize, nE, kx, ky, bendingMagnet, &isAngular, &pdistance, &xsize, &ysize, &nxp, &nyp, emin, emax);
  compute_constants(nE, nxp, nyp, nPeriod, energy, current, kx, ky, period, pdistance, emax, emin, xpc, ypc, xsize, ysize);
  
  nxp = nxp +1;
  nyp = nyp +1;
  SetupOutputFile(outputfile, &SDDSout, mode, isAngular);
  
  switch (mode_index) {
  case CLO_FLUX_DISTRIBUTION_MODE:
    /*mode = 1*/
    space_distribution(mode, bendingMagnet, nxp, nyp, nE, nPeriod, kx, ky, isAngular, emin,
                       &xpp, &ypp, &irradiance, &p1, &p2, &p3, &p4, &flux, &power);
    if (isAngular) 
      sprintf(desc, "Angular flux density distribution for %lf eV (ph/s/mr^2/0.1%bw)", emin);
    else
      sprintf(desc, "Irradiance for %lf eV @ %f m (ph/s/mm^2/0.1%bw)", emin, pdistance);
    total_rows = nxp * nyp; 
    break;
  case CLO_FLUX_SPECTRUM_MODE: /* mode =2 */
  case CLO_BRIGHTNESS_MODE: /*mode =3 */
  case CLO_BRILLIANCE_MODE: /*mode = 3 */
  case CLO_PINHOLE_SPECTRUM_MODE: /* mode =4 */
    /* isub = 2 */
    spectral_distribution(mode, nE, nxp, nyp, bendingMagnet, kx, ky, &irradiance, &p1, &p2, &p3, &p4, &flux, &power);
    if (mode==2) {
      if (isAngular)
        sprintf(desc, "Angular flux density spectrum for %f mr, %f mr (ph/s/mr^2/0.1%bw)", xpc, ypc);
      else
        sprintf(desc, "Irradiance for %f mm, %f mm at %f m W/mr^2", xpc, ypc, pdistance);
    } else if (mode==3) {
      sprintf(desc, "on axis brilliance.");
    } else if (mode==4) {
      if (isAngular)
        sprintf(desc, "Flux through %f mr, %f mr pinhole at %f mr, %f mr.", xsize, ysize, xpc, ypc);
      else
        sprintf(desc, "Flux through %f mm, %f mm pinhole at %f mm, %f mm @ %f m.", xsize, ysize, xpc, ypc, pdistance);
    }       
    total_rows = nE;
    break;
  case  CLO_INTEGRATED_SPECTRUM_MODE:
    angle_integration(2, nE, &irradiance, &p1, &p2, &p3, &p4, &flux, &power);
    total_rows = nE;
    sprintf(desc, "Agnule-integrated spectrum.");
    /*mode = 5, isub=3*/
    break;
  case CLO_POWER_DENSITY_MODE:
    if (isAngular)
      sprintf(desc, "Power density distribution (W/mr^2)");
    else
      sprintf(desc, "Power density distribution @ %f m.", pdistance);
    /* mode =6 */
    if (bendingMagnet) {
      /*isub = 5*/
      bendingMagnet_power_distribution(nxp, nyp, nE, ky, &xpp, &ypp, &irradiance, &p1, &p2, &p3, &p4, &power);
    } else {
      /* isub = 4 */
      wiggler_power_distribution(nxp, nyp, nE, &xpp, &ypp, &irradiance, &p1, &p2, &p3, &p4, &power);
    }
    total_rows = nxp*nyp;
    break;
  }
  if (!SDDS_StartPage(&SDDSout, total_rows) ||
      !SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Description", desc, "eBeamEnergy", energy, "eBeamCurrent", current,
                          "MinEnergy", emin, "MaxEnergy", emax, "EnergyStep", nE,
                          "Device", bendingMagnet ? "bending magnet" : "wiggler",
                          "Period", period, "NPeriod", nPeriod, "Kx", kx, "Ky", ky,
                          "PinholeDistance", pdistance, "PinholeXPos", xpc, "PinholeYPos", ypc, 
                          "PinholeXSize", xsize, "PinholeYSize", ysize, "PinholeXPoints", nxp, "PinholeYPoints", nyp,
                          "b0", b0, "ec0", ec0/1.0e3, "e1", e1z, "lamda1", lamda1, 
                          "TotalPowerDensity", isAngular ? totalPowerDensity : totalPowerDensity/d2, "TotalPower", totalPower, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);  
  switch (mode) {
  case 1:
  case 6:
    if (mode==6) {
      if (!SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, "IntegratedPower", power, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else {
      if (!SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, "IntegratedPower", power, "IntegratedFlux", flux, NULL))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, xpp, total_rows, 0) ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, ypp, total_rows, 1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors); 
    index = 2; /* index for next column*/
    break;
  case 2:
    if (!SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, "IntegratedPowerDensity", power, "IntegratedFluxDensity", flux, NULL) ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, photonE, total_rows, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors); 
    index = 1;
    break;
  case 3:
    if (!SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, photonE, total_rows, 0) ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, irradiance, total_rows, 1) ||
        !SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    index = 0;
    break;
  case 4:
  case 5:
    if (!SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE, "IntegratedPower", power, "IntegratedFlux", flux, NULL) ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, photonE, total_rows, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    index = 1;
    break;
  }
  if (mode!=3) {
    if ( !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, irradiance, total_rows, index) ||
         !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, p1, total_rows, index+1) ||
         !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, p2, total_rows, index+2) ||
         !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, p3, total_rows, index+3) ||
         !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_INDEX, p4, total_rows, index+4) ||
         !SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  free(photonE);
  free(xp);
  free(yp);
  free(cx);
  free(cy);
  if (xpp) free(xpp);
  if (ypp) free(ypp);
  if (irradiance) free(irradiance);
  if (p1) free(p1);
  if (p2) free(p2);
  if (p3) free(p3);
  if (p4) free(p4);
  free_scanargs(&s_arg, argc);
  return 0;
}

void SetupOutputFile(char *filename, SDDS_DATASET *SDDSout, long mode, long isAngular)
{
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(SDDSout, "Description", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "eBeamEnergy", NULL, "Gev", "Electron Beam Energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "eBeamCurrent", NULL, "mA", "Electron Beam Current", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MinEnergy", NULL, "eV", "minimum photon energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MaxEnergy", NULL, "eV", "maximum photon energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "EnergyStep", NULL, NULL, "number of points of photon energy", NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Device", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Period", NULL, "cm", "Period length of wiggler or undulator", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "NPeriod", NULL, NULL, "number of periods of wiggler", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Kx", NULL, NULL, "K value of horizontal field", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Ky", NULL, NULL, "K value of vertical field", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeDistance", NULL, "m", "pinhole distance from the source", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeXPos", NULL, isAngular ? "mrad": "mm",
                           "X-coordinate for center of pinhole (mm) or (mrad)", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeYPos", NULL, isAngular ? "mrad": "mm", 
                           "X-coordinate for center of pinhole (mm) or (mrad)", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeXSize", NULL, isAngular ? "mrad": "mm", 
                           "x-size of pinhole (full width)", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeYSize", NULL, isAngular ? "mrad": "mm", 
                           "y-size of pinhole (full width)", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeXPoints", NULL, NULL, 
                           "Number of subdivisions of pinhole in X", NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "PinholeYPoints", NULL, NULL, 
                           "Number of subdivisions of pinhole in Y", NULL, SDDS_LONG, 0)<0  ||
      SDDS_DefineParameter(SDDSout, "b0", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "ec0", NULL, "keV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "e1", NULL, "keV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "lamda1", NULL, "A", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "TotalPowerDensity",  NULL, 
                           "Watts/mm$a2$n or mrad$a2$n", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "TotalPower",  NULL, "Watts", NULL, NULL, SDDS_DOUBLE, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  switch (mode) {
  case 1:
    if (SDDS_DefineParameter(SDDSout, "IntegratedPower",  NULL, "W/eV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "IntegratedFlux",  NULL, 
                             "ph/s/0.1%bw", NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_DefineColumn(SDDSout, "X", NULL, isAngular ? "mrad" : "mm", 
                          "horizontal position of pinhole",NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "Y", NULL, isAngular ? "mrad" : "mm", 
                          "vertical position of pinhole",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (isAngular) {
      if (SDDS_DefineColumn(SDDSout, "AngularFluxDensity", NULL, "photons/s/mrad$a2$n/0.1%BW", "angular flux density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else if (SDDS_DefineColumn(SDDSout, "Irradiance", NULL,"photons/s/mm$a2$n/0.1%BW", "spatial flux density",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    break;
  case 2:
    if (SDDS_DefineParameter(SDDSout, "IntegratedPowerDensity",  NULL, isAngular ? "Watts/mrad$a2$n" : "Watts/mm$a2$n", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "IntegratedFluxDensity",  NULL,  isAngular ? "ph/s/mr$a2$n" : "ph/s/mm$a2$n", NULL, NULL, SDDS_DOUBLE, 0))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", "photon energy",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (isAngular) {
      if (SDDS_DefineColumn(SDDSout, "AngularFluxDensity", NULL, "photons/s/mrad$a2$n/0.1%BW", "angular flux density",NULL, SDDS_DOUBLE, 0)<0)
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } else if (SDDS_DefineColumn(SDDSout, "Irradiance", NULL, "photons/s/mm$a2$n/0.1%BW", "spatial flux density",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    break;
  case 3:
    if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", "photon energy",NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "Brightness", NULL, "ph/s/mr$a2$n/mm$a2$n/0.1%bw", 
                          "on axis brightness",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    break;
  case 4:
  case 5:
    if (SDDS_DefineParameter(SDDSout, "IntegratedPower",  NULL, "W", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineParameter(SDDSout, "IntegratedFlux",  NULL, 
                             "ph/s", NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", "photon energy",NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "Flux", NULL, "ph/s/0.1%bw", 
                          "flux through pinhole",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    break;
  case 6:
    if (SDDS_DefineParameter(SDDSout, "IntegratedPower",  NULL, "W", NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_DefineColumn(SDDSout, "X", NULL, isAngular ? "mrad" : "mm", 
                          "horizontal position of pinhole",NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "Y", NULL, isAngular ? "mrad" : "mm", 
                          "vertical position of pinhole",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (SDDS_DefineColumn(SDDSout, "PowerDensity", NULL, isAngular ? "Watts/mrad$a2$n" : "Watts/mm$a2$n", "power density distribution",NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    break;
  }
  if (mode!=3) {
    if (SDDS_DefineColumn(SDDSout, "P1", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "P2", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "P3", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(SDDSout, "P4", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0)
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

void checkWSInput(long mode, double *xpc, double *ypc, double xsize, double ysize, long nE, double kx, double ky, 
                  long bendingMagnet, long *isAngular, double *pdistance,
                  double *xps, double *yps, long *nxp, long *nyp, double emin, double emax)
{
  /* if (mode==3) {
    fprintf(stderr, "On-axis brightness spectrum (mode=3) not implemented in ws!\n");
    exit(1);
    } */
  if (mode<1) {
    fprintf(stderr, "Error(sddsws): the computation mode is not provided.\n");
    exit(1);
  }
  if (*xpc<0 || *ypc<0) {
    fprintf(stderr, "the pinhole position must be in the first quadrant, xposition and y position can not be less than 0.\n");
    exit(1);
  }
  if (xsize<0 || ysize<0) {
    fprintf(stderr, "the pinhole x/y size must be greater than 0.\n");
    exit(1);
  }
  if (nE>10000) {
    fprintf(stderr, "Enery array (nE) out of bounds; number of points %ld is greater than allowed (10000).\n", nE);
    exit(1);
  }
  if (kx !=0) {
    fprintf(stderr, "kx!=0 is not allowed, only kx=0 allowed.\n");
    exit(1);
  }
  if (bendingMagnet && mode==5) {
    fprintf(stderr, "Flux spectrum integrated over all angles (mode) for bending magnets not available.\n");
    exit(1);
  }
  if (*pdistance==0) {
    *isAngular=1;
    *pdistance = 1.0;
  } else {
    *isAngular=0;
  }
  if (*pdistance==0) {
    fprintf(stderr, "The pinhole distance can not be zero.\n");
    exit(1);
  }
  if (mode==1) {
    if (emin<0) {
      fprintf(stderr, "The minimum photon energy for mode=1 has to be provided.\n");
      exit(1);
    }
  }

  if (mode==2) {
    *xps = *yps = 0;
    *nxp = *nyp = 0;
  } else if (mode==3 || mode==5) {
    *xpc = *ypc = *xps = *yps = 0;
    *nxp = *nyp = 0;
  }

  
}

void compute_constants(long nE, long nxp, long nyp, double nPeriod, 
                       double energy, double current, double kx, double ky, double period, double pdistance,
                       double emax, double emin, double xpc, double ypc, double xsize, double ysize)
{
  long i;
  double xpmin, ypmin;
  
  gamma_ = energy/me_mev*1.0e3;
  g2 = gamma_ * gamma_;
  lamdar = period * C_CM_ANG/(2.0*g2); /*reduced wavelength [A]*/
  er = C_EVANG/lamdar;  /*reduced energy [eV] */
  k2 = kx * kx + ky * ky; 
  k3 = 1.0 + k2/2.0;
  lamda1 = lamdar*k3;   /*First harmonic on axis [A] */
  e1z = er / k3 ;  /* first harmonic on axis [eV] */
  d2 = pdistance * pdistance;  /* distance squared [m^2] */
  len = nPeriod * period * C_CM_M; /* length of device  [m] */
  gk = 0;
  if (kx< EPSK || ky < EPSK) {
    k_magnet = kx + ky;
    gk = k_magnet*(ipow(k_magnet, 6) + 24.0 * ipow(k_magnet,4)/7.0 + 4.0 * k_magnet * k_magnet + 16.0/7.0)/pow(1.0 + k_magnet*k_magnet, 3.5);
  } else if (abs(kx-ky)< EPSK) {
    k_magnet = kx;
    gk = 32.0/7.0*k_magnet/pow(1.0 + k_magnet*k_magnet, 3.5);
  }
  totalPower = PTOT_FAC * nPeriod * k2 * (energy*energy) * current * C_MA_A / (period * C_CM_M); /* units [W] */
  totalPowerDensity = PD_FAC * nPeriod * k_magnet * gk * energy * energy * energy * energy * current * C_MA_A/(period * C_CM_M);  /* units [W]/mrad^2 */
  b0 = ky/(CL1 * period); 
  ec0 = CL2 * energy * energy * b0; /*[eV]*/
  
  psi0 = kx / gamma_;
  facs = CL3 * g2 * BW * current/e_mks*C_MA_A/(C_RAD_MRAD*C_RAD_MRAD) * 2.0 * nPeriod / d2; /* fac spectral_distribution */
  faca = CL4 * 2.0 * k_magnet * BW * current/e_mks*C_MA_A * 2.0 * nPeriod; /* fac angle-integrated spectrum, first two is for hor. symmetry  W/mm^2*/
  facp = totalPowerDensity/d2;
  facp1 = facs/BW*e_mks; /* conversion to W/mm^2 or W/mrad^2 */
  
  dE = (emax - emin)/(nE-1);
  dxp = dyp = 0;
  if (xpc==0 && ypc==0) {
    fac = 4.0;
    xpmin = ypmin = 0;
    if (nxp>0) dxp=xsize/2.0/nxp;
    if (nyp>0) dyp=ysize/2.0/nyp;
  } else {
    fac = 1.0;
    xpmin = xpc - xsize/2.0;
    ypmin = ypc - ysize/2.0;
    if (nxp>0) dxp = xsize/nxp;
    if (nyp>0) dyp = ysize/nyp;
  }
  photonE = calloc(sizeof(*photonE), nE);
  for (i=0; i<nE; i++)
    photonE[i] = emin + i*dE;
  xp = calloc(sizeof(*xp), nxp+1);
  yp = calloc(sizeof(*yp), nyp+1);
  cx = calloc(sizeof(*cx), nxp+1);
  cy = calloc(sizeof(*cy), nyp+1);
  
  for (i=0; i<nxp+1; i++) {
    xp[i] = xpmin + i*dxp;  /* poistion [mm] */
    cx[i] = xp[i] * C_MM_M/pdistance; /*angle [rad]*/
  }
  for (i=0; i<nyp+1; i++) {
    yp[i] = ypmin + i*dyp; /* poistion [mm] */
    cy[i] = yp[i] * C_MM_M/pdistance; /*angle [rad]*/
  }
  
}

/*p_e is the photon energy */
void compute_irradiance(long nxp, long nyp, long bendingMagnet, double kx, double ky, double p_e, long mode, double *ra0, double *ra1,  double *ra3)
{
  long i, j, k;
  double  xg, yg, vp, vpmax, ecp, yg2, yg1, cc, cpi, cpis, y2, eta, y, asigma, api, fc, e1;
  
  vpmax = sqrt(1.0 - ECP_FRAC);
  for (j=0; j<nyp; j++) {
    yg = gamma_ * cy[j];
    yg2 = yg * yg;
    yg1 = 1.0 + yg2;
    cc = yg1 * yg1;
    cpi = yg2/yg1;
    if (yg>0) {
      cpis = sqrt(cpi);
    } else {
      cpis = -sqrt(cpi);
    }
    fc = facs*cc;
    e1 =  pow(yg1, 1.5);
    for (i=0; i<nxp; i++) {
      xg = gamma_ * cx[i];
      k = i * nyp + j;
      vp = fabs(xg)/ky;
      if (vp>vpmax) 
        continue;
      ecp = ec0 * sqrt(1.0-vp*vp); /* eV*/
      y = p_e/ecp;
      y2 = y*y;
      eta = 0.5* y * e1;
      /*
        asigma = dbeskv_nu(eta, 2.0/3.0);
        api = dbeskv_nu(eta, 1.0/3.0); 
        use roger's k13, and k23, seems to be faster; they were converted into c in mdbmth */
      asigma = k23(eta);
      api = k13(eta); 
      /*  fprintf(stdout, "eta=%f, asigma=%e, api=%e\n", eta, asigma, api);*/
      ra0[k] = y2*fc*(asigma*asigma + cpi*api*api);
      ra1[k] = y2*fc*(asigma*asigma - cpi*api*api);
      if (bendingMagnet) 
        ra3[k] = 2.0 * y2 * fc * (asigma * cpis * api);
      if (mode==2 || mode==3)
        break;
    }
  }
  return;
}
void space_distribution(long mode, long bendingMagnet, long nxp, long nyp, long nE, double nPeriod, double kx, double ky, long isAngular, double emin,
                        double **xpp, double **ypp, double **irradiance, 
                        double **p1, double **p2, double **p3, double **p4, double *flux, double *power)
{
  long i, j, k;
  double *ra1, *ra3, *ra0, a1, a2, a3, ra2;

  *xpp = calloc(sizeof(**xpp), nxp*nyp);
  *ypp = calloc(sizeof(**ypp), nxp*nyp);
 
  *p1 = calloc(sizeof(**p1), nxp*nyp);
  *p2 = calloc(sizeof(**p2), nxp*nyp);
  *p3 = calloc(sizeof(**p3), nxp*nyp);
  *p4 = calloc(sizeof(**p4), nxp*nyp);
  ra0 = ra1 = ra3 = NULL;
  
  ra0 = calloc(sizeof(*ra0), nxp * nyp);
  ra1 = calloc(sizeof(*ra1), nxp * nyp);
  ra3 = calloc(sizeof(*ra3), nxp * nyp);
  
  compute_irradiance(nxp, nyp, bendingMagnet,kx, ky, emin, 1, ra0, ra1, ra3);
  a2 = 0;
  ra2 = 0;
  for (i=0; i<nxp; i++) {
    for (j=0; j<nyp; j++) {
      k = i*nyp +j;
      (*xpp)[k] = xp[i];
      (*ypp)[k] = yp[j];
      if (ra0[k]>0) {
        (*p1)[k] = a1 = ra1[k]/ra0[k];
        (*p3)[k] = a3 = ra3[k]/ra0[k];
        (*p4)[k] = 1.0 - sqrt(a1*a1 + a2 * a2 + a3*a3);
      }
    }
  }
  *irradiance = ra0;
  *flux = trapz2(ra0, nxp, nyp); /*get integrated flux over observation area*/
  free(ra1);
  free(ra3);
  *flux = (*flux)*fac; /*ph/s/0.1%bw */
  *power = (*flux)/BW*e_mks; /*  W/eV */
}

double trapz2(double *ra, long nxp, long nyp)
{
  /* get integrated flux over observation  area */
  long i, j, k;
  double sum,wx,wy, area;
  
  if (!ra) {
    fprintf(stderr, "memory not allocated for irradiance.\n");
    exit(1);
  }
  sum = 0;
  for (i=0; i<nxp; i++) {
    if (i==0 || i==nxp-1)
      wx = 0.5;
    else
      wx = 1.0;
    for (j=0; j<nyp; j++) {
      k = i*nyp + j;
      if (j==0 || j==nyp-1) 
        wy = 0.5;
      else 
        wy = 1.0;
      sum += wx*wy*ra[k];
    }
  }
  area = sum * dxp * dyp;
  return area;
}

/*spatial flux density for 2<=mode<=4 */
void spectral_distribution(long mode, long nE,  long nxp, long nyp, long bendingMagnet, double kx, double ky, 
                           double **irradiance,  double **p1, double **p2, double **p3, double **p4, double *flux, double *power)
{
  /*units ph/s/mr^2.0.1%bw for angular flux density and ph/s/mm^2/0.1%bw for spatial flux density */
  double *ra0=NULL, *ra1=NULL, *ra3=NULL, area0, area1, area3=0, we=0;
  double spec0=0, spec1=0, spec3=0, a1, a3;
  long ie;
  
  
  *p1 = calloc(sizeof(**p1), nE);
  *p2 = calloc(sizeof(**p2), nE);
  *p3 = calloc(sizeof(**p3), nE);
  *p4 = calloc(sizeof(**p4), nE);
  *irradiance = calloc(sizeof(**irradiance), nE);
 
  *flux = *power =  0;
  ra0 = calloc(sizeof(*ra0), nxp * nyp);
  ra1 = calloc(sizeof(*ra1), nxp * nyp);
  ra3 = calloc(sizeof(*ra3), nxp * nyp);
  for (ie=0; ie<nE; ie++) {
    
    if (ie==0 || ie==nE-1)
      we = 0.5;
    else
      we = 1.0;
    compute_irradiance(nxp, nyp, bendingMagnet, kx, ky, photonE[ie], mode, ra0, ra1, ra3);
    if (mode==4) {
       /*fixed position */
      area0 = trapz2(ra0, nxp, nyp);
      area1 = trapz2(ra1, nxp, nyp);
      if (bendingMagnet)
        area3 = trapz2(ra3, nxp, nyp);
      spec0 = fac*area0;
      spec1 = fac*area1;
      spec3 = fac*area3;
    } else if (mode==2 || mode==3) {
      /*fixed position */
      spec0 = ra0[0];
      spec1 = ra1[0];
      spec3 = ra3[0];
    }
    if (spec0>0) {
      (*p1)[ie] = a1 = spec1/spec0;
      (*p3)[ie] = a3 = spec3/spec0;
      (*p4)[ie] = 1.0 - sqrt(a1*a1 + a3*a3);
    }
    (*irradiance)[ie] = spec0;
    *flux = *flux + we * spec0/photonE[ie];
    *power = *power + we * spec0;
  }
  free(ra0);
  free(ra1);
  free(ra3);
  *flux = *flux/BW * dE; /*  ph/s or ph/s/mm^2 or ph/s/mrad^2 */
  *power = *power/BW * dE * e_mks; /* W or W/mm^2 or W/mrad^2*/
}

void angle_integration(long lopt, long nE, 
                       double **irradiance, double **p1, double **p2, double **p3, double **p4, double *flux, double *power) 
{
  double vp, dvs, ec02, dec, *wgt, spec0, spec3, vn, ecn, ecp, *ecpa, we, dvn, sumx;
  long nxa = 100; /* Number of steps for integration over horizontal angle */
  long ie, i;
  double  g1, y1;
 
  /* c  iopt = 1 will emphasize accuracy for energies below ~ Ec/3, and iopt = 3
     c  for energies above ~ 3xEc, and iopt = 2 around Ec. In all cases, the
     c  accuracy is better than 0.5 %. */
  lopt = 2;

  vp = 0;
  dvs = 0;
  ecp = ec0;
  ec02 = ec0*ec0;
  
  ecpa = calloc(sizeof(*ecpa), nxa+1);
  wgt = calloc(sizeof(*wgt), nxa+1);
  
  ecpa[0] = ec0;
  
  switch (lopt) {
  case 1: /*Option #1: Constant stepsize in critical energy*/
    dec = ec0 * (1.0 - EC1_FRAC)/(nxa - 1);    
    ecp = ec0;
    break;
  case 2: /*Option #2: Constant stepsize in the square of the critical energy */
    dec = ec02*(1.0 -EC2_FRAC)/(nxa-1);
    ecp = ec02;
    break;
  case 3: /* Option #3: Constant stepsize in horizontal angle */
    dec = sqrt(1.0 -EC3_FRAC)/(nxa-1);
    break;
  default:
    fprintf(stderr, "Invaid option provided for anlge_integration, can only be 1, 2 or 3.\n");
    exit(1);
  }
  
  *irradiance= calloc(sizeof(**irradiance), nE);
  *p1 = calloc(sizeof(**p1), nE);
  *p2 = calloc(sizeof(**p2), nE);
  *p3 = calloc(sizeof(**p3), nE);
  *p4 = calloc(sizeof(**p4), nE);
  fprintf(stderr, "dec=%f\n", dec);
  ecpa[0] = ec0;
  for (i=0; i<nxa; i++) {
    if (lopt==1 || lopt==2) {
      if (i<nxa-1) {
        ecn = ecp - dec;
        if (lopt==1)
          vn = sqrt(1.0 - ecn*ecn/ec02);
        else
          vn = sqrt(1.0 - ecn/ec02);
        dvn = vn - vp;
        ecp = ecn;
        vp = vn;
        if (lopt==1)
          ecpa[i+1] = ecp;
        else 
          ecpa[i+1] = sqrt(ecp);
      } else {
        dvn = 0;
      }
      wgt[i] = 0.5*(dvs+dvn);
      dvs = dvn;
    } else {
      if (i==0 || i<nxa-1) {
        wgt[i] = 0.5*dec;
      } else
        wgt[i] = dec;
      vp = (i-1) * dec;
      ecp = ec0*sqrt(1.0 - vp*vp);
      ecpa[i] = ecp;
    }                                           
   
  }
  *flux = *power = 0;
  for (ie=0; ie<nE; ie++) {
    sumx = 0;
    if (ie==0 || ie==nE-1)
      we = 0.5;
    else
      we = 1;
    for (i=0; i<nxa; i++) {
      y1 = photonE[ie]/ecpa[i];
      if (y1<=40.0) {
        /* g1 = interpolate(yGy,y,rows, y1, &belowRange, &aboveRange, 1, &interpCode,1)/y1; */
        g1 = gy(1, y1);
        sumx += wgt[i]* g1;    
      } else {
        break;
      }
    }
    spec0 = faca *sumx;
    spec3 = 0;
    (*irradiance)[ie] = spec0;
    (*p4)[ie] = 1.0;
    *flux = *flux + we*spec0/photonE[ie];
    *power = *power + we * spec0;
  }
  *flux  = *flux /BW   *dE;
  *power = *power/BW* e_mks *dE;
  free(wgt);
  free(ecpa);
}

/*	subroutine power_distribution(ierror)
c  Routine for calculation of power density. The routine uses the integral
c  equation given by K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. A246, 
c  (1986) 67-70, Eq. 5, for K < KMAX (typically 100.0).
c  For K > KMAX, Eq. 10 in the same paper is used which gives the power
c  density in the limit K -> infinity. In this limit, there is no
c  radiation beyond gamma*theta/k in the horizontal plane. */
#define KMAX 100.0
void wiggler_power_distribution(long nxp, long nyp, long nE, double **xpp, double **ypp,
                                double **irradiance, double **p1, double **p2, double **p3, double **p4, double *power)
{
  double xg, yg, s0, s1, s2, s3, ra0, ra1, ra2, ra3,a1, a2, a3;
  long ix, iy, k;

  *power = 0;
  *irradiance = calloc(sizeof(**irradiance), nxp * nyp);
  *p1 = calloc(sizeof(**p1), nxp*nyp);
  *p2 = calloc(sizeof(**p2), nxp*nyp);
  *p3 = calloc(sizeof(**p3), nxp*nyp);
  *p4 = calloc(sizeof(**p4), nxp*nyp);
  *xpp = calloc(sizeof(**xpp), nxp*nyp);
  *ypp = calloc(sizeof(**ypp), nxp*nyp);
  
  for (ix=0; ix<nxp; ix++) {
    xg = gamma_ * cx[ix];
    for (iy=0; iy<nyp; iy++) {
      k = ix*nyp + iy;
      (*xpp)[k] = xp[ix];
      (*ypp)[k] = yp[iy];
      yg = gamma_ * cy[iy];
      if (k_magnet<KMAX) 
        fk(xg, yg, k_magnet, &s0, &s1, &s2, &s3);
      else
        fkl(xg, yg, k_magnet, &s0, &s1, &s2, &s3);
      ra0 = facp * s0;
      ra1 = facp * s1;
      ra2 = facp * s2;
      ra3 = facp * s3;
      (*irradiance)[k] = ra0;
      if (ra0>0) {
        (*p1)[k] = a1 = ra1/ra0;
        (*p2)[k] = a2 = ra2/ra0;
        (*p3)[k] = a3 = ra3/ra0;
        (*p4)[k] = 1.0 - sqrt(a1*a1 + a2*a2 + a3*a3);
      }
    }
  }
  *power = fac *  trapz2(*irradiance, nxp, nyp);
}

/*	subroutine power_distribution1(ierror)

c  Routine for calculation of power density. Integrates the Bessel functions over
c  the energy.  The right circular intensity appears above the horizontal plane 
c  (gamma*psi > 0.0) using the current definition. Circular polarized intensity
c  is only valid for the bending magnet and the elliptical multipole wiggler.
C  Size parameters:
	integer*4	E_SZ,P_SZ
	parameter	(E_SZ=10000,P_SZ=201)
c  NETA_SZ was determined such that the power density and integrated power
c  gives good agreement with the results from routine fkl in subroutine
c  power_distribution.  Depends also on the integration range determined
c  by ETA_MIN and ETA_MAX. */
#define NETA_SZ 500
#define ETA_MIN 1.0e-4
#define ETA_MAX 15.0
void bendingMagnet_power_distribution(long nxp, long nyp, long nE, double ky, double **xpp, double **ypp,
                                double **irradiance, double **p1, double **p2, double **p3, double **p4, double *power)
{
  long neta, ie, ix, iy, k;
  double etamin, etamax, deta, vpmax, we, *eta, *asigma, *api, ra0, ra1,ra2, ra3, xg, yg, vp, ecp, yg2, yg1, cc, cpi, cpis, c1, fc, y, y2, a1, a3;
  
  asigma = api  = NULL;
  neta   = NETA_SZ;
  etamin = ETA_MIN;
  etamax = ETA_MAX;
  deta   = (etamax -etamin)/(neta-1);
  vpmax  = sqrt(1.0 -ECP_FRAC);
  
  eta = asigma = api = 0;
  eta = malloc(sizeof(*eta)*neta);
  asigma = malloc(sizeof(*asigma)*neta);
  api = malloc(sizeof(*api) * neta);

  for (ie=0; ie<neta; ie++) {
    eta[ie] = etamin + ie*deta;
    asigma[ie] = k23(eta[ie]); /* Modified Bessel function of 2nd kind (2/3)*/
    api[ie] =  k13(eta[ie]); /* Modified Bessel function of 2nd kind (1/3) */
  }

  *irradiance = calloc(sizeof(**irradiance), nxp * nyp);
  *p1 = calloc(sizeof(**p1), nxp*nyp);
  *p2 = calloc(sizeof(**p2), nxp*nyp);
  *p3 = calloc(sizeof(**p3), nxp*nyp);
  *p4 = calloc(sizeof(**p4), nxp*nyp);
  *xpp = calloc(sizeof(**xpp), nxp*nyp);
  *ypp = calloc(sizeof(**ypp), nxp*nyp);
  *power = 0;
  for (iy=0; iy<nyp; iy++) {
    yg = gamma_ * cy[iy];
    yg2 = yg * yg;
    yg1 = 1.0 + yg2;
    cc = yg1*yg1;
    cpi = yg2/yg1;
    if (yg>0)
      cpis = sqrt(cpi);
    else
      cpis = -sqrt(cpi);
    c1 = 2.0/pow(yg1, 1.5);
    for (ix=0; ix<nxp; ix++) {
      xg = gamma_ * cx[ix];
      k = ix * nyp + iy;
      (*xpp)[k] = xp[ix];
      (*ypp)[k] = yp[iy];
      ra0 = ra1 = ra2 =ra3 =0;
      vp = fabs(xg)/ky;
      if (vp>vpmax)
        continue;
      ecp = ec0 * sqrt(1.0 - vp*vp);
      fc = facp1*cc*deta*c1*ecp;
      for (ie=0; ie<neta; ie++) {
        if (ie==0 || ie==neta-1)
          we = 0.5;
        else 
          we = 1.0;
        y = c1*eta[ie];
        y2 = y*y;
        ra0 += y2*fc*we*(asigma[ie]*asigma[ie] + cpi * api[ie] *api[ie]);
        ra1 += y2*fc*we*(asigma[ie]*asigma[ie] -cpi*api[ie]*api[ie]);
        ra3 += 2.0*y2*fc*we*(asigma[ie]*cpis*api[ie]);
      }
      (*irradiance)[k] = ra0;
      if (ra0>0) {
        (*p1)[k] = a1 = ra1/ra0;
        (*p3)[k] = a3 = ra3/ra0;
        (*p4)[k] = 1.0 - sqrt(a1*a1  + a3*a3);
      } /*end of if ra0>0 */
    } /* end of loop iy */
  } /* end of loop ix */
  *power = fac *  trapz2(*irradiance, nxp, nyp);
  free(eta);
  free(asigma);
  free(api);
  
}


/*following functions are converted from Roger Dejus's fortran codes for computing the x-ray
  intensity (magnet spectrum) */
void fkl(double xg, double yg, double k_magnet, double *s0, double *s1, double *s2, double *s3)
{
  /*
    Routine for fk in the limit K -> infinity. In reality, K > 100.0
c  is a suitable criterion.
c  Input: xg     gamma*theta
c         yg     gamma*psi
c         k_magnet     deflection parameter
c  Output:      Stokes parameters
c  Roger J. Dejus, XFD/APS, April, 1995.

  */
  double p, c, fkh, fkv;
  
  p = 1.0 + yg*yg;
  if (fabs(xg)< k_magnet)
    c = sqrt(1.0 - xg/k_magnet * xg/k_magnet);
  else
    c = 0;
  fkh = pow(c/p, 2.5);
  fkv = c*5.0*yg*yg/7.0/pow(p, 3.5);
  *s0 = fkh + fkv;
  *s1 = fkh - fkv;
  *s2 = 0.0;
  *s3 = 0.0;
  return;
}


/*c  Routine for calculation of fk. Typically used for K < 100.0. 
c  Input: x	gamma*theta
c         y	gamma*psi
c         k	deflection parameter
c  Output: 	Stokes parameters
c  Roger J. Dejus, XFD/APS, April, 1995. */
void fk(double xg, double yg, double k_magnet, double *s0, double *s1, double *s2, double *s3)
{
  double gk, c, A, B, eps, fkh, fkv;
  A=0;
  B=PI;
  eps = 1.0e-12;
  
  gk = k_magnet*(ipow(k_magnet, 6) +24.0/7.0*ipow(k_magnet, 4) +4.0*k_magnet*k_magnet +16.0/7.0)/pow(1.0 + k_magnet*k_magnet, 3.5);
  c  = 2.0*16.0/7.0*k_magnet/PI/gk;
  integ_par.xc = xg;
  integ_par.kc = k_magnet;
  integ_par.yc = yg;
  

  fkh = c*qromb(fkh_integrand, 20,A,B,eps);
  fkv = c*qromb(fkv_integrand, 20,A,B,eps);
  *s0  = fkh + fkv;
  *s1  = fkh - fkv;
  *s2  = 0.0;
  *s3  = 0.0;
  return;
}

double fkh_integrand(double alpha)
{
  /*
c  Integrand for calculation of fk for a regular planar 
c  undulator.  Typically used for K < 100.0.
c  Given by K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5.

real*8		alpha,p,d
real*8		x,y,k
common /fkc/ 	x,y,k

p = x -k*cos(alpha)
d = 1.0d0 +y**2 +p**2
fk_integrand = (1.0/d**3 -4.0d0*p**2/d**5)*sin(alpha)**2
return
END

REAL*8 FUNCTION fkh_integrand(alpha)
c  Integrand for calculation of fk in the horizontal plane for a
c  regular planar undulator.  Typically used for K < 100.0.
c  Given by K.J. Kim in "Characteristics of Synchrotron Radiation",
c  in Physics of Particle Accelerators, Vol 184, AIP Conference 
c  Proceedings, p. 565, New York (1989).
c  See also K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5, which gives the power for the 
c  horizontally and vertically polarized components added.
  */
  
  double p,d,h, result, xc, yc, kc;
  
  xc = integ_par.xc;
  yc = integ_par.yc;
  kc = integ_par.kc;
  
  p = xc -kc*cos(alpha);
  d = 1.0 + yc*yc + p*p;
  h = 1.0 +yc*yc -p*p;
  result = h*h/ipow(d, 5)*sin(alpha)*sin(alpha);
  return result;
}

double  fkv_integrand(double alpha)
{
  /*
    c  Integrand for calculation of fk in the vertical plane for a
    c  regular planar undulator.  Typically used for K < 100.0.
    c  Given by K.J. Kim in "Characteristics of Synchrotron Radiation",
c  in Physics of Particle Accelerators, Vol 184, AIP Conference 
c  Proceedings, p. 565, New York (1989).
c  See also K.J. Kim in "Angular Distribution of Undulator Power 
c  for an Arbitrary Deflection Parameter K", Nucl. Instr. Meth. 
c  A246, (1986) 67-70, Eq. 5, which gives the power for the 
c  horizontally and vertically polarized components added. */
  double p, d, v, xc, yc, kc;
  
  xc = integ_par.xc;
  yc = integ_par.yc;
  kc = integ_par.kc;
  
  p = xc -kc*cos(alpha);
  d = 1.0 +yc*yc +p*p;
  v = 2.0*yc*p;
  return (v*v/ipow(d, 5))*sin(alpha)*sin(alpha);
  
}
