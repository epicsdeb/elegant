/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/
/* program: sddsbs.c
   by Hairong Shang

$Log: not supported by cvs2svn $
Revision 1.3  2011/03/10 01:37:35  lemery
Added integrated flux computation using the obscure function gy which I found in mdb library, also used by sddsws. Verified by actual angle integration of sddsbs results for different angles.

Revision 1.2  2011/03/09 20:44:37  shang
added description

Revision 1.1  2011/03/09 02:21:25  lemery
First installation, per Shang.


*/

#include "scan.h"
#include "SDDS.h"
#include "mdb.h"
#include "oagphy.h"
#include "ws_constants.h"

#define CLO_ELECTRON_BEAM 0
#define CLO_PHOTON_ENERGY 1
#define CLO_ANGLE 2
#define CLO_MAGNET 3
#define CLO_PIPE 4
#define CLO_NO_WARNINGS 5
#define N_OPTIONS 6

static char *option[N_OPTIONS] = {
  "electronBeam", "photonEnergy", "angle", "magnetField", "pipe", "nowarnings" };

#define SPECTRA_MODES 2
static char *spectra_options[SPECTRA_MODES]={"frequency", "spatial"};

char *USAGE="sddsbs <outputFile> [-pipe[=out]] [-nowarnings] \n\
     [-electronBeam=current=<value>(mA),energy=<value>(GeV)[,xemittance=<value>][,yemittance=<value>]] \n\
     [-photonEnergy=minimum=<value>(eV),maximum=<value>(eV),points=<number>] \n\
     [-angle=start=<value>,end=<value>,steps=<value>] \n\
     [-magnetField=<B>(T)]  \n\
electronBeam     Specifies the electron beam (storage ring) parameters: \n\
                 current  electron beam current in mA. (default is 100mA). \n\
                 energy   electron energy in Gev. (default is 7.0Gev).\n\
                 for brightness calculation, the xemittance and yemittance must be provided.\n\
photonEnergy     specifies the maximum and minimum photon energy in eV, \n\
                 and the number of energy points to be computed.\n\
angle            provided the observation angle range in mrad unit, it is not need for brightness.\n\
magnetField      specifies the magnetic field of bending magnet in Tesla unit. \n\n\
sddsbs  computes bending magnet specral flux distribution for a specified \n\
        vertical angle plus the same flux integrated over all vertical angles.\n\n";

void SetupOutputFile(char *outputFile, SDDS_DATASET *SDDSout);
void compute_flux_spectra(double *photonEnergy, long nE, double angle, double cE, 
                          double energy, double current, double gamma, double *flux);
void compute_integratedFlux_spectra(double *photonE, long nE, double cE, double energy, double current, double *integratedFlux);

int main(int argc, char **argv)
{
  char  *outputFile=NULL;
  double energy=7.0, current=100.0, xemittance=0, yemittance=0, startAngle, endAngle, emax, emin, bField;
  double *photonE, *flux, *integratedFlux, deltaE, deltaA, cE=0, cLamda, ptot=0, angle, gamma;
  long nE, nAngles, i_arg, noWarnings=0, i, tmpFileUsed=0;
  SDDS_DATASET  SDDSout;
  SCANNED_ARG *s_arg;
  unsigned long dummyFlags=0, pipeFlags=0;
  char desc[2048];
  char *input="obset";
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s", USAGE);
    exit(1);
  }
  photonE = flux = NULL;
  nE  = 0;
  emax = emin = 0;
  startAngle = endAngle = 0;
  nAngles = 1;
  energy = 7; /* 7 GeV, APS */
  current = 100; /* 100mA, APS */
  bField = 1.27;
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case CLO_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        if (!(pipeFlags&USE_STDOUT))
          SDDS_Bomb("only -pipe=out syntax is valid!");
        break;
      case CLO_ELECTRON_BEAM:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -electronBeam syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "current", SDDS_DOUBLE, &current, 1, 0,
                          "energy", SDDS_DOUBLE, &energy, 1, 0,
                          "xemittance", SDDS_DOUBLE, &xemittance, 1, 0,
                          "yemittance", SDDS_DOUBLE, &yemittance, 1, 0,
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
      case CLO_ANGLE:
         if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -angle syntax.");
        s_arg[i_arg].n_items--;
        if (!scanItemList(&dummyFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "start", SDDS_DOUBLE, &startAngle, 1, 0,
                          "end", SDDS_DOUBLE, &endAngle, 1, 0,
                          "steps", SDDS_LONG, &nAngles, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -photonEnergy syntax");
        s_arg[i_arg].n_items++;
        break;
      case CLO_MAGNET:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -magnet syntax.");
        if (sscanf(s_arg[i_arg].list[1], "%lf", &bField)!=1)
          SDDS_Bomb("Invalide -magnet field value provided.");
        break;
      case CLO_NO_WARNINGS:
        noWarnings = 0;
        break;
      default:
        fprintf(stderr, "Invalid option - %s provided.\n",  s_arg[i_arg].list[0]);
        exit(1);
      }
    } else {
      if (outputFile==NULL)
        outputFile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  if (!outputFile && !pipeFlags)
    SDDS_Bomb("output file not provided.");
  if (outputFile && pipeFlags) 
    SDDS_Bomb("Too many files provided.");
  processFilenames("sddsbrightness", &input, &outputFile, pipeFlags, 0, &tmpFileUsed);
  if (nE<=1)
    SDDS_Bomb("The number of photon energy points has to be greater than 1.");
  if (emax<emin)
    SDDS_Bomb("The maximum photon energy is less than the minimum photon energy!");
  if (nAngles<1)
    SDDS_Bomb("The number of angle points has to be greater or equal to 1.");
  if (endAngle<startAngle)
    SDDS_Bomb("The end angle was smaller than the start angle.");
  if (startAngle == endAngle)
    nAngles = 1;
  
  if (nAngles>1)
    deltaA = (endAngle - startAngle)/(nAngles-1);
  else
    deltaA = 0; /*only one point */
  
  nE +=1; /* add one step to compare with online data */
  deltaE = (emax - emin)/(nE - 1);
  /*critical energy = 0.665E[GeV]^2 * B[T] */
  cE = CL2 * energy * energy * bField; /*eV */
  cLamda = C_EVANG/(cE*1.0e3); /*units is A, 1KeV=12.4A, 1eV=C_EVANG A */
  cLamda *= 10; /* change the units to nm */
  /*total radiation power = Cr *E^4/p = Cr * E^4/(3.333E/B) = Cr/3.333 * E^3 * B, Cr=8.85e-5 */
  ptot = 8.85e-5/3.333 * energy * energy * energy * bField * current * 1.0e6; /*units is Watts */
  current *= 1.0e-3; /*the input current is in mA, change it to A units */
  photonE = malloc(sizeof(*photonE)*nE);
  flux = malloc(sizeof(*flux)*nE);
  integratedFlux = malloc(sizeof(*integratedFlux)*nE);
  for (i=0; i<nE; i++)
    photonE[i] = emin + deltaE*i;
  gamma = energy/me_mev*1.0e3;
  SetupOutputFile(outputFile, &SDDSout);
  sprintf(desc, "bending magnet radiation spectra.");
  for (i=0; i<nAngles; i++) {
    angle = startAngle + deltaA * i;
    compute_flux_spectra(photonE, nE, angle, cE,  energy, current, gamma, flux);
    compute_integratedFlux_spectra(photonE, nE, cE,  energy, current, integratedFlux);

    if (!SDDS_StartPage(&SDDSout, nE) ||
        !SDDS_SetParameters(&SDDSout, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
                            "Description", desc, 
                            "eBeamEnergy", energy, 
                            "eBeamCurrent", current,
                            "MinEnergy", emin, 
                            "MaxEnergy", emax, 
                            "EnergyStep", nE,
                            "TotalPower", ptot, 
                            "MagnetField", bField,
                            "CriticalEnergy", cE/1.0e3, 
                            "CriticalWavelength", cLamda,
                            "ObservationAngle", angle, 
                            "Gamma", gamma, NULL) ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, photonE, nE, "Energy") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, flux, nE, "Flux") ||
        !SDDS_SetColumn(&SDDSout, SDDS_SET_BY_NAME, integratedFlux, nE, "IntegratedFlux") ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);  
  }
  if (!SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);  
  free(photonE);
  free(flux);
  
  free_scanargs(&s_arg, argc);
  return 0;
}

void SetupOutputFile(char *outputFile, SDDS_DATASET *SDDSout)
{
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, outputFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(SDDSout, "Description", NULL, NULL, NULL, NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "eBeamEnergy", NULL, "Gev", "Electron Beam Energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "eBeamCurrent", NULL, "A", "Electron Beam Current", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MinEnergy", NULL, "eV", "minimum photon energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MaxEnergy", NULL, "eV", "maximum photon energy", NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "EnergyStep", NULL, NULL, "number of points of photon energy", NULL, SDDS_LONG, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "MagnetField", NULL, "T", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "CriticalEnergy", NULL, "KeV", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "CriticalWavelength", NULL, "nm", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "TotalPower",  NULL, "Watts", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "ObservationAngle",  NULL, "mrad", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(SDDSout, "Gamma",  NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineColumn(SDDSout, "Energy", NULL, "eV", "photon energy",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Flux", NULL, "ph/s/0.1%bw/mrad$a2$n", NULL,NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "IntegratedFlux", NULL, "ph/s/0.1%bw/mrad", NULL,NULL, SDDS_DOUBLE, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteLayout(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}

/*bending magnet flux distribution ph/s/mrad^2/0.1%BW curent units is A in the formula 
  = Cn * E^2 * E^2/Ec^2 * I * K2/3^2(k) F(k, theta) while Cn = 1.3255e13 ph/s/mrad^2/GeV^2/A */ 
/* F(k, theta) = (1+gamma^2*theta^2) * [ 1 + gamma^2*theta^2/(1+gamma^2*theta^2) * K1/3(k)^2/K2/3(k)^2] */
void compute_flux_spectra(double *photonEnergy, long nE, double angle, double cE, 
                          double energy, double current, double gamma, double *flux)
{
  double eta, Ftheta, k1, KI13, KI23, Cr, ratio;
  long i;

  k1 = 1 + gamma*gamma*angle*angle;
  Cr = 1.3273e13;
  for (i=0; i<nE; i++) {
    ratio = photonEnergy[i]/cE;
    eta = 0.5 * ratio * k1 * sqrt(k1);
    KI13 = k13(eta);
    KI23 = k23(eta);
    Ftheta = k1 * k1 * (1 + gamma * gamma * angle * angle * KI13 * KI13 / (k1 * KI23 * KI23)); 
    flux[i] = Cr * energy * energy * current * ratio * ratio * KI23 * KI23 * Ftheta; 
  }
}

/* formula taken from Xray data booklet p2-3. */
void compute_integratedFlux_spectra(double *photonEnergy, long nE, double cE, 
	double energy, double current, double *flux)
{
  double  C, ratio;
  long i;

  C = 2.457e13;
  for (i=0; i<nE; i++) {
    ratio = photonEnergy[i]/cE;
    flux[i] = C * energy * current * gy(1, ratio);
  }
}
