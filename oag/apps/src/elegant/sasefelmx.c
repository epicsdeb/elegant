/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <math.h>
#include "mdb.h"
#include "SDDS.h"
#include "track.h"

void ComputeSASEFELParameters(
                              /* outputs */
                              double *lightWavelength,
                              double *saturationLength,
                              double *gainLength, 
                              double *noisePower,
                              double *saturationPower,
                              double *PierceParameter,
                              double *etaDiffraction,
                              double *etaEmittance,
                              double *etaEnergySpread,
                              /* inputs */
                              double charge,
                              double rmsBunchLength,
                              double undulatorPeriod,
                              double undulatorK, 
                              double beta, 
                              double emittance,
                              double sigmaDelta,
                              double pCentral, 
                              short planar)
{
  double sigma, aw, Aw, gamma, L1D, Ic, x, F;

  if (!planar)
    SDDS_Bomb("nonplanar undulators not yet supported");
  
  sigma = sqrt(beta*emittance);

  aw = undulatorK/sqrt(2.0);
  x = sqr(aw)/(2*(1+sqr(aw)));
  Aw = aw*(jn(0, x)-jn(1, x));

  gamma = sqrt(sqr(pCentral)+1);
  *lightWavelength = undulatorPeriod*(1+sqr(aw))/(2*sqr(gamma));
    
  Ic = charge/(sqrt(PIx2)*rmsBunchLength);
  
  *PierceParameter 
    = pow( (Ic/17.045e3) *
          sqr(undulatorPeriod*Aw/(PIx2*sigma)) *
          pow(0.5/gamma, 3), 1./3.);
  
  L1D = undulatorPeriod/(4*PI*sqrt(3)*(*PierceParameter));

  F = FELScalingFunction(etaDiffraction, etaEmittance, etaEnergySpread,
                         L1D, beta, emittance, *lightWavelength,
                         undulatorPeriod, sigmaDelta);
  *gainLength = L1D/F;
  *noisePower = sqr(*PierceParameter)*me_mks*pow(c_mks, 3)*gamma/(*lightWavelength);
  *saturationPower = 1.6*(*PierceParameter)*sqr(F)*(0.511e6*gamma*Ic);
  *saturationLength = (*gainLength)*log(*saturationPower/(*noisePower/9.0));
}

double FELScalingFunction(double *etaDiffraction, double *etaEmittance,
                          double *etaEnergySpread,
                          double L1D, double beta, double emittance,
                          double lightWavelength, double undulatorPeriod,
                          double sigmaDelta)
{
  double etaDiff, etaEmit, etaEnSp, LRayleigh, eta;
  double a1 = 0.45, a2=0.57, a3=0.55, a4=1.6;
  double a5=3, a6=2, a7=0.35, a8=2.9;
  double a9=2.4, a10=51, a11=0.95, a12=3;
  double a13=5.4, a14=0.7, a15=1.9, a16=1140;
  double a17=2.2, a18=2.9, a19=3.2;
  
  LRayleigh = 4*PI*emittance*beta/lightWavelength;
  etaDiff = *etaDiffraction = L1D/LRayleigh;
  etaEmit = *etaEmittance = (L1D/beta)*(4*PI*emittance/lightWavelength);
  etaEnSp = *etaEnergySpread = 4*PI*(L1D/undulatorPeriod)*sigmaDelta;

  eta = a1*pow(etaDiff, a2) +
    a3*pow(etaEmit, a4) +
      a5*pow(etaEnSp, a6) +
        a7*pow(etaEmit, a8)*pow(etaEnSp, a9) +
          a10*pow(etaDiff, a11)*pow(etaEnSp, a12) +
            a13*pow(etaDiff, a14)*pow(etaEmit, a15) +
              a16*pow(etaDiff, a17)*pow(etaEmit, a18)*pow(etaEnSp, a19);
  return 1./(1+eta);
}


