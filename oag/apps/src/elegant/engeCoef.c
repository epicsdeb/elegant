/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"

static double bnEnge[3];
static double rhoEnge, lengthEnge, angleEnge, KgEnge;
double engeProfile(double z)
{
  return 1/(1 + exp(bnEnge[0] + bnEnge[1]*z + bnEnge[2]*sqr(z)));
}

double engeProfileTheta(double theta)
{
  return engeProfile(rhoEnge*sin(theta));
}

double engeEdgeFactor(double z)
{
  double ep;
  ep = engeProfile(z);
  return ep*(1-ep);
}


double engeOptimizationFunction(double *b, long *invalid)
{
  static double F[3];
  double temp1, temp2;
  
  bnEnge[0] = b[0];
  bnEnge[1] = b[1];
  bnEnge[2] = b[2];
  *invalid = 0;
  
  /* field should go to 1 at center of arc */
  F[0] = 1 - engeProfile(-rhoEnge*tan(angleEnge/2));

  /* effective length constraint */
  if (!gaussianQuadrature(engeProfileTheta, -angleEnge/2, 0.0, 100, 1e-8, &temp1)) {
    fprintf(stderr, "GQ #1 failed\n");
    *invalid = 1;
    return 0;
  }
  temp1 *= rhoEnge;
  if (!gaussianQuadrature(engeProfile, 0.0, 10*lengthEnge, 1000, 1e-8, &temp2)) {
    fprintf(stderr, "GQ #2 failed\n");
    *invalid = 1;
    return 0;
  }
  F[1] = (lengthEnge/2 - (temp1 + temp2))/(lengthEnge/2);
  
  /* edge integral constraint */
  if (!gaussianQuadrature(engeEdgeFactor, -lengthEnge/2, 10*lengthEnge, 1000, 1e-8, &F[2])) {
    fprintf(stderr, "GQ #3 failed\n");
    *invalid = 1;
    return 0;
  }
  F[2] = (KgEnge - F[2])/KgEnge;
  
  return sqr(F[0]) + sqr(F[1]) + sqr(F[2]);
}

long computeEngeCoefficients(double *engeCoef, double rho, double length, double gap, double fint)
{
  double b[3], db[3], bMin[3], bMax[3];
  double result;
  static double lastData[4] = {-1, -1, -1, -1};
  static double lastResult[3];

  if (rho==lastData[0] && length==lastData[1] && gap==lastData[2] && fint==lastData[3]) {
    memcpy(engeCoef, lastResult, 3*sizeof(*engeCoef));
    return 1;
  }
  
  rhoEnge = rho;
  lengthEnge = length;
  angleEnge = length/rho;
  KgEnge = fint*gap;
  b[0] = -0.003183;
  b[1] = 1.911302/gap;
  b[2] = 0.0;
  db[0] = db[1] = db[2] = 1e-4;
  bMin[0] = bMin[1] = bMin[2] = -1e10;
  bMax[0] = bMax[1] = bMax[2] =  1e10;

  if (simplexMin(&result, b, db, bMin, bMax, NULL, 3, 1e-14, 1e-16, engeOptimizationFunction,
                 NULL, 500, 3, 12, 10.0, 10.0, 0)<0) {
    fprintf(stderr, "Problem finding enge coefficients, using defaults\n");
    engeCoef[0] = -0.003183;
    engeCoef[1] = 1.911302;
    engeCoef[2] = 0.0;
  }
  engeCoef[0] = b[0];
  engeCoef[1] = b[1]*gap;
  engeCoef[2] = b[2]*sqr(gap);

  lastData[0] = rho;
  lastData[1] = length;
  lastData[2] = gap;
  lastData[3] = fint;
  memcpy(lastResult, engeCoef, 3*sizeof(*lastResult));
  return 1;
}

