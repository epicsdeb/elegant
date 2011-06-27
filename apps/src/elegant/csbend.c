/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: csbend.c
 * contents:  track_through_canonical_sbend()
 *
 *
 * Michael Borland, 1991, 1992.
 */
#include "mdb.h"
#include "track.h"

#define EXSQRT(value, order) (order==0?sqrt(value):(1+0.5*((value)-1)))

static long negativeWarningsLeft = 100;

void addRadiationKick(double *Qx, double *Qy, double *dPoP, double *sigmaDelta2, long sqrtOrder,
		      double x, double h0, double Fx, double Fy,
		      double ds, double radCoef, double dsISR, double isrCoef,
                      long distributionBased, long includeOpeningAngle,
                      double meanPhotonsPerMeter,
                      double normalizedCriticalEnergy, double Po);
double pickNormalizedPhotonEnergy(double RN);

void integrate_csbend_ord2(double *Qf, double *Qi, double *sigmaDelta2, double s, long n, long sqrtOrder, double rho0, double p0);
void integrate_csbend_ord4(double *Qf, double *Qi, double *sigmaDelta2, double s, long n, long sqrtOrder, double rho0, double p0);
void exactDrift(double **part, long np, double length);
void convertFromCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode);
void convertToCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode);
void applyFilterTable(double *function, long bins, double dt, long fValues,
                      double *fFreq, double *fReal, double *fImag);

long correctDistribution(double *array, long npoints, double desiredSum);

static double **Fx_xy = NULL, **Fy_xy = NULL;
static long expansionOrder1 = 11;  /* order of expansion+1 */

static double rho0, rho_actual, rad_coef=0, isrConstant=0;
static double meanPhotonsPerRadian0, meanPhotonsPerMeter0, normalizedCriticalEnergy0;
static long distributionBasedRadiation, includeOpeningAngle;
static long photonCount = 0;
static double energyCount = 0, radiansTotal = 0;

static long particle_lost;
static double s_lost;

long inversePoissonCDF(double mu, double C);

#if !defined(PARALLEL)
/* to avoid problems with HP parallel compiler */
extern unsigned long multipoleKicksDone ;
#endif

void computeCSBENDFields(double *Fx, double *Fy, double x, double y)
{
  double xp[11], yp[11];
  double sumFx=0, sumFy=0;
  long i, j;
  static short first = 1;
  
  xp[0] = yp[0] = 1;
  for (i=1; i<expansionOrder1; i++) {
    xp[i] = xp[i-1]*x;
    yp[i] = yp[i-1]*y;
  }
  
  for (i=0; i<expansionOrder1; i++)
    for (j=1; j<expansionOrder1-i; j+=2) {
      if (Fx_xy[i][j]) {
        sumFx += Fx_xy[i][j]*xp[i]*yp[j];
/*
        if (first)
          fprintf(stdout, "Including in Fx: %e * x^%ld y^%ld\n", Fx_xy[i][j], i, j);
*/
      }
    }
  *Fx = sumFx;

  for (i=0; i<expansionOrder1; i++)
    for (j=0; j<expansionOrder1-i; j+=2) {
      if (Fy_xy[i][j]) {
        sumFy += Fy_xy[i][j]*xp[i]*yp[j];
/*
        if (first)
          fprintf(stdout, "Including in Fy: %e * x^%ld y^%ld\n", Fy_xy[i][j], i, j);
*/
      }
    }

  *Fy = sumFy;

  first = 0;
}

void computeCSBENDFieldCoefficients(double *b, double h, long nonlinear, long expansionOrder)
{
  long i;
  double h2, h3, h4, h5;

  if (expansionOrder==0) {
    /* set the order to be <highestMultipole>+2 */
    for (i=7; i>=0; i--)
      if (b[i])
        break;
    if ((expansionOrder = i+3)<4)
      /* make minimum value 4 for backward compatibility */
      expansionOrder = 4;  
  }

  expansionOrder1 = expansionOrder + 1;
  if (expansionOrder1>11)
    bombElegant("expansion order >10 for CSBEND or CSRCSBEND", NULL);
  
  if (!Fx_xy)
    Fx_xy = (double**)czarray_2d(sizeof(double), 11, 11);
  if (!Fy_xy)
    Fy_xy = (double**)czarray_2d(sizeof(double), 11, 11);
    
  for (i=0; i<expansionOrder1; i++) {
    memset(Fx_xy[i], 0, expansionOrder1*sizeof(double));
    memset(Fy_xy[i], 0, expansionOrder1*sizeof(double));
  }
  
  h2 = h*h;
  h3 = h2*h;
  h4 = h2*h2;
  h5 = h4*h;

  Fx_xy[0][1] = b[0];
  Fy_xy[0][0] = 1;
  Fy_xy[1][0] = b[0];
  
  if (nonlinear) {
    Fx_xy[1][1] = b[1];
    Fx_xy[2][1] = b[2]/2.;
    Fx_xy[3][1] = b[3]/6.;
    Fx_xy[4][1] = b[4]/24.;
    Fx_xy[5][1] = b[5]/120.;
    Fx_xy[6][1] = b[6]/720.;
    Fx_xy[7][1] = b[7]/5040.;
    Fx_xy[0][3] = (-b[2] - b[1]*h + b[0]*h2)/6.;
    Fx_xy[1][3] = (-b[3] - b[2]*h + 2*b[1]*h2 - 2*b[0]*h3)/6.;
    Fx_xy[2][3] = (-b[4] + h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))/12.;
    Fx_xy[3][3] = (-b[5] - h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2))))/36.;
    Fx_xy[4][3] = (-b[6] + h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))/ 144.;
    Fx_xy[5][3] = (-b[7] - h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2))))))/720.;
    Fx_xy[6][3] = (h*(-b[7] + 7*h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))))/ 4320.;
    Fx_xy[7][3] = (h2*(b[7] - 7*h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))))/3780.;
    Fx_xy[0][5] = (b[4] + h*(2*b[3] - 3*h*(b[2] - b[1]*h + b[0]*h2)))/120.;
    Fx_xy[1][5] = (b[5] + h*(2*b[4] + h*(-5*b[3] + 9*b[2]*h - 12*b[1]*h2 + 12*b[0]*h3)))/ 120.;
    Fx_xy[2][5] = (b[6] + h*(2*b[5] + h*(-7*b[4] + 19*b[3]*h - 39*b[2]*h2 + 60*b[1]*h3 - 60*b[0]*h4)))/240.;
    Fx_xy[3][5] = (b[7] + h*(2*b[6] + 3*h*(-3*b[5] + 11*b[4]*h - 32*b[3]*h2 + 72*b[2]*h3 - 120*b[1]*h4 + 120*b[0]*h5)))/720.;
    Fx_xy[4][5] = (h*(2*b[7] + h*(-11*b[6] - 3*h*(-17*b[5] + 5*h*(13*b[4] + 8*h*(-5*b[3] + 3*h*(4*b[2] - 7*b[1]*h + 7*b[0]*h2)))))))/2880.;
    Fx_xy[5][5] = (h2*(-13*b[7] + h*(73*b[6] + 12*h*(-29*b[5] + 5*h*(23*b[4] + 2*h*(-37*b[3] + 3*h*(31*b[2] + 56*h*(-b[1] + b[0]*h))))))))/14400.;
    Fx_xy[0][7] = (-b[6] - 3*h*(b[5] + h*(-2*b[4] + h*(4*b[3] - 9*b[2]*h + 15*b[1]*h2 - 15*b[0]*h3))))/5040.;
    Fx_xy[1][7] = (-b[7] - 3*h*(b[6] + h*(-3*b[5] + 8*b[4]*h - 21*b[3]*h2 + 51*b[2]*h3 - 90*b[1]*h4 + 90*b[0]*h5)))/5040.;
    Fx_xy[2][7] = (h*(-b[7] + h*(4*b[6] + h*(-14*b[5] + 45*b[4]*h - 135*b[3]*h2 + 345*b[2]*h3 - 630*b[1]*h4 + 630*b[0]*h5))))/ 3360.;
    Fx_xy[3][7] = (h2*(5*b[7] + h*(-22*b[6] + 3*h*(29*b[5] - 5*h*(21*b[4] + 4*h*(-17*b[3] + 45*b[2]*h - 84*b[1]*h2 + 84*b[0]*h3))))))/10080.;
    Fx_xy[0][9] = (h*(4*b[7] - 5*h*(2*b[6] + 3*h*(-2*b[5] + 7*b[4]*h - 22*b[3]*h2 + 57*b[2]*h3 - 105*b[1]*h4 + 105*b[0]*h5))))/ 362880.;
    Fx_xy[1][9] = (h2*(-14*b[7] + 5*h*(10*b[6] + 3*h*(-13*b[5] + 50*b[4]*h - 167*b[3]*h2 + 447*b[2]*h3 - 840*b[1]*h4 + 840*b[0]*h5))))/362880.;

    Fy_xy[2][0] = b[1]/2.;
    Fy_xy[3][0] = b[2]/6.;
    Fy_xy[4][0] = b[3]/24.;
    Fy_xy[5][0] = b[4]/120.;
    Fy_xy[6][0] = b[5]/720.;
    Fy_xy[7][0] = b[6]/5040.;
    Fy_xy[8][0] = b[7]/40320.;
    Fy_xy[0][2] = (-b[1] - b[0]*h)/2.;
    Fy_xy[1][2] = (-b[2] - b[1]*h + b[0]*h2)/2.;
    Fy_xy[2][2] = (-b[3] - b[2]*h + 2*b[1]*h2 - 2*b[0]*h3)/4.;
    Fy_xy[3][2] = (-b[4] + h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))/12.;
    Fy_xy[4][2] = (-b[5] - h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2))))/48.;
    Fy_xy[5][2] = (-b[6] + h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))/ 240.;
    Fy_xy[6][2] = (-b[7] - h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2))))))/1440.;
    Fy_xy[7][2] = (h*(-b[7] + 7*h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))))/ 10080.;
    Fy_xy[8][2] = (h2*(b[7] - 7*h*(b[6] + 6*h*(-b[5] + 5*h*(b[4] + 4*h*(-b[3] + 3*h*(b[2] - 2*b[1]*h + 2*b[0]*h2)))))))/10080.;
    Fy_xy[0][4] = (b[3] + h*(2*b[2] + h*(-b[1] + b[0]*h)))/24.;
    Fy_xy[1][4] = (b[4] + h*(2*b[3] - 3*h*(b[2] - b[1]*h + b[0]*h2)))/24.;
    Fy_xy[2][4] = (b[5] + h*(2*b[4] + h*(-5*b[3] + 9*b[2]*h - 12*b[1]*h2 + 12*b[0]*h3)))/48.;
    Fy_xy[3][4] = (b[6] + h*(2*b[5] + h*(-7*b[4] + 19*b[3]*h - 39*b[2]*h2 + 60*b[1]*h3 - 60*b[0]*h4)))/144.;
    Fy_xy[4][4] = (b[7] + h*(2*b[6] + 3*h*(-3*b[5] + 11*b[4]*h - 32*b[3]*h2 + 72*b[2]*h3 - 120*b[1]*h4 + 120*b[0]*h5)))/576.;
    Fy_xy[5][4] = (h*(2*b[7] + h*(-11*b[6] - 3*h*(-17*b[5] + 5*h*(13*b[4] + 8*h*(-5*b[3] + 3*h*(4*b[2] - 7*b[1]*h + 7*b[0]*h2)))))))/2880.;
    Fy_xy[6][4] = (h2*(-13*b[7] + h*(73*b[6] + 12*h*(-29*b[5] + 5*h*(23*b[4] + 2*h*(-37*b[3] + 3*h*(31*b[2] + 56*h*(-b[1] + b[0]*h))))))))/17280.;
    Fy_xy[0][6] = (-b[5] - 3*h*(b[4] + h*(-b[3] + 2*b[2]*h - 3*b[1]*h2 + 3*b[0]*h3)))/720.;
    Fy_xy[1][6] = (-b[6] - 3*h*(b[5] + h*(-2*b[4] + h*(4*b[3] - 9*b[2]*h + 15*b[1]*h2 - 15*b[0]*h3))))/720.;
    Fy_xy[2][6] = (-b[7] - 3*h*(b[6] + h*(-3*b[5] + 8*b[4]*h - 21*b[3]*h2 + 51*b[2]*h3 - 90*b[1]*h4 + 90*b[0]*h5)))/1440.;
    Fy_xy[3][6] = (h*(-b[7] + h*(4*b[6] + h*(-14*b[5] + 45*b[4]*h - 135*b[3]*h2 + 345*b[2]*h3 - 630*b[1]*h4 + 630*b[0]*h5))))/ 1440.;
    Fy_xy[4][6] = (h2*(5*b[7] + h*(-22*b[6] + 3*h*(29*b[5] - 5*h*(21*b[4] + 4*h*(-17*b[3] + 45*b[2]*h - 84*b[1]*h2 + 84*b[0]*h3))))))/5760.;
    Fy_xy[0][8] = (b[7] + h*(4*b[6] + 3*h*(-2*b[5] + 6*b[4]*h - 17*b[3]*h2 + 42*b[2]*h3 - 75*b[1]*h4 + 75*b[0]*h5)))/40320.;
    Fy_xy[1][8] = (h*(4*b[7] - 5*h*(2*b[6] + 3*h*(-2*b[5] + 7*b[4]*h - 22*b[3]*h2 + 57*b[2]*h3 - 105*b[1]*h4 + 105*b[0]*h5))))/ 40320.;
    Fy_xy[2][8] = (h2*(-14*b[7] + 5*h*(10*b[6] + 3*h*(-13*b[5] + 50*b[4]*h - 167*b[3]*h2 + 447*b[2]*h3 - 840*b[1]*h4 + 840*b[0]*h5))))/80640.;
    Fy_xy[0][10] = (h2*(2*b[7] + h*(-8*b[6] - 3*h*(-11*b[5] + h*(43*b[4] + 5*h*(-29*b[3] + 3*h*(26*b[2] + 49*h*(-b[1] + b[0]*h))))))))/725760.;
  }

  /*
  for (i=0; i<11; i++) 
    for (j=0; j<11; j++)
      fprintf(stdout, "Fx[%ld][%ld] = %e   Fy[%ld][%ld] = %e\n", i, j, Fx_xy[i][j], i, j, Fy_xy[i][j]);
      */
}


long track_through_csbend(double **part, long n_part, CSBEND *csbend, double p_error, double Po, double **accepted,
                          double z_start, double *sigmaDelta2)
{
  double h;
  long i_part, i_top;
  double rho, s, Fx, Fy;
  double x, xp, y, yp, dp, dp0;
  double n, fse, dp_prime;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double *coord;
  double angle, e1, e2, Kg;
  double psi1, psi2, he1, he2;
  double Qi[6], Qf[6];
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double delta_xp;
  double e1_kick_limit, e2_kick_limit;
  static long largeRhoWarning = 0;

  if (!csbend)
    bombElegant("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0) {
    exactDrift(part, n_part, csbend->length);
    return n_part;
  }
  
  if (!(csbend->edgeFlags&BEND_EDGE_DETERMINED)) 
    bombElegant("CSBEND element doesn't have edge flags set.", NULL);
  
  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bombElegant("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  rho0 =  csbend->length/csbend->angle;
  if (csbend->use_bn) {
    csbend->b[0] = csbend->b1;
    csbend->b[1] = csbend->b2;
    csbend->b[2] = csbend->b3;
    csbend->b[3] = csbend->b4;
    csbend->b[4] = csbend->b5;
    csbend->b[5] = csbend->b6;
    csbend->b[6] = csbend->b7;
    csbend->b[7] = csbend->b8;
  } else {
    csbend->b[0] = csbend->k1*rho0;
    csbend->b[1] = csbend->k2*rho0;
    csbend->b[2] = csbend->k3*rho0;
    csbend->b[3] = csbend->k4*rho0;
    csbend->b[4] = csbend->k5*rho0;
    csbend->b[5] = csbend->k6*rho0;
    csbend->b[6] = csbend->k7*rho0;
    csbend->b[7] = csbend->k8*rho0;
  }
  
  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    long i;
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;      /* work in rotated system */
    rho0  = csbend->length/angle;
    for (i=0; i<8; i+=2)
      csbend->b[i] *= -1;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
  }


  if (rho0>1e6) {
    if (!largeRhoWarning) {
      printf("Warning: One or more CSBENDs have radius > 1e6.  Treated as drift.\n");
      largeRhoWarning = 1;
    }
    exactDrift(part, n_part, csbend->length);
    return n_part;
  }
  
  fse = csbend->fse;
  h = 1/rho0;
  n = -csbend->b[0]/h;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  e1_kick_limit = csbend->edge1_kick_limit;
  e2_kick_limit = csbend->edge2_kick_limit;
  if (csbend->kick_limit_scaling) {
    e1_kick_limit *= rho0/rho_actual;
    e2_kick_limit *= rho0/rho_actual;
  }
  if (e1_kick_limit>0 || e2_kick_limit>0)
    fprintf(stdout, "rho0=%e  rho_a=%e fse=%e e1_kick_limit=%e e2_kick_limit=%e\n",
            rho0, rho_actual, csbend->fse, e1_kick_limit, e2_kick_limit);
    fflush(stdout);
  
  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*particleMass*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  isrConstant = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                            137.0359895/pow3(fabs(rho_actual)));
  if (!csbend->isr || (csbend->isr1Particle==0 && n_part==1))
    /* Minus sign here indicates that we accumulate ISR into sigmaDelta^2 but don't apply it to particles. */
    isrConstant *= -1; 

  if ((distributionBasedRadiation = csbend->distributionBasedRadiation)) {
    /* Sands 5.15 */
    meanPhotonsPerRadian0 = 5.0/(2.0*sqrt(3))*Po/137.0359895;  
    meanPhotonsPerMeter0 = (5*c_mks*Po*particleMass*particleRadius)/(2*sqrt(3)*hbar_mks*rho_actual);
    /* Critical energy normalized to beam energy, Sands 5.9 */
    normalizedCriticalEnergy0 = 3.0/2*hbar_mks*c_mks*pow3(Po)/fabs(rho_actual)/(Po*particleMass*sqr(c_mks));
    /* fprintf(stderr, "Mean photons per radian expected: %le   ECritical/E: %le\n", 
            meanPhotonsPerRadian0, normalizedCriticalEnergy0);
    */
    includeOpeningAngle = csbend->includeOpeningAngle;
  }
  
  computeCSBENDFieldCoefficients(csbend->b, h, csbend->nonlinear, csbend->expansionOrder);

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }

  computeEtiltCentroidOffset(dcoord_etilt, rho0, angle, etilt, tilt);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;
  
  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

  i_top = n_part-1;
#if !defined(PARALLEL)
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  if (isSlave || !notSinglePart) {
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!part) {
      fprintf(stdout, "error: null particle array found (working on particle %ld) (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (!(coord = part[i_part])) {
      fprintf(stdout, "error: null coordinate pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }

    coord[4] += dzi*EXSQRT(1 + sqr(coord[1]) + sqr(coord[3]), csbend->sqrtOrder);
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];

    x  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    y  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    s  = coord[4];
    dp = dp0 = coord[5];

    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS) {
      rho = (1+dp)*rho_actual;
      if (csbend->edge_order<2) {
        /* apply edge focusing */
        delta_xp = tan(e1)/rho*x;
        if (e1_kick_limit>0 && fabs(delta_xp)>e1_kick_limit)
          delta_xp = SIGN(delta_xp)*e1_kick_limit;
        xp += delta_xp;
        yp -= tan(e1-psi1/(1+dp))/rho*y;
      } else
        apply_edge_effects(&x, &xp, &y, &yp, rho, n, e1, he1, psi1*(1+dp), -1);
    }

    /* transform to curvilinear coordinates */
    xp *= (1+x/rho0);
    yp *= (1+x/rho0);

    /* load input coordinates into arrays */
    Qi[0] = x;  Qi[1] = xp;  Qi[2] = y;  Qi[3] = yp;  Qi[4] = 0;  Qi[5] = dp;

    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS && e1!=0 && rad_coef) {
      /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
      computeCSBENDFields(&Fx, &Fy, x, y);

      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*EXSQRT(sqr(1+x/rho0)+sqr(xp)+sqr(yp), csbend->sqrtOrder);
      Qi[5] -= dp_prime*x*tan(e1);
    }

    particle_lost = 0;
    if (!particle_lost) {
      if (csbend->integration_order==4)
        integrate_csbend_ord4(Qf, Qi, sigmaDelta2, csbend->length, csbend->n_kicks, csbend->sqrtOrder, rho0, Po);
      else
        integrate_csbend_ord2(Qf, Qi, sigmaDelta2, csbend->length, csbend->n_kicks, csbend->sqrtOrder, rho0, Po);
    }

    if (particle_lost) {
      if (!part[i_top]) {
        fprintf(stdout, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                i_part, i_top);
        fflush(stdout);
        abort();
      }
      swapParticles(part[i_part], part[i_top]);
      if (accepted) {
        if (!accepted[i_top]) {
          fprintf(stdout, 
                  "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
                  i_part, i_top);
          fflush(stdout);
          abort();
        }
        swapParticles(accepted[i_part], accepted[i_top]);
      }
      part[i_top][4] = z_start + s_lost;
      part[i_top][5] = Po*(1+part[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }

    if (csbend->edgeFlags&BEND_EDGE2_EFFECTS && e2!=0 && rad_coef) {
      /* post-adjust dp/p to correct error made by integrating over entire sector */
      x = Qf[0];
      xp = Qf[1];
      y = Qf[2];
      yp = Qf[3];
      dp = Qf[5];

      computeCSBENDFields(&Fx, &Fy, x, y);

      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+dp)*EXSQRT(sqr(1+x/rho0)+sqr(xp)+sqr(yp), csbend->sqrtOrder);
      Qf[5] -= dp_prime*x*tan(e2);
    }

    /* get final coordinates */
    if (rad_coef || isrConstant) {
      double p0, p1;
      double beta0, beta1;
      /* fix previous distance information to reflect new velocity--since distance
       * is really time-of-flight at the current velocity 
       */
      p0 = Po*(1+dp0);
      beta0 = p0/sqrt(sqr(p0)+1);
      p1 = Po*(1+Qf[5]);
      beta1 = p1/sqrt(sqr(p1)+1);
      s = beta1*s/beta0 + Qf[4];
    }
    else
      s += Qf[4];
    x = Qf[0];  xp = Qf[1];  y = Qf[2];  yp = Qf[3];  dp = Qf[5];

    /* transform to cartesian coordinates */
    xp /= (1+x/rho0);
    yp /= (1+x/rho0);

    if (csbend->edgeFlags&BEND_EDGE2_EFFECTS) {
      /* apply edge focusing */
      rho = (1+dp)*rho_actual;
      if (csbend->edge_order<2) {
        delta_xp = tan(e2)/rho*x;
        if (e2_kick_limit>0 && fabs(delta_xp)>e2_kick_limit)
          delta_xp = SIGN(delta_xp)*e2_kick_limit;
        xp += delta_xp;
        yp -= tan(e2-psi2/(1+dp))/rho*y;
      } else
        apply_edge_effects(&x, &xp, &y, &yp, rho, n, e2, he2, psi2*(1+dp), 1);
    }
    
    coord[0] =  x*cos_ttilt -  y*sin_ttilt + dcoord_etilt[0];
    coord[2] =  x*sin_ttilt +  y*cos_ttilt + dcoord_etilt[2];
    coord[1] = xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[1];
    coord[3] = xp*sin_ttilt + yp*cos_ttilt + dcoord_etilt[3];
    coord[4] = s;
    coord[5] = dp;

    coord[0] += dxf + dzf*coord[1];
    coord[2] += dyf + dzf*coord[3];
    coord[4] += dzf*EXSQRT(1+ sqr(coord[1]) + sqr(coord[3]), csbend->sqrtOrder);
  }
  }
  if (distributionBasedRadiation) {
    radiansTotal += fabs(csbend->angle);
    /*
      fprintf(stderr, "%e radians, photons/particle=%e, photons/radian = %e, mean y = %e\n",
      radiansTotal, photonCount/(1.0*i_top), photonCount/radiansTotal/(1.0*i_top), energyCount/photonCount);
    */
    distributionBasedRadiation = 0;
  }

  if (sigmaDelta2)
    /* Return average value for all particles */
    *sigmaDelta2 /= i_top+1;

  return(i_top+1);
}

void integrate_csbend_ord2(double *Qf, double *Qi, double *sigmaDelta2, double s, long n, long sqrtOrder, double rho0, double p0)
{
  long i;
  double factor, f, phi, ds, dsh, dp, dist;
  double Fx, Fy, x, y;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  
#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

  if (!Qf)
    bombElegant("NULL final coordinates pointer ()", NULL);
  if (!Qi)
    bombElegant("NULL initial coordinates pointer (integrate_csbend_ord2)", NULL);
  if (n<1)
    bombElegant("invalid number of steps (integrate_csbend_ord2)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
  dp = DPoP0;
  f = (1+dp)/EXSQRT(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0), sqrtOrder);
  QX = XP0*f;
  QY = YP0*f;

  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;

  ds = s/n;
  dsh = ds/2;
  dist = 0;

  
  for (i=0; i<n; i++) {
    if (i==0) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }

    /* calculate the scaled fields */
    x = X;
    y = Y;

    computeCSBENDFields(&Fx, &Fy, x, y);

    /* do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant)
      addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder, 
		       X, 1./rho0, Fx, Fy, 
		       ds, rad_coef, ds, isrConstant, 
                       distributionBasedRadiation, includeOpeningAngle,
                       meanPhotonsPerMeter0, normalizedCriticalEnergy0, p0);
    
    if (i==n-1) {
      /* do half-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(dsh/rho0+phi);
      if ((cosi = cos(dsh/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }
    else {
      /* do full-length drift */
      if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      f = EXSQRT(f, sqrtOrder);
      if (fabs(QX/f)>1) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      phi = asin(sin_phi=QX/f);
      sine = sin(ds/rho0+phi);
      if ((cosi = cos(ds/rho0+phi))==0) {
        particle_lost = 1;
        s_lost = dist;
        return;
      }
      tang = sine/cosi;
      cos_phi = cos(phi);
      QX = f*sine;
      Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
      dist += factor*(1+DPoP);
      f = cos_phi/cosi;
      X  = rho0*(f-1) + f*X;
    }
  }

  /* convert back to slopes */
  f = (1+X/rho0)/EXSQRT(sqr(1+DPoP)-sqr(QX)-sqr(QY), sqrtOrder);
  Qf[1] *= f;
  Qf[3] *= f;
  Qf[4] += dist;
}


void integrate_csbend_ord4(double *Qf, double *Qi, double *sigmaDelta2, double s, long n, long sqrtOrder, double rho0, double p0)
{
  long i;
  double factor, f, phi, ds, dsh, dp, dist;
  double Fx, Fy, x, y;
  double sine, cosi, tang;
  double sin_phi, cos_phi;
  
#define X0 Qi[0]
#define XP0 Qi[1]
#define Y0 Qi[2]
#define YP0 Qi[3]
#define S0 Qi[4]
#define DPoP0 Qi[5]

#define X Qf[0]
#define QX Qf[1]
#define Y Qf[2]
#define QY Qf[3]
#define S Qf[4]
#define DPoP Qf[5]

  /* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

  if (!Qf)
    bombElegant("NULL final coordinates pointer ()", NULL);
  if (!Qi)
    bombElegant("NULL initial coordinates pointer (integrate_csbend_ord4)", NULL);
  if (n<1)
    bombElegant("invalid number of steps (integrate_csbend_ord4)", NULL);

  /* calculate canonical momenta (scaled to central momentum) */
  dp = DPoP0;
  f = (1+dp)/EXSQRT(sqr(1+X0/rho0) + sqr(XP0) + sqr(YP0), sqrtOrder);
  QX = XP0*f;
  QY = YP0*f;

  X = X0;
  Y = Y0;
  S = S0;
  DPoP = DPoP0;
  
  dist = 0;

  s /= n;
  for (i=0; i<n; i++) {
    
    /* do first drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do first kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;

    computeCSBENDFields(&Fx, &Fy, x, y);
    
    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) {
      addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, rad_coef, s/3, isrConstant,
                       distributionBasedRadiation, includeOpeningAngle,
                       meanPhotonsPerMeter0, normalizedCriticalEnergy0, p0);
    }

    /* do second drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do second kick */
    ds = -s*BETA/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    computeCSBENDFields(&Fx, &Fy, x, y);

    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant)
      addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, rad_coef, s/3, isrConstant,
                       distributionBasedRadiation, includeOpeningAngle,
                       meanPhotonsPerMeter0, normalizedCriticalEnergy0, p0);

    /* do third drift */
    dsh = s*(1-BETA)/(2-BETA)/2;
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
    
    /* do third kick */
    ds = s/(2-BETA);
    /* -- calculate the scaled fields */
    x = X;
    y = Y;
    computeCSBENDFields(&Fx, &Fy, x, y);

    /* --do kicks */
    QX += -ds*(1+X/rho0)*Fy/rho_actual;
    QY += ds*(1+X/rho0)*Fx/rho_actual;
    if (rad_coef || isrConstant) 
      addRadiationKick(&QX, &QY, &DPoP, sigmaDelta2, sqrtOrder,
		       X, 1./rho0, Fx, Fy, 
		       ds, rad_coef, s/3, isrConstant,
                       distributionBasedRadiation, includeOpeningAngle,
                       meanPhotonsPerMeter0, normalizedCriticalEnergy0, p0);
    
    /* do fourth drift */
    dsh = s/2/(2-BETA);
    if ((f=sqr(1+DPoP)-sqr(QY))<=0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    f = EXSQRT(f, sqrtOrder);
    if (fabs(QX/f)>1) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    phi = asin(sin_phi=QX/f);
    sine = sin(dsh/rho0+phi);
    if ((cosi = cos(dsh/rho0+phi))==0) {
      particle_lost = 1;
      s_lost = dist;
      return;
    }
    tang = sine/cosi;
    cos_phi = cos(phi);
    QX = f*sine;
    Y += QY*(factor=(rho0+X)*cos_phi/f*(tang-sin_phi/cos_phi));
    dist += factor*(1+DPoP);
    f = cos_phi/cosi;
    X  = rho0*(f-1) + f*X;
  }

  /* convert back to slopes */
  f = (1+X/rho0)/EXSQRT(sqr(1+DPoP)-sqr(QX)-sqr(QY), sqrtOrder);
  Qf[1] *= f;
  Qf[3] *= f;
  Qf[4] += dist;
}

typedef struct {
  unsigned long lastMode;
#define CSRDRIFT_STUPAKOV          0x0001UL
#define CSRDRIFT_SALDIN54          0x0002UL
#define CSRDRIFT_OVERTAKINGLENGTH  0x0004UL
#define CSRDRIFT_ATTENUATIONLENGTH 0x0008UL
#define CSRDRIFT_SPREAD            0x0010UL
  long bins, valid;
  double dctBin, s0, ds0, zLast, z0;
  double S11, S12, S22;
  double *dGamma;
  double rho, bendingAngle, Po, perc68BunchLength, perc90BunchLength, peakToPeakWavelength, rmsBunchLength;
  /* for Saldin eq 54 (NIM A 398 (1997) 373-394) mode: */
  FILE *fpSaldin;
  long nSaldin;
  double *FdNorm;   /* Saldin's Fd(sh, x)/Fd(sh, 0), sh = bunch-length*gamma^3/rho */
  double *xSaldin;  /* distance from end of bend */
  double lastFdNorm; /* last value obtained from interpolation */
  /* for Stupakov mode */
  long SGOrder, SGHalfWidth, SGDerivHalfWidth, SGDerivOrder;
  double binRangeFactor;
  double GSConstant, MPCharge;
  char *StupakovOutput;
  SDDS_DATASET SDDS_Stupakov;
  long StupakovFileActive, StupakovOutputInterval;
  long trapazoidIntegration;
  double lowFrequencyCutoff0, lowFrequencyCutoff1;
  double highFrequencyCutoff0, highFrequencyCutoff1;
  long clipNegativeBins;
  long wffValues;
  double *wffFreqValue, *wffRealFactor, *wffImagFactor;
} CSR_LAST_WAKE;
CSR_LAST_WAKE csrWake;

#define DERBENEV_CRITERION_DISABLE 0
#define DERBENEV_CRITERION_EVAL 1
#define DERBENEV_CRITERION_ENFORCE 2
#define N_DERBENEV_CRITERION_OPTIONS 3
static char *derbenevCriterionOption[N_DERBENEV_CRITERION_OPTIONS] = {
  "disable", "evaluate", "enforce"};

void readWakeFilterFile(long *values, double **freq, double **real, double **imag, 
                        char *freqName, char *realName, char *imagName,
                        char *filename);

long track_through_csbendCSR(double **part, long n_part, CSRCSBEND *csbend, double p_error, 
                             double Po, double **accepted, double z_start, double z_end,
                             CHARGE *charge, char *rootname)
{
  double h, n, he1, he2;
  static long csrWarning = 0;
  static double *beta0=NULL, *ctHist=NULL, *ctHistDeriv=NULL;
  static double *dGamma=NULL, *T1=NULL, *T2=NULL, *denom=NULL;
  static long maxParticles = 0, maxBins = 0 ;
  static char *particleLost=NULL;
  double x=0, xp, y=0, yp, p1, beta1, p0;
  double ctLower, ctUpper, dct, slippageLength, phiBend, slippageLength13;
  long diSlippage, diSlippage4;
  long nBins, nBinned = 0;
  long i_part, i_top, kick;
  double rho=0.0, Fx, Fy;
  double fse, dp_prime;
  double tilt, etilt, cos_ttilt, sin_ttilt, ttilt;
  double *coord;
  double angle, e1, e2, Kg;
  double psi1, psi2;
  double Qi[6], Qf[6];
  double dcoord_etilt[6];
  double dxi, dyi, dzi;
  double dxf, dyf, dzf;
  double delta_xp;
  double macroParticleCharge, CSRConstant;
  long iBin, iBinBehind;
  long csrInhibit = 0, largeRhoWarning = 0;
  double derbenevRatio = 0;
  long n_partMoreThanOne = 0;
  TRACKING_CONTEXT tContext;
  VMATRIX *Msection=NULL, *Me1=NULL, *Me2=NULL;
  static double accumulatedAngle = 0;
  short accumulatingAngle = 1;
  
#if USE_MPI
  double *buffer;  
  if (notSinglePart)
    n_partMoreThanOne = 1; /* This is necessary to solve synchronization issue in parallel version*/
  else
    if (n_part > 1) n_partMoreThanOne = 1;	
#else
  if (n_part > 1) n_partMoreThanOne = 1;
#endif

  if (!(csbend->edgeFlags&SAME_BEND_PRECEDES))
    accumulatedAngle = accumulatingAngle = 0;
  
  csrWake.valid = 0;
  if (isSlave || !notSinglePart) 
    reset_driftCSR();

  getTrackingContext(&tContext);
  
  if (!csbend)
    bombElegant("null CSBEND pointer (track_through_csbend)", NULL);

  if (csbend->angle==0) {
    if (!csbend->useMatrix)
      exactDrift(part, n_part, csbend->length); 
    else {
      long i;
      if (isSlave || !notSinglePart) {
        for (i=0; i<n_part; i++) {
          part[i][0] += csbend->length*part[i][1];
          part[i][2] += csbend->length*part[i][3];
          part[i][4] += csbend->length;
        }
      }
    }
    return n_part;
  }

  if (csbend->integration_order!=2 && csbend->integration_order!=4)
    bombElegant("CSBEND integration_order is invalid--must be either 2 or 4", NULL);

  macroParticleCharge = 0;
  if (charge) {
    macroParticleCharge = charge->macroParticleCharge;
  } else if (csbend->bins && !csrWarning && csbend->csr) {
    fprintf(stdout, "Warning: you asked for CSR on CSBEND but didn't give a CHARGE element\n");
    fflush(stdout);
    csrWarning = 1;
  }
  
  if ((nBins=csbend->bins)<2)
    bombElegant("Less than 2 bins for CSR!", NULL);

  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = csbend->SGHalfWidth;
  if (csbend->SGDerivHalfWidth<=0)
    csbend->SGDerivHalfWidth = 1;

  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = csbend->SGOrder;
  if (csbend->SGDerivOrder<=0)
    csbend->SGDerivOrder = 1;
  
  if (isSlave || !notSinglePart) 
    if (n_part>maxParticles &&
	(!(beta0=SDDS_Realloc(beta0, sizeof(*beta0)*(maxParticles=n_part))) ||
	 !(particleLost=SDDS_Realloc(particleLost, sizeof(*particleLost)*n_part))))
      bombElegant("Memory allocation failure (track_through_csbendCSR)", NULL);

  rho0 = csbend->length/csbend->angle;
  if (csbend->use_bn) {
    csbend->b[0] = csbend->b1;
    csbend->b[1] = csbend->b2;
    csbend->b[2] = csbend->b3;
    csbend->b[3] = csbend->b4;
    csbend->b[4] = csbend->b5;
    csbend->b[5] = csbend->b6;
    csbend->b[6] = csbend->b7;
    csbend->b[7] = csbend->b8;
  } else {
    csbend->b[0] = csbend->k1*rho0;
    csbend->b[1] = csbend->k2*rho0;
    csbend->b[2] = csbend->k3*rho0;
    csbend->b[3] = csbend->k4*rho0;
    csbend->b[4] = csbend->k5*rho0;
    csbend->b[5] = csbend->k6*rho0;
    csbend->b[6] = csbend->k7*rho0;
    csbend->b[7] = csbend->k8*rho0;
  }

  he1 = csbend->h1;
  he2 = csbend->h2;
  if (csbend->angle<0) {
    long i;
    angle = -csbend->angle;
    e1    = -csbend->e1;
    e2    = -csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt + PI;
    rho0  = csbend->length/angle;
    for (i=0; i<8; i+=2)
      csbend->b[i] *= -1;
  }
  else {
    angle = csbend->angle;
    e1    = csbend->e1;
    e2    = csbend->e2;
    etilt = csbend->etilt;
    tilt  = csbend->tilt;
    rho0  = csbend->length/angle;
  }
  
  if (rho0>1e6) {
    if (!largeRhoWarning) {
      printf("Warning: One or more CSRCSBENDs have radius > 1e6.  Treated as drift.\n");
      largeRhoWarning = 1;
    }
    exactDrift(part, n_part, csbend->length);
    return n_part;
  }
  
  h = 1/rho0;
  n = -csbend->b[0]/h;
  fse = csbend->fse;
  if (fse>-1)
    rho_actual = 1/((1+fse)*h);
  else
    rho_actual = 1e16/h;

  /* angles for fringe-field effects */
  Kg   = 2*csbend->hgap*csbend->fint;
  psi1 = Kg/rho_actual/cos(e1)*(1+sqr(sin(e1)));
  psi2 = Kg/rho_actual/cos(e2)*(1+sqr(sin(e2)));

  /* rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of
   * the central particle.
   */
  if (csbend->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)*sqr(1+fse)/(6*PI*epsilon_o*sqr(c_mks)*particleMass*sqr(rho0));
  else
    rad_coef = 0;
  /* isrConstant is the RMS increase in dP/P per meter due to incoherent SR.  */
  if (csbend->isr && (n_part>1 || !csbend->isr1Particle)) 
    isrConstant = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*
                              137.0359895/pow3(fabs(rho_actual)));
  else
    isrConstant = 0;

  distributionBasedRadiation = 0;
  
  if (csbend->useMatrix) {
    csbend->nonlinear = 0;
    Me1 = edge_matrix(e1, 1./(rho0/(1+csbend->fse)), 0.0, n, -1, Kg, 1, 0, 0);
    Msection = bend_matrix(csbend->length/csbend->n_kicks, 
                                   angle/csbend->n_kicks, 0.0, 0.0, 
                                   0.0, 0.0, csbend->b[0]*h,  0.0,
                                   0.0, 0.0, 0.0, csbend->fse, csbend->etilt, 1, 1, 0, 0);
    Me2 = edge_matrix(e2, 1./(rho0/(1+csbend->fse)), 0.0, n, 1, Kg, 1, 0, 0);
  }
  computeCSBENDFieldCoefficients(csbend->b, h, csbend->nonlinear, csbend->expansionOrder);

  ttilt = tilt + etilt;
  if (ttilt==0) {
    cos_ttilt = 1;
    sin_ttilt = 0;
  }
  else if (fabs(fabs(ttilt)-PI)<1e-12) {
    cos_ttilt = -1;
    sin_ttilt = 0;
  }
  else if (fabs(ttilt-PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = 1;
  }
  else if (fabs(ttilt+PIo2)<1e-12) {
    cos_ttilt = 0;
    sin_ttilt = -1;
  }
  else {
    cos_ttilt = cos(ttilt);
    sin_ttilt = sin(ttilt);
  }


  if (etilt) {
    /* compute final offsets due to error-tilt of the magnet */
    /* see pages 90-93 of notebook 1 about this */
    double q1a, q2a, q3a;
    double q1b, q2b, q3b;
    double qp1, qp2, qp3; 
    double dz, tan_alpha, k;

    q1a = (1-cos(angle))*rho0*(cos(etilt)-1);
    q2a = 0;
    q3a = (1-cos(angle))*rho0*sin(etilt);
    qp1 = sin(angle)*cos(etilt);
    qp2 = cos(angle);
    k = sqrt(sqr(qp1)+sqr(qp2));
    qp1 /= k;
    qp2 /= k;
    qp3 = sin(angle)*sin(etilt)/k;
    tan_alpha = 1./tan(angle)/cos(etilt);
    q1b = q1a*tan_alpha/(tan(angle)+tan_alpha);
    q2b = -q1b*tan(angle);
    dz  = sqrt(sqr(q1b-q1a)+sqr(q2b-q2a));
    q3b = q3a + qp3*dz;

    dcoord_etilt[0] = sqrt(sqr(q1b) + sqr(q2b));
    dcoord_etilt[1] = tan(atan(tan_alpha)-(PIo2-angle));
    dcoord_etilt[2] = q3b;
    dcoord_etilt[3] = qp3;
    dcoord_etilt[4] = dz*sqrt(1+sqr(qp3));
    dcoord_etilt[5] = 0;

    /* rotate by tilt to get into same frame as bend equations. */
    rotate_coordinates(dcoord_etilt, tilt);
  }
  else
    fill_double_array(dcoord_etilt, 6L, 0.0);

  dxi = -csbend->dx;
  dzi =  csbend->dz;
  dyi = -csbend->dy;

  /* must use the original angle here because the translation is done after
   * the final rotation back
   */
  dxf =  csbend->dx*cos(csbend->angle) + csbend->dz*sin(csbend->angle);
  dzf =  csbend->dx*sin(csbend->angle) - csbend->dz*cos(csbend->angle);
  dyf = csbend->dy;

  if (isMaster) {
  if (csbend->particleOutputFile && strlen(csbend->particleOutputFile) && !csbend->particleFileActive) {
    /* set up SDDS output file for particle coordinates inside bend */
    csbend->particleFileActive = 1;
    csbend->particleOutputFile = compose_filename(csbend->particleOutputFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSpart, SDDS_BINARY, 1, 
                               NULL, NULL, csbend->particleOutputFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSpart, "Angle", NULL, SDDS_DOUBLE) ||
        (csbend->xIndex=SDDS_DefineColumn(&csbend->SDDSpart, "x", NULL, "m", 
                                          NULL, NULL, SDDS_DOUBLE, 0 ))<0 ||
        (csbend->xpIndex=SDDS_DefineColumn(&csbend->SDDSpart, "xp", NULL, NULL, 
                                           NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->tIndex=SDDS_DefineColumn(&csbend->SDDSpart, "t", NULL, "s", 
                                          NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (csbend->pIndex=SDDS_DefineColumn(&csbend->SDDSpart, "p", NULL, "m$be$nc", 
                                          NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        !SDDS_WriteLayout(&csbend->SDDSpart)) {
      SDDS_SetError("Problem setting up particle output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  }
  
  if (isMaster) { 
  if (csbend->histogramFile && strlen(csbend->histogramFile) && !csbend->wakeFileActive) {
    /* set up SDDS output file for CSR monitoring */
    csbend->wakeFileActive = 1;
    csbend->histogramFile = compose_filename(csbend->histogramFile, rootname);
    if (!SDDS_InitializeOutput(&csbend->SDDSout, SDDS_BINARY, 1, NULL, NULL, csbend->histogramFile) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Kick", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "pCentral", "m$be$nc", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "Angle", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "SlippageLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "TotalBunchLength", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "BinSize", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "dsKick", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&csbend->SDDSout, "DerbenevRatio", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "GammaDeriv", "1/m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(&csbend->SDDSout, "DeltaGammaT2", NULL, SDDS_DOUBLE) ||
        !SDDS_WriteLayout(&csbend->SDDSout)) {
      SDDS_SetError("Problem setting up wake output file for CSR");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  }
  if (csbend->wakeFilterFile && strlen(csbend->wakeFilterFile) && !csbend->wffValues) 
    readWakeFilterFile(&csbend->wffValues,
                       &csbend->wffFreqValue, &csbend->wffRealFactor, &csbend->wffImagFactor, 
                       csbend->wffFreqColumn, csbend->wffRealColumn, csbend->wffImagColumn,
                       csbend->wakeFilterFile);
  
  /*  prepare arrays for CSR integrals */
  nBins = csbend->bins;
  if (!(ctHist=SDDS_Realloc(ctHist, sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=SDDS_Realloc(ctHistDeriv, sizeof(*ctHistDeriv)*nBins)) ||
      !(denom=SDDS_Realloc(denom, sizeof(*denom)*nBins)) ||
      !(T1=SDDS_Realloc(T1, sizeof(*T1)*nBins)) ||
      !(T2=SDDS_Realloc(T2, sizeof(*T2)*nBins)) ||
      !(dGamma=SDDS_Realloc(dGamma, sizeof(*dGamma)*nBins)))
    bombElegant("memory allocation failure (track_through_csbendCSR)", NULL);

  /* prepare some data for CSRDRIFT */
  csrWake.dGamma = dGamma;
  csrWake.bins = nBins;
  csrWake.ds0 = csbend->length/csbend->n_kicks;
  csrWake.zLast = csrWake.z0 = z_end;
  csrWake.highFrequencyCutoff0 = csbend->highFrequencyCutoff0;
  csrWake.highFrequencyCutoff1 = csbend->highFrequencyCutoff1;
  csrWake.lowFrequencyCutoff0 = csbend->lowFrequencyCutoff0;
  csrWake.lowFrequencyCutoff1 = csbend->lowFrequencyCutoff1;
  csrWake.clipNegativeBins = csbend->clipNegativeBins;
  csrWake.wffValues = csbend->wffValues;
  csrWake.wffFreqValue = csbend->wffFreqValue;
  csrWake.wffRealFactor = csbend->wffRealFactor;
  csrWake.wffImagFactor = csbend->wffImagFactor;
  
#if !defined(PARALLEL)  
  multipoleKicksDone += n_part*csbend->n_kicks*(csbend->integration_order==4?4:1);
#endif

  if (isSlave || !notSinglePart) {
  /* check particle data, transform coordinates, and handle edge effects */
  for (i_part=0; i_part<n_part; i_part++) {
    if (!part) {
      fprintf(stdout, "error: null particle array found (working on particle %ld) (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (!(coord = part[i_part])) {
      fprintf(stdout, "error: null coordinate pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "error: null accepted particle pointer for particle %ld (track_through_csbend)\n", i_part);
      fflush(stdout);
      abort();
    }

    /* adjust for element offsets */
    coord[4] += dzi*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
    coord[0]  = coord[0] + dxi + dzi*coord[1];
    coord[2]  = coord[2] + dyi + dzi*coord[3];

    /* perform tilt transformations and save some data */
    x  =  coord[0]*cos_ttilt + coord[2]*sin_ttilt;
    y  = -coord[0]*sin_ttilt + coord[2]*cos_ttilt;
    coord[0] = x;
    coord[2] = y;
    xp =  coord[1]*cos_ttilt + coord[3]*sin_ttilt;
    yp = -coord[1]*sin_ttilt + coord[3]*cos_ttilt;
    coord[1] = xp;
    coord[3] = yp;
    p0 = Po*(1+coord[5]);
    beta0[i_part] = p0/sqrt(p0*p0+1);
    coord[4] /= beta0[i_part];
    particleLost[i_part] = 0;

#undef X
#undef Y
#define X coord[0]
#define Y coord[2]
#define XP coord[1]
#define YP coord[3]
#define CT coord[4]
#define DP coord[5]
    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS) {
      /* apply edge focusing */
      if (csbend->useMatrix)
        track_particles(&coord, Me1, &coord, 1);
      else {
        rho = (1+DP)*rho_actual;
        if (csbend->useMatrix || csbend->edge_order<2) {
          delta_xp = tan(e1)/rho*X;
          XP += delta_xp;
          YP -= tan(e1-psi1/(1+DP))/rho*Y;
        }
        else
          apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e1, he1, psi1*(1+DP), -1);
      }
    }

    if (!csbend->useMatrix) {
      /* transform to curvilinear coordinates */
      XP *= (1+X/rho0);
      YP *= (1+X/rho0);
    }
    if (csbend->edgeFlags&BEND_EDGE1_EFFECTS && e1!=0 && rad_coef) {
      /* pre-adjust dp/p to anticipate error made by integrating over entire sector */
      computeCSBENDFields(&Fx, &Fy, x, y);

      dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
        sqrt(sqr(1+X/rho0)+sqr(XP)+sqr(YP));
      DP -= dp_prime*X*tan(e1);
    }
  }
  }
  if (csbend->csr && n_partMoreThanOne)
    CSRConstant = 2*macroParticleCharge*particleCharge/pow(3*rho0*rho0, 1./3.)/(4*PI*epsilon_o*particleMass*sqr(c_mks));
  else
    CSRConstant = 0;
  /* Now do the body of the sector dipole */
  phiBend = accumulatedAngle; 
  for (kick=0; kick<csbend->n_kicks; kick++) {
    if (isSlave || !notSinglePart) {
      for (i_part=0; i_part<n_part; i_part++) {
	coord = part[i_part];
	if (particleLost[i_part])
	  continue;

	if (csbend->useMatrix) {
	  track_particles(&coord, Msection, &coord, 1);
	} else {
	  /* load input coordinates into arrays */
	  Qi[0] = X;
	  Qi[1] = XP;
	  Qi[2] = Y;
	  Qi[3] = YP;
	  Qi[4] = 0;  
	  Qi[5] = DP;
        
	  if (csbend->integration_order==4)
	    integrate_csbend_ord4(Qf, Qi, NULL, csbend->length/csbend->n_kicks, 1, 0, rho0, Po);
	  else
	    integrate_csbend_ord2(Qf, Qi, NULL, csbend->length/csbend->n_kicks, 1, 0, rho0, Po);
	  particleLost[i_part] = particle_lost;
      
	  /* retrieve coordinates from arrays */
	  X  = Qf[0];  
	  XP = Qf[1];  
	  Y  = Qf[2];  
	  YP = Qf[3];  
	  DP = Qf[5];

	  if (rad_coef || isrConstant) {
	    /* convert additional distance traveled to ct using mean velocity */
	    p1 = Po*(1+DP);
	    beta1 = p1/sqrt(p1*p1+1);
	    CT += Qf[4]*2/(beta0[i_part]+beta1);
	    beta0[i_part] = beta1;
	  } else {
	    CT += Qf[4]/beta0[i_part];  
	  }
	}
      }    
    }

    if (n_partMoreThanOne && csbend->derbenevCriterionMode) {
      /* evaluate Derbenev criterion from TESLA-FEL 1995-05: sigma_x/sigma_z << (R/sigma_z)^(1/3) */
      long code;
      double Sz, Sx;
      switch (code=match_string(csbend->derbenevCriterionMode, derbenevCriterionOption, N_DERBENEV_CRITERION_OPTIONS, 0)) {
      case DERBENEV_CRITERION_DISABLE:
	break;
      case DERBENEV_CRITERION_EVAL:
      case DERBENEV_CRITERION_ENFORCE:
#if !USE_MPI
	rms_emittance(part, 4, 5, n_part, &Sz, NULL, NULL);
	rms_emittance(part, 0, 1, n_part, &Sx, NULL, NULL);
#else
     if (notSinglePart) {
        /* The master will get the result from the rms_emittance routine */
	rms_emittance_p(part, 4, 5, n_part, &Sz, NULL, NULL);
	rms_emittance_p(part, 0, 1, n_part, &Sx, NULL, NULL);
     } else {
        rms_emittance(part, 4, 5, n_part, &Sz, NULL, NULL);
        rms_emittance(part, 0, 1, n_part, &Sx, NULL, NULL);
     }
#endif
	Sz = sqrt(Sz);
	Sx = sqrt(Sx);
	derbenevRatio = (Sx/Sz)/pow(rho0/Sz, 1./3.);
	if (derbenevRatio>0.1) {
	  if (code==DERBENEV_CRITERION_EVAL)
	    fprintf(stderr, "Warning: Using 1-D CSR formalism but Derbenev criterion not satisfied (%le > 0.1).\n",
		    derbenevRatio);
	  else {
	    csrInhibit = 1;
	    fprintf(stderr, "Warning: Derbenev criterion not satisfied (%le > 0.1)---not applying CSR\n",
		    derbenevRatio);
	  }
	}
	break;
      default:
	fprintf(stderr, "Error: invalid value for DERBENEV_CRITERION_MODE. Give 'disable', 'evaluate', or 'enforce'\n");
	exit(1);
	break;
      }
    }
    

#if (!USE_MPI)
    if (n_partMoreThanOne && !csrInhibit) {
#else
      if (!csrInhibit && (notSinglePart || (!notSinglePart && n_partMoreThanOne))) { /* n_part could be 0 for some processors, which could cause synchronization problem */
#endif
      /* compute CSR potential function */
      if (kick==0 || !csbend->binOnce) {
        /* - first make a density histogram */
        ctLower = ctUpper = dct = 0;
	nBinned = binParticleCoordinate(&ctHist, &maxBins,
                                   &ctLower, &ctUpper, &dct, &nBins, 
                                   csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor, 
					part, n_part, 4);
#if (!USE_MPI) 
	if (nBinned != n_part) {
          fprintf(stdout, "Only %ld of %ld particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		  nBinned, n_part, z_start, kick, csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
	  fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
		  ctLower, ctUpper, dct, nBins, maxBins);
          fflush(stdout);
        }
#else
     if (notSinglePart) {
	if (USE_MPI) {
	  long all_binned, result = 1, nBinned_total;

          if (isSlave || !notSinglePart) {
	    result = ((nBinned==n_part) ? 1 : 0);
	  }
	  MPI_Allreduce(&result, &all_binned, 1, MPI_LONG, MPI_LAND, MPI_COMM_WORLD);
	  MPI_Allreduce(&nBinned, &nBinned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
          nBinned = nBinned_total; 
	  if (!all_binned && isMaster) {
	    fprintf(stdout, "Not all particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		    z_start, kick,
		  csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
	  fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
		  ctLower, ctUpper, dct, nBins, maxBins);
          fflush(stdout);
        }
        }

	if (USE_MPI) {  /* Master needs to know the information to write the result */
	  buffer = malloc(sizeof(double) * nBins);
	  MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  memcpy(ctHist, buffer, sizeof(double)*nBins);
	  free(buffer);
	}
     }
#endif
        
        /* - smooth the histogram, normalize to get linear density, and 
           copy in preparation for taking derivative
           */
        if (csbend->highFrequencyCutoff0>0 || csbend->lowFrequencyCutoff0>=0) {
          long nz;
          nz = applyLHPassFilters(ctHist, nBins, 
                                  csbend->lowFrequencyCutoff0, csbend->lowFrequencyCutoff1,
                                  csbend->highFrequencyCutoff0, csbend->highFrequencyCutoff1,
                                  csbend->clipNegativeBins);
          if (nz && negativeWarningsLeft) {
	    fprintf(stdout, "Warning: low pass filter resulted in negative values in %ld bins\n", nz);
            if (--negativeWarningsLeft==0)
              fprintf(stdout, "         Further warnings will be suppressed for this run.\n");
            fflush(stdout);
          }
        }
        if (csbend->SGHalfWidth>0) {
          SavitzyGolaySmooth(ctHist, nBins, csbend->SGOrder, csbend->SGHalfWidth, csbend->SGHalfWidth,  0);
          correctDistribution(ctHist, nBins, 1.0*nBinned);
        }
        for (iBin=0; iBin<nBins; iBin++) {
          denom[iBin] = pow(dct*iBin, 1./3.);
          ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
        }
        /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
         * I won't scale it now as it will just fall out in the integral 
         */
        SavitzyGolaySmooth(ctHistDeriv, nBins, csbend->SGDerivOrder, 
                           csbend->SGDerivHalfWidth, csbend->SGDerivHalfWidth, 1);
      } else {
        ctLower += rho0*angle/csbend->n_kicks;
        ctUpper += rho0*angle/csbend->n_kicks;
      }
      
      
      phiBend += angle/csbend->n_kicks;
      slippageLength = rho0*ipow(phiBend, 3)/24.0;
      slippageLength13 = pow(slippageLength, 1./3.);
      diSlippage = slippageLength/dct;
      diSlippage4 = 4*slippageLength/dct;
      if (kick==0 || !csbend->binOnce) {
	for (iBin=0; iBin<nBins; iBin++) {
	  double term1, term2;
	  long count;
	  T1[iBin] = T2[iBin] = 0;
	  term1 = term2 = 0;
	  if (CSRConstant) {
	    if (csbend->steadyState) {
	      if (!csbend->trapazoidIntegration) {
		for (iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++)
		  T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
	      }
	      else {
		if ((iBinBehind=iBin+1)<nBins)
		  term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
		for (count=0, iBinBehind=iBin+1; iBinBehind<nBins; iBinBehind++, count++)
		  T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
		if ((iBin+1)<nBins)
		  T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
		if (count>1)
		  T1[iBin] -= (term1+term2)/2;
	      }
	    } else {
	      if (!csbend->trapazoidIntegration) {
		for (iBinBehind=iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; iBinBehind++)
		  T1[iBin] += ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin];
	      }
	      else {
		if ((iBinBehind = iBin+1)<nBins && iBinBehind<=(iBin+diSlippage))
		  term1 = ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]/2;
		for (count=0, iBinBehind = iBin+1; iBinBehind<=(iBin+diSlippage) && iBinBehind<nBins; 
		     count++, iBinBehind++)
		  T1[iBin] += (term2=ctHistDeriv[iBinBehind]/denom[iBinBehind-iBin]);
		if (diSlippage>0 && (iBin+1)<nBins)
		  T1[iBin] += 0.3*sqr(denom[1])*(2*ctHistDeriv[iBin+1]+3*ctHistDeriv[iBin])/dct;
		if (count>1)
		  T1[iBin] -= (term1+term2)/2;
	      }
	      if ((iBin+diSlippage)<nBins)
		T2[iBin] += ctHist[iBin+diSlippage];
	      if ((iBin+diSlippage4)<nBins)
		T2[iBin] -= ctHist[iBin+diSlippage4];
	    }
	    /* there is no negative sign here because my derivative is w.r.t. -s
	       in notation of Saldin, et. al. */
	    T1[iBin] *= CSRConstant*csbend->length/csbend->n_kicks; 
	    /* keep the negative sign on this term, which has no derivative */
	    T2[iBin] *= -CSRConstant*csbend->length/csbend->n_kicks/slippageLength13;
	  }
	  dGamma[iBin] = T1[iBin]+T2[iBin];
	}

	if (csbend->wffValues) 
	  applyFilterTable(dGamma, nBins, dct/c_mks, csbend->wffValues, csbend->wffFreqValue,
			   csbend->wffRealFactor, csbend->wffImagFactor);
      }
      if (isSlave || !notSinglePart) {
	if (CSRConstant) {
	  for (i_part=0; i_part<n_part; i_part++) {
	    long nBins1;
	    nBins1 = nBins-1;
	    coord = part[i_part];
	    if (!particleLost[i_part]) {
	      double f;
	      /* apply CSR kick */
	      iBin = (f=(CT-ctLower)/dct);
	      f -= iBin;
	      if (iBin>=0 && iBin<nBins1)
		DP += ((1-f)*dGamma[iBin]+f*dGamma[iBin+1])/Po;
	    }
	  }
	}
      }
  
      if (csbend->particleFileActive && kick%csbend->particleOutputInterval==0) {
	if (isMaster) {
        long ip;
        /* dump particle data at this location */
        if (!SDDS_StartPage(&csbend->SDDSpart, n_part) ||
            !SDDS_SetParameters(&csbend->SDDSpart, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Pass", -1, "Kick", kick, "pCentral", Po, "Angle", phiBend, 
                                NULL))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        for (ip=0; ip<n_part; ip++) {
          if (!SDDS_SetRowValues(&csbend->SDDSpart, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                                 ip, 
                                 csbend->xIndex, part[ip][0],
                                 csbend->xpIndex, part[ip][1],
                                 csbend->tIndex, part[ip][4],
                                 csbend->pIndex, Po*(1+part[ip][5]),
                                 -1)) 
            SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }
        if (!SDDS_WritePage(&csbend->SDDSpart))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!inhibitFileSync)
          SDDS_DoFSync(&csbend->SDDSpart);
	  printf ("Pelegant does not support dumping particle data inside an element now.");
	}
      }

      if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	  kick!=(csbend->n_kicks-1) &&
	  (csbend->sliceAnalysisInterval==0 ||
	   kick%csbend->sliceAnalysisInterval==0)) {
#if (!USE_MPI)
	convertFromCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
	performSliceAnalysisOutput(tContext.sliceAnalysis, part, n_part, 
				   0, tContext.step, Po, 
				   macroParticleCharge*n_part,
				   tContext.elementName, 
				   z_start + (kick*(z_end-z_start))/(csbend->n_kicks-1),
				   1);
	convertToCSBendCoords(part, n_part, rho0, cos_ttilt, sin_ttilt, 1);
#else 
      if (isMaster) 
	printf ("Pelegant does not support slice analysis output inside an element now.");
    
#endif
      }

      if (csbend->wakeFileActive && 
          ((!csbend->outputLastWakeOnly && kick%csbend->outputInterval==0) ||
           (csbend->outputLastWakeOnly && kick==(csbend->n_kicks-1)))) {
        /* scale the linear density and its derivative to get C/s and C/s^2 
         * ctHist is already normalized to dct, but ctHistDeriv requires an additional factor
         */
        for (iBin=0; iBin<nBins; iBin++) {
          ctHist[iBin] *= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] *= macroParticleCharge*sqr(c_mks)/dct;
        }
 
	if (isMaster) {
        if (!SDDS_StartPage(&csbend->SDDSout, nBins) ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, dGamma, nBins, "DeltaGamma") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "DeltaGammaT1") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "DeltaGammaT2") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "LinearDensity") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "LinearDensityDeriv") ||
            !SDDS_SetParameters(&csbend->SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Pass", -1, "Kick", kick, "dsKick", csbend->length/csbend->n_kicks,
                                "pCentral", Po, "Angle", phiBend, "SlippageLength", slippageLength,
                                "TotalBunchLength", ctUpper-ctLower,
                                "BinSize", dct, 
                                "DerbenevRatio", derbenevRatio, NULL))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	}
        if (csbend->binOnce) {
          /* fix these arrays so they can be used again */
          ctHist[iBin] /= macroParticleCharge*c_mks;
          ctHistDeriv[iBin] /= macroParticleCharge*sqr(c_mks)/dct;
        }
        /* use T1 array to output s and T2 to output dGamma/ds */
        for (iBin=0; iBin<nBins; iBin++) {
          T1[iBin] = ctLower-(ctLower+ctUpper)/2.0+dct*(iBin+0.5);
          T2[iBin] = dGamma[iBin]/(csbend->length/csbend->n_kicks);
        }
	if (isMaster){
        if (!SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T1, nBins, "s") ||
            !SDDS_SetColumn(&csbend->SDDSout, SDDS_SET_BY_NAME, T2, nBins, "GammaDeriv") ||
            !SDDS_WritePage(&csbend->SDDSout))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!inhibitFileSync)
          SDDS_DoFSync(&csbend->SDDSout);
      }
    }
  }
  }

  if (!csbend->binOnce && n_partMoreThanOne && !csrInhibit && !csbend->csrBlock) {
    /* prepare some data for use by CSRDRIFT element */
    csrWake.dctBin = dct;
    ctLower = ctUpper = dct = 0;

    nBinned =  binParticleCoordinate(&ctHist, &maxBins,
                               &ctLower, &ctUpper, &dct, &nBins, 
                               csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor, 
				     part, n_part, 4);
#if (!USE_MPI)
    if (nBinned!=n_part) {
      fprintf(stdout, "Only %ld of %ld particles binned for CSRCSBEND (z0=%le, end, BRF=%le)\n", 
	      nBinned, n_part, z_start, csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
    }
#else
    if (USE_MPI && notSinglePart) {
      long all_binned, result = 1, nBinned_total;

      if (isSlave || !notSinglePart) {
	result = ((nBinned==n_part) ? 1 : 0);
      }
      else
	nBinned = 0;
      MPI_Allreduce(&result, &all_binned, 1, MPI_LONG, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(&nBinned, &nBinned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      nBinned = nBinned_total; 
      if (!all_binned && isMaster) {
	fprintf(stdout, "Not all particles binned for CSRCSBEND (z0=%le, kick=%ld, BRF=%le)\n", 
		z_start, kick,
	      csbend->binRangeFactor<1.1?1.1:csbend->binRangeFactor);
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
      }
      if (notSinglePart) {  /* Master needs to know the information to write the result */
	buffer = malloc(sizeof(double) * nBins);
	MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	memcpy(ctHist, buffer, sizeof(double)*nBins);
	free(buffer);
      }
    }
#endif
    csrWake.s0 = ctLower + dzf;
  } else {
    ctLower = ctUpper = dct = 0;
    csrWake.dctBin = dct;
    csrWake.s0 = ctLower + dzf;
  }
  
  i_top = n_part-1;
  if (isSlave || !notSinglePart) {
    /* remove lost particles, handle edge effects, and transform coordinates */    
    for (i_part=0; i_part<=i_top; i_part++) {
      coord = part[i_part];
      if (csbend->edgeFlags&BEND_EDGE2_EFFECTS && e2!=0 && rad_coef) {
	/* post-adjust dp/p to correct error made by integrating over entire sector */
        computeCSBENDFields(&Fx, &Fy, x, y);
        
	dp_prime = -rad_coef*(sqr(Fx)+sqr(Fy))*sqr(1+DP)*
	  sqrt(sqr(1+X/rho0)+sqr(XP)+sqr(YP));
	DP -= dp_prime*X*tan(e2);
      }

      /* convert CT to distance traveled at final velocity */
      p1 = Po*(1+DP);
      beta1 = p1/sqrt(sqr(p1)+1);
      coord[4] = CT*beta1;

      if (!csbend->useMatrix) {
	/* transform to cartesian coordinates */
	XP /= (1+X/rho0);
	YP /= (1+X/rho0);
      }
    
      if (particleLost[i_part] || p1<=0) {
	if (!part[i_top]) {
	  fprintf(stdout, "error: couldn't swap particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
		  i_part, i_top);
	  fflush(stdout);
	  abort();
	}
	swapParticles(part[i_part], part[i_top]);
	if (accepted) {
	  if (!accepted[i_top]) {
	    fprintf(stdout, 
		    "error: couldn't swap acceptance data for particles %ld and %ld--latter is null pointer (track_through_csbend)\n",
		    i_part, i_top);
	    fflush(stdout);
	    abort();
	  }
	  swapParticles(accepted[i_part], accepted[i_top]);
	}
	part[i_top][4] = z_start + s_lost;
	part[i_top][5] = Po*(1+part[i_top][5]);
	i_top--;
	i_part--;
	continue;
      }

      if (csbend->edgeFlags&BEND_EDGE2_EFFECTS) {
	if (csbend->useMatrix)
	  track_particles(&coord, Me2, &coord, 1);
	else {
	  /* apply edge focusing */
	  rho = (1+DP)*rho_actual;
	  if (csbend->useMatrix || csbend->edge_order<2) {
	    delta_xp = tan(e2)/rho*X;
	    XP += delta_xp;
	    YP -= tan(e2-psi2/(1+DP))/rho*Y;
	  }
	  else 
	    apply_edge_effects(&X, &XP, &Y, &YP, rho, n, e2, he2, psi2*(1+DP), 1);
	}
      }

      coord = part[i_part];
      x  =  X*cos_ttilt -  Y*sin_ttilt + dcoord_etilt[0];
      y  =  X*sin_ttilt +  Y*cos_ttilt + dcoord_etilt[2];
      xp = XP*cos_ttilt - YP*sin_ttilt + dcoord_etilt[1];
      yp = XP*sin_ttilt + YP*cos_ttilt + dcoord_etilt[3];
      X  = x;
      Y  = y;
      XP = xp;
      YP = yp;
      coord[0] += dxf + dzf*coord[1];
      coord[2] += dyf + dzf*coord[3];
      coord[4] += dzf*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
    }
  }

  if (n_partMoreThanOne && !csbend->csrBlock) {
    /* prepare more data for CSRDRIFT */
    long imin, imax;
    double S55;

#if !USE_MPI    
    rms_emittance(part, 0, 1, i_top+1, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    rms_emittance(part, 4, 5, i_top+1, &S55, NULL, NULL);
#else
    if (notSinglePart) {	
    	rms_emittance_p(part, 0, 1, i_top+1, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    	rms_emittance_p(part, 4, 5, i_top+1, &S55, NULL, NULL);
    } else {
     	rms_emittance(part, 0, 1, i_top+1, &csrWake.S11, &csrWake.S12, &csrWake.S22);
    	rms_emittance(part, 4, 5, i_top+1, &S55, NULL, NULL);
    }
#endif

    csrWake.perc68BunchLength = approximateBeamWidth(0.6826, part, i_top+1, 4)/2;
    csrWake.perc90BunchLength = approximateBeamWidth(0.9, part, i_top+1, 4)/2;
	
    csrWake.rmsBunchLength = sqrt(S55);


#ifdef DEBUG
      fprintf(stderr, "rms bunch length = %le, percentile bunch length (68, 90) = %le, %le\n",
              csrWake.rmsBunchLength, csrWake.perc68BunchLength, csrWake.perc90BunchLength);
#endif
    if (macroParticleCharge) {
      index_min_max(&imin, &imax, csrWake.dGamma, csrWake.bins);
      csrWake.peakToPeakWavelength = 2*fabs(1.0*imax-imin)*dct;
    } else {
      csrWake.peakToPeakWavelength = csrWake.perc68BunchLength;
    }

    csrWake.valid = 1;
    csrWake.rho = rho_actual;
    csrWake.bendingAngle = accumulatingAngle ? fabs(phiBend) : fabs(angle);
    csrWake.Po = Po;
    csrWake.SGOrder = csbend->SGOrder;
    csrWake.SGDerivOrder = csbend->SGDerivOrder;
    csrWake.SGHalfWidth = csbend->SGHalfWidth;
    csrWake.SGDerivHalfWidth = csbend->SGDerivHalfWidth;
    csrWake.GSConstant = CSRConstant*pow(3*rho0*rho0, 1./3.)/2;  /* used for G. Stupakov's drift formulae */
    csrWake.MPCharge = macroParticleCharge;
    csrWake.binRangeFactor = csbend->binRangeFactor;
    csrWake.trapazoidIntegration = csbend->trapazoidIntegration;
    if (csbend->useMatrix) {
      free_matrices(Msection);
      free_matrices(Me1);
      free_matrices(Me2);
    }
  }

  if (csbend->csrBlock)
    accumulatedAngle = 0;
  else
    /* accumulate the bending angle just in case the same type of dipole follows */
    accumulatedAngle += fabs(angle);
    
#if defined(MINIMIZE_MEMORY)
  /* leave dGamma out of this because that memory is used by CSRDRIFT */
  free(beta0);
  free(ctHist);
  free(ctHistDeriv);
  free(T1);
  free(T2);
  free(denom);
  free(particleLost);
  beta0 = ctHist = ctHistDeriv = T1 = T2 = denom = NULL;
  particleLost  = NULL;
  maxBins = maxParticles = 0;
#endif

#if (!USE_MPI)
  return(i_top+1);
#else
  if (isSlave || !notSinglePart)
    return(i_top+1);
  else
    return n_part; /* i_top is not defined for master */
#endif 
}

long binParticleCoordinate(double **hist, long *maxBins,
                           double *lower, double *upper, double *binSize, long *bins,
                           double expansionFactor,
                           double **particleCoord, long nParticles, long coordinateIndex)
{
  long iBin, iParticle, nBinned;
  double value;
  
  if (*binSize<=0 && *bins<1)
    return -1;
  if (*binSize>0 && *bins>1)
    return -2;

  /* if (*lower==*upper)  This condition will be removed */ 
  if (isSlave || !notSinglePart) {
    /* find range of points */
    *upper = -(*lower = DBL_MAX);
    for (iParticle=0; iParticle<nParticles; iParticle++) {
      value = particleCoord[iParticle][coordinateIndex];
      if (value<*lower)
        *lower = value;
      if (value>*upper)
        *upper = value;
    }
  }

#if USE_MPI
  /* find the global maximum and minimum */
  if (notSinglePart) {
    if (isMaster)
      nParticles = 0;
    find_global_min_max(lower, upper, nParticles, MPI_COMM_WORLD);
  }
#endif

    if (expansionFactor>1) {
      double center, range;
      center = (*lower+*upper)/2;
      range = (*upper-*lower)*expansionFactor;
      *lower = center-range/2;
      *upper = center+range/2;
    }
  
  if (*binSize>0)
    /* bin size given, so determine the number of bins */
    *bins = (*upper-*lower)/(*binSize);
  *binSize = (*upper-*lower)/(*bins);

  /* realloc if necessary */
  if (*bins>*maxBins &&
      !(*hist=SDDS_Realloc(*hist, sizeof(**hist)*(*maxBins=*bins))))
    bombElegant("Memory allocation failure (binParticleCoordinate)", NULL);
    
  for (iBin=0; iBin<*bins; iBin++)
    (*hist)[iBin] = 0;
  nBinned = 0;
  if(isSlave || !notSinglePart) {
    for (iParticle=nBinned=0; iParticle<nParticles; iParticle++) {
      /* the coordinate of the bin center is (iBin+0.5)*(*binSize) + *lower */
      iBin = (particleCoord[iParticle][coordinateIndex] - *lower)/(*binSize);
      if (iBin<0 || iBin>(*bins-1))
        continue;
      (*hist)[iBin] += 1;
      nBinned++;
    }
  }
  return nBinned;
}

#if USE_MPI
long binParticleCoordinate_s(double **hist, long *maxBins,
                           double *lower, double *upper, double *binSize, long *bins,
                           double expansionFactor,
                           double **particleCoord, long nParticles, long coordinateIndex)
{
  long iBin, iParticle, nBinned;
  double value;
  
  if (*binSize<=0 && *bins<1)
    return -1;
  if (*binSize>0 && *bins>1)
    return -2;

  /* if (*lower==*upper)  This condition will be removed */ 
  /* find range of points */
  *upper = -(*lower = DBL_MAX);
  for (iParticle=0; iParticle<nParticles; iParticle++) {
    value = particleCoord[iParticle][coordinateIndex];
    if (value<*lower)
      *lower = value;
    if (value>*upper)
      *upper = value;
  }
  if (expansionFactor>1) {
    double center, range;
    center = (*lower+*upper)/2;
    range = (*upper-*lower)*expansionFactor;
    *lower = center-range/2;
    *upper = center+range/2;
  }
  
  if (*binSize>0)
    /* bin size given, so determine the number of bins */
    *bins = (*upper-*lower)/(*binSize);
  *binSize = (*upper-*lower)/(*bins);

  /* realloc if necessary */
  if (*bins>*maxBins &&
      !(*hist=SDDS_Realloc(*hist, sizeof(**hist)*(*maxBins=*bins))))
    bombElegant("Memory allocation failure (binParticleCoordinate)", NULL);
    
  for (iBin=0; iBin<*bins; iBin++)
    (*hist)[iBin] = 0;
  nBinned = 0;
  for (iParticle=nBinned=0; iParticle<nParticles; iParticle++) {
    /* the coordinate of the bin center is (iBin+0.5)*(*binSize) + *lower */
    iBin = (particleCoord[iParticle][coordinateIndex] - *lower)/(*binSize);
    if (iBin<0 || iBin>(*bins-1))
        continue;
    (*hist)[iBin] += 1;
    nBinned++;
  }
  return nBinned;
}
#endif

void computeSaldinFdNorm(double **FdNorm, double **x, long *n, double sMax, long ns,
                         double Po, double radius, double angle, double dx, char *normMode);
long track_through_driftCSR_Stupakov(double **part, long np, CSRDRIFT *csrDrift, 
                                     double Po, double **accepted, double zStart, char *rootname);

long track_through_driftCSR(double **part, long np, CSRDRIFT *csrDrift, 
                            double Po, double **accepted, double zStart, 
			    double revolutionLength, char *rootname)
{
  long iPart, iKick, iBin, binned=0, nKicks, iSpreadMode=0;
  double *coord, p, beta, dz, ct0=0.0, factor, dz0, dzFirst;
  double ctmin, ctmax, spreadFactor, dct;
  double zTravel, attenuationLength, thetaRad=0.0, sigmaZ, overtakingLength, criticalWavelength, wavelength=0.0;
  static char *spreadMode[3] = {"full", "simple", "radiation-only"};
  static char *wavelengthMode[3] = {"sigmaz", "bunchlength", "peak-to-peak"};
  static char *bunchlengthMode[3] = {"rms", "68-percentile", "90-percentile"};
  unsigned long mode;
  static long warned = 0, incrementWarningsLeft=100;
  long nBins1;
  TRACKING_CONTEXT tContext;
#if USE_MPI 
  long np_total=1, np_tmp=np, binned_total;
#endif
  
  getTrackingContext(&tContext);

#if (!USE_MPI)
  if (np<=1 || !csrWake.valid || !csrDrift->csr) {
#else
  if (notSinglePart){
    if (isMaster) 
      np_tmp = 0;  
    MPI_Allreduce(&np_tmp, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);   
  } else
    np_total = np;

  if (np_total<=1 || !csrWake.valid || !csrDrift->csr) 	{
    if (isSlave||!notSinglePart) {
#endif
      if (csrDrift->linearOptics) {
	long i;
	for (i=0; i<np; i++) {
	  part[i][0] += csrDrift->length*part[i][1];
	  part[i][2] += csrDrift->length*part[i][3];
	  part[i][4] += csrDrift->length;
	}
      }
      else
	exactDrift(part, np, csrDrift->length);
#if (USE_MPI)
    }
#endif
    return np;
  }	
  nBins1 = csrWake.bins - 1;

  mode = 
    (csrDrift->spread?CSRDRIFT_SPREAD:0) +
      (csrDrift->useOvertakingLength?CSRDRIFT_OVERTAKINGLENGTH:0) +
        (csrDrift->useSaldin54?CSRDRIFT_SALDIN54:0) +
          (csrDrift->attenuationLength>0?CSRDRIFT_ATTENUATIONLENGTH:0) +
            (csrDrift->useStupakov?CSRDRIFT_STUPAKOV:0) ;
  while (zStart<csrWake.zLast) {
    if (incrementWarningsLeft) {
      fprintf(stdout, "*** Warning: incrementing zStart by revolution length for CSRDRIFT (%s #%ld).\n",
              tContext.elementName, tContext.elementOccurrence);
      fprintf(stdout, "    If you are not simulating a ring, this could be a problem!\n");
      incrementWarningsLeft --;
    }
    zStart += revolutionLength;
  }
  if (bitsSet(mode)>1) {
    fprintf(stdout, "Error: Too many modes set for CSRDRIFT.\n");
    exitElegant(1);
  }
  if (csrWake.lastMode && csrWake.lastMode!=mode) {
    fprintf(stdout, "Error: CSRDRIFT mode changed between dipoles. Pick one mode following each dipole.\n");
    exitElegant(1);
  }
  csrWake.lastMode = mode;
  
  if (mode&CSRDRIFT_STUPAKOV)
    return track_through_driftCSR_Stupakov(part, np, csrDrift, Po, accepted, zStart, rootname);

  if (!warned) {
    fprintf(stdout, "Warning: USE_STUPAKOV=1 is recommended for CSRDRIFT elements.\n");
    fprintf(stdout, "This is the most physical model available at this time in elegant.\n");
    warned = 1;
  }
  
  dct = csrWake.dctBin;
  if (csrDrift->dz>0) {
    if ((nKicks = csrDrift->length/csrDrift->dz)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bombElegant("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=csrDrift->length/nKicks)/2;
  
  sigmaZ = 0;
  switch (match_string(csrDrift->bunchlengthMode, bunchlengthMode, 3, 0)) {
  case 0:
    sigmaZ = csrWake.rmsBunchLength;
    break;
  case 1:
    sigmaZ = csrWake.perc68BunchLength;
    break;
  case 2:
    sigmaZ = csrWake.perc90BunchLength;
    break;
  default:
    bombElegant("invalid bunchlength_mode for CSRDRIFT.  Use rms or percentile.", NULL);
  }
  
  overtakingLength = pow(24*sigmaZ*csrWake.rho*csrWake.rho, 1./3.);

  if (mode&CSRDRIFT_OVERTAKINGLENGTH)
    attenuationLength = overtakingLength*csrDrift->overtakingLengthMultiplier;
  else
    attenuationLength = csrDrift->attenuationLength;
  
  if (mode&CSRDRIFT_SPREAD) {
    iSpreadMode = 0;
    if (csrDrift->spreadMode && 
        (iSpreadMode=match_string(csrDrift->spreadMode, spreadMode, 3, 0))<0)
      bombElegant("invalid spread_mode for CSR DRIFT.  Use full, simple, or radiation-only", NULL);
    switch (match_string(csrDrift->wavelengthMode, wavelengthMode, 3, 0)) {
    case 0:
    case 1:
      /* bunch length */
      wavelength = sigmaZ;
      break;
    case 2:
      /* peak-to-peak */
      wavelength = csrWake.peakToPeakWavelength;
      break;
    default:
      bombElegant("invalid wavelength_mode for CSR DRIFT.  Use sigmaz or peak-to-peak", NULL);
      break;
    }
    criticalWavelength = 4.19/ipow(csrWake.Po, 3)*csrWake.rho;
    if (!particleIsElectron)
      bombElegant("CSRDRIFT spread mode is not supported for particles other than electrons", NULL);
    thetaRad = 0.5463e-3/(csrWake.Po*0.511e-3)/pow(criticalWavelength/wavelength, 1./3.);
  }

  if (mode&CSRDRIFT_SALDIN54) {
    if (csrWake.FdNorm==NULL) {
      if (csrDrift->nSaldin54Points<20) 
        csrDrift->nSaldin54Points = 20;
      computeSaldinFdNorm(&csrWake.FdNorm, &csrWake.xSaldin, &csrWake.nSaldin,
                          2*sigmaZ, csrDrift->nSaldin54Points, csrWake.Po, csrWake.rho, csrWake.bendingAngle, dz,
                          csrDrift->normMode);
      if (csrDrift->Saldin54Output)  {
        long ix;
        if (!csrDrift->fpSaldin) {
          csrDrift->Saldin54Output = compose_filename(csrDrift->Saldin54Output, rootname);
          csrDrift->fpSaldin = fopen(csrDrift->Saldin54Output, "w");
          fprintf(csrDrift->fpSaldin, "SDDS1\n&column name=z, type=double &end\n&column name=Factor, type=double &end\n");
          fprintf(csrDrift->fpSaldin, "&data mode=ascii no_row_counts=1 &end\n");
        } else
          fprintf(csrDrift->fpSaldin, "\n");
        for (ix=0; ix<csrWake.nSaldin; ix++) 
          fprintf(csrDrift->fpSaldin, "%le %le\n", csrWake.xSaldin[ix], csrWake.FdNorm[ix]);
        fflush(csrDrift->fpSaldin);
      }
    }
  }

  dzFirst = zStart - csrWake.zLast;
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */
#ifdef DEBUG
  fprintf(stdout, "CSR in drift:\n");
  fprintf(stdout, "zStart = %21.15le, zLast = %21.15le, zTravel = %21.15le\n", zStart, csrWake.zLast,
          zTravel);
  fprintf(stdout, "dzFirst = %21.15e, s0 = %21.15e\n", dzFirst, csrWake.s0);
#endif

  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;

    ctmin = DBL_MAX;
    ctmax = -DBL_MAX;

    /* propagate particles forward, converting s to c*t=s/beta */
    if (isSlave || !notSinglePart) {
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      p = Po*(1+coord[5]);
      beta = p/sqrt(p*p+1);
      if (csrDrift->linearOptics) 
        coord[4] = (coord[4]+dz)/beta;
      else 
        coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
#ifdef DEBUG
      if (coord[4]>ctmax)
        ctmax = coord[4];
      if (coord[4]<ctmin)
        ctmin = coord[4];
#endif
    }
    }

    factor = 1;
    if (csrWake.dGamma) {
      /* propagate wake forward */
      csrWake.s0 += dz+dzFirst;   /* accumulates position of back end of the radiation pulse */
      ct0 = csrWake.s0;
      
      if (attenuationLength>0) {
        /* attenuate wake */
        if ((factor = exp(-(dz+dzFirst)/attenuationLength))<1) {
          for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
      }
      /* factor to account for difference in drift lengths here and in
       * csrcsbend integration.  Use dz0 here because that is the
       * length integrated by each kick.  Add dzFirst to account for any
       * length we may have missed due to intervening non-drift elements.
       */
      factor = (dz0+dzFirst)/csrWake.ds0;
    }
    if (mode&CSRDRIFT_SPREAD) {
      /* compute loss of on-axis field due to spread of beam using a simple-minded
       * computation of beam sizes */
      switch (iSpreadMode) {
      case 0:  /* full */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + 
                                     2*zTravel*csrWake.S12 + 
                                     zTravel*zTravel*(sqr(thetaRad)+csrWake.S22))));
        break;
      case 1: /* simple */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + zTravel*zTravel*(sqr(thetaRad)+csrWake.S22))));
        break;
      case 2: /* radiation only */
        factor *= (spreadFactor =
                   sqrt(csrWake.S11/(csrWake.S11 + sqr(zTravel*thetaRad))));
        break;
      default:
        bombElegant("invalid spread code---programming error!", NULL);
        break;
      }
    }
    
    if (mode&CSRDRIFT_SALDIN54) {
      long code;
      double f0 = 0;
      if (zTravel<=csrWake.xSaldin[csrWake.nSaldin-1]) 
        factor *= (f0=interp(csrWake.FdNorm, csrWake.xSaldin, csrWake.nSaldin, zTravel, 0, 1, &code));
      else 
        factor = 0;
      csrWake.lastFdNorm = f0;
#ifdef DEBUG
      fprintf(csrWake.fpSaldin, "%le %le\n", zTravel, f0);
      fflush(csrWake.fpSaldin);
#endif
      if (!code) {
        fprintf(stderr, "Warning: interpolation failure for Saldin eq. 54\n");
        fprintf(stderr, "zTravel = %le,  csrWake available up to %le\n",
                zTravel, csrWake.xSaldin[csrWake.nSaldin-1]);
        factor = 0;
      }
    }
    
    dzFirst = 0;

    /* apply kick to each particle and convert back to normal coordinates */
    if (isSlave || !notSinglePart) {
    for (iPart=binned=0; iPart<np; iPart++) {
      coord = part[iPart];
      if (csrWake.dGamma) {
        double f;
        iBin = (f=(coord[4]-ct0)/dct);
        f -= iBin;
        if (iBin>=0 && iBin<nBins1) {
          coord[5] += ((1-f)*csrWake.dGamma[iBin]+f*csrWake.dGamma[iBin+1])/Po*factor;
          binned ++;
        }
      }
      p = (1+coord[5])*Po;
      beta = p/sqrt(p*p+1);
      coord[4] = beta*coord[4];
    }
    }
#if USE_MPI
    if (isSlave && notSinglePart) {
      MPI_Allreduce(&binned, &binned_total, 1, MPI_LONG, MPI_SUM, workers);
    }
    if ((myid==1) && (csrWake.dGamma && np_total!=binned_total)) {
      dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift %s (track_through_driftCSR)\n",
              binned_total, np_total, tContext.elementName);
#else
    if (csrWake.dGamma && np!=binned) {
      fprintf(stdout, "only %ld of %ld particles binned for CSR drift %s (track_through_driftCSR)\n",
              binned, np, tContext.elementName);
#endif
      fprintf(stdout, "beam ct min, max = %21.15e, %21.15e\n",
              ctmin, ctmax);
      fprintf(stdout, "wake ct0 = %21.15e, ct1 = %21.15e\n",
              ct0, ct0+csrWake.dctBin*csrWake.bins);
      fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif  
    }
  }
  /* do final drift of dz0/2 */
  dz = dz0/2;
  if (isSlave || !notSinglePart) {
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      if (csrDrift->linearOptics)
	coord[4] += dz;
      else
	coord[4] += dz*sqrt(1+sqr(coord[1])+sqr(coord[3]));
    }    
  }
  csrWake.zLast = zStart+csrDrift->length;
  
  if (csrWake.dGamma) {
    /* propagate wake forward */
    csrWake.s0 += dz;
    ct0 = csrWake.s0;
    
    if (attenuationLength>0) {
      /* attenuate wake */
      if ((factor = exp(-dz/attenuationLength))<1) {
        for (iBin=0; iBin<csrWake.bins; iBin++)
            csrWake.dGamma[iBin] *= factor;
        }
    }
  }

  return np;
}

/* this should be called before starting to track a beamline to make sure that
 * CSR drift elements upstream of all CSRBEND elements get treated like ordinary
 * drift spaces. */

long reset_driftCSR()
{
  csrWake.lastMode = 0;
  if (csrWake.valid && csrWake.FdNorm) {
    fprintf(stdout, "Last value of normalization factor for CSR wake was %le\n",
            csrWake.lastFdNorm);
  }
  csrWake.valid = csrWake.bins = 0;
  csrWake.dctBin = csrWake.s0 = csrWake.ds0 = csrWake.zLast =
    csrWake.z0 = csrWake.S11 = csrWake.S12 = csrWake.S22 = 0;
  csrWake.dGamma = NULL;
  csrWake.nSaldin = 0;
  if (csrWake.FdNorm) {
    free(csrWake.FdNorm);
    free(csrWake.xSaldin);
    csrWake.FdNorm = csrWake.xSaldin = NULL;
  }
  if (csrWake.StupakovFileActive) {
    if (!SDDS_Terminate(&csrWake.SDDS_Stupakov))
      bombElegant("problem terminating data file for Stupakov output from CSRDRIFT", NULL);
    csrWake.StupakovFileActive = 0;
  }
  return 1;
}

double SolveForPsiSaldin54(double xh, double sh);
double Saldin5354Factor(double xh, double sh, double phihm, double xhLowerLimit);

void computeSaldinFdNorm(double **FdNorm, double **x, long *n, double sMax, long ns,
                         double Po, double radius, double bendingAngle, double dx, 
                         char *normMode)
{
  double xEnd, sh, beta, gamma, xh, dx0;
  long ix, is;
  double phihs, phihm, xhLowerLimit, xUpperLimit, s, f, fx;
  double t1, t2, f0, fmax;
  char *allowedNormMode[2] = {"first", "peak"};

  gamma = sqrt(sqr(Po)+1);
  beta = Po/gamma;

  if ((xEnd = sMax/(1-beta))>1000 || isnan(xEnd) || isinf(xEnd)) {
    fprintf(stderr, "Warning: the extent of the CSR drift wake decay was limited at 1km\n");
    xEnd = 1000;
  }

  *n = 100;
  dx0 = xEnd/(100*(*n));
  if (dx<dx0) {
    *n = xEnd/(100*dx);
    if (*n>100000) {
      *n = 100000;
      fprintf(stderr, "Note: the CSR drift wake decay table size hit the limit of 100k points\n");
    }
  } else 
    dx = dx0;
  fx = pow(xEnd/dx, 1./(*n));

  if (!(*FdNorm = calloc(sizeof(**FdNorm), (*n))) ||
      !(*x = malloc(sizeof(**x)*(*n))))
    bombElegant("memory allocation failure (computeSaldinFdNorm)", NULL);

  for (ix=0; ix<*n; ix++)
    (*x)[ix] = ix==0 ? 0 : ipow(fx, ix-1)*dx;
  for (is=0; is<ns; is++) {
    /* don't use s=0 as it is singular */
    s = (is+1.0)*sMax/ns;
    sh = s*ipow(gamma, 3)/radius;
    phihm = bendingAngle*gamma;
    t1 = 12*sh;
    t2 = sqrt(64+144*sh*sh);
    phihs = pow(t1+t2, 1./3.) - pow(-t1 + t2, 1./3.);
    xhLowerLimit = -1;
    if (phihs>phihm)
      xhLowerLimit = sh - phihm - ipow(phihm, 3)/6 + sqrt(sqr(ipow(phihm, 3)-6*sh) + 9*ipow(phihm, 4))/6;
    xUpperLimit = 0.999*s/(1-beta);
    for (ix=0; ix<*n; ix++) {
      if ((*x)[ix]>=xUpperLimit)
        break;
      xh = (*x)[ix]*gamma/radius;
      (*FdNorm)[ix] += Saldin5354Factor(xh, sh, phihm, xhLowerLimit);
    }
  }

  /* average over s */
  for (ix=0; ix<*n; ix++)
    (*FdNorm)[ix] /= ns;
  
  /* get the first nonzero and also the maximum value of Fd */
  for (ix=f0=fmax=0; ix<*n; ix++) {
    f= (*FdNorm)[ix];
    if (f0==0 && f>0)
      f0 = f;
    if (fmax<f)
      fmax = f;
  }
  if (fmax>f0/0.99) {
    fprintf(stderr, "Warning: possible problem with SALDIN54 drift mode: too few (%ld) points. Max/start-1 is %le\n",
            ns,
            fmax/f0-1);
  }
  switch (match_string(normMode, allowedNormMode, 2, 0)) {
  case 0:
    /* first */
    f = f0;
    break;
  case 1:
    /* peak */
    f = fmax;
    break;
  default:
    fprintf(stderr, "Error: unknown Saldin-54 normalization mode: %s\n", normMode);
    exitElegant(1);
    break;
  }
  if (f)
    for (ix=0; ix<*n; ix++)
      (*FdNorm)[ix] /= f;
  else
    for (ix=0; ix<*n; ix++)
      (*FdNorm)[ix] = 0;
}

double SolveForPsiSaldin54(double xh, double sh)
{
  double s_sum, s_diff2, bestSol;
  double solList[4] = {-1, -1, -1, -1};
  long nSols=0, sol;

  s_sum = (-2*xh - sqrt(-8 + 4*pow(xh,2) - 
                         (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                               3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                    pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                          3*pow(xh,4),2)),0.3333333333333333) + 
                         2*pow(2,0.6666666666666666)*
                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                               3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                    pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                          3*pow(xh,4),2)),0.3333333333333333)))/2.;
  if (!isnan(s_sum)) {
    s_diff2 = (-16 + 8*pow(xh,2) + (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              2*pow(2,0.6666666666666666)*
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) + 
              (16*(-3*sh + pow(xh,3)))/
              sqrt(-8 + 4*pow(xh,2) - 
                   (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
            3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                 pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                       3*pow(xh,4),2)),0.3333333333333333) + 
                   2*pow(2,0.6666666666666666)*
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333)))/4.;
    if (s_diff2>=0) {
      solList[0] = s_sum+sqrt(s_diff2);
      solList[1] = s_sum+sqrt(s_diff2);
      nSols = 2;
    }
  }
  
  s_sum =    (-2*xh + sqrt(-8 + 4*pow(xh,2) - 
                            (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                            pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                  3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                       pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                             3*pow(xh,4),2)),0.3333333333333333) + 
                            2*pow(2,0.6666666666666666)*
                            pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                  3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                                       pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                             3*pow(xh,4),2)),0.3333333333333333)))/2.;
  if (!isnan(s_sum)) {
    s_diff2 = (-16 + 8*pow(xh,2) + (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              2*pow(2,0.6666666666666666)*
              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                    3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                         pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                               3*pow(xh,4),2)),0.3333333333333333) - 
              (16*(-3*sh + pow(xh,3)))/
              sqrt(-8 + 4*pow(xh,2) - 
                   (4*pow(2,0.3333333333333333)*(-1 + pow(xh,2)))/
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333) + 
                   2*pow(2,0.6666666666666666)*
                   pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                         3*pow(xh,4) + sqrt(4*pow(-1 + pow(xh,2),3) + 
                                              pow(2 + 9*pow(sh,2) - 3*pow(xh,2) - 6*sh*pow(xh,3) + 
                                                    3*pow(xh,4),2)),0.3333333333333333)))/4.;
  
    if (s_diff2>=0) {
      solList[nSols] = s_sum+sqrt(s_diff2);
      solList[nSols+1] = s_sum-sqrt(s_diff2);
      nSols += 2;
    }
  }
  bestSol = solList[0];
  for (sol=0; sol<nSols; sol++) {
    if (solList[sol]>bestSol) {
      bestSol = solList[sol];
    }     
  }
  return bestSol;
}

double Saldin5354Factor(double xh, double sh, double phihm, double xhLowerLimit)
{
  double t1, t2, f, psi, psi2;
  if (xh<xhLowerLimit) {
    /* use Saldin 53 */
    t1 = (ipow(phihm, 3) + 3*xh*sqr(phihm) - 6*sh);
    t2 = 3*(phihm+2*xh);
    f = 2/(phihm+2*xh)*(1 + (t1 + t2)/sqrt(t1*t1+sqr(phihm*t2))) - 1/sh;
  } else {
    if ((psi = SolveForPsiSaldin54(xh, sh))>=0) {
      psi2 = psi*psi;
      f =  4*(2*xh*(psi2+1)+psi*(psi2+2))/
        (4*xh*xh*(psi2+1)+4*xh*psi*(psi2+2)+psi2*(psi2+4)) - 1/sh;
    } else
      return 0;
  }
  if (isnan(f) || isinf(f))
    f = 0;
  return f;
}

void exactDrift(double **part, long np, double length)
{
  long i;
  double *coord;
  for (i=0; i<np; i++) {
    coord = part[i];
    coord[0] += coord[1]*length;
    coord[2] += coord[3]*length;
    coord[4] += length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
  }
}


double SolveForPhiStupakov(double x, double ds, double phim);
void DumpStupakovOutput(char *filename, SDDS_DATASET *SDDSout, long *active,
                        double zTravel, double *ctHist, double *ctHistDeriv,
                        double *dGamma, long nBins, double dct, 
                        double MPCharge, double dz,
                        long nCaseC, long nCaseD1,long nCaseD2,
                        double x, double dsMax, double phi0, double phi1) ;

static double SolveForPhiStupakovDiffSum = 0;
static long SolveForPhiStupakovDiffCount = 0;

long track_through_driftCSR_Stupakov(double **part, long np, CSRDRIFT *csrDrift, 
                            double Po, double **accepted, double zStart, char *rootname)
{
  long iPart, iKick, iBin, binned=0, nKicks;
  long nCaseC, nCaseD1, nCaseD2;
  double ctLower, ctUpper, ds;
  long nBins, maxBins, nBinned, diBin;
  double *coord, p, beta, dz, factor, dz0, dzFirst;
  double zTravel, dct, zOutput;
  double *ctHist=NULL, *ctHistDeriv=NULL, *phiSoln=NULL;
  double length;
  long nBins1;
  double dsMax, x;
  TRACKING_CONTEXT tContext;
#if USE_MPI
  long binned_total=1, np_total=1;
  double *buffer;
#endif

  getTrackingContext(&tContext);

  SolveForPhiStupakovDiffCount = 0;
  SolveForPhiStupakovDiffSum = 0;
  
  length = csrDrift->length;
  if (zStart!=csrWake.zLast) {
    length += (dzFirst = zStart-csrWake.zLast);
    /* propagate beam back so we can tranverse the missing length including CSR
     */
    if (isSlave || !notSinglePart)
      for (iPart=0; iPart<np; iPart++) {
	coord = part[iPart];
	coord[0] -= dzFirst*coord[1];
	coord[2] -= dzFirst*coord[3];
	if (csrDrift->linearOptics)
	  coord[4] -= dzFirst;
	else
	  coord[4] -= dzFirst*sqrt(1+sqr(coord[1])+sqr(coord[3]));
      }
    zStart = csrWake.zLast;
  }
  zOutput = zStart;  /* absolute coordinate used for output of data vs z or s */
  
  if (csrDrift->dz>0) {
    if ((nKicks = length/csrDrift->dz+0.5)<1)
      nKicks = 1;
  } else 
    nKicks = csrDrift->nKicks;
  if (nKicks<=0)
    bombElegant("nKicks=0 in CSR drift.", NULL);
  dz = (dz0=length/nKicks)/2;
  
  zTravel = zStart-csrWake.z0;  /* total distance traveled by radiation to reach this point */

  maxBins = nBins = csrWake.bins;
  nBins1 = nBins-1;
  if (!(ctHist=SDDS_Malloc(sizeof(*ctHist)*nBins)) ||
      !(ctHistDeriv=SDDS_Malloc(sizeof(*ctHistDeriv)*nBins)) ||
      !(phiSoln=SDDS_Malloc(sizeof(*phiSoln)*nBins)))
    bombElegant("memory allocation failure (track_through_driftCSR)", NULL);
  
  for (iKick=0; iKick<nKicks; iKick++) {
    /* first drift is dz=dz0/2, others are dz0 */
    if (iKick==1)
      dz = dz0;
    zTravel += dz;
    zOutput += dz;
    
    x = zTravel/csrWake.rho;
    dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
      *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
    /* propagate particles forward, converting s to c*t=s/beta */
    if (isSlave || !notSinglePart) {
      for (iPart=0; iPart<np; iPart++) {
	coord = part[iPart];
	coord[0] += coord[1]*dz;
	coord[2] += coord[3]*dz;
	p = Po*(1+coord[5]);
	beta = p/sqrt(p*p+1);
	if (csrDrift->linearOptics)
	  coord[4] = (coord[4]+dz)/beta;
	else
	  coord[4] = (coord[4]+dz*sqrt(1+sqr(coord[1])+sqr(coord[3])))/beta;
      }
    }
    /* bin the particle distribution */
    ctLower = ctUpper = dct = 0;
    nBinned = binParticleCoordinate(&ctHist, &maxBins,
				    &ctLower, &ctUpper, &dct, &nBins, 
				    csrWake.binRangeFactor<1.1?1.1:csrWake.binRangeFactor,
				    part, np, 4);
#if USE_MPI
  if (notSinglePart) {
    if (isSlave)
      MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);    
    MPI_Allreduce(&nBinned, &binned_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  }

  if (notSinglePart) {  /* Master needs to know the information to write the result */
    buffer = malloc(sizeof(double) * nBins);
    MPI_Allreduce(ctHist, buffer, nBins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(ctHist, buffer, sizeof(double)*nBins);
    free(buffer);
  }
  if ((myid==1) && (np_total!=binned_total)) {
    dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
    fprintf(stdout, "Only %ld of %ld particles binned for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	    binned_total, np_total,
	    tContext.elementName, csrWake.binRangeFactor);
    fflush(stdout);
#else
  if (nBinned!=np) {
    fprintf(stdout, "Only %ld of %ld particles binned for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	    nBinned, np,
	    tContext.elementName, csrWake.binRangeFactor);
#endif
    fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	    ctLower, ctUpper, dct, nBins, maxBins);
    fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif 
  }
      
    /* - smooth the histogram, normalize to get linear density, and 
       copy in preparation for taking derivative
       */
    if (csrWake.highFrequencyCutoff0>0 || csrWake.lowFrequencyCutoff0>=0) {
      long nz;
      nz = applyLHPassFilters(ctHist, nBins, 
                              csrWake.lowFrequencyCutoff0, csrWake.lowFrequencyCutoff1,
                              csrWake.highFrequencyCutoff0, csrWake.highFrequencyCutoff1,
                              csrWake.clipNegativeBins);
      if (nz && negativeWarningsLeft) {
        fprintf(stdout, "Warning: low pass filter resulted in negative values in %ld bins\n",
                nz);
        if (--negativeWarningsLeft==0)
          fprintf(stdout, "         Further warnings will be suppressed for this run.\n");
        fflush(stdout);
      }
    }

    if (csrWake.SGHalfWidth>0) {
      SavitzkyGolaySmooth(ctHist, nBins, csrWake.SGOrder, csrWake.SGHalfWidth, csrWake.SGHalfWidth,  0);
#if (!USE_MPI)
      correctDistribution(ctHist, nBins, 1.0*nBinned);
#else
      if (notSinglePart)
	correctDistribution(ctHist, nBins, 1.0*binned_total);
      else
	correctDistribution(ctHist, nBins, 1.0*nBinned);
#endif
    }
    for (iBin=0; iBin<nBins; iBin++)
      ctHistDeriv[iBin] = (ctHist[iBin] /= dct);
    /* - compute derivative with smoothing.  The deriv is w.r.t. index number and
     * I won't scale it now as it will just fall out in the integral 
     */
    SavitzkyGolaySmooth(ctHistDeriv, nBins, csrWake.SGDerivOrder, 
                       csrWake.SGDerivHalfWidth, csrWake.SGDerivHalfWidth, 1);

    /* Case C */ 
    nCaseC = 0;
    nCaseD1 = 0;
    nCaseD2 = 0;
    for (iBin=0; iBin<nBins; iBin++) {
      double f;
      ds = csrWake.rho/6*sqr(csrWake.bendingAngle)*(csrWake.bendingAngle + 3*x);
      diBin = ds/dct;
      if (iBin+diBin<nBins) {
        f = -1/(csrWake.bendingAngle+2*x); 
        csrWake.dGamma[iBin] = f*ctHist[iBin+diBin];
        nCaseC++;
      } else
        csrWake.dGamma[iBin] = 0;
    }
    /* Case D */
    for (iBin=0; iBin<nBins; iBin++) {
      phiSoln[iBin] = -1;
      if ((ds = iBin*dct)>dsMax)
        break;
      phiSoln[iBin] = SolveForPhiStupakov(x, iBin*dct/csrWake.rho, csrWake.bendingAngle);
    }
    for (iBin=0; iBin<nBins; iBin++) {
      long jBin, first, count;
      double term1=0, term2=0;
      diBin = dsMax/dct;
      if (iBin+diBin<nBins) {
        nCaseD1 ++;
        csrWake.dGamma[iBin] += ctHist[iBin+diBin]/(csrWake.bendingAngle+2*x);
      }
      first = 1;
      count = 0;
      for (jBin=iBin; jBin<nBins; jBin++) {
        double phi;
        if ((phi = phiSoln[jBin-iBin])>=0) {
          /* I put in a negative sign here because my s is opposite in direction to 
           * Saldin et al. and Stupakov, so my derivative has the opposite sign.
           * Note lack of ds factor here as I use the same one in my unnormalized derivative.
           */
          if (phi>0) {
            /* ^^^ If I test phi+2*x here, I get noisy, unphysical results very close
             * to the dipole exit 
             */
            term2 = ctHistDeriv[jBin]/(phi+2*x);
            csrWake.dGamma[iBin] -= term2;
            if (first) {
              term1 = term2;
              first = 0;
            }
            count++;
            nCaseD2++;
          }
        } else
          break;
      }
      if (count>1 && csrWake.trapazoidIntegration)
        /* trapazoid rule correction for ends */
        csrWake.dGamma[iBin] += (term1+term2)/2;
    }
    /* the minus sign adjusts for Stupakov using wake<0 to indicate energy gain
     */
    factor = -4/csrWake.rho*csrWake.GSConstant*dz0;
    for (iBin=0; iBin<nBins; iBin++)
      csrWake.dGamma[iBin] *= factor;

    if (csrWake.wffValues) 
      applyFilterTable(csrWake.dGamma, nBins, dct/c_mks, csrWake.wffValues, csrWake.wffFreqValue,
                       csrWake.wffRealFactor, csrWake.wffImagFactor);

    if ((csrDrift->StupakovOutput || csrWake.StupakovFileActive) && 
        (csrDrift->StupakovOutputInterval<2 || iKick%csrDrift->StupakovOutputInterval==0)) {
      double x, dsMax, phi0, phi1;
      if (!csrWake.StupakovFileActive) {
        if (!SDDS_CopyString(&csrWake.StupakovOutput, csrDrift->StupakovOutput))
          bombElegant("string copying problem preparing Stupakov output for CSRDRIFT", NULL);
        csrWake.StupakovOutput = compose_filename(csrWake.StupakovOutput, rootname);
      }
      x = zTravel/csrWake.rho;
      dsMax = csrWake.rho/24*pow(csrWake.bendingAngle, 3)
        *(csrWake.bendingAngle+4*x)/(csrWake.bendingAngle+x);
      phi0 = SolveForPhiStupakov(x, 0.0, csrWake.bendingAngle);
      phi1 = SolveForPhiStupakov(x, dsMax/csrWake.rho*0.999, csrWake.bendingAngle);
      
      /* note that the contents of ctHist and ctHistDeriv are corrupted by this operation */
      DumpStupakovOutput(csrWake.StupakovOutput, &csrWake.SDDS_Stupakov, 
                         &csrWake.StupakovFileActive, zTravel,
                         ctHist, ctHistDeriv, csrWake.dGamma, nBins, dct, csrWake.MPCharge,
                         dz0, nCaseC, nCaseD1, nCaseD2,
                         x, dsMax/csrWake.rho, phi0, phi1);
    }
    
    /* apply kick to each particle and convert back to normal coordinates */
    if (isSlave || !notSinglePart) {
      for (iPart=binned=0; iPart<np; iPart++) {
	double f;
	coord = part[iPart];
	iBin = (f=(coord[4]-ctLower)/dct);
	f -= iBin;
	if (iBin>=0 && iBin<nBins1) {
	  coord[5] += ((1-f)*csrWake.dGamma[iBin] + f*csrWake.dGamma[iBin+1])/Po;
	  binned ++;
	} else {
	  fprintf(stdout, "Particle out of bin range---not kicked: ct-ctLower=%21.15e, dct=%21.15e, iBin=%ld\n",
		  coord[4]-ctLower, dct, iBin);
	}
	p = (1+coord[5])*Po;
	beta = p/sqrt(p*p+1);
	coord[4] = beta*coord[4];
      }
    }

    if (tContext.sliceAnalysis && tContext.sliceAnalysis->active &&
	(csrDrift->sliceAnalysisInterval==0 ||
	 iKick%csrDrift->sliceAnalysisInterval==0)) {
#if USE_MPI
      /* This function will be parallelized in the future */
      fprintf(stdout, "performSliceAnalysisOutput is not supported in parallel mode currently.\n");
      MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
      MPI_Abort(MPI_COMM_WORLD, 1); 
#endif
	performSliceAnalysisOutput(tContext.sliceAnalysis, part, np, 
				   0, tContext.step, Po, 
				   csrWake.MPCharge*np,
				   tContext.elementName, 
				   zOutput, 0);
    }
#if USE_MPI
    if (isSlave && notSinglePart) {
      MPI_Allreduce(&binned, &binned_total, 1, MPI_LONG, MPI_SUM, workers);
    }
    if ((myid==1) && (np_total!=binned_total)) {
      dup2(fd,fileno(stdout)); /* Let the first slave processor write the output */
      fprintf(stdout, "Only %ld of %ld particles kicked for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	      binned_total, np_total,
	      tContext.elementName, csrWake.binRangeFactor);
#else
    if (np!=binned) {
      fprintf(stdout, "Only %ld of %ld particles kicked for CSRDRIFT (%s, BRF=%le, Stupakov)\n", 
	      binned, np,
	      tContext.elementName, csrWake.binRangeFactor);
#endif
      fprintf(stdout, "ct min, max = %21.15e, %21.15e, dct = %21.15e, nBins=%ld, maxBins=%ld\n",
	      ctLower, ctUpper, dct, nBins, maxBins);
      fflush(stdout);
#if USE_MPI
#if defined(_WIN32)
    freopen("NUL","w",stdout); 
#else
    freopen("/dev/null","w",stdout); 
#endif
#endif 
    }
  }
  
  /* do final drift of dz0/2 */
  dz = dz0/2;
  if (isSlave || !notSinglePart)
    for (iPart=0; iPart<np; iPart++) {
      coord = part[iPart];
      coord[0] += coord[1]*dz;
      coord[2] += coord[3]*dz;
      if (csrDrift->linearOptics)
	coord[4] += dz;
      else
	coord[4] += dz*sqrt(1+sqr(coord[1])+sqr(coord[3]));
    }    

  csrWake.zLast = zStart + length;
  free(ctHist);
  free(ctHistDeriv);
  free(phiSoln);
#if DEBUG
  if (SolveForPhiStupakovDiffCount)
    fprintf(stdout, "Phi solution accuracy for %ld solutions: %le\n",
            SolveForPhiStupakovDiffCount, SolveForPhiStupakovDiffSum/SolveForPhiStupakovDiffCount);
#endif
  return np;
}

static double SolveForPhiStupakov_x, SolveForPhiStupakov_4x;

double SolveForPhiStupakovFn(double phi)
{
  return phi*phi*phi*(phi+SolveForPhiStupakov_4x)/(phi+SolveForPhiStupakov_x);
}

/* solve for phi:  ds=phi^3/24*(phi+4*x)/(phi+x), where ds = (s-s')/rho */

double SolveForPhiStupakov(double x, double ds, double phim)
{
  double phi;
  static double phiLast = -1;
  
  if (ds<0)
    return -1;
  if (ds==0)
    return 0;
  
  ds *= 24;
  SolveForPhiStupakov_x = x;
  SolveForPhiStupakov_4x = 4*x;

  if (phiLast==-1)
    phiLast = phim/2;

  /* try phim first */
  if (fabs(ds-SolveForPhiStupakovFn(phim))<ds/1e4) {
    phiLast = phim;
    return phim;
  }
  
  /* try a solution with Newton's method */
  phi = zeroNewton(SolveForPhiStupakovFn, ds, phiLast, phim/1000, 3, ds/1e4);
  if (phi<0 || phi>phim || fabs(ds - SolveForPhiStupakovFn(phi))>ds/1e4) 
    /* try a more plodding method */
    phi = zeroInterp(SolveForPhiStupakovFn, ds, 0, phim*1.01, phim/100, ds/1e4);
  if (phi<0 || phi>phim)
    return -1;
  phiLast = phi;
  SolveForPhiStupakovDiffCount ++;
  SolveForPhiStupakovDiffSum += fabs(ds - SolveForPhiStupakovFn(phi));
  return phi;
}


/* this procedure destroys the contents of ctHist and ctHistDeriv ! */

void DumpStupakovOutput(char *filename, SDDS_DATASET *SDDSout, long *active,
                        double zTravel, double *ctHist, double *ctHistDeriv,
                        double *dGamma, long nBins, double dct, 
                        double MPCharge, double dz,
                        long nCaseC, long nCaseD1, long nCaseD2,
                        double x, double dsMax, double phi0, double phi1) 
{
  long i;
  if (!*active) {
    if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, filename) ||
        !SDDS_DefineSimpleParameter(SDDSout, "z", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseC", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseD1", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "CaseD2", "#", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDSout, "x", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "dsMax", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "phi0", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(SDDSout, "phi1", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "s", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "LinearDensity", "C/s", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "LinearDensityDeriv", "C/s$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "DeltaGamma", NULL, SDDS_DOUBLE) ||
        !SDDS_DefineSimpleColumn(SDDSout, "GammaDeriv", "1/m", SDDS_DOUBLE) ||
        !SDDS_WriteLayout(SDDSout)) {
      SDDS_SetError("Problem setting up output file for CSRDRIFT (Stupakov mode)");
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    *active = 1;
  }
  for (i=0; i<nBins; i++) {
    ctHist[i] *= MPCharge*c_mks;
    ctHistDeriv[i] *= MPCharge*sqr(c_mks)/dct;
  }
  if (!SDDS_StartPage(SDDSout, nBins) ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, dGamma, nBins, "DeltaGamma")  ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "LinearDensity")  ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "LinearDensityDeriv")) {
    SDDS_SetError("Problem writing to output file for CSRDRIFT (Stupakov mode)");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  /* use ctHist array for output of s and ctHistDeriv for dGamma/ds */
  for (i=0; i<nBins; i++) {
    ctHist[i] = dct*(i+0.5-nBins/2);
    ctHistDeriv[i] = dGamma[i]/dz;
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHist, nBins, "s") ||
      !SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, ctHistDeriv, nBins, "GammaDeriv") ||
      !SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "z", zTravel, "CaseC", nCaseC,
                          "CaseD1", nCaseD1, "CaseD2", nCaseD2, 
                          "x", x, "dsMax", dsMax, "phi0", phi0, "phi1", phi1,
                          NULL) ||
      !SDDS_WritePage(SDDSout)) {
    SDDS_SetError("Problem writing to output file for CSRDRIFT (Stupakov mode)");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(SDDSout);
}


void apply_edge_effects(
                        double *x, double *xp, double *y, double *yp, 
                        double rho, double n, double beta, double he, double psi, long which_edge
                        )
{
  double h, tan_beta, tan2_beta, sec_beta, sec2_beta, h2;
  double R21, R43;
  double T111, T133, T211, T441, T331, T221, T233, T243, T431, T432;
  double x0, xp0, y0, yp0;

  h = 1/rho;
  R21 = h*(tan_beta=tan(beta));
  R43 = -h*tan(beta-psi);

  h2 = sqr(h);
  T111 = which_edge*h/2*(tan2_beta=sqr(tan_beta));
  T133 = -which_edge*h/2*(sec2_beta=sqr(sec_beta=1./cos(beta)));
  T211 = which_edge==-1?
    -n*h2*tan_beta:
    -h2*(n+tan2_beta/2)*tan_beta;
  T441 = -(T331 = T221 = -which_edge*h*tan2_beta);
  T233 =  which_edge==-1?
    h2*(n+.5+tan2_beta)*tan_beta:
    h2*(n-tan2_beta/2)*tan_beta;
  T243 = which_edge*h*tan2_beta;
  T431 = h2*(2*n+(which_edge==1?sec2_beta:0))*tan_beta;
  T432 = which_edge*h*sec2_beta;
  if (he!=0) {
    double term;
    term = h/2*he*sec2_beta*sec_beta;
    T211 += term;
    T233 -= term;
    T431 -= 2*term;
  }

  x0 = *x;  xp0 = *xp;  y0 = *y;  yp0 = *yp;
  *x  = x0  + T111*sqr(x0) + T133*sqr(y0);
  *xp = xp0 + R21*x0 + T211*sqr(x0) + T221*x0*xp0 + T233*sqr(y0) + T243*y0*yp0;
  *y  = y0  + T331*x0*y0;
  *yp = yp0 + R43*y0 + T441*yp0*x0 + T431*x0*y0 + T432*xp0*y0;
}

/* this is used solely to convert coordinates inside the element for
 * the purpose of generating output.  It ignores misalignments.
 */

void convertFromCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, 
			     long ctMode)
{
  long ip;
  double x, y, xp, yp, *coord;

  for (ip=0; ip<np; ip++) {
    coord = part[ip];

    XP /= (1+X/rho0);
    YP /= (1+X/rho0);
    
    x  =  X*cos_ttilt -  Y*sin_ttilt;
    y  =  X*sin_ttilt +  Y*cos_ttilt;
    xp = XP*cos_ttilt - YP*sin_ttilt;
    yp = XP*sin_ttilt + YP*cos_ttilt;

    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;

    if (ctMode)
      coord[4] /= c_mks;
  }
}


/* this is used solely to undo the transformation done by 
 * convertFromCSBendCoords
 */

void convertToCSBendCoords(double **part, long np, double rho0, 
			     double cos_ttilt, double sin_ttilt, long ctMode)
{
  long ip;
  double x, y, xp, yp, *coord;

  for (ip=0; ip<np; ip++) {
    coord = part[ip];

    x  =   X*cos_ttilt +  Y*sin_ttilt;
    y  =  -X*sin_ttilt +  Y*cos_ttilt;
    xp =  XP*cos_ttilt + YP*sin_ttilt;
    yp = -XP*sin_ttilt + YP*cos_ttilt;

    X  = x;
    Y  = y;
    XP = xp;
    YP = yp;

    XP *= (1+X/rho0);
    YP *= (1+X/rho0);

    if (ctMode)
      coord[4] *= c_mks;
  }
}

#include "fftpackC.h"
long applyLowPassFilter(double *histogram, long bins, 
			double start,   /* in units of Nyquist frequency */
			double end      /* in units of Nyquist frequency */
			)
{
  long i, i1, i2;
  double fraction, dfraction, sum;
  double *realimag;
  long frequencies;

  if (!(realimag = (double*)malloc(sizeof(*realimag)*(bins+2))))
    SDDS_Bomb("allocation failure");

  if (end<start)
    end = start;

  frequencies = bins/2 + 1;
  realFFT2(realimag, histogram, bins, 0);

  i1 = start*frequencies;
  if (i1<0) 
    i1=0;
  if (i1>frequencies-1)
    i1 = frequencies-1;
  
  i2 = end*frequencies;
  if (i2<0) 
    i2=0;
  if (i2>frequencies-1)
    i2 = frequencies-1;
  
  dfraction = i1==i2? 0 : 1./(i2-i1);
  fraction = 1;
  for (i=i1; i<=i2; i++) {
    realimag[2*i  ] *= fraction;
    realimag[2*i+1] *= fraction;
    if ((fraction -= dfraction)<0)
      fraction = 0;
  }
  for (; i<frequencies; i++) {
    realimag[2*i  ] = 0;
    realimag[2*i+1] = 0;
  }

  realFFT2(realimag, realimag, bins, INVERSE_FFT);

  /* copy data to input buffer.
   * normalize to keep the sum constant
   * don't allow negative values 
   */
  for (i=sum=0; i<bins; i++) {
    sum += histogram[i];
    histogram[i] = realimag[i];
  }
  free(realimag);
  return correctDistribution(histogram, bins, sum);
}

long applyLHPassFilters(double *histogram, long bins, 
			double startHP,   /* in units of Nyquist frequency */
			double endHP,     /* in units of Nyquist frequency */
			double startLP,   /* in units of Nyquist frequency */
			double endLP,     /* in units of Nyquist frequency */
                        long clipNegative
			)
{
  long i, i1, i2;
  double fraction, dfraction, sum;
  double *realimag;
  long frequencies;

  if (!(realimag = (double*)malloc(sizeof(*realimag)*(bins+2))))
    SDDS_Bomb("allocation failure");

  if (endLP<startLP)
    endLP = startLP;
  if (endHP<startHP)
    endHP = startHP;
  
  frequencies = bins/2 + 1;
  realFFT2(realimag, histogram, bins, 0);

  if (startLP>0) {
    i1 = startLP*frequencies;
    if (i1<0) 
      i1=0;
    if (i1>frequencies-1)
      i1 = frequencies-1;
    
    i2 = endLP*frequencies;
    if (i2<0) 
      i2=0;
    if (i2>frequencies-1)
      i2 = frequencies-1;
    
    dfraction = i1==i2? 0 : 1./(i2-i1);
    fraction = 1;
    for (i=i1; i<=i2; i++) {
      realimag[2*i  ] *= fraction;
      realimag[2*i+1] *= fraction;
      if ((fraction -= dfraction)<0)
        fraction = 0;
    }
    for (; i<frequencies; i++) {
      realimag[2*i  ] = 0;
      realimag[2*i+1] = 0;
    }
  }
  
  if (startHP>0) {
    i1 = startHP*frequencies;
    if (i1<0) 
      i1=0;
    if (i1>frequencies-1)
      i1 = frequencies-1;
    
    i2 = endHP*frequencies;
    if (i2<0) 
      i2=0;
    if (i2>frequencies-1)
      i2 = frequencies-1;
    
    dfraction = i1==i2? 0 : 1./(i2-i1);
    fraction = 0;
    for (i=0; i<i1; i++) {
      realimag[2*i  ] = 0;
      realimag[2*i+1] = 0;
    }
    for (i=i1; i<=i2; i++) {
      realimag[2*i  ] *= fraction;
      realimag[2*i+1] *= fraction;
      if ((fraction += dfraction)>1)
        fraction = 1;
    }
  }
  
  realFFT2(realimag, realimag, bins, INVERSE_FFT);

  /* copy data to input buffer  */
  for (i=sum=0; i<bins; i++) {
    sum += histogram[i];
    histogram[i] = realimag[i];
  }
  free(realimag);

  if (clipNegative)
    /* normalize to keep the sum constant
     * don't allow negative values 
     */
    return correctDistribution(histogram, bins, sum);
  else
    return 0;
}

long correctDistribution(double *array, long npoints, double desiredSum)
{
  double sum, factor;
  long nz, i;
  for (i=nz=sum=0; i<npoints; i++) {
    if (array[i]<0) {
      nz ++;
      array[i] = 0;
    }
    sum += array[i];
  }
  if (!sum)
    return nz;
  factor = desiredSum/sum;
  for (i=0; i<npoints; i++)
    array[i] *= factor;
  return nz;
}

void computeEtiltCentroidOffset(double *dcoord_etilt, double rho0, double angle, double etilt, double tilt)
{
  /* compute final offsets due to error-tilt of the magnet */
  /* see pages 90-93 of notebook 1 about this */
  double q1a, q2a, q3a;
  double q1b, q2b, q3b;
  double qp1, qp2, qp3; 
  double dz, tan_alpha, k;
  
  if (!etilt) {
    fill_double_array(dcoord_etilt, 6L, 0.0);
    return;
  }

  q1a = (1-cos(angle))*rho0*(cos(etilt)-1);
  q2a = 0;
  q3a = (1-cos(angle))*rho0*sin(etilt);
  qp1 = sin(angle)*cos(etilt);
  qp2 = cos(angle);
  k = sqrt(sqr(qp1)+sqr(qp2));
  qp1 /= k;
  qp2 /= k;
  qp3 = sin(angle)*sin(etilt)/k;
  tan_alpha = 1./tan(angle)/cos(etilt);
  q1b = q1a*tan_alpha/(tan(angle)+tan_alpha);
  q2b = -q1b*tan(angle);
  dz  = sqrt(sqr(q1b-q1a)+sqr(q2b-q2a));
  q3b = q3a + qp3*dz;

  dcoord_etilt[0] = sqrt(sqr(q1b) + sqr(q2b));
  dcoord_etilt[1] = tan(atan(tan_alpha)-(PIo2-angle));
  dcoord_etilt[2] = q3b;
  dcoord_etilt[3] = qp3;
  dcoord_etilt[4] = dz*sqrt(1+sqr(qp3));
  dcoord_etilt[5] = 0;
#ifdef DEBUG
  fprintf(stdout, "pre-tilt offsets due to ETILT=%le:  %le %le %le %le %le\n",
          etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
          dcoord_etilt[3], dcoord_etilt[4]);
  fflush(stdout);
#endif

  /* rotate by tilt to get into same frame as bend equations. */
  rotate_coordinates(dcoord_etilt, tilt);
#ifdef DEBUG
  fprintf(stdout, "offsets due to ETILT=%le:  %le %le %le %le %le\n",
          etilt, dcoord_etilt[0], dcoord_etilt[1], dcoord_etilt[2],
          dcoord_etilt[3], dcoord_etilt[4]);
  fflush(stdout);
#endif
}

void readWakeFilterFile(long *values, 
                        double **freq, double **real, double **imag, 
                        char *freqName, char *realName, char *imagName,
                        char *filename)
{
  SDDS_DATASET SDDSin;
  long i;
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename) || !SDDS_ReadPage(&SDDSin)) {
    fprintf(stderr, "Error: unable to open or read CSRCSBEND wake filter file %s\n", filename);
    exitElegant(1);
  }
  if ((*values = SDDS_RowCount(&SDDSin))<2) {
    fprintf(stderr, "Error: too little data in CSRCSBEND wake filter file %s\n", filename);
    exitElegant(1);
  }
  if (!freqName || !strlen(freqName))
    SDDS_Bomb("WFF_FREQ_COLUMN is blank in CSRCSBEND");
  if (SDDS_CheckColumn(&SDDSin, freqName, "Hz", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: column %s invalid in CSRCSBEND wake filter file %s---check existence, type, and units (Hz).\n", 
            freqName, filename);
    exitElegant(1);
  }
  if (!realName || !strlen(realName))
    SDDS_Bomb("WFF_REAL_COLUMN is blank in CSRCSBEND");
  if (SDDS_CheckColumn(&SDDSin, realName, NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: column %s invalid in CSRCSBEND wake filter file %s---check existence and type.\n", 
            realName, filename);
    exitElegant(1);
  }
  if (!imagName || !strlen(imagName))
    SDDS_Bomb("WFF_IMAG_COLUMN is blank in CSRCSBEND");
  if (SDDS_CheckColumn(&SDDSin, imagName, NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stderr, "Error: column %s invalid in CSRCSBEND wake filter file %s---check existence and type.\n", 
            imagName, filename);
    exitElegant(1);
  }
  if (!(*freq=SDDS_GetColumnInDoubles(&SDDSin, freqName)) ||
      !(*real=SDDS_GetColumnInDoubles(&SDDSin, realName)) ||
      !(*imag=SDDS_GetColumnInDoubles(&SDDSin, imagName))) {
    fprintf(stderr, "Problem getting data from CSRCSBEND wake filter file %s.\n", filename);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  for (i=1; i<*values; i++) {
    if ((*freq)[i-1]>=(*freq)[i]) {
      fprintf(stderr, "Error: frequency data is not monotonically increasing in CSRCSBEND wake filter file %s.\n", filename);
      exitElegant(1);
    }
  }
}

void applyFilterTable(double *function, long bins, double dx, long fValues,
                     double *fFreq, double *fReal, double *fImag)
{
  long i, i1, i2;
  double f;
  double *realimag, dfrequency, length;
  long frequencies;
  double sum;
  
  if (!(realimag = (double*)malloc(sizeof(*realimag)*(bins+2))))
    SDDS_Bomb("allocation failure");

  frequencies = bins/2 + 1;
  length = dx*(bins-1);
  dfrequency = 1.0/length;
  realFFT2(realimag, function, bins, 0);
  
  for (i=0; i<frequencies; i++) {
    long code;
    i1 = 2*i+0;
    i2 = 2*i+1;
    f = i*dfrequency;
    realimag[i1] *= interp(fReal, fFreq, fValues, f, 0, 1, &code);
    realimag[i2] *= interp(fImag, fFreq, fValues, f, 0, 1, &code);
  }
  realFFT2(realimag, realimag, bins, INVERSE_FFT);

  /* copy data to input buffer.
   */
  for (i=sum=0; i<bins; i++) 
    function[i] = realimag[i];
  free(realimag);
}

void addRadiationKick(double *Qx, double *Qy, double *dPoP, double *sigmaDelta2, long sqrtOrder,
		      double x, double h0, double Fx, double Fy,
		      double ds, double radCoef, double dsISR, double isrCoef,
                      long distributionBased, long includeOpeningAngle, double meanPhotonsPerMeter,
                      double normalizedCriticalEnergy0, double Po)
{
  double f, xp, yp, F2, F, deltaFactor, dsFactor;
  double nMean, dDelta, thetaRms;
  long i, nEmitted;
  double y, logy;
  double normalizedCriticalEnergy;
  
  f = (1+x*h0)/EXSQRT(sqr(1+*dPoP)-sqr(*Qx)-sqr(*Qy), sqrtOrder);
  xp = *Qx*f;
  yp = *Qy*f;
  dsFactor = EXSQRT(sqr(1+x*h0)+sqr(xp)+sqr(yp), sqrtOrder);
  F2 = sqr(Fx)+sqr(Fy);

  if (!distributionBased) {
    deltaFactor = sqr(1 + *dPoP);
    *Qx /= (1 + *dPoP);
    *Qy /= (1 + *dPoP);
    if (radCoef)
      *dPoP -= radCoef*deltaFactor*F2*ds*dsFactor;
    if (isrCoef>0)
      /* The minus sign is for consistency with the previous version. */
      *dPoP -= isrCoef*deltaFactor*pow(F2,0.75)*sqrt(dsISR*dsFactor)*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
    if (sigmaDelta2)
      *sigmaDelta2 += sqr(isrCoef*deltaFactor)*pow(F2,1.5)*dsISR*dsFactor;
    *Qx *= (1 + *dPoP);
    *Qy *= (1 + *dPoP);
  } else {
    F = sqrt(F2);
    /* Compute the mean number of photons emitted = meanPhotonsPerMeter*meters */
    nMean = meanPhotonsPerMeter*dsISR*dsFactor*F;
    /* Pick the actual number of photons emitted from Poisson distribution */
    nEmitted = inversePoissonCDF(nMean, random_2(1));
    /* Adjust normalized critical energy to local field strength (FSE is already included via rho_actual) */
    /* Variation with energy offset included below. */
    normalizedCriticalEnergy = normalizedCriticalEnergy0*(1+*dPoP)*F;
    /* For each photon, pick its energy and emission angles */
    for (i=0; i<nEmitted; i++) {
      /* Pick photon energy normalized to critical energy */
      y=pickNormalizedPhotonEnergy(random_2(1));
      /* Multiply by critical energy normalized to beam energy, adjusting for variation with
       * individual electron energy offset. */
      dDelta = normalizedCriticalEnergy*(1 + *dPoP)*y;
      photonCount ++;
      energyCount += y;
      /* Change the total electron momentum */
      *dPoP -= dDelta;
      if (includeOpeningAngle) {
        /* Compute rms spread in electron angle = (rms photon angle)*dDelta */
        logy = log10(y);
        thetaRms = dDelta*pow(10,
                              -2.418673276661232e-01
                              + logy*(-4.472680955382907e-01+logy*(-4.535350424882360e-02
                                                                   -logy*6.181818621278201e-03)))/Po;
        /* Compute change in electron angle due to photon angle */
        xp += thetaRms*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
        yp += thetaRms*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      }
    }
    f = (1 + *dPoP)/EXSQRT(sqr(1+x*h0)+sqr(xp)+sqr(yp), sqrtOrder);
    *Qx = xp*f;
    *Qy = yp*f;
  }
  
}

long inversePoissonCDF(double mu, double C)
{
  double sum, expMinusMu, term;
  long r, rMax;
  
  r = 0;
  if ((rMax = 50*mu)<10)
    rMax = 10;
  expMinusMu = exp(-mu);
  term = sum = expMinusMu;
  while (r<=rMax && C>=sum) {
    term *= mu/(++r);
    sum += term;
  }
  /* fprintf(stderr, "inversePoissonCDF: r=%ld for mu=%e, C=%e\n", r, mu, C); */
  return r;
}

/* Return randomly-chosen photon energy normalized to the critical energy */
double pickNormalizedPhotonEnergy(double RN) 
{
  long interpCode;
  double value;
  static double ksiTable[200] = {
1.000000000000000e-07, 1.103351074554523e-07, 1.217383310646075e-07, 1.343200559586096e-07, 1.482020747927429e-07,
1.635189113578160e-07, 1.804187271717404e-07, 1.990651067932036e-07, 2.196385495378513e-07, 2.423384085535426e-07,
2.673842889192374e-07, 2.950186137528527e-07, 3.255088868939687e-07, 3.591505332405233e-07, 3.962690515521040e-07,
4.372236992560780e-07, 4.824109227537678e-07, 5.322685173365927e-07, 5.872789374272963e-07, 6.479745823257918e-07,
7.149429940622619e-07, 7.888329490825732e-07, 8.703595415651299e-07, 9.603117573208774e-07, 1.059560346180786e-06,
1.169066741218290e-06, 1.289890855361064e-06, 1.423201921186706e-06, 1.570290408531488e-06, 1.732581081736759e-06,
1.911644944229382e-06, 2.109214732888093e-06, 2.327202951554652e-06, 2.567720983189001e-06, 2.833097369770374e-06,
3.125899929822278e-06, 3.448963034857157e-06, 3.805415579829938e-06, 4.198708930670648e-06, 4.632648449870741e-06,
5.111434736502657e-06, 5.639704554113462e-06, 6.222573518736223e-06, 6.865680963315421e-06, 7.575252265320169e-06,
8.358158724774669e-06, 9.221982709850737e-06, 1.017508139438384e-05, 1.122668092062936e-05, 1.238696398672945e-05,
1.366716918178818e-05, 1.507968136669955e-05, 1.663817388115531e-05, 1.835773664119261e-05, 2.025502748788774e-05,
2.234840008736815e-05, 2.465811862080412e-05, 2.720654509486478e-05, 3.001836976776838e-05, 3.312079173855694e-05,
3.654384295607505e-05, 4.032066206311618e-05, 4.448784496925588e-05, 4.908569920644159e-05, 5.415873276357582e-05,
5.975605436563109e-05, 6.593190648243325e-05, 7.274602260436053e-05, 8.026436452046312e-05, 8.855971070475770e-05,
9.771244913354602e-05, 1.078111117026042e-04, 1.189534517802125e-04, 1.312473286188584e-04, 1.448118705606887e-04,
1.597782996280689e-04, 1.762914808166189e-04, 1.945112638414439e-04, 2.146141870658145e-04, 2.367947478384840e-04,
2.612676279365202e-04, 2.882697280504828e-04, 3.180626634504074e-04, 3.509347182037930e-04, 3.872040386295810e-04,
4.272217166086854e-04, 4.713754443547046e-04, 5.200925166529275e-04, 5.738444084509028e-04, 6.331514454110419e-04,
6.985881553542330e-04, 7.707878748140519e-04, 8.504493031356108e-04, 9.383435740575402e-04, 1.035322090321881e-03,
1.142323583631410e-03, 1.260383487096962e-03, 1.390644637378946e-03, 1.534368742954572e-03, 1.692947192375836e-03,
1.867914439234830e-03, 2.060964191854564e-03, 2.273966197872862e-03, 2.508982774009835e-03, 2.768287894108917e-03,
3.054391670052556e-03, 3.370064913211098e-03, 3.718364390303022e-03, 4.102660002531562e-03, 4.526671789275104e-03,
4.994505863796759e-03, 5.510692966986859e-03, 6.080227102658105e-03, 6.708621448554406e-03, 7.401960924834442e-03,
8.166960997457158e-03, 9.011022498794720e-03, 9.942316075907171e-03, 1.096985907697749e-02, 1.210360516825099e-02,
1.335452196639184e-02, 1.473471854468308e-02, 1.625755786686624e-02, 1.793779326078131e-02, 1.979167811637735e-02,
2.183715833862031e-02, 2.409403673367596e-02, 2.658418068941497e-02, 2.933167681381903e-02, 3.236312132092089e-02,
3.570786033002121e-02, 3.939830569107346e-02, 4.347015239952341e-02, 4.796281661196050e-02, 5.291978786613488e-02,
5.838910396032759e-02, 6.442366592223334e-02, 7.108188831787103e-02, 7.842822357081021e-02, 8.653385973120035e-02,
9.547720673026440e-02, 1.053448316706813e-01, 1.162322544203413e-01, 1.282449697899047e-01, 1.414991970199129e-01,
1.561232236887442e-01, 1.722586122799903e-01, 1.900616973883389e-01, 2.097047392306722e-01, 2.313778527830923e-01,
2.552908366685073e-01, 2.816753658394816e-01, 3.107867644741113e-01, 3.429067722728065e-01, 3.783463152734567e-01,
4.174487166108054e-01, 4.605924179038830e-01, 5.081949415377639e-01, 5.607170865377987e-01, 6.186676294134220e-01,
6.826074965088431e-01, 7.531554325044472e-01, 8.309943513916955e-01, 9.168782178988778e-01, 1.011638437307961e+00,
1.116191954411967e+00, 1.231550862322391e+00, 1.358832474428039e+00, 1.499269100004806e+00, 1.654219599620752e+00,
1.825183916868001e+00, 2.013817817791925e+00, 2.221947831729780e+00, 2.451587713539253e+00, 2.704960411634972e+00,
2.984519634347505e+00, 3.292972664221137e+00, 3.633303772560254e+00, 4.008807416353689e+00, 4.423119788364888e+00,
4.880253623874576e+00, 5.384631444881934e+00, 5.941135706944927e+00, 6.555154946882635e+00, 7.232636842499024e+00,
7.980135322277263e+00, 8.804886289535018e+00, 9.714875109915180e+00, 1.071891743371295e+01, 1.182672581369469e+01,
1.304902401296616e+01, 1.439764568247248e+01, 1.588565738231238e+01, 1.752745256838863e+01, 1.933892408641842e+01,
2.133760842747432e+01, 2.354287285156119e+01, 2.597604764912193e+01, 2.866068635656761e+01, 3.162277660168377e+01,
  };
  static double FTable[200] = {
0.000000000000000e+00, 1.916076787477782e-04, 3.896006996482199e-04, 5.941918318862451e-04, 8.056009324383097e-04,
1.024055848381587e-03, 1.249790750550654e-03, 1.483048166648730e-03, 1.724078746036354e-03, 1.973142196708657e-03,
2.230505581886648e-03, 2.496445345396121e-03, 2.771247236692068e-03, 3.055207274452791e-03, 3.348630028390361e-03,
3.651830594751359e-03, 3.965134731031564e-03, 4.288879835176022e-03, 4.623413241840414e-03, 4.969094094184835e-03,
5.326293748409966e-03, 5.695396753910589e-03, 6.076799201662367e-03, 6.470910424341261e-03, 6.878153743490802e-03,
7.298967431515415e-03, 7.733803160081558e-03, 8.183127443537616e-03, 8.647422816544402e-03, 9.127188753348749e-03,
9.622940276393589e-03, 1.013520903749105e-02, 1.066454503150837e-02, 1.121151744173274e-02, 1.177671348365064e-02,
1.236073899369794e-02, 1.296422080554008e-02, 1.358780748228750e-02, 1.423216848805338e-02, 1.489799412460975e-02,
1.558599872144761e-02, 1.629692120583610e-02, 1.703152471578071e-02, 1.779059568966145e-02, 1.857494805017956e-02,
1.938542355142036e-02, 2.022289197174982e-02, 2.108824911944652e-02, 2.198242222447622e-02, 2.290636998738409e-02,
2.386108352008151e-02, 2.484758298036615e-02, 2.586692442491193e-02, 2.692019946744418e-02, 2.800853716992024e-02,
2.913309895920114e-02, 3.029508724043934e-02, 3.149574455091229e-02, 3.273635664445640e-02, 3.401824527268329e-02,
3.534277891084329e-02, 3.671137126806664e-02, 3.812548585471972e-02, 3.958662612641901e-02, 4.109634873515033e-02,
4.265626122542306e-02, 4.426802842804741e-02, 4.593335936539433e-02, 4.765402350963371e-02, 4.943184769110832e-02,
5.126872354805775e-02, 5.316659282366442e-02, 5.512746484065590e-02, 5.715341374017695e-02, 5.924658623150298e-02,
6.140918640941097e-02, 6.364349310252398e-02, 6.595185823803101e-02, 6.833671466591608e-02, 7.080056062143326e-02,
7.334597653614160e-02, 7.597562485484224e-02, 7.869225776043576e-02, 8.149870143512994e-02, 8.439787179667740e-02,
8.739277614110098e-02, 9.048652052440445e-02, 9.368229400251835e-02, 9.698338259579789e-02, 1.003931731738744e-01,
1.039151601553213e-01, 1.075529299247917e-01, 1.113101721273651e-01, 1.151906862423360e-01, 1.191983871317649e-01,
1.233372898347266e-01, 1.276115170379167e-01, 1.320253087906072e-01, 1.365830262538971e-01, 1.412891370418868e-01,
1.461482173909104e-01, 1.511649654280279e-01, 1.563442022251273e-01, 1.616908577721921e-01, 1.672099659551390e-01,
1.729066816736924e-01, 1.787862779608226e-01, 1.848541324986933e-01, 1.911157129897795e-01, 1.975765981791920e-01,
2.042424693223850e-01, 2.111190968366426e-01, 2.182123129823302e-01, 2.255280364093208e-01, 2.330722555939911e-01,
2.408510146807283e-01, 2.488703695479055e-01, 2.571364147677684e-01, 2.656552557248533e-01, 2.744329918358574e-01,
2.834756510338721e-01, 2.927892168723024e-01, 3.023795848158047e-01, 3.122525396990814e-01, 3.224136624350130e-01,
3.328683532505147e-01, 3.436217660758924e-01, 3.546787751835536e-01, 3.660438466342482e-01, 3.777210511874600e-01,
3.897139689097260e-01, 4.020256374046105e-01, 4.146583795759362e-01, 4.276137967068266e-01, 4.408926345108801e-01,
4.544946994027770e-01, 4.684186411623745e-01, 4.826619073829413e-01, 4.972205662043195e-01, 5.120891723276200e-01,
5.272605134477985e-01, 5.427255018050126e-01, 5.584729557729631e-01, 5.744894096219003e-01, 5.907588309168045e-01,
6.072624513895929e-01, 6.239785205673239e-01, 6.408820784452591e-01, 6.579446823895981e-01, 6.751342192611512e-01,
6.924146905220633e-01, 7.097460264751286e-01, 7.270839321656986e-01, 7.443798099419868e-01, 7.615807315005799e-01,
7.786294876113157e-01, 7.954647879517510e-01, 8.120215670459850e-01, 8.282314411609915e-01, 8.440233375507971e-01,
8.593244148619452e-01, 8.740611430503441e-01, 8.881606055288961e-01, 9.015520524341847e-01, 9.141687939214495e-01,
9.259502059677074e-01, 9.368438193470870e-01, 9.468075155760183e-01, 9.558117659236037e-01, 9.638415906288208e-01,
9.708980284014210e-01, 9.769991556010101e-01, 9.821804314564566e-01, 9.864941142754962e-01, 9.900075473722704e-01,
9.928006021230259e-01, 9.949622418923878e-01, 9.965864066620354e-01, 9.977674409043544e-01, 9.985957390441925e-01,
9.991538996011060e-01, 9.995138052403904e-01, 9.997348433078885e-01, 9.998634876235176e-01, 9.999340435105117e-01,
9.999702897374164e-01, 9.999876125809349e-01, 9.999952573365697e-01, 9.999983471293699e-01, 9.999994807411241e-01,
9.999998545219742e-01, 9.999999640793696e-01, 9.999999922833868e-01, 9.999999985785893e-01, 9.999999997790499e-01,
9.999999999715266e-01, 9.999999999970198e-01, 9.999999999997560e-01, 9.999999999999872e-01, 1.000000000000000e+00,
  };
  value = interp(ksiTable, FTable, 200, RN, 0, 2, &interpCode);
  if (!interpCode)
    return ksiTable[0];
  return value;
}

void addCorrectorRadiationKick(double **coord, long np, ELEMENT_LIST *elem, long type, double Po, double *sigmaDelta2, long disableISR)
{
  double F2;
  double kick, length;
  double isrCoef, radCoef, dp, p, beta0, beta1, deltaFactor;
  short isr, sr;
  long i;

  if (!np)
    return;

  isr = sr = 0;

  switch (type) {
  case T_HCOR:
    kick = ((HCOR*)elem->p_elem)->kick;
    if ((length = ((HCOR*)elem->p_elem)->length)==0) 
      length = ((HCOR*)elem->p_elem)->lEffRad;
    if (((HCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((HCOR*)elem->p_elem)->isr) 
	isr = 1;
    }
    break;
  case T_VCOR:
    kick = ((VCOR*)elem->p_elem)->kick;
    if ((length = ((VCOR*)elem->p_elem)->length)==0) 
      length = ((VCOR*)elem->p_elem)->lEffRad;
    if (((VCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((VCOR*)elem->p_elem)->isr) 
	isr = 1;
    }
    break;
  case T_HVCOR:
    kick = sqrt(sqr(((HVCOR*)elem->p_elem)->xkick)+sqr(((HVCOR*)elem->p_elem)->ykick));
    if ((length = ((HVCOR*)elem->p_elem)->length)==0) 
      length = ((HVCOR*)elem->p_elem)->lEffRad;
    if (((HVCOR*)elem->p_elem)->synchRad) {
      sr = 1;
      if (((HVCOR*)elem->p_elem)->isr) 
	isr = 1;
    }
    break;
  }
  if (sr==0 || length==0) 
    return ;
  if (disableISR)
    isr = 0;
  radCoef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
  isrCoef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);

  F2 = sqr(kick/length);
  for (i=0; i<np; i++) {
    dp = coord[i][5];
    p = Po*(1+dp);
    beta0 = p/sqrt(sqr(p)+1);
    deltaFactor = sqr(1+dp);
    dp -= radCoef*deltaFactor*F2*length;
    if (isr)
      dp += isrCoef*deltaFactor*pow(F2, 0.75)*sqrt(length)*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
    if (sigmaDelta2)
      *sigmaDelta2 += sqr(isrCoef*deltaFactor)*pow(F2, 1.5)*length;
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[i][5] = dp;
    coord[i][4] = beta1*coord[i][4]/beta0;
  }
  if (sigmaDelta2)
    *sigmaDelta2 /= np;
}
  
