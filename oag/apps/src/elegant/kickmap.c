/*************************************************************************\
* Copyright (c) 2008 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2008 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: kickmap.c
 * purpose: track particles through a kick map (function of x and y)
 * 
 * Michael Borland, Weiming Guo, 2008
 */
#include "mdb.h"
#include "track.h"

void initializeUndulatorKickMap(UKICKMAP *map);
long interpolateUndulatorKickMap(double *xpFactor, double *ypFactor, UKICKMAP *map, double x, double y);

long trackUndulatorKickMap(
               double **particle,    /* array of particles */
               double **accepted,    /* acceptance array */
               long nParticles,      /* number of particles */
               double pRef,          /* central momentum */
               UKICKMAP *map,
               double zStart
               )
{
  long ip, iTop, ik, nKicks;
  double *coord;
  double eomc, H;
  double dxpFactor, dypFactor;
  double length, fieldFactor;
  double radCoef, isrCoef, sxpCoef, beta, alpha, sqrtBeta, deltaFactor, delta;
  double I1, I2, I3, I4, I5;
  double sqrtI3, sqrtI5;
  
  length = map->length;
  fieldFactor = map->fieldFactor;
  if (!map->initialized)
    initializeUndulatorKickMap(map);

  if ((nKicks=map->nKicks)<1)
    bomb("N_KICKS must be >=1 for UKICKMAP", NULL);

  radCoef = isrCoef = sxpCoef = 0;
  if (map->synchRad)
    /* radCoef*I2 is d((P-Po)/Po) per step for the on-axis, on-momentum particle */
    radCoef = 2./3*particleRadius*ipow(pRef, 3);
  if (map->isr) {
    /* isrCoef*sqrt(I3) is the RMS increase in dP/P per step due to incoherent SR.  */
    isrCoef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(pRef)*137.0359895);
    /* sxpCoef*sqrt(I5) is related to the increase in the RMS divergence per step due to incoherent SR */
    sxpCoef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(pRef)*137.0359895);
  }
  
  length /= nKicks;
  
  eomc = particleCharge/particleMass/c_mks; 

  if (map->dx || map->dy || map->dz)
    offsetBeamCoordinates(particle, nParticles, map->dx, map->dy, map->dz);
  if (map->tilt)
    rotateBeamCoordinates(particle, nParticles, map->tilt);
  
  iTop = nParticles-1;
  for (ik=0; ik<nKicks; ik++) {
    sqrtI3 = sqrtI5 = 0;
    if (radCoef || sxpCoef) {
      double S11, S12, S22, S16, S26, S66;
      double S11beta, S12beta, S22beta, emit;

      if (nKicks!=map->periods)
        bomb("Number of kicks must equal number of periods for UKICKMAP radiation tracking", NULL);
      
      if (sxpCoef) {
        /* These calls are inefficient (repeat some calculations), but better than other routines we have now */
#if !USE_MPI
        rms_emittance(particle, 0, 1, iTop+1, &S11, &S12, &S22);
        rms_emittance(particle, 0, 5, iTop+1, &S11, &S16, &S66);
        rms_emittance(particle, 1, 5, iTop+1, &S22, &S26, &S66);
#else
        rms_emittance_p(particle, 0, 1, iTop+1, &S11, &S12, &S22);
        rms_emittance_p(particle, 0, 5, iTop+1, &S11, &S16, &S66);
        rms_emittance_p(particle, 1, 5, iTop+1, &S22, &S26, &S66);
#endif

        if (S66) {
          S11beta = S11 - sqr(S16)/S66;
          S12beta = S12 - S16*S26/S66;
          S22beta = S22 - sqr(S26)/S66;
        } else {
          S11beta = S11;
          S12beta = S12;
          S22beta = S22;
        }
        beta = alpha = sqrtBeta = 0;
        if ((emit = S11beta*S22beta-sqr(S12beta))>0) {
          emit = sqrt(emit);
          beta = sqrt(S11beta/emit);
          sqrtBeta = sqrt(beta);
          alpha = -S12beta/emit;
        }
      } else {
        emit = 0;
        sqrtBeta = beta = 1;
        alpha = 0;
      }
      I1 = I2 = I3 = I4 = I5 = 0;
      AddWigglerRadiationIntegrals(map->length/map->periods, 2, map->radiusInternal,
                                   0.0, 0.0, beta, alpha,
                                   &I1, &I2, &I3, &I4, &I5);
      if (sxpCoef) {
        if (I3<0 || I5<0)
          bomb("I3 or I5 is negative in UKICKMAP", NULL);
        sqrtI3 = sqrt(I3);
        sqrtI5 = sqrt(I5);
      }
    }
    
    if (isSlave || !notSinglePart) {
      for (ip=0; ip<=iTop; ip++) {
        coord = particle[ip];
        
        /* 1. go through half length */
        coord[0] += coord[1]*length/2.0;
        coord[2] += coord[3]*length/2.0;
        coord[4] += length/2.0*sqrt(1+sqr(coord[1])+sqr(coord[3]));
        
        /* 2. apply the kicks 
         * use interpolation to get dxpFactor and dypFactor 
         */
        if (!interpolateUndulatorKickMap(&dxpFactor, &dypFactor, map, coord[0], coord[2])) {
          /* particle is lost */
          swapParticles(particle[ip], particle[iTop]); 
          if (accepted)
            swapParticles(accepted[ip], accepted[iTop]);
          particle[iTop][4] = zStart;
          particle[iTop][5] = pRef*(1+particle[iTop][5]);
          iTop--;
          ip--;
        } else {
          H = pRef*(1+coord[5])/eomc;
          coord[1] += dxpFactor*sqr(fieldFactor/H)/nKicks;
          coord[3] += dypFactor*sqr(fieldFactor/H)/nKicks;
          
          /* 3. go through another half length */
          coord[0] += coord[1]*length/2.0;
          coord[2] += coord[3]*length/2.0;
          coord[4] += length/2.0*sqrt(1+sqr(coord[1])+sqr(coord[3]));
        }
        
        /* 3. Optionally apply synchrotron radiation kicks */
        if (radCoef || isrCoef) {
          delta = coord[5];
          deltaFactor = ipow(1+delta, 2);
          if (radCoef) 
            coord[5] -= radCoef*I2*deltaFactor;
          if (isrCoef)
            coord[5] += isrCoef*sqrtI3*deltaFactor*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
          if (sxpCoef && sqrtBeta)
            coord[1] += sxpCoef*sqrtI5*(1+delta)/sqrtBeta*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
        }
      }
    }
  }
  

  if (map->tilt)
    rotateBeamCoordinates(particle, nParticles, -map->tilt);
  if (map->dx || map->dy || map->dz)
    offsetBeamCoordinates(particle, nParticles, -map->dx, -map->dy, -map->dz);

  return iTop+1;
}

void initializeUndulatorKickMap(UKICKMAP *map)
{
  SDDS_DATASET SDDSin;
  double *x=NULL, *y=NULL, *xpFactor=NULL, *ypFactor=NULL;
  long nx;

  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, map->inputFile) ||
      SDDS_ReadPage(&SDDSin)<=0 ||
      !(x=SDDS_GetColumnInDoubles(&SDDSin, "x")) || !(y=SDDS_GetColumnInDoubles(&SDDSin, "y"))) {
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!check_sdds_column(&SDDSin, "x", "m") ||
      !check_sdds_column(&SDDSin, "y", "m")) {
    fprintf(stderr, "UKICKMAP input file must have x and y in m (meters)\n");
    exit(1);
  }

  if (!(xpFactor=SDDS_GetColumnInDoubles(&SDDSin, "xpFactor")) || !(ypFactor=SDDS_GetColumnInDoubles(&SDDSin, "ypFactor"))) {
    fprintf(stderr, "UKICKMAP input file must have both (xpFactor, ypFactor)\n");
    exit(1);
  }
  if (!check_sdds_column(&SDDSin, "xpFactor", "(T*m)$a2$n") ||
      !check_sdds_column(&SDDSin, "ypFactor", "(T*m)$a2$n")) {
    fprintf(stderr, "UKICKMAP input file must have xpFactor and ypFactor with units of (T*m)$a2$n\n");
    exit(1);
  }
  
  if (!(map->points=SDDS_CountRowsOfInterest(&SDDSin)) || map->points<2) {
    fprintf(stdout, "file %s for UKICKMAP element has insufficient data\n", map->inputFile);
    fflush(stdout);
    exit(1);
  }
  SDDS_Terminate(&SDDSin);

  if (map->xyFactor!=1) {
    long i;
    for (i=0; i<map->points; i++) {
      x[i] *= map->xyFactor;
      y[i] *= map->xyFactor;
    }
  }
  
  /* It is assumed that the data is ordered so that x changes fastest.
   * This can be accomplished with sddssort -column=y,incr -column=x,incr
   * The points are assumed to be equipspaced.
   */
  nx = 1;
  map->xmin = x[0];
  while (nx<map->points) {
    if (x[nx-1]>x[nx])
      break;
    nx ++;
  }
  if ((nx>=map->points) || nx<=1 || y[0]>y[nx] || (map->ny = map->points/nx)<=1) {
    fprintf(stdout, "file %s for UKICKMAP element doesn't have correct structure or amount of data\n",
            map->inputFile);
    fflush(stdout);
    fprintf(stdout, "nx = %ld, ny=%ld\n", map->nx, map->ny);
    fflush(stdout);
    exit(1);
  }
  map->nx = nx;
  map->xmax = x[nx-1];
  map->dxg = (map->xmax-map->xmin)/(nx-1);
  map->ymin = y[0];
  map->ymax = y[map->points-1];
  map->dyg = (map->ymax-map->ymin)/(map->ny-1);
  fprintf(stdout, "UKICKMAP element from file %s: nx=%ld, ny=%ld, dxg=%e, dyg=%e, x:[%e, %e], y:[%e, %e]\n",
          map->inputFile, map->nx, map->ny, map->dxg, map->dyg, 
          map->xmin, map->xmax, 
          map->ymin, map->ymax);
  free(x);
  free(y);
  map->xpFactor = xpFactor;
  map->ypFactor = ypFactor;
  map->initialized = 1;
}

long interpolateUndulatorKickMap(double *xpFactor, double *ypFactor, UKICKMAP *map, double x, double y)
{
  double Fa, Fb, fx, fy;
  long ix, iy;
  
  if (isnan(x) || isnan(y) || isinf(x) || isinf(y))
    return 0;
  
  ix = (x-map->xmin)/map->dxg ;
  iy = (y-map->ymin)/map->dyg ;
  if (ix<0 || iy<0 || ix>map->nx-1 || iy>map->ny-1)
    return 0;
  
  fx = (x-(ix*map->dxg+map->xmin))/map->dxg;
  fy = (y-(iy*map->dyg+map->ymin))/map->dyg;

  Fa = (1-fy)*map->xpFactor[ix+iy*map->nx] + fy*map->xpFactor[ix+(iy+1)*map->nx];
  Fb = (1-fy)*map->xpFactor[ix+1+iy*map->nx] + fy*map->xpFactor[ix+1+(iy+1)*map->nx];
  *xpFactor = (1-fx)*Fa+fx*Fb;

  Fa = (1-fy)*map->ypFactor[ix+iy*map->nx] + fy*map->ypFactor[ix+(iy+1)*map->nx];
  Fb = (1-fy)*map->ypFactor[ix+1+iy*map->nx] + fy*map->ypFactor[ix+1+(iy+1)*map->nx];
  *ypFactor = (1-fx)*Fa+fx*Fb;

  return 1;
}


