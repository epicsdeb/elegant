/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <ctype.h>
#include "mdb.h"
#include "track.h"
#include "sort.h"

long evaluateLostWithOpenSides(long code, double dx, double dy, double xsize, double ysize);

/* routine: rectangular_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          rectangular collimator
 *
 * Michael Borland, 1989
 */

#ifdef VAX_VMS
#define isnan(x) 0
#define isinf(x) 0
#endif

long rectangular_collimator(
                            double **initial, RCOL *rcol, long np, double **accepted, double z,
                            double Po
                            )
{
  double length, *ini;
  long ip, itop, is_out, lost, openCode;
  double xsize, ysize;
  double x_center, y_center;
  double x1, y1, zx, zy, dx, dy;

  xsize  = rcol->x_max;
  ysize  = rcol->y_max;
  x_center = rcol->dx;
  y_center = rcol->dy;

  if (xsize<=0 && ysize<=0) {
    exactDrift(initial, np, rcol->length);
    return(np);
  }
  openCode = determineOpenSideCode(rcol->openSide);
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    dx = ini[0] - x_center;
    dy = ini[2] - y_center;
    lost = 0;
    if ((xsize>0 && fabs(dx) > xsize) ||
        (ysize>0 && fabs(dy) > ysize)) {
      lost = openCode ? evaluateLostWithOpenSides(openCode, dx, dy, xsize, ysize) : 1;
    } else if (isinf(ini[0]) || isinf(ini[2]) ||
               isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
     if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }
  if (np==0 || (length=rcol->length)<=0) {
    return(np);
  }

  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    x1 = ini[0] + length*ini[1];
    y1 = ini[2] + length*ini[3];
    dx = x1 - x_center;
    dy = y1 - y_center;
    is_out = 0;
    if (xsize>0 && fabs(dx)>xsize)
      is_out += 1*(openCode?evaluateLostWithOpenSides(openCode, dx, 0, xsize, ysize):1);
    if (ysize>0 && fabs(dy)>ysize)
      is_out += 2*(openCode?evaluateLostWithOpenSides(openCode, 0, dy, xsize, ysize):1);
    if (isinf(x1) || isinf(y1) || isnan(x1) || isnan(y1) )
      is_out += 4;
    if (is_out&4) {
      ini[4] = z+length;
      ini[0] = x1;
      ini[2] = y1;
      ini[5] = Po*(1+ini[5]);
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      --np;
    } else if (is_out) {
      if (!openCode) {
        zx = zy = DBL_MAX;
        if (is_out&1 && ini[1]!=0)
          zx = (SIGN(ini[1])*xsize-(ini[0]-x_center))/ini[1];
        if (is_out&2 && ini[3]!=0)
          zy = (SIGN(ini[3])*ysize-(ini[2]-y_center))/ini[3];
        if (zx<zy) {
          ini[0] += ini[1]*zx;
          ini[2] += ini[3]*zx;
          ini[4] = z+zx; 
        }
        else {
          ini[0] += ini[1]*zy;
          ini[2] += ini[3]*zy;
          ini[4] = z+zy;
        }
      }
      ini[5] = Po*(1+ini[5]);
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      --np;
    }
    else {
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
      ini[0] = x1;
      ini[2] = y1;
    }
  }
  return(np);
}


/* routine: limit_amplitudes()
 * purpose: eliminate particles with (x,y) larger than given values
 *
 * Michael Borland, 1989
 */

long limit_amplitudes(
    double **coord, double xmax, double ymax, long np, double **accepted,
    double z, double Po, long extrapolate_z, long openCode)
{
    long ip, itop, is_out;
    double *part;
    double dz, dzx, dzy;

    if (xmax<0 && ymax<0) {
      return(np);
    }

    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        is_out = 0;
        if (xmax>0 && fabs(part[0])>xmax)
            is_out += 1;
        if (ymax>0 && fabs(part[2])>ymax)
            is_out += 2;
        if (openCode)
          is_out *= evaluateLostWithOpenSides(openCode, part[0], part[2], xmax, ymax);
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) )
            is_out += 4;
        dz = 0;
        if (is_out && !(is_out&4) && !openCode && extrapolate_z) {
            /* find the actual position of loss, assuming a drift preceded with 
             * the same aperture 
             */
            dzx = dzy = -DBL_MAX;
            if (is_out&1 && part[1]!=0)
                dzx = (part[0]-SIGN(part[1])*xmax)/part[1];
            if (is_out&2 && part[3]!=0)
                dzy = (part[2]-SIGN(part[3])*ymax)/part[3];
            if (dzx>dzy) 
                dz = -dzx;
            else
                dz = -dzy;
            if (dz==-DBL_MAX)
                dz = 0;
            part[0] += dz*part[1];
            part[2] += dz*part[3];
            }
        if (is_out) {
          if (ip!=itop)
            swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z+dz;  /* record position of loss */
          coord[itop][5] = Po*(1+coord[itop][5]);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    return(np);
  }


long removeInvalidParticles(
    double **coord, long np, double **accepted,
    double z, double Po)
{
    long ip, itop, is_out, ic;
    double *part;

    itop = np-1;
    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        is_out = 0;
        for (ic=0; ic<6; ic++)
          if (isnan(part[ic])) {
            is_out = 1;
            break;
          }
        if (part[5]<=-1)
          is_out = 1;
        if (is_out) {
          swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z; 
          coord[itop][5] = Po*(1+coord[itop][5]);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    return(np);
  }

            
/* routine: elliptical_collimator()
 * purpose: eliminate particles that hit the walls of a non-zero length
 *          elliptical collimator
 *
 * Michael Borland, 1989
 */

long elliptical_collimator(
                           double **initial, ECOL *ecol, long np, double **accepted, double z,
                           double Po)
{
  double length, *ini;
  long ip, itop, lost, openCode;
  double a2, b2;
  double dx, dy, xo, yo, xsize, ysize;
  TRACKING_CONTEXT context;

  xsize = ecol->x_max;
  ysize = ecol->y_max;
  if (ecol->exponent<2 || ecol->exponent%2) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
            context.elementName, ecol->exponent);
    exit(1);
  }
  a2 = ipow(ecol->x_max, ecol->exponent);
  b2 = ipow(ecol->y_max, ecol->exponent);
  dx = ecol->dx;
  dy = ecol->dy;

  if (ecol->x_max<=0 || ecol->y_max<=0) {
    /* At least one of x_max or y_max is non-positive */
    if (ecol->x_max>0 || ecol->y_max>0) {
      /* One of x_max or y_max is positive. Use rectangular collimator routine to implement this. */
      RCOL rcol;
      rcol.length = ecol->length;
      rcol.x_max = ecol->x_max;
      rcol.y_max = ecol->y_max;
      rcol.dx = ecol->dx;
      rcol.dy = ecol->dy;
      rcol.openSide = ecol->openSide;
      return rectangular_collimator(initial, &rcol, np, accepted, z, Po);
    } 
    exactDrift(initial, np, ecol->length);
    return(np);
  }
  openCode = determineOpenSideCode(ecol->openSide);

  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    lost = 0;
    xo = ini[0] - dx;
    yo = ini[2] - dy;
    if ((ipow(xo, ecol->exponent)/a2 + ipow(yo, ecol->exponent)/b2)>1)
      lost = openCode ? evaluateLostWithOpenSides(openCode, xo, yo, xsize, ysize) : 1;
    else if (isinf(ini[0]) || isinf(ini[2]) ||
             isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
      initial[itop][4] = z;
      initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      np--;
    }
  }

  if (np==0 || (length=ecol->length)<=0)
    return(np);
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    lost = 0;
    ini[0] += length*ini[1];
    ini[2] += length*ini[3];
    xo = ini[0]/xsize;
    yo = ini[2]/ysize;
    if ((ipow(xo, ecol->exponent) + ipow(yo, ecol->exponent))>1)
      lost = openCode ? evaluateLostWithOpenSides(openCode, xo, yo, 1, 1) : 1;
    else if (isinf(ini[0]) || isinf(ini[2]) ||
             isnan(ini[0]) || isnan(ini[2]) )
      lost = 1;
    if (lost) {
      swapParticles(initial[ip], initial[itop]);
      initial[itop][4] = z + length;
      initial[itop][5] = sqrt(sqr(Po*(1+initial[itop][5]))+1);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      np--;
    }
    else 
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
  }

  return(np);
}



/* routine: elimit_amplitudes()
 * purpose: eliminate particles outside an ellipse with given semi-major
 *          and semi-minor axes.
 *
 * Michael Borland, 1989
 */

long elimit_amplitudes(
    double **coord,
    double xmax,       /* half-axis in x direction */ 
    double ymax,       /* half-axis in y direction */
    long np,
    double **accepted,
    double z,
    double Po,
    long extrapolate_z,
    long openCode,
    long exponent                       
    )
{
    long ip, itop, lost;
    double *part;
    double a2, b2, c1, c2, c0, dz, det;
    TRACKING_CONTEXT context;

    if (exponent<2 || exponent%2) {
      getTrackingContext(&context);
      fprintf(stderr, "Error for %s: exponent=%ld is not valid.  Give even integer >=2\n",
              context.elementName, exponent);
      exit(1);
    }
    if (xmax<=0 || ymax<=0) {
      /* At least one of the dimensions is non-positive and therefore ignored */
      if (xmax>0 || ymax>0) 
        return limit_amplitudes(coord, xmax, ymax, np, accepted, z, Po, extrapolate_z, openCode);
      return(np);
    }

    a2 = ipow(xmax, exponent);
    b2 = ipow(ymax, exponent);
    itop = np-1;

    for (ip=0; ip<np; ip++) {
        part = coord[ip];
        if (isinf(part[0]) || isinf(part[2]) || isnan(part[0]) || isnan(part[2]) ) {
            swapParticles(coord[ip], coord[itop]);
            coord[itop][4] = z;
            coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
            if (accepted)
                swapParticles(accepted[ip], accepted[itop]);
            --itop;
            --ip;
            np--;
            continue;
            }
        lost = 0;
        if ((ipow(part[0], exponent)/a2 + ipow(part[2], exponent)/b2) > 1)
          lost = 1;
        else
          continue;
        if (openCode)
          lost *= evaluateLostWithOpenSides(openCode, part[0], part[2], xmax, ymax);
        if (lost) {
          dz = 0;
          if (extrapolate_z && !openCode && exponent==2) {
            c0 = sqr(part[0])/a2 + sqr(part[2])/b2 - 1;
            c1 = 2*(part[0]*part[1]/a2 + part[2]*part[3]/b2);
            c2 = sqr(part[1])/a2 + sqr(part[3])/b2;
            det = sqr(c1)-4*c0*c2;
            if (z>0 && c2 && det>=0) {
              if ((dz = (-c1+sqrt(det))/(2*c2))>0)
                dz = (-c1-sqrt(det))/(2*c2);
              if ((z+dz)<0) 
                dz = -z;
              part[0] += dz*part[1];
              part[2] += dz*part[3];
            }
          }
          swapParticles(coord[ip], coord[itop]);
          coord[itop][4] = z+dz;
          coord[itop][5] = sqrt(sqr(Po*(1+coord[itop][5]))+1);
          if (accepted)
            swapParticles(accepted[ip], accepted[itop]);
          --itop;
          --ip;
          np--;
        }
      }
    log_exit("elimit_amplitudes");
    return(np);
  }

#define SQRT_3 1.7320508075688772
            
long beam_scraper(
                  double **initial, SCRAPER *scraper, long np, double **accepted, double z,
                  double Po
                  )
{
  double length, *ini;
  long do_x, do_y, ip, itop;
  double limit;

  log_entry("beam_scraper");

  if (scraper->direction<0 || scraper->direction>3)
    return np;
  
  if (scraper->direction==0 || scraper->direction==2) {
    do_x = scraper->direction==0 ? 1 : -1;
    do_y = 0;
    limit = scraper->position*do_x;
  }
  else {
    do_x = 0;
    do_y = scraper->direction==1 ? 1 : -1;
    limit = scraper->position*do_y;
  }    

  if (scraper->length && (scraper->Xo || scraper->Z)) {
    /* scraper has material properties that scatter beam and
     * absorb energy
     */
    MATTER matter;
    matter.length = scraper->length;
    matter.Xo = scraper->Xo;
    matter.elastic = scraper->elastic;
    matter.energyStraggle = scraper->energyStraggle;
    matter.Z = scraper->Z;
    matter.A = scraper->A;
    matter.rho = scraper->rho;
    matter.pLimit = scraper->pLimit;
    
    for (ip=0; ip<np; ip++) {
      ini = initial[ip];
      if ((do_x && do_x*(ini[0]-scraper->dx)>limit) ||
          (do_y && do_y*(ini[2]-scraper->dy)>limit)) {
        /* scatter and/or absorb energy */
        track_through_matter(&ini, 1, &matter, Po);
      } else {
        ini[0] = ini[0] + ini[1]*scraper->length;
        ini[2] = ini[2] + ini[3]*scraper->length;
        ini[4] += scraper->length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
      }
    }
    log_exit("beam_scraper");
    return(np);
  }

  /* come here for idealized scraper that just absorbs particles */
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    if ((do_x && do_x*(ini[0]-scraper->dx) > limit) ||
        (do_y && do_y*(ini[2]-scraper->dy) > limit) ) {
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }
  if (np==0 || (length=scraper->length)<=0) {
    log_exit("beam_scraper");
    return(np);
  }

  z += length;
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    ini[0] += length*ini[1];
    ini[2] += length*ini[3];
    if ((do_x && do_x*(ini[0]-scraper->dx) > limit) ||
        (do_y && do_y*(ini[2]-scraper->dy) > limit) ) {
      swapParticles(initial[ip], initial[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      --itop;
      --ip;
      --np;
    }
    else
      ini[4] += length*sqrt(1+sqr(ini[1])+sqr(ini[3]));
  }
  log_exit("beam_scraper");
  return(np);
}

long track_through_pfilter(
    double **initial, PFILTER *pfilter, long np, double **accepted, double z,
    double Po
    )
{
  long ip, itop;
  static double *deltaBuffer=NULL;
  static long maxBuffer = 0;
  double reference;
#ifdef USE_KAHAN
  double error = 0.0; 
#endif
  
  itop = np-1;

  if ((pfilter->lowerFraction || pfilter->upperFraction) && !pfilter->limitsFixed) {
    double level[2]={-1,-1}, limit[2], upper[2];
    long count = 0, i;
    if (maxBuffer<np &&
        !(deltaBuffer=SDDS_Realloc(deltaBuffer, sizeof(*deltaBuffer)*(maxBuffer=np))))
      SDDS_Bomb("memory allocation failure");  
    if (isSlave)
    for (ip=0; ip<np; ip++)
      deltaBuffer[ip] = initial[ip][5];
    /* eliminate lowest lowerfraction of particles and highest
       upperfraction */
    if (pfilter->lowerFraction>0 && pfilter->lowerFraction<1) {
      upper[count] = 0;
      level[count++] = pfilter->lowerFraction*100;
    }
    if (pfilter->upperFraction>0 && pfilter->upperFraction<1) {
      upper[count] = 1;
      level[count++] = 100-pfilter->upperFraction*100;
    }
#if SDDS_MPI_IO
    approximate_percentiles_p(limit, level, count, deltaBuffer, np, pfilter->bins); 
#else     
    compute_percentiles(limit, level, count, deltaBuffer, np);
#endif
    itop = np-1;
    for (i=0; i<2; i++) {
      if (level[i]<0)
        break;
      if (upper[i]) {
        pfilter->pUpper = (1+limit[i])*Po;
        pfilter->hasUpper = 1;
      }
      else {
        pfilter->pLower = (1+limit[i])*Po;
        pfilter->hasLower = 1;
      }
      if (!pfilter->fixPLimits) {
	/* filter in next block so there are no discrepancies due to
	 * small numerical differences
	 */
	if (isSlave)
	for (ip=0; ip<=itop; ip++) {
	  if ((upper[i] && initial[ip][5]>limit[i]) ||
	      (!upper[i] && initial[ip][5]<limit[i])) {
	    swapParticles(initial[ip], initial[itop]);
	    initial[itop][4] = z;  /* record position of particle loss */
	    initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
	    if (accepted)
	      swapParticles(accepted[ip], accepted[itop]);
	    --itop;
	    --ip;
	  }
	}
      }
    }
    if (pfilter->fixPLimits)
      pfilter->limitsFixed = 1;
  }
  
  if (pfilter->limitsFixed) {
    double p;
    if(isSlave)
    for (ip=0; ip<=itop; ip++) {
      p = (1+initial[ip][5])*Po;
      if ((pfilter->hasUpper && p>pfilter->pUpper) ||
	  (pfilter->hasLower && p<pfilter->pLower)) {
	swapParticles(initial[ip], initial[itop]);
	initial[itop][4] = z;  /* record position of particle loss */
	initial[itop][5] = p;  /* momentum at loss */
	if (accepted)
	  swapParticles(accepted[ip], accepted[itop]);
	--itop;
	--ip;
      }
    }
  }
  
  if (pfilter->deltaLimit<0) {
#if defined(MINIMIZE_MEMORY)
    free(deltaBuffer);
    deltaBuffer = NULL;
    maxBuffer= 0;
#endif
    return itop+1;
  }
  reference = 0.0;
  if (pfilter->beamCentered) {
    if (isSlave)
    for (ip=0; ip<=itop; ip++) {
#ifndef USE_KAHAN
      reference += initial[ip][5];
#else
      reference = KahanPlus(reference, initial[ip][5], &error);
#endif
    }
#if USE_MPI
    if (isMaster)
      itop = 0; 
    if (USE_MPI) {
      long itop_total;
      MPI_Allreduce(&itop, &itop_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      reference /= (itop_total+1);
    }
#else
    reference /= (itop+1);
#endif
  }
  if (isSlave)
  for (ip=0; ip<=itop; ip++) {
    if (fabs(initial[ip][5]-reference)<pfilter->deltaLimit)
      continue;
    swapParticles(initial[ip], initial[itop]);
    initial[itop][4] = z;  /* record position of particle loss */
    initial[itop][5] = Po*(1+initial[itop][5]);  /* momentum at loss */
    if (accepted)
      swapParticles(accepted[ip], accepted[itop]);
    --itop;
    --ip;
  }
#if defined(MINIMIZE_MEMORY)
  free(deltaBuffer);
  deltaBuffer = NULL;
  maxBuffer= 0;
#endif
  return itop+1;
}

long remove_outlier_particles(
                              double **initial, CLEAN *clean, long np, double **accepted, double z,
                              double Po
                              )
{
  double *ini, beta, p, *sSave;
  long ip, itop, is_out, j, mode;
  double limit[6], centroid[6], stDev[6];
  long count;
#define CLEAN_STDEV 0
#define CLEAN_ABSDEV 1
#define CLEAN_ABSVAL 2
#define CLEAN_MODES 3
  static char *modeName[CLEAN_MODES] = {"stdeviation", "absdeviation", "absvalue"};

  switch ((mode=match_string(clean->mode, modeName, 3, 0))) {
  case CLEAN_STDEV:
  case CLEAN_ABSDEV:
  case CLEAN_ABSVAL:
    break;
  default:
    fprintf(stderr, "Error: mode for CLEAN element must be one of the following:\n");
    for (j=0; j<CLEAN_MODES; j++)
      fprintf(stderr, "%s  ", modeName[j]);
    fprintf(stderr, "\n");
    exit(1);
    break;
  }
  
  if (!(clean->xLimit>0 || clean->xpLimit>0 ||
        clean->yLimit>0 || clean->ypLimit>0 ||
        clean->tLimit>0 || clean->deltaLimit>0)) {
    return(np);
  }

  /* copy limits into array */
  limit[0] = clean->xLimit;
  limit[1] = clean->xpLimit;
  limit[2] = clean->yLimit;
  limit[3] = clean->ypLimit;
  limit[4] = clean->tLimit;
  limit[5] = clean->deltaLimit;

  if (!(sSave = malloc(sizeof(*sSave)*np)))
    bomb("memory allocation failure", NULL);

  /* compute centroids for each coordinate */
  for (j=0; j<6; j++) {
    centroid[j] = stDev[j] = 0;
  }
  for (ip=count=0; ip<np; ip++) {
    ini = initial[ip];
    sSave[ip] = ini[4];
    for (j=0; j<6; j++)
      if (isnan(ini[j]) || isinf(ini[j]))
        break;
    if (j!=6)
      continue;
    /* compute time of flight and store in place of s */
    p = Po*(1+ini[5]);
    beta = p/sqrt(p*p+1);
    ini[4] /= beta*c_mks;
    count++;
    for (j=0; j<6; j++) {
      centroid[j] += ini[j];
    }
  }

#if USE_MPI
  if (isSlave) {
    double sum[6];
    long count_total;

    MPI_Allreduce(centroid,sum,6,MPI_DOUBLE,MPI_SUM,workers);
    memcpy(centroid,sum,sizeof(double)*6);
    MPI_Allreduce(&count,&count_total,1,MPI_LONG,MPI_SUM,workers);
    count = count_total;
  }
#endif

  if (!count) {
    for (ip=0; ip<np; ip++) {
      initial[ip][4] = z; /* record position of particle loss */
      initial[ip][5] = Po*(1+initial[ip][5]);
    }
    free(sSave);
    return 0;
  }
  for (j=0; j<6; j++)
    centroid[j] /= count;

  /* compute standard deviation of coordinates, if needed */
  if (mode==CLEAN_STDEV) {
    for (ip=0; ip<np; ip++) {
      ini = initial[ip];
      for (j=0; j<6; j++) {
        if (!isnan(ini[j]) && !isinf(ini[j]))
          stDev[j] += sqr(centroid[j] - ini[j]);
      }
    }
#if USE_MPI
  if (isSlave) {
    double sum[6];

    MPI_Allreduce(stDev,sum,6,MPI_DOUBLE,MPI_SUM,workers);
    memcpy(stDev,sum,sizeof(double)*6);
  }
#endif
    for (j=0; j<6; j++)
      stDev[j] = sqrt(stDev[j]/count);
  }
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    for (j=is_out=0; j<6; j++) {
      if (limit[j]<=0)
        continue;
      switch (mode) {
      case CLEAN_STDEV:
        if (fabs(ini[j]-centroid[j])/stDev[j]>limit[j])
          is_out = 1;
        break;
      case CLEAN_ABSDEV:
        if (fabs(ini[j]-centroid[j])>limit[j])
          is_out = 2;
        break;
      case CLEAN_ABSVAL:
        if (fabs(ini[j])>limit[j])
          is_out = 3;
        break;
      default:
        fprintf(stderr, "invalid mode in remove_outlier_particles---programming error!\n");
        exit(1);
        break;
      }
      if (is_out)
        break;
    }
    if (is_out) {
      swapParticles(initial[ip], initial[itop]);
      SWAP_DOUBLE(sSave[ip], sSave[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }

  /* restore s data to particle array */
  for (ip=0; ip<np; ip++)
    initial[ip][4] = sSave[ip];
  free(sSave);

  return(np);
}


long evaluateLostWithOpenSides(long code, double dx, double dy, double xsize, double ysize)
{
  long lost = 1;
  switch (code) {
  case OPEN_PLUS_X:
    if (dx>0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_MINUS_X:
    if (dx<0 && fabs(dy)<ysize)
      lost = 0;
    break;
  case OPEN_PLUS_Y:
    if (dy>0 && fabs(dx)<xsize)
      lost = 0;
    break;
          case OPEN_MINUS_Y:
    if (dy<0 && fabs(dx)<xsize)
      lost = 0;
    break;
  default:
    break;
  }
  return lost;
}

long determineOpenSideCode(char *openSide)
{
  TRACKING_CONTEXT context;
  long value;
  
  if (!openSide)
    return 0;
  if (strlen(openSide)==2 &&
      (openSide[0]=='+' || openSide[0]=='-')) {
    switch (openSide[1]) {
    case 'x':
    case 'X':
    case 'h':
    case 'H':
      value = (openSide[0]=='+' ? OPEN_PLUS_X : OPEN_MINUS_X);
      break;
    case 'y':
    case 'Y':
    case 'v':
    case 'V':
      value = (openSide[0]=='+' ? OPEN_PLUS_Y : OPEN_MINUS_Y);
      break;
    default:
      value = -1;
      break;
    }
  }
  if (value==-1) {
    getTrackingContext(&context);
    fprintf(stderr, "Error for %s: open_side=%s is not valid\n",
            context.elementName, openSide);
    exit(1);
  }
  return value;
}


/* routine: imposeApertureData()
 * purpose: eliminate particles with (x,y) larger than values given in aperture input file.
 *
 * Michael Borland, 2005
 */

#define DEBUG_APERTURE 0

long interpolateApertureData(double z, APERTURE_DATA *apData,
                             double *xCenter, double *yCenter, double *xSize, double *ySize)
{
  double z0, period;
  long iz;
  
  z0 = z;
  if (apData->s[apData->points-1]<z) {
    if (apData->periodic) {
      period = apData->s[apData->points-1]-apData->s[0];
      z0 -= period*((long) ((z0 - apData->s[0])/period));
    } else
      return 0;
  }
  if ((iz = binaryArraySearch(apData->s, sizeof(apData->s[0]), apData->points, &z0, 
                              double_cmpasc, 1))<0)
    return 0;
  if (iz==apData->points-1)
    iz -= 1;

  if (apData->s[iz]==apData->s[iz+1]) {
    *xCenter = apData->dx[iz];
    *yCenter = apData->dy[iz];
    *xSize = apData->xMax[iz];
    *ySize = apData->yMax[iz];
  }
  else {
    *xCenter = INTERPOLATE(apData->dx[iz], apData->dx[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *yCenter = INTERPOLATE(apData->dy[iz], apData->dy[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *xSize = INTERPOLATE(apData->xMax[iz], apData->xMax[iz+1], apData->s[iz], apData->s[iz+1], z0);
    *ySize = INTERPOLATE(apData->yMax[iz], apData->yMax[iz+1], apData->s[iz], apData->s[iz+1], z0);
  }
  return 1;
}


long imposeApertureData(
                        double **initial, long np, double **accepted,
                        double z, double Po, APERTURE_DATA *apData)
{
  long ip, itop, iz, lost;
  double *ini;
  double z0, period;
  double xSize, ySize;
  double xCenter, yCenter;
  double dx, dy;

#if DEBUG_APERTURE
  static FILE *fp = NULL;
  if (!fp) {
    TRACKING_CONTEXT tcontext;
    char s[1000];
    getTrackingContext(&tcontext);
    sprintf(s, "%s.aplos", tcontext.rootname);
    fflush(stdout);
    if (!(fp = fopen(s, "w")))
      bomb("unable to open debug file for aperture losses", NULL);
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&column name=z type=double units=m &end\n");
    fprintf(fp, "&column name=x type=double units=m &end\n");
    fprintf(fp, "&column name=y type=double units=m &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  }  
#endif

  if (!interpolateApertureData(z, apData, 
                               &xCenter, &yCenter, &xSize, &ySize))
    return np;
  
  itop = np-1;
  for (ip=0; ip<np; ip++) {
    ini = initial[ip];
    dx = ini[0] - xCenter;
    dy = ini[2] - yCenter;
    lost = 0;
    if ((xSize && fabs(dx) > xSize) ||
        (ySize && fabs(dy) > ySize))
      lost = 1;
    if (lost) {
#if DEBUG_APERTURE
      fprintf(fp, "%e %e %e\n", z, initial[ip][0], initial[ip][2]);
#endif
      swapParticles(initial[ip], initial[itop]);
      if (accepted)
        swapParticles(accepted[ip], accepted[itop]);
      initial[itop][4] = z; /* record position of particle loss */
      initial[itop][5] = Po*(1+initial[itop][5]);
      --itop;
      --ip;
      --np;
    }
  }
  
#if DEBUG_APERTURE
  fflush(fp);
#endif

  return(np);
}

