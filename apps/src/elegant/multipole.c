/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: multipole()
 * purpose: apply kicks due to an arbitrary multipole
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "track.h"

double *expansion_coefficients(long n);
void apply_canonical_multipole_kicks(double *qx, double *qy,
                                   double *sum_Fx, double *sum_Fy,
                                   double x, double y,
                                   long order, double KnL, long skew);
void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder);
void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult);

#define EXSQRT(value, order) (order==0?sqrt(value):(1+0.5*((value)-1)))

unsigned long multipoleKicksDone = 0;

#define ODD(j) ((j)%2)

void readErrorMultipoleData(MULTIPOLE_DATA *multData,
                               char *multFile, long steering)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  if (!multFile || !strlen(multFile)) {
    multData->orders = 0;
    multData->initialized = 0;
    return;
  }
  if (multData->initialized)
    return;
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multFile)) {
    fprintf(stdout, "Problem opening file %s\n", multFile);
    fflush(stdout);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "an", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      (!steering && SDDS_CheckColumn(&SDDSin, "bn", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) ||
      SDDS_CheckParameter(&SDDSin, "referenceRadius", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK) {
    fprintf(stdout, "Problems with data in multipole file %s\n", multFile);
    fflush(stdout);
    exitElegant(1);
  }
  if (steering && SDDS_CheckColumn(&SDDSin, "bn", NULL, SDDS_ANY_FLOATING_TYPE, NULL)==SDDS_CHECK_OK) {
    fprintf(stdout, "Warning: Steering multipole file %s should not have bn.\n",
            multFile);
    fprintf(stdout, "Use an to specify multipole content for a horizontal steerer.\n");
    fprintf(stdout, "Multipole content for vertical steerer is deduced from this.\n");
    fflush(stdout);
  }
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading multipole file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    fprintf(stdout, "Warning: no data in multipole file %s\n", multFile);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
  }
  if (!SDDS_GetParameterAsDouble(&SDDSin, "referenceRadius", &multData->referenceRadius) ||
      !(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->an=SDDS_GetColumnInDoubles(&SDDSin, "an")) || 
      (!steering && !(multData->bn=SDDS_GetColumnInDoubles(&SDDSin, "bn")))) {
    sprintf(buffer, "Unable to read multipole data for file %s\n", multFile);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }    
  if (steering &&
      !(multData->bn=SDDS_Malloc(sizeof(*(multData->bn))*multData->orders))) {
    fprintf(stdout, "Memory allocation failure (readErrorMultipoleData)\n");
    exitElegant(1);
  }
  if (SDDS_ReadPage(&SDDSin)==2) {
    fprintf(stdout, "Warning: multipole file %s has multiple pages, which are ignored\n",
            multFile);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  if (steering) {
    long i, j;
    /* check for disallowed multipoles */
    for (i=0; i<multData->orders; i++) {
      if (ODD(multData->order[i])) {
        fprintf(stdout, "Error: steering multipole file %s has disallowed odd orders.\n",
                multFile);
        exitElegant(1);
      }
    }
    /* normalize to n=0 if present */
    /* find i such that order[i] is 0 (dipole) */
    for (i=0; i<multData->orders; i++) {
      if (multData->order[i]==0)
        break;
    }
    if (multData->orders>1 && i!=multData->orders) {
      /* dipole present */
      if (!multData->an[i]) {
        fprintf(stdout, "Steering multipole data in %s is invalid: an is zero for order=0\n",
                multFile);
        exitElegant(1);
      }
      /* normalize to dipole for normal and skew separately */
      for (j=0; j<multData->orders; j++)
        if (j!=i)
          multData->an[j] /= multData->an[i];
      /* remove the dipole data */
      for (j=i+1; j<multData->orders; j++) {
        multData->an[j-1] = multData->an[j];
        multData->order[j-1] = multData->order[j];
      }
      multData->orders -= 1;
      for (i=0; i<multData->orders; i++)
        multData->bn[i] = multData->an[i]*ipow(-1.0, multData->order[i]/2);
#ifdef DEBUG
      fprintf(stdout, "Steering multipole data: \n");
      for (i=0; i<multData->orders; i++)
        fprintf(stdout, "%ld: %e %e\n", multData->order[i],
                multData->an[i], multData->bn[i]);
#endif
    }
  }
  multData->initialized = 1;
}

void initialize_fmultipole(FMULT *multipole)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  MULTIPOLE_DATA *multData;
  
  multData = &(multipole->multData);
  if (multData->initialized)
    return;
  if (!multipole->filename)
    bombElegant("FMULT element doesn't have filename", NULL);
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, multipole->filename)) {
    sprintf(buffer, "Problem opening file %s (FMULT)\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if (SDDS_CheckColumn(&SDDSin, "order", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "KnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "JnL", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK)
    bombElegant("problems with data in FMULT input file", NULL);
  if (SDDS_ReadPage(&SDDSin)!=1)  {
    sprintf(buffer, "Problem reading FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if ((multData->orders = SDDS_RowCount(&SDDSin))<=0) {
    fprintf(stdout, "Warning: no data in FMULT file %s\n", multipole->filename);
    fflush(stdout);
    SDDS_Terminate(&SDDSin);
    return;
  }
  multData->JnL = NULL;
  if (!(multData->order=SDDS_GetColumnInLong(&SDDSin, "order")) ||
      !(multData->KnL=SDDS_GetColumnInDoubles(&SDDSin, "KnL")) ||
      !(multData->JnL=SDDS_GetColumnInDoubles(&SDDSin, "JnL"))) {
    sprintf(buffer, "Unable to read data for FMULT file %s\n", multipole->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }    
  if (SDDS_ReadPage(&SDDSin)==2) {
    fprintf(stdout, "Warning: FMULT file %s has multiple pages, which are ignored\n",
            multipole->filename);
    fflush(stdout);
  }
  SDDS_Terminate(&SDDSin);
  multData->initialized = 1;
}


long fmultipole_tracking(
                         double **particle,  /* initial/final phase-space coordinates */
                         long n_part,        /* number of particles */
                         FMULT *multipole,   /* multipole structure */
                         double p_error,     /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start
                         )
{
  double dummy;
  long n_kicks;       /* number of kicks to split multipole into */
  long i_part, i_top, is_lost=0, i_order;
  double *coord;
  double drift;
  double x=0.0, xp=0.0, y=0.0, yp=0.0;
  double rad_coef;
  MULTIPOLE_DATA multData;
  
  if (!particle)
    bombElegant("particle array is null (fmultipole_tracking)", NULL);

  if (!multipole)
    bombElegant("null MULT pointer (fmultipole_tracking)", NULL);
  
  if (!multipole->multData.initialized)
    initialize_fmultipole(multipole);
  
  if ((n_kicks=multipole->n_kicks)<=0)
    bombElegant("n_kicks<=0 in multipole()", NULL);

  drift = multipole->length;

  if (multipole->synch_rad)
    rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
  else
    rad_coef = 0;

  multData = multipole->multData;
  for (i_order=0; i_order<multData.orders; i_order++) {
    multData.KnL[i_order] *= (1+multipole->fse);
    if (multData.JnL)
      multData.JnL[i_order] *= (1+multipole->fse);
  }
  
  if (multipole->dx || multipole->dy || multipole->dz)
    offsetBeamCoordinates(particle, n_part, multipole->dx, multipole->dy, multipole->dz);
  if (multipole->tilt)
    rotateBeamCoordinates(particle, n_part, multipole->tilt);

  i_top = n_part-1;
  multipoleKicksDone += (i_top+1)*multData.orders*n_kicks*4;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      fprintf(stdout, "null coordinate pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "null accepted coordinates pointer for particle %ld (fmultipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    if (!integrate_kick_multipole_ord4(coord, multipole->dx, multipole->dy, 0.0, 0.0, Po, rad_coef, 0.0,
                                       1, multipole->sqrtOrder, 0.0, n_kicks, drift, &multData, NULL, NULL,
                                       &dummy, NULL)) {
      is_lost = 1;
      break;
    }
    
    if (!is_lost) {
      x = coord[0];
      y = coord[2];
      xp = coord[1];
      yp = coord[3];

#if defined(IEEE_MATH)
      if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
        swapParticles(particle[i_part], particle[i_top]);
        if (accepted)
          swapParticles(accepted[i_part], accepted[i_top]);
        particle[i_top][4] = z_start;
        particle[i_top][5] = Po*(1+particle[i_top][5]);
        i_top--;
        i_part--;
        continue;
      }
#endif
    }
    if (is_lost || FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
        FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
  }

  if (multipole->tilt)
    rotateBeamCoordinates(particle, n_part, -multipole->tilt);
  if (multipole->dx || multipole->dy || multipole->dz)
    offsetBeamCoordinates(particle, n_part, -multipole->dx, -multipole->dy, -multipole->dz);

  log_exit("fmultipole_tracking");
  return(i_top+1);
}


long multipole_tracking(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    MULT *multipole,    /* multipole structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
    double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
    long order;         /* order (n) */
    long n_kicks;       /* number of kicks to split multipole into */
    long i_part, i_kick, i, i_top, is_lost;
    double sum_Fx, sum_Fy, xypow, denom, qx, qy;
    double *coord;
    double drift;
    double *coef;
    double x, xp, y, yp, s, dp;
    double ratio, rad_coef;
    double beta0, beta1, p;

    log_entry("multipole_tracking");

    if (!particle)
        bombElegant("particle array is null (multipole_tracking)", NULL);

    if (!multipole)
        bombElegant("null MULT pointer (multipole_tracking)", NULL);

    if ((n_kicks=multipole->n_kicks)<=0) {
      TRACKING_CONTEXT tc;
      getTrackingContext(&tc);
      fprintf(stderr, "error: n_kicks<=0 in multipole() for %s #%ld\n", tc.elementName, tc.elementOccurrence);
      exitElegant(1);
    }

    if ((order=multipole->order)<0)
        bombElegant("order < 0 in multipole()", NULL);

    if (!(coef = expansion_coefficients(order)))
        bombElegant("expansion_coefficients() returned null pointer (multipole_tracking)", NULL);

    drift = multipole->length/n_kicks/2;
    if (multipole->bore)
        /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) */
        KnL = dfactorial(multipole->order)*multipole->BTipL/ipow(multipole->bore, multipole->order)*
              (particleCharge/(particleMass*c_mks*Po))*multipole->factor;
    else
      KnL = multipole->KnL*multipole->factor/n_kicks;

    if (multipole->synch_rad)
        rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    else
        rad_coef = 0;

    if (multipole->dx || multipole->dy || multipole->dz)
      offsetBeamCoordinates(particle, n_part, multipole->dx, multipole->dy, multipole->dz);
    if (multipole->tilt)
      rotateBeamCoordinates(particle, n_part, multipole->tilt);

    i_top = n_part-1;
    multipoleKicksDone += (i_top+1)*n_kicks*4;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stdout, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stdout, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
            fflush(stdout);
            abort();
            }
        if (KnL==0) {
            coord[4] += multipole->length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            coord[0] += multipole->length*coord[1];
            coord[2] += multipole->length*coord[3];
            continue;
            }

        x = coord[0];
        xp = coord[1];
        y = coord[2];
        yp = coord[3];
        s  = 0;
        dp = coord[5];
        p = Po*(1+dp);
        beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

        /* calculate initial canonical momenta */
        qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
        qy = (1+dp)*yp/denom;
        is_lost = 0;
        for (i_kick=0; i_kick<n_kicks; i_kick++) {
            if (drift) {
                x += xp*drift*(i_kick?2:1);
                y += yp*drift*(i_kick?2:1);
                s += (i_kick?2:1)*drift*sqrt(1+sqr(xp)+sqr(yp));
                }
            if (x==0) {
                xypow = ipow(y, order);
                ratio = 0;
                i = order;
                }
            else {
                xypow = ipow(x, order);
                ratio = y/x;
                i = 0;
                }
            /* now sum up the terms for the multipole expansion */
            for (sum_Fx=sum_Fy=0; i<=order; i++) {
                if (ODD(i))
                    sum_Fx += coef[i]*xypow;
                else
                    sum_Fy += coef[i]*xypow;
                xypow *= ratio;
                }
            /* apply kicks canonically */
            qx -= KnL*sum_Fy;
            qy += KnL*sum_Fx;
            if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
                is_lost = 1;
                break;
                }
            xp = qx/(denom=sqrt(denom));
            yp = qy/denom;
            if (rad_coef && drift) {
                qx /= (1+dp);
                qy /= (1+dp);
                dp -= rad_coef*sqr(KnL*(1+dp))*(sqr(sum_Fy)+sqr(sum_Fx))*sqrt(1+sqr(xp)+sqr(yp))/(2*drift);
                qx *= (1+dp);
                qy *= (1+dp);
                }
            }
        if (drift && !is_lost) {
            /* go through final drift */
            x += xp*drift;
            y += yp*drift;
            s += drift*sqrt(1 + sqr(xp) + sqr(yp));
            }

        coord[0] = x;
        coord[1] = xp;
        coord[2] = y;
        coord[3] = yp;
        coord[5] = dp;

        if (rad_coef) {
            p = Po*(1+dp);
            beta1 = p/sqrt(sqr(p)+1);
            coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
            }
        else 
            coord[4] += s;

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            swapParticles(particle[i_part], particle[i_top]);
            if (accepted)
                swapParticles(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT || is_lost) {
          swapParticles(particle[i_part], particle[i_top]);
          if (accepted)
            swapParticles(accepted[i_part], accepted[i_top]);
          particle[i_top][4] = z_start;
          particle[i_top][5] = Po*(1+particle[i_top][5]);
          i_top--;
          i_part--;
          continue;
        }

      }
    
    if (multipole->tilt)
      rotateBeamCoordinates(particle, n_part, -multipole->tilt);
    if (multipole->dx || multipole->dy || multipole->dz)
      offsetBeamCoordinates(particle, n_part, -multipole->dx, -multipole->dy, -multipole->dz);
    
    log_exit("multipole_tracking");
    return(i_top+1);
    }



double *expansion_coefficients(long n)
{
  static double **expansion_coef=NULL;
  static long *orderDone=NULL;
  static long maxOrder = -1;
  long i;
  
  if (n<=maxOrder && orderDone[n]) 
    return(expansion_coef[n]);

  if (n>maxOrder) {
    expansion_coef = trealloc(expansion_coef, sizeof(*expansion_coef)*(n+1));
    orderDone      = trealloc(orderDone, sizeof(*orderDone)*(n+1));
    for (i=maxOrder+1; i<=n; i++)
      orderDone[i] = 0;
    maxOrder = n;
  }

  expansion_coef[n] = tmalloc(sizeof(**expansion_coef)*(n+1));
  
  /* calculate expansion coefficients with signs for (x+iy)^n/n! */
  for (i=0; i<=n; i++) {
    expansion_coef[n][i] = (ODD(i/2)?-1.0:1.0)/(dfactorial(i)*dfactorial(n-i));
  }             
  orderDone[n] = 1;
  
  return(expansion_coef[n]);
}

long multipole_tracking2(
                         double **particle,   /* initial/final phase-space coordinates */
                         long n_part,         /* number of particles */
                         ELEMENT_LIST *elem,  /* element pointer */
                         double p_error,      /* p_nominal/p_central */
                         double Po,
                         double **accepted,
                         double z_start,
                         /* from MAXAMP element */
                         double x_max,
                         double y_max,
                         long elliptical,
                         /* from aperture_data command */
                         APERTURE_DATA *apFileData,
                         /* For return of accumulated change in sigmaDelta^2 */
                         double *sigmaDelta2
                         )
{
  double KnL;         /* integrated strength = L/(B.rho)*(Dx^n(By))_o for central momentum */
  double dx, dy, dz;  /* offsets of the multipole center */
  long order;         /* order (n) */
  long n_kicks, integ_order;
  long i_part, i_top, n_parts;
  double *coef, *coord;
  double drift;
  double tilt, rad_coef, isr_coef, xkick, ykick, dzLoss;
  KQUAD *kquad;
  KSEXT *ksext;
  KQUSE *kquse;
  MULTIPOLE_DATA *multData = NULL, *steeringMultData = NULL;
  long sqrtOrder;
  MULT_APERTURE_DATA apertureData;
  double K2L;
  
  log_entry("multipole_tracking2");

  if (!particle)
    bombElegant("particle array is null (multipole_tracking)", NULL);

  if (!elem)
    bombElegant("null element pointer (multipole_tracking2)", NULL);
  if (!elem->p_elem)
    bombElegant("null p_elem pointer (multipole_tracking2)", NULL);

  rad_coef = xkick = ykick = isr_coef = 0;
  sqrtOrder = 0;

  switch (elem->type) {
  case T_KQUAD:
    kquad = ((KQUAD*)elem->p_elem);
    n_kicks = kquad->n_kicks;
    order = 1;
    if (kquad->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = kquad->B/kquad->bore*(particleCharge/(particleMass*c_mks*Po))*kquad->length*(1+kquad->fse);
    else
      KnL = kquad->k1*kquad->length*(1+kquad->fse);
    drift = kquad->length;
    tilt = kquad->tilt;
    dx = kquad->dx;
    dy = kquad->dy;
    dz = kquad->dz;
    xkick = kquad->xkick;
    ykick = kquad->ykick;
    integ_order = kquad->integration_order;
    sqrtOrder = kquad->sqrtOrder?1:0;
    if (kquad->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!kquad->isr || (kquad->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (!kquad->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(kquad->systematicMultipoleData),
                             kquad->systematic_multipoles, 0);
      readErrorMultipoleData(&(kquad->randomMultipoleData),
                             kquad->random_multipoles, 0);
      readErrorMultipoleData(&(kquad->steeringMultipoleData), 
                             kquad->steering_multipoles, 1);
      kquad->multipolesInitialized = 1;
    }
    computeTotalErrorMultipoleFields(&(kquad->totalMultipoleData),
                                     &(kquad->systematicMultipoleData),
                                     &(kquad->randomMultipoleData),
                                     &(kquad->steeringMultipoleData),
                                     KnL, 1);
    multData = &(kquad->totalMultipoleData);
    steeringMultData = &(kquad->steeringMultipoleData);
    break;
  case T_KSEXT:
    ksext = ((KSEXT*)elem->p_elem);
    n_kicks = ksext->n_kicks;
    order = 2;
    if (ksext->bore)
      /* KnL = d^nB/dx^n * L/(B.rho) = n! B(a)/a^n * L/(B.rho) * (1+FSE) */
      KnL = 2*ksext->B/sqr(ksext->bore)*(particleCharge/(particleMass*c_mks*Po))*ksext->length*(1+ksext->fse);
    else
      KnL = ksext->k2*ksext->length*(1+ksext->fse);
    drift = ksext->length;
    tilt = ksext->tilt;
    dx = ksext->dx;
    dy = ksext->dy;
    dz = ksext->dz;
    integ_order = ksext->integration_order;
    sqrtOrder = ksext->sqrtOrder?1:0;
    if (ksext->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass);
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!ksext->isr || (ksext->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    if (!ksext->multipolesInitialized) {
      /* read the data files for the error multipoles */
      readErrorMultipoleData(&(ksext->systematicMultipoleData),
                             ksext->systematic_multipoles, 0);
      readErrorMultipoleData(&(ksext->randomMultipoleData),
                             ksext->random_multipoles, 0);
      ksext->multipolesInitialized = 1;
    }
    computeTotalErrorMultipoleFields(&(ksext->totalMultipoleData),
                                     &(ksext->systematicMultipoleData),
                                     &(ksext->randomMultipoleData),
                                     NULL,
                                     KnL, 2);
    multData = &(ksext->totalMultipoleData);
    break;
  case T_KQUSE:
    /* Implemented as a quadrupole with sextupole as a secondary multipole */
    kquse = ((KQUSE*)elem->p_elem);
    n_kicks = kquse->n_kicks;
    order = 1;
    KnL = kquse->k1*kquse->length*(1+kquse->fse1);
    drift = kquse->length;
    tilt = kquse->tilt;
    dx = kquse->dx;
    dy = kquse->dy;
    dz = kquse->dz;
    integ_order = kquse->integration_order;
    sqrtOrder = 0;
    if (kquse->synch_rad)
      rad_coef = sqr(particleCharge)*pow3(Po)/(6*PI*epsilon_o*sqr(c_mks)*particleMass); 
    isr_coef = particleRadius*sqrt(55.0/(24*sqrt(3))*pow5(Po)*137.0359895);
    if (!kquse->isr || (kquse->isr1Particle==0 && n_part==1))
      /* Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles */
      isr_coef *= -1;
    K2L = kquse->k2*kquse->length*(1+kquse->fse2);
    if (K2L) {
      multData = tmalloc(sizeof(*multData));
      multData->orders = multData->initialized = 1;
      multData->randomized = 0;
      multData->order = tmalloc(sizeof(*(multData->order))*1);
      multData->order[0] = 2;
      multData->KnL = tmalloc(sizeof(*(multData->KnL))*1);
      multData->KnL[0] = K2L;
      multData->JnL = NULL;
    }
    break;
  default:
    fprintf(stdout, "error: multipole_tracking2() called for element %s--not supported!\n", elem->name);
    fflush(stdout);
    exitElegant(1);
    break;
  }
  if (multData && !multData->initialized)
    multData = NULL;
  
  if (n_kicks<=0)
    bombElegant("n_kicks<=0 in multipole()", NULL);
  if (order<=0)
    bombElegant("order <= 0 in multipole()", NULL);
  if (integ_order!=2 && integ_order!=4) 
    bombElegant("multipole integration_order must be 2 or 4", NULL);
  
  if (!(coef = expansion_coefficients(order)))
    bombElegant("expansion_coefficients() returned null pointer (multipole_tracking)", NULL);

  i_top = n_part-1;
  if (integ_order==4) {
    if ((n_parts = ceil(n_kicks/4.0))<1)
      n_parts = 1;
    n_kicks = n_parts*4;
  }
  else
    n_parts = n_kicks;
  
  multipoleKicksDone += (i_top+1)*n_kicks;
  if (multData)
    multipoleKicksDone += (i_top+1)*n_kicks*multData->orders;

  apertureData.xCen = apertureData.yCen = 0;
  apertureData.xMax = x_max;
  apertureData.yMax = y_max;
  apertureData.elliptical = elliptical;
  apertureData.present = x_max>0 || y_max>0;
  if (apFileData && apFileData->initialized) {
    /* If there is file-based aperture data, it may override MAXAMP data. */
    double xCenF, yCenF, xMaxF, yMaxF;
    apertureData.present = 1;
    if (interpolateApertureData(z_start+drift/2, apFileData, 
                                &xCenF, &yCenF, &xMaxF, &yMaxF)) {
      if (x_max<=0 || (x_max>fabs(xCenF+xMaxF) && x_max>fabs(xCenF-xMaxF))) {
        apertureData.xMax = xMaxF;
        apertureData.xCen = xCenF;
        apertureData.elliptical = 0;
      }
      if (y_max<=0 || (y_max>fabs(yCenF+yMaxF) && y_max>fabs(yCenF-yMaxF))) {
        apertureData.yMax = yMaxF;
        apertureData.yCen = yCenF;
        apertureData.elliptical = 0;
      }
    }
  }
  
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, dx, dy, dz);
  if (tilt)
    rotateBeamCoordinates(particle, n_part, tilt);

  if (sigmaDelta2)
    *sigmaDelta2 = 0;
  for (i_part=0; i_part<=i_top; i_part++) {
    if (!(coord = particle[i_part])) {
      fprintf(stdout, "null coordinate pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }
    if (accepted && !accepted[i_part]) {
      fprintf(stdout, "null accepted coordinates pointer for particle %ld (multipole_tracking)", i_part);
      fflush(stdout);
      abort();
    }

    if ((integ_order==4 &&
         !integrate_kick_multipole_ord4(coord, dx, dy, xkick, ykick,
                                        Po, rad_coef, isr_coef, order, sqrtOrder, KnL,
                                        n_parts, drift, 
                                        multData, steeringMultData,
                                        &apertureData, &dzLoss, sigmaDelta2)) ||
        (integ_order==2 &&
         !integrate_kick_multipole_ord2(coord, dx, dy, xkick, ykick,
                                        Po, rad_coef, isr_coef, order, sqrtOrder, KnL, 
                                        n_parts, drift,
                                        multData, steeringMultData,
                                        &apertureData, &dzLoss, sigmaDelta2))) {
      swapParticles(particle[i_part], particle[i_top]);
      if (accepted)
        swapParticles(accepted[i_part], accepted[i_top]);
      particle[i_top][4] = z_start+dzLoss;
      particle[i_top][5] = Po*(1+particle[i_top][5]);
      i_top--;
      i_part--;
      continue;
    }
  }
  if (sigmaDelta2)
    *sigmaDelta2 /= i_top+1;
  
  if (tilt)
    rotateBeamCoordinates(particle, n_part, -tilt);
  if (dx || dy || dz)
    offsetBeamCoordinates(particle, n_part, -dx, -dy, -dz);

  log_exit("multipole_tracking2");
  return(i_top+1);
}

int integrate_kick_multipole_ord2(double *coord, double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef, double isr_coef,
                                  long order, long sqrtOrder, double KnL, long n_kicks, double drift,
                                  MULTIPOLE_DATA *multData, 
                                  MULTIPOLE_DATA *steeringMultData,
                                  MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, imult;
  
  drift = drift/n_kicks/2.0;
  KnL = KnL/n_kicks;
  
  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  denom = 1+sqr(xp)+sqr(yp);
  denom = EXSQRT(denom, sqrtOrder);
  qx = (1+dp)*xp/denom;
  qy = (1+dp)*yp/denom;

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
  }

  *dzLoss = 0;
  for (i_kick=0; i_kick<n_kicks; i_kick++) {
    if (drift) {
      x += xp*drift*(i_kick?2:1);
      y += yp*drift*(i_kick?2:1);
      s += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
      *dzLoss += drift*(i_kick?2:1)*sqrt(1 + sqr(xp) + sqr(yp));
    }
    if (apData && apData->present &&
        ((apData->xMax && fabs(x + dx - apData->xCen)>apData->xMax) ||
         (apData->yMax && fabs(y + dy - apData->yCen)>apData->yMax) )) {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    
    apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, order, KnL, 0);

    /* do kicks for spurious multipoles */
    if (multData) {
      for (imult=0; imult<multData->orders; multData++) {
        if (multData->KnL && multData->KnL[imult]) 
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->KnL[imult]/n_kicks, 0);
        if (multData->JnL && multData->JnL[imult]) 
          apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                          multData->order[imult], 
                                          multData->JnL[imult]/n_kicks, 1);
      }
    }

    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
    if ((rad_coef || isr_coef) && drift) {
      double deltaFactor, F2, dsFactor;
      deltaFactor = sqr(1+dp);
      F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/(2*drift));
      dsFactor = EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder)*2*drift;
      qx /= (1+dp);
      qy /= (1+dp);
      if (rad_coef)
	dp -= rad_coef*deltaFactor*F2*dsFactor;
      if (isr_coef>0)
	dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsFactor)*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
      if (sigmaDelta2)
        *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
      qx *= (1+dp);
      qy *= (1+dp);
    }
  }
  if (drift) {
    /* go through final drift */
    x += xp*drift;
    y += yp*drift;
    s += drift*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
    *dzLoss += drift*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
  }
  if (apData && apData->present &&
      ((apData->xMax && fabs(x + dx - apData->xCen)>apData->xMax) ||
       (apData->yMax && fabs(y + dy - apData->yCen)>apData->yMax) )) {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->KnL[imult]*xkick/2, 0);
      apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                      steeringMultData->order[imult], 
                                      steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    denom = EXSQRT(denom, sqrtOrder);
    xp = qx/denom;
    yp = qy/denom;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef || isr_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}


/* BETA is 2^(1/3) */
#define BETA 1.25992104989487316477

int integrate_kick_multipole_ord4(double *coord, double dx, double dy, double xkick, double ykick,
                                  double Po, double rad_coef, double isr_coef,
                                  long order, long sqrtOrder, double KnL, long n_parts, double drift,
                                  MULTIPOLE_DATA *multData, MULTIPOLE_DATA *steeringMultData,
                                  MULT_APERTURE_DATA *apData, double *dzLoss, double *sigmaDelta2) 
{
  double p, qx, qy, denom, beta0, beta1, dp, s;
  double x, y, xp, yp, sum_Fx, sum_Fy;
  long i_kick, step, imult;
  double dsh;
  static double driftFrac[4] = {
    0.5/(2-BETA),  (1-BETA)/(2-BETA)/2,  (1-BETA)/(2-BETA)/2,  0.5/(2-BETA)
    } ;
  static double kickFrac[4] = {
    1./(2-BETA),  -BETA/(2-BETA),  1/(2-BETA),  0
    } ;
  
  drift = drift/n_parts;
  KnL = KnL/n_parts;

  x = coord[0];
  xp = coord[1];
  y = coord[2];
  yp = coord[3];
  s  = 0;
  dp = coord[5];
  p = Po*(1+dp);
  beta0 = p/sqrt(sqr(p)+1);

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  /* calculate initial canonical momenta */
  qx = (1+dp)*xp/(denom=EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder));
  qy = (1+dp)*yp/denom;

  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult])
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      return 0;
    }
    xp = qx/(denom=EXSQRT(denom, sqrtOrder));
    yp = qy/denom;
  }

  *dzLoss = 0;
  for (i_kick=0; i_kick<n_parts; i_kick++) {
    for (step=0; step<4; step++) {
      if (drift) {
        dsh = drift*driftFrac[step];
        x += xp*dsh;
        y += yp*dsh;
        s += dsh*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
        *dzLoss += dsh*EXSQRT(1 + sqr(xp) + sqr(yp), sqrtOrder);
      }
      if (apData && apData->present &&
          ((apData->xMax>=0 && fabs(x + dx - apData->xCen)>apData->xMax) ||
           (apData->yMax>=0 && fabs(y + dy - apData->yCen)>apData->yMax) )) {
        coord[0] = x;
        coord[2] = y;
        return 0;
      }

      if (!kickFrac[step])
        break;
      apply_canonical_multipole_kicks(&qx, &qy, &sum_Fx, &sum_Fy, x, y, 
                                      order, KnL*kickFrac[step], 0);
      if (multData) {
        /* do kicks for spurious multipoles */
        for (imult=0; imult<multData->orders; imult++) {
          if (multData->KnL && multData->KnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->KnL[imult]*kickFrac[step]/n_parts,
                                            0);
          }
          if (multData->JnL && multData->JnL[imult]) {
            apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                            multData->order[imult], 
                                            multData->JnL[imult]*kickFrac[step]/n_parts,
                                            1);
          }
        }
      }
      if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
        return 0;
      }
      xp = qx/(denom=EXSQRT(denom, sqrtOrder));
      yp = qy/denom;
      if ((rad_coef || isr_coef) && drift) {
	double deltaFactor, F2, dsFactor, dsISRFactor;
        qx /= (1+dp);
        qy /= (1+dp);
	deltaFactor = sqr(1+dp);
	F2 = (sqr(sum_Fy)+sqr(sum_Fx))*sqr(KnL/drift);
	dsFactor = EXSQRT(1+sqr(xp)+sqr(yp), sqrtOrder);
	dsISRFactor = dsFactor*drift/3;   /* recall that kickFrac may be negative */
	dsFactor *= drift*kickFrac[step]; /* that's ok here, since we don't take sqrt */
	if (rad_coef)
	  dp -= rad_coef*deltaFactor*F2*dsFactor;
	if (isr_coef>0)
	  dp -= isr_coef*deltaFactor*pow(F2, 0.75)*sqrt(dsISRFactor)*gauss_rn_lim(0.0, 1.0, 3.0, random_2);
        if (sigmaDelta2)
          *sigmaDelta2 += sqr(isr_coef*deltaFactor)*pow(F2, 1.5)*dsFactor;
        qx *= (1+dp);
        qy *= (1+dp);
      }
    }
  }
  
  if (apData && apData->present &&
      ((apData->xMax && fabs(x + dx - apData->xCen)>apData->xMax) ||
       (apData->yMax && fabs(y + dy - apData->yCen)>apData->yMax) ))  {
    coord[0] = x;
    coord[2] = y;
    return 0;
  }
  
  if (steeringMultData) {
    /* apply steering corrector multipoles */
    for (imult=0; imult<steeringMultData->orders; imult++) {
      if (steeringMultData->KnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->KnL[imult]*xkick/2, 0);
      if (steeringMultData->JnL[imult]) 
        apply_canonical_multipole_kicks(&qx, &qy, NULL, NULL, x, y, 
                                        steeringMultData->order[imult], 
                                        steeringMultData->JnL[imult]*ykick/2, 1);
    }
    if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
      coord[0] = x;
      coord[2] = y;
      return 0;
    }
    xp = qx/(denom=EXSQRT(denom, sqrtOrder));
    yp = qy/denom;
  }

  /* apply steering corrector kick */
  xp += xkick/(1+dp)/2;
  yp += ykick/(1+dp)/2;

  coord[0] = x;
  coord[1] = xp;
  coord[2] = y;
  coord[3] = yp;
  if (rad_coef) {
    p = Po*(1+dp);
    beta1 = p/sqrt(sqr(p)+1);
    coord[4] = beta1*(coord[4]/beta0 + 2*s/(beta0+beta1));
  }
  else 
    coord[4] += s;
  coord[5] = dp;

#if defined(IEEE_MATH)
  if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
    return 0;
  }
#endif
  if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
      FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
    return 0;
  }
  return 1;
}

void apply_canonical_multipole_kicks(double *qx, double *qy, 
                                     double *sum_Fx_return, double *sum_Fy_return,
                                     double x, double y,
                                     long order, double KnL, long skew)
{
  long i;
  double sum_Fx, sum_Fy, xypow, ratio;
  double *coef;
  if (sum_Fx_return)
    *sum_Fx_return = 0;
  if (sum_Fy_return)
    *sum_Fy_return = 0;
  coef = expansion_coefficients(order);
  if (x==0) {
    if (y==0)
      return;
    xypow = ipow(y, order);
    i = order;
    ratio = 0;
  }
  else {
    xypow = ipow(x, order);
    ratio = y/x;
    i = 0;
  }
  /* now sum up the terms for the multipole expansion */
  for (sum_Fx=sum_Fy=0; i<=order; i++) {
    if (ODD(i))
      sum_Fx += coef[i]*xypow;
    else
      sum_Fy += coef[i]*xypow;
    xypow *= ratio;
  }
  if (skew) {
    SWAP_DOUBLE(sum_Fx, sum_Fy);
    sum_Fx = -sum_Fx;
  }
  /* add the kicks */
  *qx -= KnL*sum_Fy;
  *qy += KnL*sum_Fx;
  if (sum_Fx_return)
    *sum_Fx_return = sum_Fx;
  if (sum_Fy_return)
    *sum_Fy_return = sum_Fy;
}

void randomizeErrorMultipoleFields(MULTIPOLE_DATA *randomMult)
{
  long i;
  double nFactorial, rpow, rn1, rn2;

  if (!randomMult || randomMult->randomized)
    return;
  for (i=0; i<randomMult->orders; i++) {
    nFactorial = dfactorial(randomMult->order[i]);
    rpow = ipow(randomMult->referenceRadius, randomMult->order[i]);
    rn1 = gauss_rn_lim(0.0, 1.0, 2.0, random_1_elegant);
    rn2 = gauss_rn_lim(0.0, 1.0, 2.0, random_1_elegant);
    randomMult->anMod[i] = randomMult->an[i]*nFactorial*rn1/rpow;
    randomMult->bnMod[i] = randomMult->bn[i]*nFactorial*rn2/rpow;
  }
  randomMult->randomized = 1;
}

void computeTotalErrorMultipoleFields(MULTIPOLE_DATA *totalMult,
                                      MULTIPOLE_DATA *systematicMult,
                                      MULTIPOLE_DATA *randomMult,
                                      MULTIPOLE_DATA *steeringMult,
                                      double KmL, long rootOrder)
{
  long i;
  double sFactor=0.0, rFactor=0.0;

  if (!totalMult->initialized) {
    totalMult->initialized = 1;
    if (steeringMult && steeringMult->orders) {
      if (!(steeringMult->KnL = SDDS_Malloc(sizeof(*steeringMult->KnL)*steeringMult->orders)) ||
          !(steeringMult->JnL = SDDS_Malloc(sizeof(*steeringMult->JnL)*steeringMult->orders)) )
        bombElegant("memory allocation failure (computeTotalMultipoleFields)", NULL);
    }
    /* make a list of unique orders for random and systematic multipoles */
    if (systematicMult->orders && randomMult->orders &&
        systematicMult->orders!=randomMult->orders)
      bombElegant("The number of systematic and random multipole error orders must be the same for any given element", NULL);
    if (systematicMult->orders)
      totalMult->orders = systematicMult->orders;
    else
      totalMult->orders = randomMult->orders;
    if (totalMult->orders) {
      if (!(totalMult->order=SDDS_Malloc(sizeof(*totalMult->order)*(totalMult->orders))))
        bombElegant("memory allocation failure (computeTotalMultipoleFields)", NULL);
      if (systematicMult->orders &&
          (!(systematicMult->anMod=SDDS_Malloc(sizeof(*systematicMult->anMod)*systematicMult->orders)) ||
           !(systematicMult->bnMod=SDDS_Malloc(sizeof(*systematicMult->bnMod)*systematicMult->orders)) ||
           !(systematicMult->KnL=SDDS_Malloc(sizeof(*systematicMult->KnL)*systematicMult->orders)) ||
           !(systematicMult->JnL=SDDS_Malloc(sizeof(*systematicMult->JnL)*systematicMult->orders))))
        bombElegant("memory allocation failure (computeTotalMultipoleFields)", NULL);
      if (randomMult->orders &&
          (!(randomMult->anMod=SDDS_Malloc(sizeof(*randomMult->anMod)*randomMult->orders)) ||
           !(randomMult->bnMod=SDDS_Malloc(sizeof(*randomMult->bnMod)*randomMult->orders)) ||
           !(randomMult->KnL=SDDS_Malloc(sizeof(*randomMult->KnL)*randomMult->orders)) ||
           !(randomMult->JnL=SDDS_Malloc(sizeof(*randomMult->JnL)*randomMult->orders))))
        bombElegant("memory allocation failure (computeTotalMultipoleFields", NULL);
      if (!(totalMult->KnL = SDDS_Malloc(sizeof(*totalMult->KnL)*totalMult->orders)) ||
          !(totalMult->JnL = SDDS_Malloc(sizeof(*totalMult->JnL)*totalMult->orders)) )
        bombElegant("memory allocation failure (computeTotalMultipoleFields)", NULL);
      for (i=0; i<totalMult->orders; i++) {
        if (systematicMult->orders && randomMult->orders &&
            systematicMult->order[i]!=randomMult->order[i])
          bombElegant("multipole orders in systematic and random lists must match up for any given element.",
               NULL);
        if (systematicMult->orders) {
          totalMult->order[i] = systematicMult->order[i] ;
          systematicMult->anMod[i] = systematicMult->an[i]*dfactorial(systematicMult->order[i])/
            ipow(systematicMult->referenceRadius, systematicMult->order[i]);
          systematicMult->bnMod[i] = systematicMult->bn[i]*dfactorial(systematicMult->order[i])/
            ipow(systematicMult->referenceRadius, systematicMult->order[i]);
        } else {
          totalMult->order[i] = randomMult->order[i];
          /* anMod and bnMod will be computed later for randomized multipoles */
        }
      }
    }
  }
  
  if (randomMult->orders)
    randomizeErrorMultipoleFields(randomMult);
  
  /* compute normal (KnL) and skew (JnL) from an and bn
   * KnL = an*n!/r^n*(KmL*r^m/m!), 
   * JnL = bn*n!/r^n*(KmL*r^m/m!), where m is the root order 
   * of the magnet with strength KmL
   */
  if (systematicMult->orders)
    sFactor = KmL/dfactorial(rootOrder)*ipow(systematicMult->referenceRadius, rootOrder);
  if (randomMult->orders)
    rFactor = KmL/dfactorial(rootOrder)*ipow(randomMult->referenceRadius, rootOrder);
  for (i=0; i<totalMult->orders; i++) {
    totalMult->KnL[i] = totalMult->JnL[i] = 0;
    if (systematicMult->orders) {
      totalMult->KnL[i] += sFactor*systematicMult->anMod[i];
      totalMult->JnL[i] += sFactor*systematicMult->bnMod[i];
    }
    if (randomMult->orders) {
      totalMult->KnL[i] += rFactor*randomMult->anMod[i];
      totalMult->JnL[i] += rFactor*randomMult->bnMod[i];
    }
  }
  if (steeringMult) {
    /* same for steering multipoles, but compute KnL/theta and JnL/theta (in this case m=0) */
    for (i=0; i<steeringMult->orders; i++) {
      steeringMult->KnL[i] = 
        -1*steeringMult->an[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
      steeringMult->JnL[i] = 
        -1*steeringMult->bn[i]*dfactorial(steeringMult->order[i])/ipow(steeringMult->referenceRadius, steeringMult->order[i]);
    }
  }
}

