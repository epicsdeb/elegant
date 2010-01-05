/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: polynomial_kicks()
 * purpose: apply kicks due to a polynomial dependence on x or y
 * 
 * Michael Borland, 1992 
 */
#include "mdb.h"
#include "track.h"

long polynomial_kicks(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    KPOLY *kpoly,       /* kick-polynomial structure */
    double p_error,     /* p_nominal/p_central */
    double Po,
    double **accepted,
    double z_start
    )
{
    double KL0;         /* on-momentum integrated strength at x=1 or y=1 */
    double KL;
    double dx, dy, dz;  /* offsets of the multipole center */
    long order;         /* order (n) */
    long i_part, i_top, yplane;
    double *coord;
    double cos_tilt, sin_tilt;
    double x, xp, y, yp;

    log_entry("polynomial_kicks");

    if (!particle)
        bomb("particle array is null (polynomial_kicks)", NULL);

    if (!kpoly)
        bomb("null KPOLY pointer (polynomial_kicks)", NULL);

    if ((order=kpoly->order)<0)
        bomb("order < 0 for KPOLY element (polynomial_kicks)", NULL);

    if (kpoly->plane && (kpoly->plane[0]=='y' || kpoly->plane[0]=='Y'))
        yplane = 1;
    else {
        if (kpoly->plane && !(kpoly->plane[0]=='x' || kpoly->plane[0]=='X')) {
            fputs("warning: KPOLY plane not recognized--x plane assumed.", stdout);
            cp_str(&kpoly->plane, "x");
            }
        yplane = 0;
        }

    KL0 = kpoly->coefficient*kpoly->factor;

    cos_tilt = cos(kpoly->tilt);
    sin_tilt = sin(kpoly->tilt);
    dx = kpoly->dx;
    dy = kpoly->dy;
    dz = kpoly->dz;

    i_top = n_part-1;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stdout, "null coordinate pointer for particle %ld (polynomial_kicks)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stdout, "null accepted coordinates pointer for particle %ld (polynomial_kicks)", i_part);
            fflush(stdout);
            abort();
            }

        /* adjust strength for momentum offset--good to all orders */
        KL = KL0/(1+coord[5]);

        /* calculate coordinates in rotated and offset frame */
        coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
        coord[0] += -dx + dz*coord[1];
        coord[2] += -dy + dz*coord[3];

        x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
        y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
        xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
        yp = - sin_tilt*coord[1] + cos_tilt*coord[3];

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

        if (yplane)
            yp += KL*ipow(y, order);
        else
            xp += KL*ipow(x, order);
                
        /* undo the rotation and store in place of initial coordinates */
        /* don't need to change coord[0] or coord[2] since x and y are unchanged */
        coord[1] = cos_tilt*xp - sin_tilt*yp;
        coord[3] = sin_tilt*xp + cos_tilt*yp;

        /* remove the coordinate offsets */
        coord[0] += dx - coord[1]*dz;
        coord[2] += dy - coord[3]*dz;
        coord[4] -= dz*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
        }
    log_exit("polynomial_kicks");
    return(i_top+1);
    }
