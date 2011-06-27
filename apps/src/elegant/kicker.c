/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: kicker.c
 * contents: track_through_kicker()
 * 
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "table.h"

void set_up_kicker(KICKER *kicker);
void set_up_mkicker(MKICKER *kicker);

void track_through_kicker(
    double **part, long np, KICKER *kicker, double p_central, long pass, long default_order
    )
{
    long i, n, ip;
    double time, time_offset, angle, t0, *coord, sum_amp, amplitude, ds;
    double x, xp, y, yp, dp, s, curv, dx;
    double theta_i, alpha_i, alpha_f;

    if (np<=0)
        return;

    log_entry("track_through_kicker");

    if (!kicker->t_wf)
        set_up_kicker(kicker);

    if (kicker->fire_on_pass>pass) {
        drift_beam(part, np, kicker->length, default_order); 
        log_exit("track_through_kicker");
        return;
        }

    if (!kicker->t_wf || !kicker->amp_wf || !kicker->n_wf ||
        kicker->tmin>=kicker->tmax)        
        bombElegant("no (valid) waveform data for kicker", NULL);

    if (kicker->phase_reference==0) 
        kicker->phase_reference = unused_phase_reference();

    switch (get_phase_reference(&time_offset, kicker->phase_reference)) {
        case REF_PHASE_RETURNED:
            break;
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
            if (!kicker->fiducial_seen) {
                /* set reference phase so that the center of this bunch goes through
                 * at the desired phase.
                 */
                t0 = 0;
                for (ip=t0=0; ip<np; ip++)
                    t0 += part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5]));
#if (!USE_MPI)
                t0 /= np;
#else
		if (notSinglePart) {
		  double t0_sum; 
		  long np_total;
		  
		  MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers); 
		  MPI_Allreduce(&t0, &t0_sum, 1, MPI_DOUBLE, MPI_SUM, workers);
		  t0 = t0_sum/np_total;
		} else
		  t0 /= np;
#endif
                kicker->t_fiducial = t0;
                kicker->fiducial_seen = 1;
                }
            set_phase_reference(kicker->phase_reference, time_offset = kicker->t_fiducial);
            break;
        default:
            bombElegant("unknown return value from get_phase_reference()", NULL);
            break;
        }

    if (kicker->length<=0) {
        log_exit("track_through_kicker");
        return;
        }

    time_offset += kicker->time_offset;

    if (kicker->dx || kicker->dy || kicker->dz)
      offsetBeamCoordinates(part, np, kicker->dx, kicker->dy, kicker->dz);
    if (kicker->tilt)
      rotateBeamCoordinates(part, np, kicker->tilt);

    sum_amp = 0;
    for (ip=0; ip<np; ip++) {
        angle  = kicker->angle;
        time = part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5])) - time_offset;

        if (time>kicker->tmax || time<kicker->tmin) {
            if (kicker->periodic) {
                n = (time-kicker->tmin)/(kicker->tmax-kicker->tmin);
                time -= n*(kicker->tmax-kicker->tmin);
                }
            }
  
        amplitude = 0;
        if (time>kicker->tmax || time<kicker->tmin)
            angle *= 0;
        else {
            /* should to a binary search here ! */
            for (i=1; i<kicker->n_wf; i++)
                if (kicker->t_wf[i]>=time)
                    break;
            if (i==kicker->n_wf) {
                fprintf(stdout, "error: waveform interpolation problem in track_through_kicker()\n");
                fflush(stdout);
                fprintf(stdout, "particle time is %21.15e\nwaveform limits are %21.15e to %21.15e\n",
                        time, kicker->tmin, kicker->tmax);
                fflush(stdout);
                exitElegant(1);
                }
            i--;
            angle *= (amplitude=INTERPOLATE(kicker->amp_wf[i], kicker->amp_wf[i+1],
                               kicker->t_wf[i],   kicker->t_wf[i+1], time));
            }
        sum_amp += amplitude;

        coord = part[ip];

        if (!angle) {
            coord[0] += kicker->length*coord[1];
            coord[2] += kicker->length*coord[3];
            coord[4] += kicker->length*sqrt(1+sqr(coord[1])+sqr(coord[3]));
            continue;
            }

        x  = coord[0];
        xp = coord[1];
        y  = coord[2];
        yp = coord[3];
        s  = coord[4];
        dp = coord[5];

        if (kicker->n_kicks<=0) {
          curv = sin(-angle)/kicker->length/(1+dp);
          theta_i = atan(xp);
          alpha_i = -theta_i;
          alpha_f = asin(kicker->length*curv + sin(alpha_i));
          if (fabs(curv*kicker->length)<1e-12) {
            x += (dx=kicker->length*xp);
            s += (ds=fabs(kicker->length/cos(alpha_i)));
          }
          else {
            x += (dx=(cos(alpha_f)  - cos(alpha_i))/curv);
            s += (ds=fabs((alpha_f-alpha_i)/curv));
          }
          xp = tan(theta_i - (alpha_f - alpha_i));
          y += kicker->length*yp;
        }
        else {
          double l1, x0, y0;
          l1 = kicker->length/kicker->n_kicks;
          for (i=0; i<kicker->n_kicks; i++) {
            curv = sin(-angle)/kicker->length/(1+dp)*(1+kicker->b2*(x*x-y*y));
            theta_i = atan(xp);
            alpha_i = -theta_i;
            alpha_f = asin(l1*curv + sin(alpha_i));
            x0 = x;
            if (fabs(curv*l1)<1e-12) {
              x += (dx=l1*xp);
              s += (ds=fabs(l1/cos(alpha_i)));
            }
            else {
              x += (dx=(cos(alpha_f)  - cos(alpha_i))/curv);
              s += (ds=fabs((alpha_f-alpha_i)/curv));
            }
            xp = tan(theta_i - (alpha_f - alpha_i));
            x0 = (x+x0)/2;
            y0 = y;
            y += l1*yp ;
            y0 = (y+y0)/2;
            yp += -2*sin(angle)/kicker->length/(1+dp)*ds*x0*y0*kicker->b2;
          }
        }
    
        coord[0] = x;
        coord[1] = xp;
        coord[2] = y;
        coord[3] = yp;
        coord[4] = s;
      }
    
/*
    if (np)
        fprintf(stdout, "average kicker amplitude = %f\n", sum_amp/np);
        fflush(stdout);
 */

    if (kicker->tilt)
      rotateBeamCoordinates(part, np, -kicker->tilt);
    if (kicker->dx || kicker->dy || kicker->dz)
      offsetBeamCoordinates(part, np, -kicker->dx, -kicker->dy, -kicker->dz);

    log_exit("track_through_kicker");
    }

void set_up_kicker(KICKER *kicker)
{
    TABLE data;
    long i;

    log_entry("set_up_kicker");

    if (!kicker->waveform)
        bombElegant("no waveform filename given for kicker", NULL);

    if (!getTableFromSearchPath(&data, kicker->waveform, 1, 0))
        bombElegant("unable to read waveform for kicker", NULL);

    if (data.n_data<=1)
        bombElegant("kicker waveform contains less than 2 points", NULL);

    kicker->t_wf   = data.c1;
    kicker->amp_wf = data.c2;
    kicker->n_wf   = data.n_data;
    for (i=0; i<kicker->n_wf-1; i++)
        if (kicker->t_wf[i]>kicker->t_wf[i+1])
            bombElegant("time values not monotonically increasing in kicker waveform", NULL);
    kicker->tmin = kicker->t_wf[0];
    kicker->tmax = kicker->t_wf[kicker->n_wf-1];
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    log_exit("set_up_kicker");
    }


void track_through_mkicker(
                          double **part, long np, MKICKER *kicker, double p_central, long pass, long default_order
                          )
{
  long i, n, ip;
  double time, time_offset, strength, t0, amplitude;
  double dummy;
  long n_parts;
  
  if (np<=0)
    return;

  if (!kicker->t_wf)
    set_up_mkicker(kicker);

  if (kicker->order<1) 
    bombElegant("order<1 for MKICKER", NULL);
  
  if (kicker->fire_on_pass>pass) {
    drift_beam(part, np, kicker->length, default_order); 
    return;
  }

  if (!kicker->t_wf || !kicker->amp_wf || !kicker->n_wf ||
      kicker->tmin>=kicker->tmax)        
    bombElegant("no (valid) waveform data for MBUMPER", NULL);

  if (kicker->phase_reference==0) 
    kicker->phase_reference = unused_phase_reference();

  switch (get_phase_reference(&time_offset, kicker->phase_reference)) {
  case REF_PHASE_RETURNED:
    break;
  case REF_PHASE_NOT_SET:
  case REF_PHASE_NONEXISTENT:
    if (!kicker->fiducial_seen) {
      /* set reference phase so that the center of this bunch goes through
       * at the desired phase.
       */
      t0 = 0;
      for (ip=t0=0; ip<np; ip++)
        t0 += part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5]));
#if (!USE_MPI)
      t0 /= np;
#else
      if (notSinglePart) {
        double t0_sum; 
        long np_total;
        
        MPI_Allreduce(&np, &np_total, 1, MPI_LONG, MPI_SUM, workers); 
        MPI_Allreduce(&t0, &t0_sum, 1, MPI_DOUBLE, MPI_SUM, workers);
        t0 = t0_sum/np_total;
      } else
        t0 /= np;
#endif
      kicker->t_fiducial = t0;
      kicker->fiducial_seen = 1;
    }
    set_phase_reference(kicker->phase_reference, time_offset = kicker->t_fiducial);
    break;
  default:
    bombElegant("unknown return value from get_phase_reference()", NULL);
    break;
  }

  if (kicker->length<=0)
    return;

  time_offset += kicker->time_offset;

  if ((n_parts = ceil(kicker->n_kicks/4.0))<1)
    n_parts = 1;

  if (kicker->dx || kicker->dy || kicker->dz)
    offsetBeamCoordinates(part, np, kicker->dx, kicker->dy, kicker->dz);
  if (kicker->tilt)
    rotateBeamCoordinates(part, np, kicker->tilt);

  for (ip=0; ip<np; ip++) {
    strength = kicker->strength;
    time = part[ip][4]/(c_mks*beta_from_delta(p_central, part[ip][5])) - time_offset;
    
    if (time>kicker->tmax || time<kicker->tmin) {
      if (kicker->periodic) {
        n = (time-kicker->tmin)/(kicker->tmax-kicker->tmin);
        time -= n*(kicker->tmax-kicker->tmin);
      }
    }
    
    amplitude = 0;
    if (time>kicker->tmax || time<kicker->tmin)
      strength *= 0;
    else {
      /* should to a binary search here ! */
      for (i=1; i<kicker->n_wf; i++)
        if (kicker->t_wf[i]>=time)
          break;
      if (i==kicker->n_wf) {
        fprintf(stdout, "error: waveform interpolation problem in track_through_mkicker()\n");
        fflush(stdout);
        fprintf(stdout, "particle time is %21.15e\nwaveform limits are %21.15e to %21.15e\n",
                time, kicker->tmin, kicker->tmax);
        fflush(stdout);
        exitElegant(1);
      }
      i--;
      strength *= (amplitude=INTERPOLATE(kicker->amp_wf[i], kicker->amp_wf[i+1],
                                      kicker->t_wf[i],   kicker->t_wf[i+1], time));
    }

    integrate_kick_multipole_ord4(part[ip], kicker->dx, kicker->dy, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  kicker->order, 0, strength*kicker->length, n_parts,
                                  kicker->length, NULL, NULL, NULL, &dummy, NULL);
  }

  if (kicker->tilt)
    rotateBeamCoordinates(part, np, -kicker->tilt);
  if (kicker->dx || kicker->dy || kicker->dz)
    offsetBeamCoordinates(part, np, -kicker->dx, -kicker->dy, -kicker->dz);
}

void set_up_mkicker(MKICKER *kicker)
{
  TABLE data;
  long i;

  if (!kicker->waveform)
    bombElegant("no waveform filename given for MBUMPER", NULL);

  if (!getTableFromSearchPath(&data, kicker->waveform, 1, 0))
    bombElegant("unable to read waveform for MBUMPER", NULL);

  if (data.n_data<=1)
    bombElegant("kicker waveform contains less than 2 points", NULL);

  kicker->t_wf   = data.c1;
  kicker->amp_wf = data.c2;
  kicker->n_wf   = data.n_data;
  for (i=0; i<kicker->n_wf-1; i++)
    if (kicker->t_wf[i]>kicker->t_wf[i+1])
      bombElegant("time values not monotonically increasing in MBUMPER waveform", NULL);
  kicker->tmin = kicker->t_wf[0];
  kicker->tmax = kicker->t_wf[kicker->n_wf-1];
  tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
  data.xlab = data.ylab = data.title = data.topline = NULL;
}

