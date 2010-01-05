/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: ramp_momentum()
 * purpose: do momentum ramping
 *
 * Michael Borland, 1992, 1993
 */
#include "mdb.h"
#include "track.h"
#include "table.h"

void set_up_ramp_p(RAMPP *rampp);
long find_nearby_array_entry(double *entry, long n, double key);
double linear_interpolation(double *y, double *t, long n, double t0, long i);

void ramp_momentum(
    double **coord,
    long np,
    RAMPP *rampp,
    double *P_central,    /* current beta*gamma on input, changed on output */
    long pass
    )
{
    long ip, i_time;
    double P_new, t, t0;

    log_entry("ramp_momentum");

    if (!rampp->t_Pf)
        set_up_ramp_p(rampp);

    if (!rampp->t_Pf || !rampp->Pfactor || rampp->n_pts<2)
        bomb("NULL data for rampp element", NULL);

    if (!rampp->Po)
        rampp->Po = *P_central;

#if defined(DEBUG)
    fprintf(stdout, "*P_central = %15.8e\n", *P_central);
    fflush(stdout);
#endif

    for (ip=t0=0; ip<np; ip++) {
        t0 += (t=coord[ip][4]/(c_mks*beta_from_delta(*P_central, coord[ip][5])));
#if defined(IEEE_MATH)
        if (isnan(t) || isinf(t)) {
            long i;
            fprintf(stdout, "error: bad time coordinate for particle %ld\n", ip);
            fflush(stdout);
            for (i=0; i<6; i++)
                fprintf(stdout, "%15.8e ", coord[ip][i]);
                fflush(stdout);
            fputc('\n', stdout);
            abort();
            }
#endif
        }
    t0 /= np;
    i_time = find_nearby_array_entry(rampp->t_Pf, rampp->n_pts, t0);
    P_new = rampp->Po*linear_interpolation(rampp->Pfactor, rampp->t_Pf, rampp->n_pts, t0, i_time);

#if defined(DEBUG)
    fprintf(stdout, "new momentum for pass %ld, <t> = %.15e s:  %.15e\n", pass, t0, P_new);
    fflush(stdout);
#endif

    for (ip=0; ip<np; ip++)
        coord[ip][5] = (1+coord[ip][5])*(*P_central)/ P_new-1;
        
    *P_central =  P_new;
    log_exit("ramp_momentum");
    }
    
void set_up_ramp_p(RAMPP *rampp)
{
    TABLE data;
    long i;

    log_entry("set_up_rampp");

    if (!rampp->waveform)
        bomb("no waveform filename given for rampp", NULL);

    if (!getTableFromSearchPath(&data, rampp->waveform, 1, 0))
        bomb("unable to read waveform for rampp", NULL);

    if (data.n_data<=1)
        bomb("rampp waveform contains less than 2 points", NULL);

    rampp->Po      = 0;
    rampp->t_Pf    = data.c1;
    rampp->Pfactor = data.c2;
    rampp->n_pts   = data.n_data;
    for (i=0; i<rampp->n_pts-1; i++)
        if (rampp->t_Pf[i]>rampp->t_Pf[i+1])
            bomb("time values are not monotonically increasing in rampp waveform", NULL);
    tfree(data.xlab); tfree(data.ylab); tfree(data.title); tfree(data.topline);
    data.xlab = data.ylab = data.title = data.topline = NULL;
    log_exit("set_up_rampp");
    }
