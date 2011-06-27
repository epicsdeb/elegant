/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: copy_particles()
 * purpose: copy particle data from one array to another
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

void copy_particles(
    double **copy,
    double **original,
    long n_particles
    )
{
    register long ip, ic;
    double *cptr, *optr;

    log_entry("copy_particles");
    if (!copy)
        bombElegant("can't copy particles--target is NULL pointer (copy_particles)", NULL);
    if (!original)
        bombElegant("can't copy particles--source is NULL pointer (copy_particles)", NULL);

    for (ip=n_particles-1; ip>=0; ip--) {
        cptr = *copy++;
        optr = *original++;
        if (!cptr)
            bombElegant("element of target array is NULL pointer (copy_particles)", NULL);
        if (!optr)
            bombElegant("element of source array is NULL pointer (copy_particles)", NULL);
        for (ic=6; ic>=0; ic--)
            *cptr++ = *optr++;
        }
    log_exit("copy_particles");
    }
 
