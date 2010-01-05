/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: misalign_matrix()
 * purpose: alter a matrix to reflect misalignments
 * N.B.:  to be used correctly, this routine should only be called AFTER 
 *        the tilt or rotation of the element is applied.
 *
 * Michael Borland, 1991.
 */
#include "mdb.h"
#include "track.h"

void misalign_matrix(VMATRIX *M, double dx, double dy, double dz, double bend_angle)
{
    VMATRIX *M1, *M2, *tmp;
    static MALIGN mal;

    log_entry("misalign_matrix");

    /* entrance misalignment corresponds to opposite tranverse offset of the beam
     * plus propagation of beam forward/backward to the new entrance of the magnet 
     */
    mal.dx = -dx; mal.dxp = 0;
    mal.dy = -dy; mal.dyp = 0;
    mal.dz = dz;
    mal.dp = mal.dt = 0;
    M1 = misalignment_matrix(&mal, M->order);

    /* exit misalignment corresponds to undoing the entrance misalignment of the
     * beam, with consideration given to the fact that the curvilinear coordinate
     * system at the exit may be rotated relative to that at the exit.  This
     * code is incorrect if the bending magnet is tilted!
     */
    mal.dx =   dx*cos(bend_angle) + dz*sin(bend_angle);
    mal.dz =   dx*sin(bend_angle) - dz*cos(bend_angle);
    mal.dy = dy;
    M2 = misalignment_matrix(&mal, M->order);
    
    tmp = tmalloc(sizeof(*tmp));
    initialize_matrices(tmp, tmp->order=M->order);

    concat_matrices(tmp, M, M1, 0);
    concat_matrices(M, M2, tmp, 0);

    free_matrices(tmp); tfree(tmp); tmp = NULL;
    free_matrices(M1); tfree(M1); M1 = NULL;
    free_matrices(M2); tfree(M2); M2 = NULL;

    log_exit("misalign_matrix");
    }

/* routine: misalignment_matrix()
 * purpose: provide a matrix for a misalignment
 *
 * Michael Borland, 1991.
 */
    
VMATRIX *misalignment_matrix(MALIGN *malign, long order)
{
    VMATRIX *M;
    static VMATRIX *M1 = NULL;
    double *C, **R;

    log_entry("misalignment_matrix");

    M = tmalloc(sizeof(*M));
    M->order = order;
    initialize_matrices(M, M->order);
    R = M->R;
    C = M->C;

    C[0] = malign->dx;
    C[1] = malign->dxp;
    C[2] = malign->dy;
    C[3] = malign->dyp;
    C[4] = 0;
    /* ignore malign->dt.  Not sure what it really means. */
    C[5] = malign->dp + malign->de;            /* missing 1/beta^2 on malign->de term here ... */

    R[0][0] = R[1][1] = R[2][2] = R[3][3] = R[4][4] = R[5][5] = 1;

    if (malign->dz) {
        /* do longitudinal misalignment */
        VMATRIX *tmp;

        if (!M1) {
            M1 = tmalloc(sizeof(*M1));
            initialize_matrices(M1, M1->order=M->order);
            }
        tmp = drift_matrix(malign->dz, M->order);
        concat_matrices(M1, M, tmp, 0);
        free_matrices(tmp); tfree(tmp); tmp = NULL;
        copy_matrices1(M, M1);
        }
    log_exit("misalignment_matrix");
    return(M);
    }


void offset_matrix(VMATRIX *M, double dx, double dxp, double dy, double dyp)
{
    VMATRIX *M1, M2;
    static MALIGN mal;

    log_entry("offset_matrix");

    mal.dx = dx; mal.dxp = dxp;
    mal.dy = dy; mal.dyp = dyp;
    mal.dz = mal.dp = mal.dt = 0;
    M1 = misalignment_matrix(&mal, M->order);

    copy_matrices(&M2, M);
    concat_matrices(M, &M2, M1, 0);

    free_matrices(M1); tfree(M1); M1 = NULL;
    free_matrices(&M2); 

    log_exit("offset_matrix");
    }

/* The name is a misnomer: we offset the coordinates of a beam
 * to emulate offsetting of the upcoming element by the values
 */

void offsetBeamCoordinates(double **coord, long np, double dx, double dy, double dz)
{
  long ip;
  double *part;
  
  for (ip=np-1; ip>=0; ip--) {
    part = coord[ip];
    part[4] += dz*sqrt(1 + sqr(part[1]) + sqr(part[3]));
    part[0]  = part[0] - dx + dz*part[1];
    part[2]  = part[2] - dy + dz*part[3];
  }
}

