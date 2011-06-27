/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "matlib.h"

/* structure for tune correction information */
typedef struct {
    double tunex, tuney;    /* desired tunes */
    char **name;            /* names of quadrupole families */
    long n_families;        /* number of families */
    long n_iterations;      /* number of times to repeat correction */
    double gain;            /* gain for correction */
    long use_perturbed_matrix;
    double tolerance;       /* how close to get to desired tunes */
    long step_up_interval;  /* number of steps to take before increasing gain */
    double delta_gain;      /* amount by which to change the gain */
    double maximum_gain;    /* maximum gain to use */
    MATRIX *T;              /* Nfx2 matrix to give quadrupole strength changes to change 
                               chromaticities by given amount */
    MATRIX *dK1;           /* Nfx1 matrix of quadrupole strength changes */
    MATRIX *dtune;         /* 2x1 matrix of desired tune changes */
    } TUNE_CORRECTION;


/* prototypes for tune.c */
void setup_tune_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune);
long do_tune_correction(TUNE_CORRECTION *tune, RUN *run, LINE_LIST *beamline, double *clorb, long step, long last_iteration);
void computeTuneCorrectionMatrix(RUN *run, LINE_LIST *beamline, TUNE_CORRECTION *tune, long printout);
