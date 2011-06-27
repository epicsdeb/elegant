/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "matlib.h"
#ifdef __cplusplus
extern "C" {
#endif

/* structure for chromaticity correction information */
typedef struct {
    double chromx, chromy;    /* desired chromaticities */
    double strengthLimit;     /* maximum absolute value of strength */
    char **name;              /* names of sextupole families */
    long n_families;          /* number of families */
    long n_iterations;        /* number of times to repeat correction */
    double correction_fraction;  /* to prevent unstable correction */
    long use_perturbed_matrix;
    double sextupole_tweek;
    double tolerance;         /* how close to get to desired chromaticities */
    long exit_on_failure;     /* exit if fails to converge */
    MATRIX *T;                /* Nfx2 matrix to give sextupole strength changes to change 
                                 chromaticities by given amount */
    MATRIX *dK2;              /* Nfx1 matrix of sextupole strength changes */
    MATRIX *dchrom;           /* 2x1 matrix of desired chromaticity changes */
    } CHROM_CORRECTION;


/* prototypes for chrom.c */
void setup_chromaticity_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom);
long do_chromaticity_correction(CHROM_CORRECTION *chrom, RUN *run, LINE_LIST *beamline, double *clorb,
        long step, long last_iteration);
void computeChromaticities(double *chromx, double *chromy, 
                           double *dbetax, double *dbetay,
                           double *dalphax, double *dalphay,
                           TWISS *twiss0, TWISS *twiss1, VMATRIX *M);
void computeHigherOrderChromaticities(LINE_LIST *beamline, double *clorb, RUN *run,
				      long concatOrder, double deltaStep, long deltaPoints, long quickMode);
void computeChromCorrectionMatrix(RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom);
void computeChromaticTuneLimits(LINE_LIST *beamline);
#ifdef __cplusplus
}
#endif

