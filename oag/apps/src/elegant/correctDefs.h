/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "matrixOp.h"

/* see correction.c for additional explanation of the next three structures */

typedef struct {
    /* arrays for information on individual correcting elements */
    char **corr_name;                      /* names of groups of correcting elements */
    long *corr_type;                       /* type numbers */
    char **corr_param;                     /* parameter names */
    long *param_offset;                    /* offset of correcting parameter in element structure */
    long *param_index;                     /* index of correcting parameter in entity description */
    double *corr_tweek;                    /* tweek values--amount to change parameter by to get dkick/dparam */
    double *corr_limit;                    /* limiting absolute value of the parameter */
    long n_corr_types;
    } STEERING_LIST;

typedef struct {
    /* information on useful monitors and correctors */
    long *mon_index;                       /* index of monitor in trajectory array */
    int32_t nmon, ncor;                       /* numbers of monitors and correctors */
    ELEMENT_LIST **umoni, **ucorr;         /* arrays of pointers to monitor and corrector elements */
    double *kick_coef;                     /* dkick/dparam (==1 for hkick, vkick, and hvkick elements) */
    long *sl_index;                        /* index of steering list entry for each corrector */

    /* arrays for holding corrector information for output */
    double **kick, **posi;
    /* copies of input specifications for correction */
    double corr_fraction, corr_accuracy, corr_limit, bpm_noise, default_tweek, bpm_noise_cutoff;
    long fixed_length, bpm_noise_distribution, default_threading_divisor;
    long remove_smallest_SVs, keep_largest_SVs, auto_limit_SVs;
    double minimum_SV_ratio;
    /* correction matrix and inverse, respectively: */
    /* Mij(C, i, j) = dX(monitor i)/dK(corrector j) */
    MAT *C, *T; 
    /* information about last correction */
    long n_cycles_done;
    } CORMON_DATA;

typedef struct {
    long mode;
#define TRAJECTORY_CORRECTION 0
#define ORBIT_CORRECTION 1
    long method, verbose, track_before_and_after, n_iterations, n_xy_cycles, minimum_cycles;
    long prezero_correctors, start_from_centroid, use_actual_beam, response_only;
    double clorb_accuracy;
    double clorb_iterations;
    double clorb_iter_fraction;
    STEERING_LIST SLx, SLy;
    CORMON_DATA *CMx, *CMy;
    TRAJECTORY **traj;
    } CORRECTION ;
    

extern void compute_trajcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert);
extern void compute_orbcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert, long fixed_length, long verbose);

extern void setup_corrector_output(char *filename, RUN *run);
extern void dump_corrector_data(CORMON_DATA *CM, STEERING_LIST *SL, long index, char *plane, long step);
extern void setup_cormon_stats(char *filename, RUN *run);
extern void dump_cormon_stats(long verbose, long plane, double **kick, long n_kicks, 
    double **position, long n_positions, double *Cdp, long n_iterations, long cycle,
    long final_cycle, long step, long textOnly);
extern void setup_orb_traj_output(char *filename, char *mode, RUN *run);
extern void dump_orb_traj(TRAJECTORY *traj, long n_elems, char *description, long step);

extern void setup_correction_matrix_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, 
                                           CORRECTION *correct,
                                           long *do_response, long BnLUnitsOk);
extern void run_response_output(RUN *run, LINE_LIST *beamline, CORRECTION *correct, long tune_corrected);
extern void correction_setup(CORRECTION *_correct, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline);
extern long do_correction(CORRECTION *correct, RUN *run, LINE_LIST *beamline, double *starting_coords, 
        BEAM *beam, long sim_step, long initial_correction);
extern void add_steering_element(CORRECTION *correct, LINE_LIST *beamline, RUN *run, NAMELIST_TEXT *nltext);
void compute_amplification_factors(NAMELIST_TEXT *nltext, RUN *run, CORRECTION *correct,
    long closed_orbit, LINE_LIST *beamline);

