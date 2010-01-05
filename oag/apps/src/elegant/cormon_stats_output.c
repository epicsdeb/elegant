/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: cormon_stats_output.c
 * purpose: output of corrector/monitor statistics for orbit/trajectory correction
 * 
 * Michael Borland, 1993 
 */
#include "mdb.h"
#include "track.h"
#include "SDDS.h"

static long SDDS_cormon_initialized = 0;
static SDDS_TABLE SDDS_cormon;

#define IC_ITERATION 0
#define IC_CYCLE 1
#define IC_KRMS 2
#define IC_PRMS 3
#define IC_KMAX 4
#define IC_PMAX 5
#define IC_CDP  6
#define IC_STAGE 7
#define N_COLUMNS 8
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"Iteration", "&column name=Iteration, type=long, description=\"Correction iteration\" &end"},
    {"Cycle", "&column name=Cycle, type=long, description=\"Correction cycle\" &end"},
    {"Krms", "&column name=Krms, units=rad, type=double, symbol=\"K$brms$n\", description=\"RMS corrrection kick\" &end"},
    {"Prms", "&column name=Prms, units=m, type=double, symbol=\"P$brms$n\", description=\"RMS monitor position\" &end"},
    {"Kmax", "&column name=Kmax, units=rad, type=double, symbol=\"K$bmax$n\", description=\"Maximum corrrection kick\" &end"},
    {"Pmax", "&column name=Pmax, units=m, type=double, symbol=\"P$bmax$n\", description=\"Maximum monitor position\" &end"},
    {"Cdelta", "&column name=Cdelta, type=double, symbol=\"<$gd$r>\", description=\"Momentum centroid offset\" &end"},
    {"Stage", "&column name=Stage, type=string, description=\"Stage of correction\" &end"},
    } ;

#define IP_STEP 0
#define IP_PLANE 1
#define IP_FINAL 2
#define N_PARAMETERS 3
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"Plane", "&parameter name=Plane, type=string, description=\"Transverse plane\" &end"},
    {"Final", "&parameter name=Final, type=character, description=\"Final results?\" &end"},
    } ;

void setup_cormon_stats(char *filename, RUN *run)
{

    log_entry("setup_cormon_stats");

    if (!SDDS_cormon_initialized) {
        SDDS_cormon_initialized = 1;
        zero_memory(&SDDS_cormon, sizeof(SDDS_cormon));
        }
    if (!filename)
        bomb("NULL filename passed (setup_cormon_stats)", NULL);
    if (!run)
        bomb("NULL RUN pointer passed (setup_cormon_stats)", NULL);

    SDDS_ElegantOutputSetup(&SDDS_cormon, filename, SDDS_BINARY, 1, "corrector/monitor statistics", run->runfile,
                            run->lattice, parameter_definition, N_PARAMETERS,
                            column_definition, N_COLUMNS, "setup_cormon_stats",
                            SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    log_exit("setup_cormon_stats");
    }

void dump_cormon_stats(long verbose, long plane, double **kick, long n_kicks, 
    double **position, long n_positions, double *Cdp, long n_iterations, long cycle,
    long final_cycle, long sim_step, long textOnly)
{
    long i, j;
    static double **data = NULL;
    static char **stage = NULL;
    static long n_iter = 0;

    if (!verbose && !SDDS_cormon_initialized)
        return;

    log_entry("dump_cormon_stats");

    if (!data || n_iterations>n_iter) {
        if (data)
            free_czarray_2d((void**)data, n_iter+1, N_COLUMNS);
        if (stage)
            free(stage);
        data = (double**)czarray_2d(sizeof(**data), n_iterations+1, N_COLUMNS);    
        stage = tmalloc(sizeof(*stage)*(n_iterations+1));
        n_iter = n_iterations;
        }

    for (j=0; j<=n_iterations; j++) {
        if (j==0)
            stage[j] = "uncorrected";
        else if (j==n_iterations)
            stage[j] = "corrected";
        else
            stage[j] = "intermediate";
        data[j][IC_ITERATION] = j;
        data[j][IC_CYCLE] = cycle;
        data[j][IC_KRMS] = data[j][IC_KMAX] = data[j][IC_PRMS] = 
            data[j][IC_PMAX] = 0;
        for (i=0; i<n_kicks; i++) {
            if (FABS(kick[j][i]) > data[j][IC_KMAX])
                data[j][IC_KMAX] = FABS(kick[j][i]);
            data[j][IC_KRMS] += sqr(kick[j][i]);
            }
        if (n_kicks) 
            data[j][IC_KRMS] = sqrt(data[j][IC_KRMS]/n_kicks);
        for (i=0; i<n_positions; i++) {
            if (FABS(position[j][i]) > data[j][IC_PMAX])
                data[j][IC_PMAX] = FABS(position[j][i]);
            data[j][IC_PRMS] += sqr(position[j][i]);
            }
        if (n_positions)
            data[j][IC_PRMS] = sqrt(data[j][IC_PRMS]/n_positions);
        data[j][IC_CDP] = (Cdp?Cdp[j]:0);
        }

    if (SDDS_cormon_initialized && !textOnly) {
        if (!SDDS_StartTable(&SDDS_cormon, n_iterations+1)) {
            SDDS_SetError("Unable to start SDDS table (dump_cormon_stats)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!SDDS_SetParameters(&SDDS_cormon, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 
                                IP_STEP, sim_step, IP_PLANE, plane?"vertical":"horizontal", 
                                IP_FINAL, final_cycle?'y':'n', -1)) {
            SDDS_SetError("Unable to set SDDS parameter values (dump_cormon_stats)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        for (j=0; j<=n_iterations; j++) {
            if (!SDDS_SetRowValues(&SDDS_cormon, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, j,
                                   IC_ITERATION, j,
                                   IC_CYCLE, cycle,
                                   IC_KRMS, data[j][IC_KRMS],
                                   IC_PRMS, data[j][IC_PRMS],
                                   IC_KMAX, data[j][IC_KMAX],
                                   IC_PMAX, data[j][IC_PMAX],
                                   IC_CDP, data[j][IC_CDP],
                                   IC_STAGE, stage[j], -1)) {
                fprintf(stdout, "Unable to set row %ld values (dump_cormon_stats)\n", j);
                fflush(stdout);
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                exit(1);
                }
            }
        if (!SDDS_WriteTable(&SDDS_cormon)) {
            SDDS_SetError("Unable to write corrector data (dump_cormon_stats)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!inhibitFileSync)
          SDDS_DoFSync(&SDDS_cormon);
        if (!SDDS_EraseData(&SDDS_cormon)) {
            SDDS_SetError("Unable to erase corrector data (dump_cormon_stats)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }

    if (verbose) {
        for (j=0; j<=n_iterations; j++) 
          if (Cdp) {
                fprintf(stdout, "   %c     %4ld     %8.3f    %8.3f     %8.3f     %8.3f   %8.3f %s\n",
                    'x'+(plane?1:0), j,
                    1e3*data[j][IC_KRMS], 1e3*data[j][IC_PRMS],
                    1e3*data[j][IC_KMAX], 1e3*data[j][IC_PMAX],
                    1e2*data[j][IC_CDP],
                    (j==n_iterations?"**":"") );
                fflush(stdout);
              }
          else {
                fprintf(stdout, "   %c     %4ld     %8.3f    %8.3f     %8.3f     %8.3f %s\n",
                    'x'+(plane?1:0), j,
                    1e3*data[j][IC_KRMS], 1e3*data[j][IC_PRMS],
                    1e3*data[j][IC_KMAX], 1e3*data[j][IC_PMAX],
                    (j==n_iterations?"**":"") );
                fflush(stdout);
              }
        }

    log_exit("dump_cormon_stats");    
    }

void finish_cormon_stats()
{
    if (!SDDS_cormon_initialized)
        return;
    if (!SDDS_Terminate(&SDDS_cormon)) {
        fprintf(stdout, "Unable to terminate SDDS output for correctors (finish_cormon_stats)\n");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }
    SDDS_cormon_initialized = 0;
    }
