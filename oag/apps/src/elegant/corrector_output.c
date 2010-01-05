/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: corrector_output.c
 * purpose: trajectory/orbit correction corrector output for elegant
 * 
 * Michael Borland, 1993.
 */
#include "mdb.h"
#include "track.h"
#include "SDDS.h"
#include "correctDefs.h"

/* Variables for output of corrector data */
static long SDDS_cor_initialized = 0;
static SDDS_TABLE SDDS_cor;

#define IC_S 0
#define IC_KICK 1
#define IC_ELEMENT 2
#define IC_OCCURENCE 3
#define IC_PARAMETER 4
#define IC_VALUE 5
#define IC_PARUNITS 6
#define N_COLUMNS 7
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double &end"},
    {"Kick", "&column name=Kick, units=rad, type=double &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementParameter", "&column name=ElementParameter, type=string, description=\"Parameter name\" &end"},
    {"ParameterValue", "&column name=ParameterValue, type=double, description=\"Parameter value\" &end"},
    {"ParameterUnits", "&column name=ParameterUnits, type=string &end"},
    } ;

#define IP_STEP 0
#define IP_PLANE 1
#define N_PARAMETERS 2
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    {"Plane", "&parameter name=Plane, type=string, description=\"Transverse plane\" &end"},
    } ;

void setup_corrector_output(char *filename, RUN *run)
{
    log_entry("setup_corrector_output");
    if (!SDDS_cor_initialized) {
        SDDS_cor_initialized = 1;
        zero_memory(&SDDS_cor, sizeof(SDDS_cor));
        }

    if (!filename)
        bomb("NULL filename passed (setup_corrector_output)", NULL);
    if (!run)
        bomb("NULL RUN pointer passed (setup_corrector_output)", NULL);

    SDDS_ElegantOutputSetup(&SDDS_cor, filename, SDDS_BINARY, 1, "corrector data", run->runfile,
                            run->lattice, parameter_definition, N_PARAMETERS,
                            column_definition, N_COLUMNS, "setup_corrector_output",
                            SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    log_exit("setup_corrector_output");
    }

void dump_corrector_data(CORMON_DATA *CM, STEERING_LIST *SL, long index, char *plane, long step)
{
    long i, sl_index;
    double value;
    ELEMENT_LIST *eptr;

    if (!SDDS_cor_initialized)
        return;

    log_entry("dump_corrector_data");

    if (!SDDS_StartTable(&SDDS_cor, CM->ncor)) {
        SDDS_SetError("Unable to start SDDS table (dump_corrector_data)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (!SDDS_SetParameters(&SDDS_cor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            IP_STEP, step, IP_PLANE, plane, -1)) {
        SDDS_SetError("Unable to set SDDS parameters (dump_corrector_data)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    for (i=0; i<CM->ncor; i++) {
        eptr = CM->ucorr[i];
        sl_index = CM->sl_index[i];
        value = *(double*)((char*)eptr->p_elem+SL->param_offset[sl_index]);
        if (!SDDS_SetRowValues(&SDDS_cor, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               IC_S, eptr->end_pos,
                               IC_KICK, CM->kick[index][i],
                               IC_ELEMENT, eptr->name,
                               IC_OCCURENCE, eptr->occurence,
                               IC_PARAMETER, SL->corr_param[sl_index],
                               IC_VALUE, value,
                               IC_PARUNITS, entity_description[eptr->type].parameter[SL->param_index[sl_index]].unit,
                               -1)) {
            fprintf(stdout, "Unable to set row %ld values (dump_corrector_data)\n", i);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        }

    if (!SDDS_WriteTable(&SDDS_cor)) {
        SDDS_SetError("Unable to write corrector data (dump_corrector_data)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_EraseData(&SDDS_cor)) {
        SDDS_SetError("Unable to erase corrector data (dump_corrector_data)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    log_exit("dump_corrector_data");
    }

void finish_corrector_output()
{
    if (!SDDS_cor_initialized)
        return;
    if (!SDDS_Terminate(&SDDS_cor)) {
        fprintf(stdout, "Unable to terminate SDDS output for correctors (finished_corrector_output)\n");
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
        }
    SDDS_cor_initialized = 0;
    }

