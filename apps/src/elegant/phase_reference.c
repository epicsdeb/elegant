/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: get_phase_reference()
 * purpose: store/retrieve reference phases for time-dependent elements
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"

/* flag bit values */
#define FL_REF_PHASE_SET 1

struct phase_reference {
    long ref_number;
    long flags;
    double phase;
    } *reference = NULL;
static long n_references=0;


long get_phase_reference(
    double *phase,
    long phase_ref_number
    )
{
    long i;
#ifdef DEBUG    
    fprintf(stdout, "get_phase_reference(%le, %ld)\n", *phase, phase_ref_number);
    fflush(stdout);
#endif

    log_entry("get_phase_reference");

    if (phase_ref_number==0) {
        log_exit("get_phase_reference");
        return(0);
        }

    for (i=0; i<n_references; i++) {
        if (reference[i].ref_number==phase_ref_number) {
            if (reference[i].flags&FL_REF_PHASE_SET) {
                *phase = reference[i].phase;
#ifdef DEBUG
                fprintf(stdout, "returning REF_PHASE_RETURNED\n");
                fflush(stdout);
#endif
                log_exit("get_phase_reference");
                return(REF_PHASE_RETURNED);
                }
#ifdef DEBUG
                fprintf(stdout, "returning REF_PHASE_NOT_SET\n");
                fflush(stdout);
#endif
            log_exit("get_phase_reference");
            return(REF_PHASE_NOT_SET);
            }
        }
#ifdef DEBUG
                fprintf(stdout, "returning REF_PHASE_NONEXISTENT\n");
                fflush(stdout);
#endif
    log_exit("get_phase_reference");
    return(REF_PHASE_NONEXISTENT);
    }

long set_phase_reference(
    long phase_ref_number,  /* number of the phase reference group */
    double phase            /* phase to assert for fiducial particle */
    )
{
    long i;

    log_entry("set_phase_reference");
    
#ifdef DEBUG
    fprintf(stdout, "set_phase_reference(%ld, %le)\n", phase_ref_number, phase);
    fflush(stdout);
#endif

    if (phase_ref_number==0) {
        log_exit("set_phase_reference");
        return(0);
        }

    for (i=0; i<n_references; i++) {
        if (reference[i].ref_number==phase_ref_number) {
            reference[i].phase = phase;
            reference[i].flags = FL_REF_PHASE_SET;
#ifdef DEBUG
            fprintf(stdout, "existing phase reference set\n");
            fflush(stdout);
#endif
            log_exit("set_phase_reference");
            return(1);
            }
        }
    if (phase_ref_number>LONG_MAX/2)
        bombElegant("please use a small integer for the phase_reference number", NULL);
    reference = trealloc(reference, sizeof(*reference)*(++n_references));
    reference[i].ref_number = phase_ref_number;
    reference[i].phase = phase;
    reference[i].flags = FL_REF_PHASE_SET;
#ifdef DEBUG
    fprintf(stdout, "new phase reference set\n");
    fflush(stdout);
#endif
    log_exit("set_phase_reference");
    return(1);
    }

void delete_phase_references(void)
{
    long i;

    log_entry("delete_phase_references");

#ifdef DEBUG
    fprintf(stdout, "phase references deleted\n");
    fflush(stdout);
#endif
    for (i=0; i<n_references; i++)
        reference[i].flags = 0;
    log_exit("delete_phase_references");
    }

long unused_phase_reference()
{
    static long big=LONG_MAX;

    log_entry("unused_phase_reference");

    reference = trealloc(reference, sizeof(*reference)*(n_references+1));
    reference[n_references].ref_number = big;
    reference[n_references].phase = 0;
    reference[n_references].flags = 0;
    n_references++;
    log_exit("unused_phase_reference");
    return(big--);
    }

double get_reference_phase(long phase_ref, double phase0)
    /* routine to maintain emulate old get_reference_phase(),
     * which is obsolete and should be phased out (pun intended)
     */
{
    double phase;

    log_entry("get_reference_phase");
#ifdef DEBUG
    fprintf(stdout, "obsolete routine get_reference_phase called\n");
    fflush(stdout);
#endif
    switch (get_phase_reference(&phase, phase_ref)) {
        case REF_PHASE_RETURNED:
            log_exit("get_reference_phase");
            return(phase);
        case REF_PHASE_NOT_SET:
        case REF_PHASE_NONEXISTENT:
        default:
            set_phase_reference(phase_ref, phase0);
            log_exit("get_reference_phase");
            return(phase0);
        }
    }





