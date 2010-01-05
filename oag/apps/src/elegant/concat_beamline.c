/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: concat_beamline
 * purpose: concatenate a beamline into a series of matrices and non-matrix elements
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"

#define DEBUG 0

void copy_matrices1(VMATRIX *M1,  VMATRIX *M0);
void free_elements1(ELEMENT_LIST *elemlist);

void concatenate_beamline(LINE_LIST *beamline, RUN *run)
{
    ELEMENT_LIST *lastelem, *elem, *ecat;
    VMATRIX  *M1, *M2, *tmp;
    long new_seq, n_seqs, in_seq;
    long n_matrices, n_nonmatrices;
    char s[100];
    double z, last_z;
    ELEMENT_LIST *pred, *succ;
    
    log_entry("concatenate_beamline");

    if (beamline->flags&BEAMLINE_CONCAT_CURRENT) {
        if (!(beamline->flags&BEAMLINE_CONCAT_DONE))
            bomb("logic problem with beamline concatenation", NULL);
        log_exit("concatenate_beamline");
        return;
        }

    elem = &(beamline->elem);
    ecat = &(beamline->ecat);
    if (beamline->flags&BEAMLINE_CONCAT_DONE)
        free_elements1(&(beamline->ecat));
    ecat->pred = NULL;

#if DEBUG
    fprintf(stdout, "run->concat_order = %ld\n", run->concat_order);
    fflush(stdout);
#endif

    M1 = tmalloc(sizeof(*M1));    initialize_matrices(M1, run->concat_order);
    M2 = tmalloc(sizeof(*M2));    initialize_matrices(M2, run->concat_order);

    new_seq = 1;
    n_seqs = in_seq = 0;
    n_matrices = n_nonmatrices = 0;
    z = 0;
    beamline->ncat_elems = 0;
    do {
#if DEBUG
        fprintf(stdout, "working on %s\n", (elem->name?elem->name:"NULL"));
        fflush(stdout);
#endif        
        if (entity_description[elem->type].flags&HAS_MATRIX && elem->matrix==NULL) {
            compute_matrix(elem, run, NULL);
            }
        last_z = z;
        if (entity_description[elem->type].flags&HAS_LENGTH)
            z += ((DRIFT*)elem->p_elem)->length;
        if (entity_description[elem->type].flags&HAS_MATRIX && elem->matrix->order<=run->concat_order &&
            !(entity_description[elem->type].flags&DONT_CONCAT) ) {
#if DEBUG
            fprintf(stdout, "element has matrix of order %ld\n", elem->matrix->order);
            fflush(stdout);
#endif        
            if (new_seq) {
                /* start concatenating new sequence of matrices */
#if DEBUG
                fprintf(stdout, "starting new sequence of concatenated matrices\n");
                fflush(stdout);
#endif
                copy_matrices1(M1, elem->matrix);
                new_seq = 0;
                }
            else {
                concat_matrices(M2, elem->matrix, M1,
                                entity_description[elem->type].flags&HAS_RF_MATRIX?
                                CONCAT_EXCLUDE_S0:0);
                tmp = M1;  M1 = M2;  M2 = tmp;
                }
            in_seq = 1;
            }
        else {
            if (in_seq) {
                /*  end of sequence--copy concatenated matrix into ecat list */
#if DEBUG
                fprintf(stdout, "ending sequence of concatenated matrices\n");
                fflush(stdout);
#endif
                ecat->matrix = tmalloc(sizeof(*(ecat->matrix)));
                copy_matrices(ecat->matrix, M1);
                sprintf(s, "M%ld", n_seqs++);
                cp_str(&ecat->name, s);
                ecat->type = T_MATR;
                ecat->end_pos = last_z;
                ecat->flags = 0;
                ecat->p_elem = NULL;
#if DEBUG
                fprintf(stdout, "concatenated matrix %s has order %ld\n", ecat->name, ecat->matrix->order);
                fflush(stdout);
#endif
                extend_elem_list(&ecat);
                in_seq = 0;
                n_matrices++;
                beamline->ncat_elems++;
                }
            /* non-matrix element--just copy everything and extend the list */
#if DEBUG
            fprintf(stdout, "copying non-matrix element %s\n", elem->name);
            fflush(stdout);
#endif
            pred = ecat->pred;
            succ = ecat->succ;
            memcpy((char*)ecat, (char*)elem, sizeof(*ecat));
            ecat->pred = pred;
            ecat->succ = succ;
            ecat->flags   = 0;
            if (ecat->type==T_RECIRC)
                beamline->ecat_recirc = ecat;
            if (entity_description[ecat->type].flags&HAS_MATRIX) {
                ecat->matrix = tmalloc(sizeof(*(ecat->matrix)));
                copy_matrices(ecat->matrix, elem->matrix);
                if (entity_description[ecat->type].flags&DONT_CONCAT)
                    n_nonmatrices++;
                else
                    n_matrices++;
                }
            else {
                ecat->matrix = NULL;
                n_nonmatrices++;
                }
            beamline->ncat_elems++;
            extend_elem_list(&ecat);
            new_seq = 1;
            }
        lastelem = elem;
        } while ((elem=elem->succ));

    if (in_seq) {
        /*  end of sequence--copy concatenated matrix into ecat list */
        elem = lastelem;
        ecat->matrix = tmalloc(sizeof(*(ecat->matrix)));
        copy_matrices(ecat->matrix, M1);
        sprintf(s, "M%ld", n_seqs++);
        cp_str(&ecat->name, s);
        ecat->type = T_MATR;
        ecat->end_pos = z;
        ecat->flags = 0;
        ecat->p_elem = NULL;
        ecat->twiss = lastelem->twiss;
        if (ecat->succ) {
            ecat = ecat->succ;
            if (ecat->matrix) {
                free_matrices(ecat->matrix);
                tfree(ecat->matrix);
                ecat->matrix = NULL;
                }
            if (ecat->name) {
                tfree(ecat->name);
                ecat->name = NULL;
                }
            }
        n_matrices++;
        beamline->ncat_elems++;
        ecat->succ = NULL;
        }
    else
        ecat->pred->succ = NULL;

#if DEBUG
    fprintf(stdout, "concatenated element list:\n");
    fflush(stdout);
    print_elem_list(stdout, &(beamline->ecat));
#endif

    beamline->flags |= BEAMLINE_CONCAT_CURRENT+BEAMLINE_CONCAT_DONE;
    log_exit("concatenate_beamline");
    }

void copy_matrices1(VMATRIX *M1,  VMATRIX *M0)
{
    register long i, j, k, l;

    log_entry("copy_matrices1");

    if (!M0)
        bomb("null VMATRIX pointer for source matrix (copy_matrices1)", NULL);
    if (!M1)
        bomb("null VMATRIX pointer for target matrix (copy_matrices1)", NULL);
    if (M0->order<1)
        bomb("invalid VMATRIX order for source matrix (copy_matrices1)", NULL);
    if (M1->order<1)
        bomb("invalid VMATRIX order for target matrix (copy_matrices1)", NULL);
    if (!M0->C || !M0->R)
        bomb("null C and/or R for source matrix (copy_matrices1)", NULL);
    if (!M1->C || !M1->R)
        bomb("null C and/or R for target matrix (copy_matrices1)", NULL);

    for (i=0; i<6; i++) {
        M1->C[i] = M0->C[i];
        for (j=0; j<6; j++)
            M1->R[i][j] = M0->R[i][j];
        }


    if (M1->order>=2) {
        if (!M1->T)
            bomb("null T for target matrix (copy_matrices1)", NULL);
        if (M0->order>=2) {
            if (!M0->T)
                bomb("null T for source matrix (copy_matrices1)", NULL);
            for (i=0; i<6; i++)
                for (j=0; j<6; j++)
                    for (k=0; k<=j; k++)
                        M1->T[i][j][k] = M0->T[i][j][k];
            }
        else {
            for (i=0; i<6; i++)
                for (j=0; j<6; j++)
                    for (k=0; k<=j; k++)
                        M1->T[i][j][k] = 0;
            }
        }

    if (M1->order>=3) {
        if (!M1->Q)
            bomb("null Q for target matrix (copy_matrices1)", NULL);
        if (M0->order>=3) {
            if (!M0->Q)
                bomb("null Q for source matrix (copy_matrices1)", NULL);
            for (i=0; i<6; i++)
                for (j=0; j<6; j++)
                    for (k=0; k<=j; k++)
                        for (l=0; l<=k; l++)
                            M1->Q[i][j][k][l] = M0->Q[i][j][k][l];
            }
        else {
            for (i=0; i<6; i++)
                for (j=0; j<6; j++)
                    for (k=0; k<=j; k++)
                        for (l=0; l<=k; l++)
                            M1->Q[i][j][k][l] = 0;
            }
        }
    log_exit("copy_matrices1");
    }

void free_elements1(ELEMENT_LIST *eptr)
/* similar to free_elements() only I don't free shared memory in p_elem and name
 * in concatenated element lists formed by concatenate_beamline()
 */
{
    long is_root = 1;

    log_entry("free_elements1");

    while (eptr) {
        if (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix ) {
            if (eptr->matrix->order<1 || eptr->matrix->order>3)
                bomb("invalid matrix order (free_elements1)", NULL);
/*            tfree(eptr->name); */
            free_matrices(eptr->matrix);
            tfree(eptr->matrix);
            eptr->matrix = NULL;
            eptr->name = NULL;
            }
        if (eptr->succ) {
            eptr = eptr->succ;
            if (!is_root)
                tfree(eptr->pred);
            else
                is_root = 0;
            eptr->pred = NULL;
            }
        else {
            if (!is_root) {
                tfree(eptr);
                eptr = NULL;
                }
            else
                is_root = 0;
            break;
            }
        }
    log_exit("free_elements1");
    }
