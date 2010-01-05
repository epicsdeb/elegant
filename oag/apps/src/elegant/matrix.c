/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: matrix.c--3rd-order matrix routines for tracking.
 * contents:  print_matrices(), initialize_matrices(), null_matrices(),
 *            free_matrices(), track_particles(), 
 *
 * It is assumed that Tijk and Tikj for k<j have been combined into 
 * one element, namely, Tijk, which has been multiplied by 2.
 * Similarly:
 *    -- since Qijkk==Qikjk==Qikkj, for j>k, these elements have been combined
 *       as 3*Qijkk.
 *    -- since Qijjk==Qijkj==Qikjj, for j>k, these elements have been combined
 *       as 3*Qijjk.
 *    -- since Qijkl==Qijlk, for j>k and k>l, these elements have been combined
 *       as 6*Qijkl.
 *
 * Michael Borland, 1989.
 */
#include "mdb.h"
#include "track.h"

void print_matrices(FILE *fp, char *string, VMATRIX *M)
{
    register long i, j, k, l;
    double *C;
    double **R;
    double ***T;
    double ****Q;

    log_entry("print_matrices");

    set_matrix_pointers(&C, &R, &T, &Q, M);

    fprintf(fp, "%s\nC:   ", string);
    for (i=0; i<6; i++) 
        fprintf(fp, "%22.15e ", C[i]);
    fputc('\n', fp);

    for (i=0; i<6; i++) {
        fprintf(fp, "R%ld: ", i+1);
        for (j=0; j<6; j++)
            fprintf(fp, "%22.15e ", R[i][j]);
        fputc('\n', fp);
        }

    if (M->order>=2) {
        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                fprintf(fp, "T%ld%ld: ", i+1, j+1);
                for (k=0; k<=j; k++) 
                    fprintf(fp, "%22.15e ", T[i][j][k]);
                fputc('\n', fp);
                }
            }
        }

    if (M->order>=3) {
        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                for (k=0; k<=j; k++) {
                    fprintf(fp, "Q%ld%ld%ld: ", i+1, j+1, k+1);
                    for (l=0; l<=k; l++)
                        fprintf(fp, "%22.15e ", Q[i][j][k][l]);
                    fputc('\n', fp);
                    }
                }
            }
        }
    log_exit("print_matrices");
    }


void initialize_matrices(VMATRIX *M, long order)
{
    register long i, j, k, l;
    double *Tij, **Qij, *Qijk;
    double *C, **R;
    double ***T;
    double ****Q;

    log_entry("initialize_matrices");

    M->eptr = NULL;
    
    switch (M->order = order) {
        case 3:
            M->C = C = tmalloc(sizeof(*C)*6);
            M->R = R = tmalloc(sizeof(*R)*6);
            M->T = T = tmalloc(sizeof(*T)*6);
            M->Q = Q = tmalloc(sizeof(*Q)*6);
            for (i=0; i<6; i++) {
                C[i] = 0;
                R[i] = tmalloc(sizeof(**R)*6);
                T[i] = tmalloc(sizeof(**T)*6);
                Q[i] = tmalloc(sizeof(**Q)*6);
                for (j=0; j<6; j++) {
                    R[i][j] = 0;
                    Tij = T[i][j] = tmalloc(sizeof(***T)*(j+1));
                    Qij = Q[i][j] = tmalloc(sizeof(***Q)*(j+1));
                    for (k=0; k<=j; k++) {
                        *Tij++ = 0;
                        Qijk = *Qij++ = tmalloc(sizeof(****Q)*(k+1));
                        for (l=0; l<=k; l++)
                            *Qijk++ = 0;
                        }
                    }
                }
            break;
        case 2:
            M->C = C = tmalloc(sizeof(*C)*6);
            M->R = R = tmalloc(sizeof(*R)*6);
            M->T = T = tmalloc(sizeof(*T)*6);
            M->Q = NULL;
            for (i=0; i<6; i++) {
                C[i] = 0;
                R[i] = tmalloc(sizeof(**R)*6);
                T[i] = tmalloc(sizeof(**T)*6);
                for (j=0; j<6; j++) {
                    R[i][j] = 0;
                    Tij = T[i][j] = tmalloc(sizeof(***T)*(j+1));
                    for (k=0; k<=j; k++) {
                        *Tij++ = 0;
                        }
                    }
                }
            break;
        case 1:
            M->C = C = tmalloc(sizeof(*C)*6);
            M->R = R = tmalloc(sizeof(*R)*6);
            M->T = NULL;
            M->Q = NULL;
            for (i=0; i<6; i++) {
                C[i] = 0;
                R[i] = tmalloc(sizeof(**R)*6);
                for (j=0; j<6; j++) {
                    R[i][j] = 0;
                    }
                }
            break;
        default:
            fprintf(stdout, "invalid order: %ld  (initialize_matrices)\n", 
                order);
            fflush(stdout);
            exit(1);
            break;
        }
    log_exit("initialize_matrices");
    }


void null_matrices(VMATRIX *M, unsigned long flags)
{
    register long i, j, k, l;
    double *Tij, **Qij, *Qijk;
    double *C, **R, ***T, ****Q;

    M->eptr = NULL;
    
    set_matrix_pointers(&C, &R, &T, &Q, M);

    if (M->order==3 && !(flags&EXCLUDE_Q) && Q)
        for (j=0; j<6; j++) {
          Qij = Q[i][j];
          for (k=0; k<=j; k++) {
            Qijk   = *Qij++;
            for (l=0; l<=k; l++) 
              *Qijk++ = 0;
          }
        }
    if (M->order>=2 && !(flags&EXCLUDE_T) && T)
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          Tij = T[i][j];
          for (k=0; k<=j; k++)
            *Tij++ = 0;
        }
      }
    if (M->order>=1 && !(flags&EXCLUDE_R) && R)
      for (i=0; i<6; i++) 
        for (j=0; j<6; j++) 
          R[i][j] = (flags&SET_UNIT_R) && (i==j) ? 1 : 0;

    if (!(flags&EXCLUDE_C) && C)
      for (i=0; i<6; i++) 
        C[i] = 0;
  }


void track_particles(double **final, VMATRIX *M, double  **initial, long n_part)
{
    double sum1;
    double *ini_k, *Tij, **Ti, sum, coord_j;
    long k, j;
    long i, l, i_part;
    double coord_k, coord_jk;
    double **Qij, *Qijk;
    double *C, **R, ***T, ****Q;
    double *fin, *ini;
    double *Ri;
    static double temp[6];
    
    log_entry("track_particles");
    
    if (!M)
        bomb("NULL VMATRIX pointer in track_particles", NULL);
#if !SDDS_MPI_IO
    if (!final)
        bomb("NULL final coordinates pointer in track_particles", NULL);
    if (!initial)
        bomb("NULL initial coordinates pointer in track_particles", NULL);
#endif

    set_matrix_pointers(&C, &R, &T, &Q, M);

    switch (M->order) {
      case 3:
        if (!C)
            bomb("NULL C pointer (track_particles)", NULL);
        if (!R)
            bomb("NULL R pointer (track_particles)", NULL);
        if (!T)
            bomb("NULL T pointer (track_particles)", NULL);
        if (!Q)
            bomb("NULL Q pointer (track_particles)", NULL);
        for (i_part=n_part-1; i_part>=0; i_part--) {
            if (!(fin = final[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            if (!(ini = initial[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            fin[6] = ini[6];       /* copy particle ID # */
            for (i=5; i>=0; i--) {
                sum = C[i];
                for (j=5; j>=0; j--) {
                    if ((coord_j=ini[j])!=0) {
                        sum += R[i][j]*coord_j;
                        Tij  = T[i][j]+j;
                        Qij  = Q[i][j]+j;
                        for (k=j; k>=0; k--, Tij--) {
                            Qijk = *(Qij--)+k;
                            if ((coord_k=ini[k])!=0) {
                                coord_jk = coord_j*coord_k;
                                sum += *Tij*coord_jk;
                                for (l=k; l>=0; l--) 
                                    sum += *(Qijk--)*coord_jk*ini[l];
                                }
                            }
                        }
                    }
                temp[i] = sum; /* to prevent changing initial values in 
                                  the event initial and final are same*/
                }
            for (i=5; i>=0; i--)
              fin[i] = temp[i];
            }
        break;
      case 2:
        if (!C)
            bomb("NULL C pointer (track_particles)", NULL);
        if (!R)
            bomb("NULL R pointer (track_particles)", NULL);
        if (!T)
            bomb("NULL T pointer (track_particles)", NULL);
        for (i_part=n_part-1; i_part>=0; i_part--) {
            if (!(fin = final[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            if (!(ini = initial[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            fin[6] = ini[6]; /* copy particle ID # */
            for (i=5; i>=0; i--) {
                sum = C[i];
                if (!(Ri  = R[i]+5))
                    bomb("NULL R[i] pointer (track_particles)", NULL);
                if (!(Ti  = T[i]+5))
                    bomb("NULL T[i] pointer (track_particles)", NULL);
                for (j=5; j>=0; j--, Ri--, Ti--) {
                    if ((coord_j= *(ini_k=ini+j))) {
                        sum1 = *Ri;
                        if (!(Tij  = *Ti+j))
                            bomb("NULL T[i][j] pointer (tracking_particles)", NULL);
                        for (k=j; k>=0; k--, Tij--, ini_k--)
                            sum1 += *Tij**ini_k;
                        sum += sum1*coord_j;
                        }
                    }
                temp[i] = sum;
                }
            for (i=5; i>=0; i--)
                fin[i] = temp[i];
            }
        break;
      case 1:
        if (!C)
            bomb("NULL C pointer (track_particles)", NULL);
        if (!R)
            bomb("NULL R pointer (track_particles)", NULL);
        for (i_part=n_part-1; i_part>=0; i_part--) {
            if (!(fin = final[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            if (!(ini = initial[i_part])) {
                fprintf(stdout, "error: final coordinate pointer is NULL for particle %ld (track_particles)\n", i_part);
                fflush(stdout);
                abort();
                }
            fin[6] = ini[6]; /* copy particle ID # */
            for (i=5; i>=0; i--) {
                sum = C[i];
                if (!(Ri=R[i]+5))
                    bomb("NULL R[i] pointer (track_particles)", NULL);
                for (j=5; j>=0; Ri--, j--)
                    sum += *Ri * ini[j];
                temp[i] = sum;
                }
            for (i=5; i>=0; i--)
                fin[i] = temp[i];
            }
        break;
      default:
        fprintf(stdout, "invalid order: %ld  (track_particle)\n", 
               M->order);
        fflush(stdout);
        exit(1);
        break;
        }
    
    log_exit("track_particles");
    }

void free_matrices(VMATRIX *M)
{
    register long i, j, k;
    double **Qij;
    double *C, **R;
    double ***T;
    double ****Q;

    log_entry("free_matrices");
    if (!M)
        bomb("NULL matrix passed to free_matrices", NULL);
    
    set_matrix_pointers(&C, &R, &T, &Q, M);
    switch (M->order) {
        case 3:
            if (!Q || !T || !R || !C)
                bomb("NULL Q, T, R, or C entry for matrix (free_matrices)", NULL);
            for (i=0; i<6; i++) {
                if (!R[i])
                    bomb("NULL R[i] entry for matrix (free_matrices)", NULL);
                if (!T[i])
                    bomb("NULL T[i] entry for matrix (free_matrices)", NULL);
                if (!Q[i])
                    bomb("NULL Q[i] entry for matrix (free_matrices)", NULL);
                for (j=0; j<6; j++) {
                    if (!(Qij = Q[i][j]))
                        bomb("NULL Q[i][j] entry for matrix (free_matrices)", NULL);
                    for (k=0; k<=j; k++) {
                        if (!*Qij)
                            bomb("NULL Q[i][j][k] entry for matrix (free_matrices)", NULL);
                        tfree(*Qij++);
                        }
                    if (!T[i][j])
                        bomb("NULL T[i][j] entry for matrix (free_matrices)", NULL);
                    if (!Q[i][j])
                        bomb("NULL Q[i][j] entry for matrix (free_matrices)", NULL);
                    tfree(T[i][j]); T[i][j] = NULL;
                    tfree(Q[i][j]); Q[i][j] = NULL;
                    }
                tfree(R[i]); R[i] = NULL;
                tfree(T[i]); T[i] = NULL;
                tfree(Q[i]); Q[i] = NULL;
                }
            tfree(C);
            tfree(R);
            tfree(T);
            tfree(Q);
            break;
        case 2:
            if (!T || !R || !C)
                bomb("NULL T, R, or C entry for matrix (free_matrices)", NULL);
            for (i=0; i<6; i++) {
                if (!R[i])
                    bomb("NULL R[i] entry for matrix (free_matrices)", NULL);
                if (!T[i])
                    bomb("NULL T[i] entry for matrix (free_matrices)", NULL);
                for (j=0; j<6; j++) {
                    if (!T[i][j])
                        bomb("NULL T[i][j] entry for matrix (free_matrices)", NULL);
                    tfree(T[i][j]); T[i][j] = NULL;
                    }
                tfree(R[i]); R[i] = NULL;
                tfree(T[i]); T[i] = NULL;
                }
            tfree(C);
            tfree(R);
            tfree(T);
            break;
        case 1:
            if (!R || !C)
                bomb("NULL R or C entry for matrix (free_matrices)", NULL);
            for (i=0; i<6; i++) {
                if (!R[i])
                    bomb("NULL R[i] entry for matrix (free_matrices)", NULL);
                tfree(R[i]); R[i] = NULL;
                }
            tfree(C);
            tfree(R);
            break;
        default:
            fprintf(stdout, "invalid order: %ld  (free_matrices)\n", 
                M->order);
            fflush(stdout);
            exit(1);
            break;
        }
    M->C = NULL;
    M->R = NULL;
    M->T = NULL;
    M->Q = NULL;
    log_exit("free_matrices");
    }

void free_nonlinear_matrices(VMATRIX *M)
{
    register long i, j, k;
    double **Qij;
    double *C, **R;
    double ***T;
    double ****Q;

    log_entry("free_nonlinear_matrices");
    
    set_matrix_pointers(&C, &R, &T, &Q, M);
    switch (M->order) {
        case 3:
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    Qij = Q[i][j];
                    for (k=0; k<=j; k++) {
                        tfree(*Qij++);
                        }
                    tfree(T[i][j]); T[i][j] = NULL;
                    tfree(Q[i][j]); Q[i][j] = NULL;
                    }
                tfree(T[i]); T[i] = NULL;
                tfree(Q[i]); Q[i] = NULL;
                }
            tfree(T);
            tfree(Q);
            break;
        case 2:
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    tfree(T[i][j]); T[i][j] = NULL;
                    }
                tfree(T[i]); T[i] = NULL;
                }
            tfree(T);
            break;
        case 1:
            break;
        default:
            fprintf(stdout, "invalid order: %ld  (free_matrices)\n", 
                M->order);
            fflush(stdout);
            exit(1);
            break;
        }
    M->T = NULL;
    M->Q = NULL;
    M->order = 1;
    log_exit("free_nonlinear_matrices");
    }

void set_matrix_pointers(double **C, double ***R, double ****T, double *****Q, VMATRIX *M)
{
    log_entry("set_matrix_pointers");
    if (!M) {
        fprintf(stdout, "error: NULL VMATRIX pointer\n");
        fflush(stdout);
        abort();
        }
    *C = M->C;
    *R = M->R;
    *T = M->T;
    *Q = M->Q;
    log_exit("set_matrix_pointers");
    }

long read_matrices(VMATRIX *M, FILE *fp)
{
    long order, i, j, k, l;
    char s[256], *ptr;
    double *C, **R, ***T, ****Q;

    log_entry("read_matrices");

    set_matrix_pointers(&C, &R, &T, &Q, M);
    order = M->order;

    if (order>=1) {
        if (!fgets(s, 256, fp) || !(ptr=strchr(s, ':'))) {
            log_exit("read_matrices");
            return(0);
            }
        for (i=0; i<6; i++)
            if (!get_double(C+i, ptr)) {
                log_exit("read_matrices");
                return(0);
                }
        for (i=0; i<6; i++) {
            if (!fgets(s, 256, fp) || !(ptr=strchr(s, ':'))) {
                log_exit("read_matrices");
                return(0);
                }
            for (j=0; j<6; j++)
                if (!get_double(R[i]+j, ptr)) {
                    log_exit("read_matrices");
                    return(0);
                    }
            }
        }

    if (order>=2) {
        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                if (!fgets(s, 256, fp) || !(ptr=strchr(s, ':'))) {
                    log_exit("read_matrices");
                    return(1);
                    }
                for (k=0; k<=j; k++) 
                    if (!get_double(T[i][j]+k, ptr)) {
                        log_exit("read_matrices");
                        return(1);
                        }
                }
            }
        }

    if (order>=3) {
        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                for (k=0; k<=j; k++) {
                    if (!fgets(s, 256, fp) || !(ptr=strchr(s, ':'))) {
                        log_exit("read_matrices");
                        return(2);
                        }
                    for (l=0; l<=k; l++) 
                        if (!get_double(Q[i][j][k]+l, ptr)) {
                            log_exit("read_matrices");
                            return(2);
                            }
                    }
                }
            }
        }

    log_exit("read_matrices");
    return(order);
    }

void filter_matrices(VMATRIX *M, double threshold)
{
    register long i, j, k, l;
    double *Tij, **Qij, *Qijk;
    double *C, **R, ***T, ****Q;

    log_entry("filter_matrices");

    set_matrix_pointers(&C, &R, &T, &Q, M);

    switch (M->order) {
        case 3:
            for (i=0; i<6; i++) {
                if (fabs(C[i])<threshold)
                    C[i] = 0;
                for (j=0; j<6; j++) {
                    if (fabs(R[i][j])<threshold)
                        R[i][j] = 0;
                    Tij = T[i][j];
                    Qij = Q[i][j];
                    for (k=0; k<=j; k++) {
                        if (fabs(Tij[k])<threshold)
                            Tij[k] = 0;
                        Qijk   = *Qij++;
                        for (l=0; l<=k; l++) 
                            if (fabs(Qijk[l])<0)
                                Qijk[l] = 0;
                        }
                    }
                }
            break;
        case 2:
            for (i=0; i<6; i++) {
                if (fabs(C[i])<threshold)
                    C[i] = 0;
                for (j=0; j<6; j++) {
                    if (fabs(R[i][j])<threshold)
                        R[i][j] = 0;
                    Tij = T[i][j];
                    for (k=0; k<=j; k++) {
                        if (fabs(Tij[k])<threshold)
                            Tij[k] = 0;
                        }
                    }
                }
            break;
        case 1:
            for (i=0; i<6; i++) {
                if (fabs(C[i])<threshold)
                    C[i] = 0;
                for (j=0; j<6; j++) {
                    if (fabs(R[i][j])<threshold)
                        R[i][j] = 0;
                    }
                }
            break;
        default:
            fprintf(stdout, "invalid order: %ld  (filter_matrices)\n", 
                M->order);
            fflush(stdout);
            exit(1);
            break;
        }
    log_exit("filter_matrices");
    }

void random_matrices(VMATRIX *M, double C0, double R0, double T0, double Q0)
{
    register long i, j, k, l;
    double *Tij, **Qij, *Qijk;
    double *C, **R, ***T, ****Q;

    log_entry("random_matrices");

    set_matrix_pointers(&C, &R, &T, &Q, M);

    switch (M->order) {
        case 3:
            for (i=0; i<6; i++) {
                C[i] = C0*(2*random_1_elegant(1)-1);
                for (j=0; j<6; j++) {
                    R[i][j] = R0*(2*random_1_elegant(1)-1);
                    Tij = T[i][j];
                    Qij = Q[i][j];
                    for (k=0; k<=j; k++) {
                        *Tij++ = T0*(2*random_1_elegant(1)-1);
                        Qijk   = *Qij++;
                        for (l=0; l<=k; l++) 
                            *Qijk++ = Q0*(2*random_1_elegant(1)-1);
                        }
                    }
                }
            break;
        case 2:
            for (i=0; i<6; i++) {
                C[i] = C0*(2*random_1_elegant(1)-1);
                for (j=0; j<6; j++) {
                    R[i][j] = R0*(2*random_1_elegant(1)-1);
                    Tij = T[i][j];
                    for (k=0; k<=j; k++) {
                        *Tij++ = T0*(2*random_1_elegant(1)-1);
                        }
                    }
                }
            break;
        case 1:
            for (i=0; i<6; i++) {
                C[i] = C0*(2*random_1_elegant(1)-1);
                for (j=0; j<6; j++) {
                    R[i][j] = R0*(2*random_1_elegant(1)-1);
                    }
                }
            break;
        default:
            fprintf(stdout, "invalid order: %ld  (null_matrices)\n", 
                M->order);
            fflush(stdout);
            exit(1);
            break;
        }
    log_exit("random_matrices");
    }

void copy_matrices(VMATRIX *M1, VMATRIX *M0)
{
    long i, j, k, l;

    log_entry("copy_matrices");

    initialize_matrices(M1, M1->order=M0->order);

    for (i=0; i<6; i++) {
        M1->C[i] = M0->C[i];
        for (j=0; j<6; j++)
            M1->R[i][j] = M0->R[i][j];
        }

    if (M1->order>=2) {
        for (i=0; i<6; i++)
            for (j=0; j<6; j++)
                for (k=0; k<=j; k++)
                    M1->T[i][j][k] = M0->T[i][j][k];
        }

    if (M1->order>=3) {
        for (i=0; i<6; i++)
            for (j=0; j<6; j++)
                for (k=0; k<=j; k++)
                    for (l=0; l<=k; l++)
                        M1->Q[i][j][k][l] = M0->Q[i][j][k][l];
        }

    log_exit("copy_matrices");
    }

long check_matrix(VMATRIX *M, char *comment)
{
    long i, j, k;

    log_entry("check_matrix");

    if (M==NULL) {
        fprintf(stdout, "error: NULL matrix pointer---%s\n", comment);
        fflush(stdout);
        abort();
        }
    if (M->order<=0 || M->order>3) {
        fprintf(stdout, "error: matrix order out of range---%s\n", comment);
        fflush(stdout);
        abort();
        }
    if (M->R==NULL) {
        fprintf(stdout, "error: NULL R matrix---%s\n", comment);
        fflush(stdout);
        abort();
        }
    for (i=0; i<6; i++) {
        if (M->R[i]==NULL)
            fprintf(stdout, "error: NULL R[%ld] row---%s\n", i, comment);
            fflush(stdout);
        }
    if (M->order==1) {
        log_exit("check_matrix");
        return(1);
        }
    if (M->T==NULL) {
        fprintf(stdout, "error: NULL Tmatrix---%s\n", comment);
        fflush(stdout);
        abort();
        }
    for (i=0; i<6; i++) {
        if (M->T[i]==NULL) {
            fprintf(stdout, "error: NULL T[%ld] row---%s\n", i, comment);
            fflush(stdout);
            abort();
            }
        for (j=0; j<6; j++) {
            if (M->T[i][j]==NULL) {
                fprintf(stdout, "error: NULL T[%ld][%ld] row---%s\n", i, j, comment);
                fflush(stdout);
                abort();
                }
            }
        }
    if (M->order==2) {
        log_exit("check_matrix");
        return(2);
        }
    if (M->Q==NULL) {
        fprintf(stdout, "error: NULL Q matrix---%s\n", comment);
        fflush(stdout);
        abort();
        }
    for (i=0; i<6; i++) {
        if (M->Q[i]==NULL) {
            fprintf(stdout, "error: NULL T[%ld] row---%s\n", i, comment);
            fflush(stdout);
            abort();
            }
        for (j=0; j<6; j++) {
            if (M->Q[i][j]==NULL) {
                fprintf(stdout, "error: NULL T[%ld][%ld] row---%s\n", i, j, comment);
                fflush(stdout);
                abort();
                }
            for (k=0; k<=j; k++) {
                if (M->Q[i][j][k]==NULL) {
                    fprintf(stdout, "error: NULL Q[%ld][%ld][%ld] row---%s\n", i, j, k, comment);
                    fflush(stdout);
                    abort();
                    }
                }
            }
        }
    log_exit("check_matrix");
    return(3);
    }


