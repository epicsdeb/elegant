/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: concat_matrices()
 * 
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"


#define MULT_T(j,k) (j==k?2:1)
#define MULT_Q(m, n, p) (m==n ? (n==p?6:2) : (m==p ? 2: (n==p?2:1)))
#define SYMM_R_T(R, T, j, k, m, n, p)   (R[j][m]*T[k][n][p] + R[k][m]*T[j][n][p])
#define SYMM_C_Q(C, Q, j, k, m, n, p)   (C[j]*Q[k][m][n][p] + C[k]*Q[j][m][n][p])
#define SYMM_C_C_Q(C, Q, j, k, l, n, m, p)  (n>=m && m>=p? \
    C[j]*C[k]*Q[l][n][m][p] + C[j]*C[l]*Q[k][n][m][p] +\
    C[l]*C[k]*Q[j][n][m][p]  : 0)
#define SYMM_C_C_T(C, T, j, k, l, n, m) (n>=m ? \
    C[j]*C[k]*T[l][n][m] + C[j]*C[l]*T[k][n][m] +\
    C[k]*C[l]*T[j][n][m] : 0)
#define SYMM_C_R_T(C, R, T, j, k, l, n, m, p) (m>=p? \
    C[j]*R[k][n]*T[l][m][p] + C[k]*R[l][n]*T[j][m][p] + C[l]*R[j][n]*T[k][m][p] + \
    C[j]*R[l][n]*T[k][m][p] + C[k]*R[j][n]*T[l][m][p] + C[l]*R[k][n]*T[j][m][p] : 0)

void concat_matrices(VMATRIX *M2, VMATRIX *M1, VMATRIX *M0, unsigned long mode)
/* M2 = M1*M0.  M0 is applied to particles first. */
{
  double *C2, **R2, ***T2, ****Q2;
  double *C1, **R1, ***T1, ****Q1;
  double *C0, **R0, ***T0, ****Q0;
  double M0C4;
  long order;
  long i, j, k, l, m, n, p;
  double sum;

  log_entry("concat_matrices");

  set_matrix_pointers(&C2, &R2, &T2, &Q2, M2);
  set_matrix_pointers(&C1, &R1, &T1, &Q1, M1);
  set_matrix_pointers(&C0, &R0, &T0, &Q0, M0);

  check_matrix(M0, "M0 matrix passed to concat_matrices");
  check_matrix(M1, "M1 matrix passed to concat_matrices");
  check_matrix(M2, "M2 matrix passed to concat_matrices");
  
  order = M2->order;
  if (M1->order<1 || M0->order<1 || order<1)
    bombElegant("order<1 in concat_matrices", NULL);

  M0C4 = C0[4];
  if (mode&CONCAT_EXCLUDE_S0) {
    /* necessary with rf elements because path length is not differential quantity */
    C0[4] = 0; 
  }
  
  /* calculate new zero-th order terms */
  for (i=0; i<6; i++) {
    /* sum up contributions to C[i]: */
    for (j=sum=0; j<6; j++) {
      sum += R1[i][j]*C0[j];
    }
    C2[i] = sum + C1[i];
  }

  /* calculate new zero-th order terms */
  if (M1->order>1) {
    for (i=0; i<6; i++) {
      /* sum up contributions to C[i]: */
      for (j=sum=0; j<6; j++) {
        for (k=0; k<=j; k++) {
          sum += T1[i][j][k]*C0[j]*C0[k];
          if (M1->order>2) {
            for (l=0; l<=k; l++)
              sum += Q1[i][j][k][l]*C0[j]*C0[k]*C0[l];
          }
        }
      }
      C2[i] += sum;
    }
  }

  if (mode&CONCAT_EXCLUDE_S0)
    C2[4] += M0C4;
  
  /* calculate new first order terms */
  for (i=0; i<6; i++) {
    for (m=0; m<6; m++) {
      /* sum up contributions to R[i][m]: */
      for (j=sum=0; j<6; j++) 
        sum += R1[i][j]*R0[j][m];
      R2[i][m] = sum;
    }
  }
  if (M1->order>1) {
    for (i=0; i<6; i++) {
      for (m=0; m<6; m++) {
        /* sum up contributions to R[i][m]: */
        for (j=sum=0; j<6; j++) {
          for (k=0; k<=j; k++) {
            sum += T1[i][j][k]*(C0[j]*R0[k][m]+C0[k]*R0[j][m]);
            if (M1->order>2) {
              for (l=0; l<=k; l++) {
                sum += Q1[i][j][k][l]*(
                                       C0[j]*C0[k]*R0[l][m] +
                                       C0[k]*C0[l]*R0[j][m] +
                                       C0[l]*C0[j]*R0[k][m] );
              }
            }
          }
        }
        R2[i][m] += sum;
      }
    }
  }

  /* calculate new second order terms */
  if (order>1) {
    for (i=0; i<6; i++) {
      for (m=0; m<6; m++) {
        for (n=0; n<=m; n++) {
          /* sum up contributions to T[i][m][n]: */
          for (j=sum=0; j<6; j++) {
            if (M0->order>1)
              sum += R1[i][j]*T0[j][m][n]*MULT_T(m, n);
            if (M1->order>1) {
              for (k=0; k<=j; k++) {
                sum += T1[i][j][k]*(R0[j][m]*R0[k][n] + R0[k][m]*R0[j][n]);
                if (M0->order>1) 
                  sum += T1[i][j][k]*(C0[j]*T0[k][m][n] + C0[k]*T0[j][m][n])*MULT_T(m, n);
                if (M1->order>2) {
                  for (l=0; l<=k; l++) {
                    sum += Q1[i][j][k][l]*(
                                           C0[j]*R0[k][m]*R0[l][n] + 
                                           C0[j]*R0[k][n]*R0[l][m] + 
                                           C0[k]*R0[l][m]*R0[j][n] + 
                                           C0[k]*R0[l][n]*R0[j][m] + 
                                           C0[l]*R0[j][m]*R0[k][n] + 
                                           C0[l]*R0[j][n]*R0[k][m] );
                    if (M0->order>1) 
                      sum += Q1[i][j][k][l]*(
                                             SYMM_C_C_T(C0, T0, j, k, l, n, m) +
                                             SYMM_C_C_T(C0, T0, j, k, l, m, n) );
                  }
                }
              }
            }
          }
          T2[i][m][n] = sum/MULT_T(m, n);
        }
      }
    }
  }

  /* calculate new third order terms */
  if (order>2) {
    for (i=0; i<6; i++) {
      for (n=0; n<6; n++) {
        for (m=0; m<=n; m++) {
          for (p=0; p<=m; p++) {
            /* sum up contributions to Q[i][m][n][p]: */
            sum = 0;
            if (M0->order>1 && M1->order>1) {
              for (j=0; j<6; j++) {
                for (k=0; k<=j; k++) {
                  sum += T1[i][j][k]*(
                                      (n>=p?SYMM_R_T(R0, T0, j, k, m, n, p):0) +
                                      (p>=n?SYMM_R_T(R0, T0, j, k, m, p, n):0) +
                                      (m>=p?SYMM_R_T(R0, T0, j, k, n, m, p):0) +
                                      (p>=m?SYMM_R_T(R0, T0, j, k, n, p, m):0) +
                                      (m>=n?SYMM_R_T(R0, T0, j, k, p, m, n):0) +
                                      (n>=m?SYMM_R_T(R0, T0, j, k, p, n, m):0)
                                      );
                }
              }
            }
            if (M1->order>2) {
              for (j=0; j<6; j++) {
                for (k=0; k<=j; k++) {
                  for (l=0; l<=k; l++) {
                    sum += Q1[i][j][k][l]*(
                                           R0[j][m]*R0[k][n]*R0[l][p] +
                                           R0[j][m]*R0[k][p]*R0[l][n] +
                                           R0[j][n]*R0[k][m]*R0[l][p] +
                                           R0[j][n]*R0[k][p]*R0[l][m] +
                                           R0[j][p]*R0[k][m]*R0[l][n] +
                                           R0[j][p]*R0[k][n]*R0[l][m] );
                    if (M0->order>1) {
                      sum += Q1[i][j][k][l]*(
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, n, m, p) +
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, n, p, m) +
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, m, p, n) +
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, m, n, p) +
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, p, n, m) +
                                             SYMM_C_R_T(C0, R0, T0, j, k, l, p, m, n) );
                    }
                    if (M0->order>2) {
                      sum += Q1[i][j][k][l]*(
                                             SYMM_C_C_Q(C0, Q0, j, k, l, n, m, p) +
                                             SYMM_C_C_Q(C0, Q0, j, k, l, n, p, m) +
                                             SYMM_C_C_Q(C0, Q0, j, k, l, m, n, p) +
                                             SYMM_C_C_Q(C0, Q0, j, k, l, m, p, n) +
                                             SYMM_C_C_Q(C0, Q0, j, k, l, p, m, n) +
                                             SYMM_C_C_Q(C0, Q0, j, k, l, p, n, m) );
                    }
                  }
                }
              }
            }
            if (M0->order>2 && M1->order>1) {
              for (j=0; j<6; j++) {
                for (k=0; k<=j; k++) {
                  sum += T1[i][j][k]*(
                                      (n>=m && m>=p ? SYMM_C_Q(C0, Q0, j, k, n, m, p) : 0) +
                                      (n>=p && p>=m ? SYMM_C_Q(C0, Q0, j, k, n, p, m) : 0) +
                                      (m>=n && n>=p ? SYMM_C_Q(C0, Q0, j, k, m, n, p) : 0) +
                                      (m>=p && p>=n ? SYMM_C_Q(C0, Q0, j, k, m, p, n) : 0) +
                                      (p>=n && n>=m ? SYMM_C_Q(C0, Q0, j, k, p, n, m) : 0) +
                                      (p>=m && m>=n ? SYMM_C_Q(C0, Q0, j, k, p, m, n) : 0) );
                }
              }
            }
            sum /= MULT_Q(n, m, p);
            if (M0->order>2) {
              for (j=0; j<6; j++)
                sum += R1[i][j]*Q0[j][n][m][p];
            }
            Q2[i][n][m][p] = sum;
          }
        }
      }
    }
  }

  if (mode&CONCAT_EXCLUDE_S0)
    C0[4] = M0C4;

  log_exit("concat_matrices");
}

