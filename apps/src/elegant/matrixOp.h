/* redefine Matrix so that meschach lib is no longer needed 
 * These versions assume that matrices are stored in column-major order, to be compatible with 
 * LAPACK SVD routines.
 * Example:
 * A = matrix_get(m, n) creates a matrix with m rows and n columns
 * However, A->me[i][j] is the ith column and jth row, 0<=i<=n, 0<=j<=m 
 * Using the macro Mij allows one to avoid confusion.
 * Mij(A, i, j) = A->me[j][i]
 */

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**                          Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


#include "SDDS.h"

#ifndef TRMATRIX_H
#define TRMATRIX_H

/* if nothing is defined, choose double precision */
#ifndef REAL_DBL
#ifndef REAL_FLT
#define REAL_DBL 1
#endif
#endif

/* single precision */
#ifdef REAL_FLT
#define  Real float
#define  LongReal float
#define REAL FLOAT
#define LONGREAL FLOAT
#endif

/* double precision */
#ifdef REAL_DBL
#define Real double
#define LongReal double
#define REAL DOUBLE
#define LONGREAL DOUBLE
#endif

typedef	struct	{
  unsigned int	dim, max_dim;
  Real	*ve;
} VEC;

/*matrix is in column major order */
typedef	struct	{
  unsigned int	m, n;
  unsigned int	max_m, max_n, max_size;
  Real	**me,*base;	/* base is base of alloc'd mem */
} MAT;
 
/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),(size_t)sizeof(type)))

/* get  vector/matrix to given dimension */
extern	VEC *vec_get(int);
extern	MAT *matrix_get(int,int); 
extern	int matrix_free(MAT *),vec_free(VEC *);
extern  MAT *matrix_copy(MAT *);
/* A is the input*/
extern  MAT *matrix_transpose(MAT *A); 

/*in following routines, mat1 and mat2 are input, no change */
extern  MAT *matrix_add(MAT *mat1, MAT *mat2);
extern  MAT *matrix_sub(MAT *mat1, MAT *mat2);
extern  MAT *matrix_h_mult(MAT *mat1, MAT *mat2);
extern  MAT *matrix_h_divide(MAT *mat1, MAT *mat2);
extern  MAT *matrix_mult(MAT *mat1, MAT *mat2);

/*in following routines, mat1 is input and output, changed; mat2 is input, no change*/
extern  int32_t matrix_add_sm(MAT *mat1, MAT *mat2); 
extern  int32_t matrix_sub_sm(MAT *mat1, MAT *mat2);
extern  int32_t matrix_h_mult_sm(MAT *mat1, MAT *mat2);
extern  int32_t matrix_h_divide_sm(MAT *mat1, MAT *mat2);
extern  void matrix_scmul(MAT *mat1, double scalar);

/*invert routine*/
extern  MAT *matrix_invert(MAT *A, double *weight, int32_t largestSValue, int32_t smallestSValue, double minRatio,
                   int32_t deleteVectors, int32_t *deleteVector, char **deletedVector,
                   VEC **S_Vec, int32_t *sValues,
                   VEC **S_Vec_used, int32_t *usedSValues,
                   MAT **U_matrix, MAT **Vt_matrix,
                   double *conditionNum);

/*generate identity matrix*/
MAT *matrix_identity(int32_t m, int32_t n);

/*compute matrix determinant */
double matrix_det(MAT *D);

/*read column major ordered matrix array*/
void *SDDS_GetCastMatrixOfRowsByColumn(SDDS_DATASET *SDDS_dataset, int32_t *n_rows, long sddsType);

#include <stdio.h>
extern void matrix_show(MAT *A, char *format, char *label, FILE *fp);
extern int matrix_check(MAT *A);


/* access matrices with indices in the expected order, i.e., a[row][column] */

#define Mij(A, i, j) A->me[j][i]
 
#endif
