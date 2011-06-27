/* GSL
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jorma Olavi Tähtinen, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 */

#include "mdb.h"

typedef struct
{
  double dat[2];
}
gsl_complex;

struct gsl_sf_result_struct {
  double val;
  double err;
};
struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};


typedef gsl_complex fcomplex; 
typedef unsigned int gsl_mode_t;
typedef struct gsl_sf_result_struct gsl_sf_result;
typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);
typedef void gsl_stream_handler_t (const char * label, const char * file,
                                   int line, const char * reason);
typedef struct cheb_series_struct cheb_series;


#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif
#define GSL_SUCCESS  0
#define GSL_EDOM  0
#define GSL_PREC_DOUBLE  0
#define GSL_EUNDRFLW  15
#define GSL_EOVRFLW  16
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_REAL(z)     ((z).dat[0])
#define GSL_IMAG(z)     ((z).dat[1])
#define GSL_SET_COMPLEX(zp,x,y) do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)
#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)
#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;
#define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);
#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)

#if !defined(M_PI)
#define M_PI   3.14159265358979323846
#define M_PI_4 0.78539816339744830962
#endif


gsl_complex gsl_complex_mul_real (gsl_complex a, double x);
gsl_complex gsl_complex_rect (double x, double y);
gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);
double gsl_sf_airy_Ai(const double x, gsl_mode_t mode);
double gsl_sf_airy_Bi(const double x, gsl_mode_t mode);
double gsl_sf_airy_Ai_deriv(const double x, gsl_mode_t mode);
double gsl_sf_airy_Bi_deriv(const double x, gsl_mode_t mode);
int gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, 
                     gsl_sf_result * result);
int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);

void gsl_error (const char * reason, const char * file, int line, int gsl_errno);
void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);
int gsl_sf_cos_e(double x, gsl_sf_result * result);
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);
int gsl_sf_sin_e(double x, gsl_sf_result * result);
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);

