/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsTouschekInteg.c
 * purpose: Evaluates the Touschek lifetime Integral C(zeta) as
 *  per Murphy's Synchrotron Light Source Data Book, page 16.
 *
 */
#include "mdb.h"
#if defined(_WIN32)
#include <float.h>
#define isnan(x) _isnan(x)
#else
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif
#endif
double eacc;
double f1;

long gaussianQuadrature(double(*fn)(), double a, double b, long n, double err, double *result);

int main()
{
  double a, b, err, factor, eaccFinal;
  long i, n, limit, evals, decades, decade;
  double fn(), result=0, lastResult;
  double fn1(double x);
  char filename[500];
  FILE *fp;
  
  b = PI/2;

  err = query_double("Error for adaptive integration: ", 1e-10);
  n = query_long("Number of initial evaluations: ", 100);
  eacc = query_double("Starting value of zeta: ", 1e-9);
  eaccFinal = query_double("Ending value of zeta: ", 500);
  limit = query_long("Number of values per decade: ", 10);
  queryn("Output file: ", filename, 500);
  chop_nl(filename);
  factor = exp(log(result/eacc)/limit);
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stdout, "Couldn't open: %s\n", filename);
    fflush(stdout);
    exit(1);
  }
  fprintf(fp, "SDDS1\n");
  fprintf(fp, "&parameter name=IntegError, type=double, fixed_value=%e &end\n",
          err);
  fprintf(fp, "&parameter name=NInitial, type=long, fixed_value=%ld &end\n",
          n);
  fprintf(fp, "&parameter name=ZetaStart, type=double, fixed_value=%e &end\n",
          eacc);
  fprintf(fp, "&parameter name=ZetaEnd, type=double, fixed_value=%e &end\n",
          eaccFinal);
  fputs("&column name=Evaluations type=long &end\n", fp);
  fputs("&column name=zeta type=double &end\n", fp);
  fputs("&column name=C type=double &end\n", fp);
  fputs("&column name=Evaluations1 type=long &end\n", fp);
  fputs("&column name=F type=double &end\n", fp);
  fputs("&data mode=ascii no_row_counts=1 &end\n", fp);
  fflush(fp);
  
  lastResult = 0;
  decades = log10(eaccFinal/eacc)+0.5;
  for (decade=0; decade<decades; decade++) {
    factor = pow(10.0, 1./limit);
    for (i=0; i<limit; i++) {
      f1 = 3*eacc-eacc*log(eacc)+2;
      a = atan(eacc);
      evals = gaussianQuadrature(fn, a, b, n, err, &result);
      if (result==0)
        break;
      fprintf(fp, "%ld %.16e %.16e ", evals, eacc, 0.5*result-1.5*exp(-eacc));
      evals = gaussianQuadrature(fn1, 0.0, PI/2, n, err, &result);
      fprintf(fp, "%ld %.16e \n", evals, result);
      fflush(fp);
      eacc *= factor;
      if (eacc>eaccFinal)
        break;
      if (lastResult) {
        /* err *= result/lastResult;
           fprintf(stdout, "err -> %le\n", err);
           fflush(stdout);
           */
      }
      lastResult = result;
    }
    if (eacc>eaccFinal)
      break;
  }
  fclose(fp);
  return(0);
}
  
double fn(double x)
{
  double sinx, cosx, tanx, prod, result;
  
  sinx = sin(x);
  cosx = cos(x);
  if ((prod=sinx*cosx)==0)
    return 0;
  tanx = tan(x);
  if (isnan(tanx))
    return 0;
  result = (eacc*log(tanx)+f1)*exp(-tanx)/prod;
  if (isnan(result))
    return 0;
  return result;
}

double fn1(double x)
{
  double tanx, cosx, tanxp1;
  double result;
  tanx = tan(x);
  cosx = cos(x);
  if (isnan(tanx) || isinf(tanx))
    return 0;
  if (cosx==0)
    return 0;
  tanxp1 = tanx+1;
  result = (tanx/tanxp1-0.5*log(tanxp1)/tanxp1)*exp(-eacc*tanxp1)/tanxp1;
  if (result==0)
    return 0;
  result = result/(cosx*cosx);
  if (isnan(result) || isinf(result))
    return 0;
  return result;
}

