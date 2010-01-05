/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#if defined(_WIN32)
#include <float.h>
#define isnan(x) _isnan(x)
#endif


int FindLineCircleIntersections2(double *x, double *y, 
                                double slope, double yi,
                                double xc, double yc, 
                                double r);

/* Find the intersections of a circle with a line, when the
 * line is specified by an angle and a point.
 * Returns the number of intersections found.
 * x[2] and y[2] have the intersection coordinates
 */

int FindLineCircleIntersections1(double *x, double *y,
                                 double x0, double y0, double theta, 
                                 double xc, double yc, double r)
{
  double slope=0.0, intercept, tmp;
  long rotate = 0, solutions;

  if (theta==PIo2 || theta==-PIo2)
    rotate = 1;
  else {
    slope = tan(theta);
    if (isnan(slope) || isinf(slope) || slope>1e10 || slope<-1e10)
      rotate = 1;
  }
  
  if (rotate) {
    theta += PIo2;
    slope = tan(theta);
    
    tmp = x0;
    x0 = -y0;
    y0 = tmp;
    tmp = xc;
    xc = -yc;
    yc = tmp;
    intercept = y0 - x0*slope;
    
    solutions = FindLineCircleIntersections2(y, x, slope, intercept,
                                             xc, yc, r);
    if (solutions>0)
      y[0] = -y[0];
    if (solutions>1)
      y[1] = -y[1];
  
    return solutions;
  } else {
    intercept = y0 - x0*slope;
    return FindLineCircleIntersections2(x, y, slope, intercept,
                                             xc, yc, r);
  }
}


/* Find the intersections of a circle with a line, when the
 * line is specified by slope and y intercept.
 * Returns the number of intersections found.
 * x[2] and y[2] have the intersection coordinates
 */
int FindLineCircleIntersections2(double *x, double *y, 
                                double slope, double yi,
                                double xc, double yc, 
                                double r)
{
  double a, b, c, d;

/*
  fprintf(stdout,"Finding intersection:\nLine: slope=%le, yi=%le\nCircle: c=%le, %le, r=%le\n",
          slope, yi,xc,yc,r);
  fflush(stdout);
*/
  
  a = 1+sqr(slope);
  b = 2*(slope*yi-xc-slope*yc);
  c = sqr(yc-yi) - sqr(r) + sqr(xc);

  d = sqr(b)-4*a*c;
  if (d<0)
    return 0;
  if (d==0) {
    x[0] = x[1] = -b/(2*a);
    y[0] = y[1] = x[0]*slope+yi;
    return 1;
  }
  else {
    d = sqrt(d);
    x[0] = (-b+d)/(2*a);
    x[1] = (-b-d)/(2*a);
    y[0] = x[0]*slope+yi;
    y[1] = x[1]*slope+yi;
    return 2;
  }
}





