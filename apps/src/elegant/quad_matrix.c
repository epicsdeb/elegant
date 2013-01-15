/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* contents: quadrupole_matrix(), qfringe_matrix(), quad_fringe(),
 *           qfringe_R_matrix(), qfringe_T_matrix()
 *
 * Michael Borland, 1989.
 */
#include "mdb.h"
#include "track.h"

static double swap_tmp;
#define swap_double(x, y) (swap_tmp=(x),(x)=(y),(y)=swap_tmp)

#define INSET_FRINGE 0
#define FIXED_STRENGTH_FRINGE 1
#define INTEGRALS_FRINGE 2
static char *fringeTypeOpt[3] = {"inset", "fixed-strength", "integrals"};

VMATRIX *quadrupole_matrix(double K1, double lHC, long maximum_order,
                           double tilt, double fse,
                           double xkick, double ykick,
                           double edge1_effects, double edge2_effects,
                           char *fringeType, double ffringe,
                           double *fringeIntM, double *fringeIntP,
			   long radial
                           )
{
    VMATRIX *M;
    VMATRIX *Mfringe, *Mtot, *Md, *tmp;
    double *C, **R, ***T, ****U;
    double kl, k, sin_kl, cos_kl, cosh_kl, sinh_kl;
    double lNominal, lEdge=0;
    long fringeCode = 0;
    
    K1 *= (1+fse);

    fringeCode = -1;
    if (fringeType && (fringeCode = match_string(fringeType, fringeTypeOpt, 3, 0))<0)
      bombElegant("Unrecognized fringe type for QUAD or KQUAD", NULL);
    
    if (K1==0 || lHC==0) {
      M = drift_matrix(lHC, maximum_order);
    }
    else {
      M = tmalloc(sizeof(*M));
      initialize_matrices(M, M->order = MIN(3,maximum_order));
      R = M->R;
      C = M->C;

      /* lHC is the "hard core" length (wherein K1 is constant)
       * lNominal is the effective length.
       * If fringe effects are off, these are the same.
       */
      lNominal = lHC;
      if (ffringe && fringeCode!=INTEGRALS_FRINGE) {
        /* This is the old default behavior, triggered just by having FFRINGE non-zero */
        if (edge1_effects==0 || edge2_effects==0)
          bombElegant("EDGE1_EFFECTS and EDGE2_EFFECTS must both be non-zero for FFRINGE-based quadrupole fringe effects", NULL);
        
        /* If mode is fixedStrength, then the sloped area is symmetric about the nominal
         * entrance and exit points.  This means that the integrated strength is not
         * changed.   FFRINGE=(2*fringeLengthOnOneSide)/HardEdgeLength
         * If the mode is "inset", then the sloped areas end at the nominal ends of
         * the quad.  The integrated strength changes as the fringe fraction changes.
         */
        /* length of each edge */
        lEdge = lNominal*ffringe/2;
        switch (fringeCode) {
        case 1:
          /* only half the total edge-field length is inside the nominal length */
          lHC = lNominal-lEdge;
          break;
        case 0:
          lHC = lNominal-2*lEdge;
          break;
        default:
          break;
        }
      }
      

      kl = (k=sqrt(fabs(K1)))*lHC;
      sin_kl  = sin(kl);
      cos_kl  = cos(kl);
      cosh_kl = cosh(kl);
      sinh_kl = sinh(kl);

      R[4][4] = R[5][5] = 1;
      C[4] = lHC;
      if (K1>0) {
        /* focussing in horizontal plane */
        R[0][0] = R[1][1] = cos_kl;
        R[0][1] = sin_kl/k;
        R[1][0] = -k*sin_kl;
        R[2][2] = R[3][3] = cosh_kl;
        R[2][3] = sinh_kl/k;
        R[3][2] = k*sinh_kl;
        
        if (M->order>=2) {
          double k2,k3,k4,l2,l3,l4,cos_kl,cos_2kl,cos_3kl,sin_kl,sin_2kl,sin_3kl,cosh_kl,cosh_2kl,cosh_3kl,sinh_kl,sinh_2kl,sinh_3kl ;
          k2 = pow(k,2) ;
          k3 = pow(k,3) ;
          k4 = pow(k,4) ;
          l2 = pow(lHC,2) ;
          l3 = pow(lHC,3) ;
          l4 = pow(lHC,4) ;
          cos_kl = cos(k*lHC) ;
          cos_2kl = cos(2*k*lHC) ;
          cos_3kl = cos(3*k*lHC) ;
          sin_kl = sin(k*lHC) ;
          sin_2kl = sin(2*k*lHC) ;
          sin_3kl = sin(3*k*lHC) ;
          cosh_kl = cosh(k*lHC) ;
          cosh_2kl = cosh(2*k*lHC) ;
          cosh_3kl = cosh(3*k*lHC) ;
          sinh_kl = sinh(k*lHC) ;
          sinh_2kl = sinh(2*k*lHC) ;
          sinh_3kl = sinh(3*k*lHC) ;

          T = M->T;
          T[0][5][0] = T[1][5][1] = kl*sin_kl/2;
          T[0][5][1] = sin_kl/(2*k) - lHC*cos_kl/2;
          T[1][5][0] = k/2*(kl*cos_kl + sin_kl);
          T[2][5][2] = -kl/2*sinh_kl;
          T[2][5][3] = (sinh_kl/k - lHC*cosh_kl)/2;
          T[3][5][2] = -k/2*(kl*cosh_kl + sinh_kl);
          T[3][5][3] = -kl/2*sinh_kl;
          T[4][0][0] = sqr(k)*(lHC - sin_kl/k*cos_kl)/4;
          T[4][1][0] = -sqr(sin_kl)/2;
          T[4][1][1] = (lHC + sin_kl/k*cos_kl)/4;
          T[4][2][2] = -sqr(k)*(lHC - sinh_kl/k*cosh_kl)/4;
          T[4][3][2] = sqr(sinh_kl)/2;
          T[4][3][3] = (lHC + sinh_kl/k*cosh_kl)/4;
          
          if (M->order>=3) {
            U = M->Q;
            U[0][0][0][0] = (-3*k3*lHC*sin_kl)/16. + (3*k2*sin_kl*sin_2kl)/32. ;
            U[0][1][0][0] = (-3*k2*lHC*cos_kl)/16. + (21*k*sin_kl)/64. - (3*k*sin_3kl)/64. ;
            U[0][1][1][0] = (-9*k*lHC*sin_kl)/16. - (3*sin_kl*sin_2kl)/32. ;
            U[0][1][1][1] = (3*lHC*cos_kl)/16. - (21*sin_kl)/(64.*k) + (3*sin_3kl)/(64.*k) ;
            U[0][2][2][0] = (k3*lHC*sin_kl)/8. + (k2*cos_kl*pow(sinh_kl,2))/16. - (3*k2*sin_kl*sinh_2kl)/32. ;
            U[0][2][2][1] = -(k2*lHC*cos_kl)/8. - (3*k*sin_kl)/32. + (k*cosh_2kl*sin_kl)/32. + (3*k*cos_kl*sinh_2kl)/32. ;
            U[0][3][2][0] = (k2*lHC*cos_kl)/4. - (7*k*sin_kl)/32. - (3*k*cosh_2kl*sin_kl)/32. + (k*cos_kl*sinh_2kl)/32. ;
            U[0][3][2][1] = (k*lHC*sin_kl)/4. + (3*cos_kl*pow(sinh_kl,2))/16. + (sin_kl*sinh_2kl)/32. ;
            U[0][3][3][0] = -(k*lHC*sin_kl)/8. + (cos_kl*pow(sinh_kl,2))/16. - (3*sin_kl*sinh_2kl)/32. ;
            U[0][3][3][1] = (lHC*cos_kl)/8. - (11*sin_kl)/(32.*k) + (cosh_2kl*sin_kl)/(32.*k) + (3*cos_kl*sinh_2kl)/(32.*k) ;
            U[0][5][5][0] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8. ;
            U[0][5][5][1] = (lHC*cos_kl)/8. - sin_kl/(8.*k) - (k*l2*sin_kl)/8. ;
            U[1][0][0][0] = (-3*k4*lHC*cos_kl)/16. - (3*k3*sin_kl)/16. + (3*k3*cos_2kl*sin_kl)/16. + (3*k3*cos_kl*sin_2kl)/32. ;
            U[1][1][0][0] = (9*k2*cos_kl)/64. - (9*k2*cos_3kl)/64. + (3*k3*lHC*sin_kl)/16. ;
            U[1][1][1][0] = (-9*k2*lHC*cos_kl)/16. - (9*k*sin_kl)/16. - (3*k*cos_2kl*sin_kl)/16. - (3*k*cos_kl*sin_2kl)/32. ;
            U[1][1][1][1] = (-9*cos_kl)/64. + (9*cos_3kl)/64. - (3*k*lHC*sin_kl)/16. ;
            U[1][2][2][0] = (k4*lHC*cos_kl)/8. + (k3*sin_kl)/8. - (3*k3*cosh_2kl*sin_kl)/16. + (k3*cos_kl*cosh_kl*sinh_kl)/8. - (k3*sin_kl*pow(sinh_kl,2))/16. - (3*k3*cos_kl*sinh_2kl)/32. ;
            U[1][2][2][1] = (-7*k2*cos_kl)/32. + (7*k2*cos_kl*cosh_2kl)/32. + (k3*lHC*sin_kl)/8. - (k2*sin_kl*sinh_2kl)/32. ;
            U[1][3][2][0] = (k2*cos_kl)/32. - (k2*cos_kl*cosh_2kl)/32. - (k3*lHC*sin_kl)/4. - (7*k2*sin_kl*sinh_2kl)/32. ;
            U[1][3][2][1] = (k2*lHC*cos_kl)/4. + (k*sin_kl)/4. + (k*cosh_2kl*sin_kl)/16. + (3*k*cos_kl*cosh_kl*sinh_kl)/8. - (3*k*sin_kl*pow(sinh_kl,2))/16. + (k*cos_kl*sinh_2kl)/32. ;
            U[1][3][3][0] = -(k2*lHC*cos_kl)/8. - (k*sin_kl)/8. - (3*k*cosh_2kl*sin_kl)/16. + (k*cos_kl*cosh_kl*sinh_kl)/8. - (k*sin_kl*pow(sinh_kl,2))/16. - (3*k*cos_kl*sinh_2kl)/32. ;
            U[1][3][3][1] = (-7*cos_kl)/32. + (7*cos_kl*cosh_2kl)/32. - (k*lHC*sin_kl)/8. - (sin_kl*sinh_2kl)/32. ;
            U[1][5][5][0] = (-5*k2*lHC*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
            U[1][5][5][1] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8. ;
            U[2][2][0][0] = (k2*cosh_kl*pow(sin_kl,2))/16. + (k3*lHC*sinh_kl)/8. - (3*k2*sin_2kl*sinh_kl)/32. ;
            U[2][2][1][0] = -(k2*lHC*cosh_kl)/4. - (k*cosh_kl*sin_2kl)/32. + (7*k*sinh_kl)/32. + (3*k*cos_2kl*sinh_kl)/32. ;
            U[2][2][1][1] = -(cosh_kl*pow(sin_kl,2))/16. + (k*lHC*sinh_kl)/8. + (3*sin_2kl*sinh_kl)/32. ;
            U[2][2][2][2] = (-3*k3*lHC*sinh_kl)/16. + (3*k2*sinh_kl*sinh_2kl)/32. ;
            U[2][3][0][0] = (k2*lHC*cosh_kl)/8. - (3*k*cosh_kl*sin_2kl)/32. + (3*k*sinh_kl)/32. - (k*cos_2kl*sinh_kl)/32. ;
            U[2][3][1][0] = (-3*cosh_kl)/32. + (3*cos_2kl*cosh_kl)/32. - (k*lHC*sinh_kl)/4. - (sin_2kl*sinh_kl)/32. ;
            U[2][3][1][1] = (lHC*cosh_kl)/8. + (3*cosh_kl*sin_2kl)/(32.*k) - (11*sinh_kl)/(32.*k) + (cos_2kl*sinh_kl)/(32.*k) ;
            U[2][3][2][2] = (3*k2*lHC*cosh_kl)/16. - (21*k*sinh_kl)/64. + (3*k*sinh_3kl)/64. ;
            U[2][3][3][2] = (9*k*lHC*sinh_kl)/16. + (3*sinh_kl*sinh_2kl)/32. ;
            U[2][3][3][3] = (3*lHC*cosh_kl)/16. - (21*sinh_kl)/(64.*k) + (3*sinh_3kl)/(64.*k) ;
            U[2][5][5][2] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8. ;
            U[2][5][5][3] = (lHC*cosh_kl)/8. - sinh_kl/(8.*k) + (k*l2*sinh_kl)/8. ;
            U[3][2][0][0] = (k4*lHC*cosh_kl)/8. + (k3*cos_kl*cosh_kl*sin_kl)/8. - (3*k3*cosh_kl*sin_2kl)/32. + (k3*sinh_kl)/8. - (3*k3*cos_2kl*sinh_kl)/16. + (k3*pow(sin_kl,2)*sinh_kl)/16. ;
            U[3][2][1][0] = -(k2*cosh_kl)/32. + (k2*cos_2kl*cosh_kl)/32. - (k3*lHC*sinh_kl)/4. - (7*k2*sin_2kl*sinh_kl)/32. ;
            U[3][2][1][1] = (k2*lHC*cosh_kl)/8. - (k*cos_kl*cosh_kl*sin_kl)/8. + (3*k*cosh_kl*sin_2kl)/32. + (k*sinh_kl)/8. + (3*k*cos_2kl*sinh_kl)/16. - (k*pow(sin_kl,2)*sinh_kl)/16. ;
            U[3][2][2][2] = (-3*k4*lHC*cosh_kl)/16. - (3*k3*sinh_kl)/16. + (3*k3*cosh_2kl*sinh_kl)/16. + (3*k3*cosh_kl*sinh_2kl)/32. ;
            U[3][3][0][0] = (7*k2*cosh_kl)/32. - (7*k2*cos_2kl*cosh_kl)/32. + (k3*lHC*sinh_kl)/8. - (k2*sin_2kl*sinh_kl)/32. ;
            U[3][3][1][0] = -(k2*lHC*cosh_kl)/4. - (7*k*cosh_kl*sin_2kl)/32. - (11*k*sinh_kl)/32. + (k*cos_2kl*sinh_kl)/32. ;
            U[3][3][1][1] = (-7*cosh_kl)/32. + (7*cos_2kl*cosh_kl)/32. + (k*lHC*sinh_kl)/8. + (sin_2kl*sinh_kl)/32. ;
            U[3][3][2][2] = (-9*k2*cosh_kl)/64. + (9*k2*cosh_3kl)/64. + (3*k3*lHC*sinh_kl)/16. ;
            U[3][3][3][2] = (9*k2*lHC*cosh_kl)/16. + (9*k*sinh_kl)/16. + (3*k*cosh_2kl*sinh_kl)/16. + (3*k*cosh_kl*sinh_2kl)/32. ;
            U[3][3][3][3] = (-9*cosh_kl)/64. + (9*cosh_3kl)/64. + (3*k*lHC*sinh_kl)/16. ;
            U[3][5][5][2] = (5*k2*lHC*cosh_kl)/8. + (3*k*sinh_kl)/8. + (k3*l2*sinh_kl)/8. ;
            U[3][5][5][3] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8.  ;
          }
        }
      } else {
        /* defocussing in horizontal plane */
        R[2][2] = R[3][3] = cos_kl;
        R[2][3] = sin_kl/k;
        R[3][2] = -k*sin_kl;
        R[0][0] = R[1][1] = cosh_kl;
        R[0][1] = sinh_kl/k;
        R[1][0] = k*sinh_kl;
        
        if (M->order>=2) {
          double k2,k3,k4,l2,l3,l4,cos_kl,cos_2kl,cos_3kl,sin_kl,sin_2kl,sin_3kl,cosh_kl,cosh_2kl,cosh_3kl,sinh_kl,sinh_2kl,sinh_3kl ;
          k2 = pow(k,2) ;
          k3 = pow(k,3) ;
          k4 = pow(k,4) ;
          l2 = pow(lHC,2) ;
          l3 = pow(lHC,3) ;
          l4 = pow(lHC,4) ;
          cos_kl = cos(k*lHC) ;
          cos_2kl = cos(2*k*lHC) ;
          cos_3kl = cos(3*k*lHC) ;
          sin_kl = sin(k*lHC) ;
          sin_2kl = sin(2*k*lHC) ;
          sin_3kl = sin(3*k*lHC) ;
          cosh_kl = cosh(k*lHC) ;
          cosh_2kl = cosh(2*k*lHC) ;
          cosh_3kl = cosh(3*k*lHC) ;
          sinh_kl = sinh(k*lHC) ;
          sinh_2kl = sinh(2*k*lHC) ;
          sinh_3kl = sinh(3*k*lHC) ;

          T = M->T;
          T[2][5][2] = T[3][5][3] = kl*sin_kl/2;
          T[2][5][3] = sin_kl/(2*k) - lHC*cos_kl/2;
          T[3][5][2] = k/2*(kl*cos_kl + sin_kl);
          T[0][5][0] = T[1][5][1] = -kl/2*sinh_kl;
          T[0][5][1] = (sinh_kl/k - lHC*cosh_kl)/2;
          T[1][5][0] = -k/2*(kl*cosh_kl + sinh_kl);
          T[4][0][0] = -sqr(k)*(lHC - sinh_kl/k*cosh_kl)/4;
          T[4][1][0] = sqr(sinh_kl)/2;
          T[4][1][1] = (lHC + sinh_kl/k*cosh_kl)/4;
          T[4][2][2] = sqr(k)*(lHC - sin_kl/k*cos_kl)/4;
          T[4][3][2] = -sqr(sin_kl)/2;
          T[4][3][3] = (lHC + sin_kl/k*cos_kl)/4;

          if (M->order>=3) {
            U = M->Q;
            U[0][0][0][0] = (-3*k3*lHC*sinh_kl)/16. + (3*k2*sinh_kl*sinh_2kl)/32. ;
            U[0][1][0][0] = (3*k2*lHC*cosh_kl)/16. - (21*k*sinh_kl)/64. + (3*k*sinh_3kl)/64. ;
            U[0][1][1][0] = (9*k*lHC*sinh_kl)/16. + (3*sinh_kl*sinh_2kl)/32. ;
            U[0][1][1][1] = (3*lHC*cosh_kl)/16. - (21*sinh_kl)/(64.*k) + (3*sinh_3kl)/(64.*k) ;
            U[0][2][2][0] = (k2*cosh_kl*pow(sin_kl,2))/16. + (k3*lHC*sinh_kl)/8. - (3*k2*sin_2kl*sinh_kl)/32. ;
            U[0][2][2][1] = (k2*lHC*cosh_kl)/8. - (3*k*cosh_kl*sin_2kl)/32. + (3*k*sinh_kl)/32. - (k*cos_2kl*sinh_kl)/32. ;
            U[0][3][2][0] = -(k2*lHC*cosh_kl)/4. - (k*cosh_kl*sin_2kl)/32. + (7*k*sinh_kl)/32. + (3*k*cos_2kl*sinh_kl)/32. ;
            U[0][3][2][1] = (-3*cosh_kl)/32. + (3*cos_2kl*cosh_kl)/32. - (k*lHC*sinh_kl)/4. - (sin_2kl*sinh_kl)/32. ;
            U[0][3][3][0] = -(cosh_kl*pow(sin_kl,2))/16. + (k*lHC*sinh_kl)/8. + (3*sin_2kl*sinh_kl)/32. ;
            U[0][3][3][1] = (lHC*cosh_kl)/8. + (3*cosh_kl*sin_2kl)/(32.*k) - (11*sinh_kl)/(32.*k) + (cos_2kl*sinh_kl)/(32.*k) ;
            U[0][5][5][0] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8. ;
            U[0][5][5][1] = (lHC*cosh_kl)/8. - sinh_kl/(8.*k) + (k*l2*sinh_kl)/8. ;
            U[1][0][0][0] = (-3*k4*lHC*cosh_kl)/16. - (3*k3*sinh_kl)/16. + (3*k3*cosh_2kl*sinh_kl)/16. + (3*k3*cosh_kl*sinh_2kl)/32. ;
            U[1][1][0][0] = (-9*k2*cosh_kl)/64. + (9*k2*cosh_3kl)/64. + (3*k3*lHC*sinh_kl)/16. ;
            U[1][1][1][0] = (9*k2*lHC*cosh_kl)/16. + (9*k*sinh_kl)/16. + (3*k*cosh_2kl*sinh_kl)/16. + (3*k*cosh_kl*sinh_2kl)/32. ;
            U[1][1][1][1] = (-9*cosh_kl)/64. + (9*cosh_3kl)/64. + (3*k*lHC*sinh_kl)/16. ;
            U[1][2][2][0] = (k4*lHC*cosh_kl)/8. + (k3*cos_kl*cosh_kl*sin_kl)/8. - (3*k3*cosh_kl*sin_2kl)/32. + (k3*sinh_kl)/8. - (3*k3*cos_2kl*sinh_kl)/16. + (k3*pow(sin_kl,2)*sinh_kl)/16. ;
            U[1][2][2][1] = (7*k2*cosh_kl)/32. - (7*k2*cos_2kl*cosh_kl)/32. + (k3*lHC*sinh_kl)/8. - (k2*sin_2kl*sinh_kl)/32. ;
            U[1][3][2][0] = -(k2*cosh_kl)/32. + (k2*cos_2kl*cosh_kl)/32. - (k3*lHC*sinh_kl)/4. - (7*k2*sin_2kl*sinh_kl)/32. ;
            U[1][3][2][1] = -(k2*lHC*cosh_kl)/4. - (7*k*cosh_kl*sin_2kl)/32. - (11*k*sinh_kl)/32. + (k*cos_2kl*sinh_kl)/32. ;
            U[1][3][3][0] = (k2*lHC*cosh_kl)/8. - (k*cos_kl*cosh_kl*sin_kl)/8. + (3*k*cosh_kl*sin_2kl)/32. + (k*sinh_kl)/8. + (3*k*cos_2kl*sinh_kl)/16. - (k*pow(sin_kl,2)*sinh_kl)/16. ;
            U[1][3][3][1] = (-7*cosh_kl)/32. + (7*cos_2kl*cosh_kl)/32. + (k*lHC*sinh_kl)/8. + (sin_2kl*sinh_kl)/32. ;
            U[1][5][5][0] = (5*k2*lHC*cosh_kl)/8. + (3*k*sinh_kl)/8. + (k3*l2*sinh_kl)/8. ;
            U[1][5][5][1] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8. ;
            U[2][2][0][0] = (k3*lHC*sin_kl)/8. + (k2*cos_kl*pow(sinh_kl,2))/16. - (3*k2*sin_kl*sinh_2kl)/32. ;
            U[2][2][1][0] = (k2*lHC*cos_kl)/4. - (7*k*sin_kl)/32. - (3*k*cosh_2kl*sin_kl)/32. + (k*cos_kl*sinh_2kl)/32. ;
            U[2][2][1][1] = -(k*lHC*sin_kl)/8. + (cos_kl*pow(sinh_kl,2))/16. - (3*sin_kl*sinh_2kl)/32. ;
            U[2][2][2][2] = (-3*k3*lHC*sin_kl)/16. + (3*k2*sin_kl*sin_2kl)/32. ;
            U[2][3][0][0] = -(k2*lHC*cos_kl)/8. - (3*k*sin_kl)/32. + (k*cosh_2kl*sin_kl)/32. + (3*k*cos_kl*sinh_2kl)/32. ;
            U[2][3][1][0] = (k*lHC*sin_kl)/4. + (3*cos_kl*pow(sinh_kl,2))/16. + (sin_kl*sinh_2kl)/32. ;
            U[2][3][1][1] = (lHC*cos_kl)/8. - (11*sin_kl)/(32.*k) + (cosh_2kl*sin_kl)/(32.*k) + (3*cos_kl*sinh_2kl)/(32.*k) ;
            U[2][3][2][2] = (-3*k2*lHC*cos_kl)/16. + (21*k*sin_kl)/64. - (3*k*sin_3kl)/64. ;
            U[2][3][3][2] = (-9*k*lHC*sin_kl)/16. - (3*sin_kl*sin_2kl)/32. ;
            U[2][3][3][3] = (3*lHC*cos_kl)/16. - (21*sin_kl)/(64.*k) + (3*sin_3kl)/(64.*k) ;
            U[2][5][5][2] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8. ;
            U[2][5][5][3] = (lHC*cos_kl)/8. - sin_kl/(8.*k) - (k*l2*sin_kl)/8. ;
            U[3][2][0][0] = (k4*lHC*cos_kl)/8. + (k3*sin_kl)/8. - (3*k3*cosh_2kl*sin_kl)/16. + (k3*cos_kl*cosh_kl*sinh_kl)/8. - (k3*sin_kl*pow(sinh_kl,2))/16. - (3*k3*cos_kl*sinh_2kl)/32. ;
            U[3][2][1][0] = (k2*cos_kl)/32. - (k2*cos_kl*cosh_2kl)/32. - (k3*lHC*sin_kl)/4. - (7*k2*sin_kl*sinh_2kl)/32. ;
            U[3][2][1][1] = -(k2*lHC*cos_kl)/8. - (k*sin_kl)/8. - (3*k*cosh_2kl*sin_kl)/16. + (k*cos_kl*cosh_kl*sinh_kl)/8. - (k*sin_kl*pow(sinh_kl,2))/16. - (3*k*cos_kl*sinh_2kl)/32. ;
            U[3][2][2][2] = (-3*k4*lHC*cos_kl)/16. - (3*k3*sin_kl)/16. + (3*k3*cos_2kl*sin_kl)/16. + (3*k3*cos_kl*sin_2kl)/32. ;
            U[3][3][0][0] = (-7*k2*cos_kl)/32. + (7*k2*cos_kl*cosh_2kl)/32. + (k3*lHC*sin_kl)/8. - (k2*sin_kl*sinh_2kl)/32. ;
            U[3][3][1][0] = (k2*lHC*cos_kl)/4. + (k*sin_kl)/4. + (k*cosh_2kl*sin_kl)/16. + (3*k*cos_kl*cosh_kl*sinh_kl)/8. - (3*k*sin_kl*pow(sinh_kl,2))/16. + (k*cos_kl*sinh_2kl)/32. ;
            U[3][3][1][1] = (-7*cos_kl)/32. + (7*cos_kl*cosh_2kl)/32. - (k*lHC*sin_kl)/8. - (sin_kl*sinh_2kl)/32. ;
            U[3][3][2][2] = (9*k2*cos_kl)/64. - (9*k2*cos_3kl)/64. + (3*k3*lHC*sin_kl)/16. ;
            U[3][3][3][2] = (-9*k2*lHC*cos_kl)/16. - (9*k*sin_kl)/16. - (3*k*cos_2kl*sin_kl)/16. - (3*k*cos_kl*sin_2kl)/32. ;
            U[3][3][3][3] = (-9*cos_kl)/64. + (9*cos_3kl)/64. - (3*k*lHC*sin_kl)/16. ;
            U[3][5][5][2] = (-5*k2*lHC*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
            U[3][5][5][3] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8.  ;
          }
        }
      }
      

      if (fringeCode==INSET_FRINGE || fringeCode==FIXED_STRENGTH_FRINGE) {
        if (lEdge && K1) {
          Md = NULL;
          Mtot = tmalloc(sizeof(*Mtot));
          initialize_matrices(Mtot, M->order);

          /* entrance fringe fields */
          Mfringe = quad_fringe(lEdge, K1, M->order, 0, 0.0);
        
          if (fringeCode==FIXED_STRENGTH_FRINGE) {
            /* drift back to fringe entrance */
            Md = drift_matrix(-lEdge/2, M->order);
            concat_matrices(Mtot, Mfringe, Md, 0);
            tmp = Mfringe;
            Mfringe = Mtot;
            Mtot = tmp;
          }
          
          concat_matrices(Mtot, M, Mfringe, 0);
          tmp  = Mtot;
          Mtot = M;
          M    = tmp;
          free_matrices(Mfringe); tfree(Mfringe); Mfringe = NULL;
          
          /* exit fringe fields */
          Mfringe = quad_fringe(lEdge, K1, M->order, 1, 0.0);
          concat_matrices(Mtot, Mfringe, M, 0);
          tmp  = Mtot;
          Mtot = M;
          M    = tmp;
          
          if (fringeCode==FIXED_STRENGTH_FRINGE) {
            /* drift back to quad exit plane */
            concat_matrices(Mtot, Md, M, 0);
            tmp = M;
            M = Mtot;
            Mtot = tmp;
            free_matrices(Md); tfree(Md); Md = NULL;
          }
          free_matrices(Mfringe); tfree(Mfringe); Mfringe = NULL;
          free_matrices(Mtot); tfree(Mtot); Mtot = NULL;
        }
      } else if (fringeCode==INTEGRALS_FRINGE && (edge1_effects || edge2_effects)) {
        short hasFringeIntegrals = 0, i;
        for (i=0; i<5; i++) 
          if (fringeIntM[i] || fringeIntP[i]) {
            hasFringeIntegrals = 1;
            break;
          }
        if (hasFringeIntegrals) {
          Mtot = tmalloc(sizeof(*Mtot));
          initialize_matrices(Mtot, M->order);
          Mfringe = NULL;
          if (edge1_effects) {
            Mfringe = quadFringeMatrix(NULL, K1, -1, fringeIntM, fringeIntP);
            concat_matrices(Mtot, M, Mfringe, 0);
            tmp = M; M = Mtot; Mtot = tmp;
          }
          if (edge2_effects) {
            Mfringe = quadFringeMatrix(Mfringe, K1, 1, fringeIntM, fringeIntP);
            concat_matrices(Mtot, Mfringe, M, 0);
            tmp = M; M = Mtot; Mtot = tmp;
          }
          free_matrices(Mtot); free(Mtot); Mtot = NULL;
          free_matrices(Mfringe); free(Mfringe); Mfringe = NULL;
        }
      }
    }

    if (radial) {
      long i, j, k, l;
      /* make a radially-focusing lens by copying the horizontal matrix to the vertical */ 
      for (i=2; i<4; i++)
	for (j=2; j<4; j++)
	  M->R[i][j] = M->R[i-2][j-2];
      if (M->order>1) {
	for (i=0; i<2; i++) {
	  for (j=0; j<6; j++) {
	    if (j==2 || j==3) continue;
	    for (k=0; k<=j; k++) {
	      if (j<4 && k<4) { 
		M->T[i+2][j+2][k+2] = M->T[i][j][k];
		/* printf("T[%ld][%ld][%ld] = T[%ld][%ld][%ld] = %le\n", 
		   i+2, j+2, k+2, i, j, k, M->T[i][j][k]); */
	      }
	      else if (k<4 && (k+2)<=j) {
		M->T[i+2][j][k+2] = M->T[i][j][k];
		/* printf("T[%ld][%ld][%ld] = T[%ld][%ld][%ld] = %le\n", 
		   i+2, j, k+2, i, j, k, M->T[i][j][k]); */
	      }
	      else {
		M->T[i+2][j][k] = M->T[i][j][k];
		/* printf("T[%ld][%ld][%ld] = T[%ld][%ld][%ld] = %le\n", 
		   i+2, j, k, i, j, k, M->T[i][j][k]); */
	      }
	    }
	  }
	}
	for (i=4; i<6; i++) {
	  for (j=0; j<6; j++) {
	    if (j==2 || j==3) continue;
	    for (k=0; k<=j; k++) {
	      if (j<4 && k<4) {
		M->T[i][j+2][k+2] = M->T[i][j][k];
		/* printf("T[%ld][%ld][%ld] = T[%ld][%ld][%ld] = %le\n", 
		   i, j+2, k+2, i, j, k, M->T[i][j][k]); */
	      }
	      else if (k<4 && (k+2)<=j) {
		M->T[i][j][k+2] = M->T[i][j][k];
		/* printf("T[%ld][%ld][%ld] = T[%ld][%ld][%ld] = %le\n", 
		   i, j, k+2, i, j, k, M->T[i][j][k]); */
	      }
	    }
	  }
	}
	if (M->order>2) {
	  for (i=0; i<2; i++) {
	    for (j=0; j<6; j++) {
	      if (j==2 || j==3) continue;
	      for (k=0; k<=j; k++) {
		if (j<4 && k<4) { 
		  for (l=0; l<=k; l++) {
		    if (l<4 && (l+2)<=k) 
		      M->Q[i+2][j+2][k+2][l+2] = M->Q[i][j][k][l];
		    else
		      M->Q[i+2][j+2][k+2][l] = M->Q[i][j][k][l];
		  }
		}
		else if (k<4 && (k+2)<=j) {
		  for (l=0; l<=k; l++) {
		    if (l<4 && (l+2)<=k)
		      M->Q[i+2][j][k+2][l+2] = M->Q[i][j][k][l];
		    else
		      M->Q[i+2][j][k+2][l] = M->Q[i][j][k][l];
		  }
		}
		else {
		  for (l=0; l<=k; l++) {
		    if (l<4 && (l+2)<=k)
		      M->Q[i+2][j][k][l+2] = M->Q[i][j][k][l];
		    else
		      M->Q[i+2][j][k][l] = M->Q[i][j][k][l];
		  }
		}
	      }
	    }
	  }
	  for (i=4; i<6; i++) {
	    for (j=0; j<6; j++) {
	      if (j==2 || j==3) continue;
	      for (k=0; k<=j; k++) {
		if (j<4 && k<4) {
		  for (l=0; l<=k; l++) {
		    if (l<4 && (l+2)<=k)
		      M->Q[i][j+2][k+2][l+4] = M->Q[i][j][k][l];
		    else
		      M->Q[i][j+2][k+2][l+4] = M->Q[i][j][k][l];
		  }
		}
		else if (k<4 && (k+2)<=j) {
		  for (l=0; l<=k; l++) {
		    if (l<4 && (l+2)<=k)
		      M->Q[i][j][k+2][l+2] = M->Q[i][j][k][l];
		    else
		      M->Q[i][j][k+2][l] = M->Q[i][j][k][l];
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    if (xkick || ykick) {
      /* put identical kicks at the entrance and exit */
      Mtot = tmalloc(sizeof(*Mtot));
      initialize_matrices(Mtot, M->order);
      Mfringe = hvcorrector_matrix(0, xkick/2, ykick/2, 0.0, 0.0, 1.0, 1.0, 0, M->order);
      concat_matrices(Mtot, Mfringe, M, 0);
      concat_matrices(M, Mtot, Mfringe, 0);
      free_matrices(Mfringe); tfree(Mfringe); Mfringe = NULL;
      free_matrices(Mtot); tfree(Mtot); Mtot = NULL;
    }
    
    tilt_matrices(M, tilt);
    return(M);
    }


VMATRIX *quad_fringe(double l, double ko, long order, long reverse, double fse)
{
    VMATRIX *M;

    log_entry("quad_fringe");

    ko *= (1+fse);
    M  = tmalloc(sizeof(*M));
    initialize_matrices(M, order);

    M->C[4] = l;
    M->R[5][5] = M->R[4][4] = 1;

    qfringe_R_matrix(
        &M->R[0][0], &M->R[0][1], &M->R[1][0], &M->R[1][1],  ko/l, l);
    qfringe_R_matrix(
        &M->R[2][2], &M->R[2][3], &M->R[3][2], &M->R[3][3], -ko/l, l);
    if (reverse) {
        swap_double(M->R[0][0], M->R[1][1]);
        swap_double(M->R[2][2], M->R[3][3]);
        }

    if (order>=2) {
        qfringe_T_matrix(
            &M->T[0][5][0], &M->T[0][5][1], &M->T[1][5][0], &M->T[1][5][1],
            &M->T[4][0][0], &M->T[4][1][0], &M->T[4][1][1], 
            ko/l, l, reverse);
        qfringe_T_matrix(
            &M->T[2][5][2], &M->T[2][5][3], &M->T[3][5][2], &M->T[3][5][3], 
            &M->T[4][2][2], &M->T[4][3][2], &M->T[4][3][3],
            -ko/l, l, reverse);
        }

    log_exit("quad_fringe");
    return(M);
    }

void qfringe_R_matrix(
    double *R11, double *R12,
    double *R21, double *R22,
    double dk_dz, double l
    )
{
    double term, l3;
    long n;

    log_entry("qfringe_R_matrix");

    if (!l) {
        *R11 = *R22 = 1;
        *R21 = *R12 = 0;
        log_exit("qfringe_R_matrix");
        return;
        }
    if (!dk_dz) {
        *R11 = *R22 = 1;
        *R21 = 0;
        *R12 = l;
        log_exit("qfringe_R_matrix");
        return;
        }

    l3 = pow3(l);

    /* compute R11, R21 */
    *R11 = *R21 = 0;
    term = 1;
    n = 0;
    do {
        *R11 += term;
        *R21 += n*term;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*R11))>1e-16);
    *R21 /= l;

    /* compute R12, R22 */
    *R12 = 0;
    *R22 = 0;
    term = l;
    n = 1;
    do {
        *R12 += term;
        *R22 += n*term;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*R12))>1e-16);
    *R22 /= l;
    log_exit("qfringe_R_matrix");
    }

void qfringe_T_matrix(
    double *T116, double *T126, double *T216, double *T226, 
    double *T511, double *T512, double *T522,
    double dk_dz, double l, long reverse
    )
{
    double term, l3, ko, mult;
    long n;

    log_entry("qfringe_T_matrix");

    if (!l || !dk_dz) {
        *T116 = *T226 = *T216 = *T126 = *T511 = *T512 = *T522 = 0;
        log_exit("qfringe_T_matrix");
        return;
        }

    l3 = pow3(l);

    /* compute T116, T216 */
    *T116 = *T216 = 0;
    term = 1;
    n = 0;
    do {
        *T116 += -n*term/3;
        *T216 += -sqr(n)*term/3;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*T116?*T116:1))>1e-16);
    *T216 /= l;

    /* compute T126, T226 */
    *T126 = 0;
    *T226 = 0;
    term = l;
    n = 1;
    do {
        *T126 += -(n-1)*term/3;
        *T226 += -n*(n-1)*term/3;
        term *= -l3*dk_dz/((n+3)*(n+2));
        n += 3;
        } while (FABS(term/(*T126?*T126:1))>1e-16);
    *T226 /= l;

    if (reverse) 
        swap_double(*T116, *T226);

    /* compute path-length terms using truncated series (truncated at z^12) */
    if (!reverse) {
        /* entrance fringe field */
        ko = dk_dz*l;
        /* T511 terms */
        term = sqr(ko)*pow3(l);
        if ((mult = ko*sqr(l))>1)
            fprintf(stdout, "warning: path-length terms for qfringe may be inaccurate: ko*lf^2>1\n");
            fflush(stdout);
        *T511  = 1./20.*term;         term *= mult;
        *T511 += -1./240.*term;       term *= mult; 
        *T511 +=  13./79200.*term;    term *= mult;
        *T511 += -19./4989600.*term;  term *= mult;
        *T511 += 19./325721088.*term;
        /* T522 terms */
        *T522  = term = l;            term *= mult;
        *T522 += -1./6.*term;         term *= mult;
        *T522 += 5./252.*term;        term *= mult;
        *T522 += -11./11340.*term;    term *= mult;
        *T522 += 187./7076160.*term;   term *= mult;
        *T522 += -391./849139200.*term; 
        *T522 /= 2;
        /* T512 terms */
        term = mult;
        *T512  = -1./6.*term;             term *= mult;
        *T512 += 1./30.*term;             term *= mult;
        *T512 += -1./480.*term;           term *= mult;
        *T512 += 1./14784.*term;          term *= mult;
        *T512 += -1./739200.*term; term *= mult;
        *T512 += 1./54454400.*term;
        }
    else {
        /* exit fringe field */
        ko = dk_dz*l;
        /* T511 terms */
        term = sqr(ko)*pow3(l);
        if ((mult = ko*sqr(l))>1)
            fprintf(stdout, "warning: path-length terms for qfringe may be inaccurate: ko*lf^2>1\n");
            fflush(stdout);
        *T511  = 2./15.*term;        term *= mult;
        *T511 += -1./80.*term;       term *= mult; 
        *T511 +=  1./1848.*term;     term *= mult;
        *T511 += -1./73920*term;     term *= mult;
        *T511 += 3./13613600.*term;
        /* T522 terms */
        *T522  = term = l;            term *= mult;
        *T522 += -1./6.*term;          term *= mult;
        *T522 += 1./70.*term;          term *= mult;
        *T522 += -1./1680.*term;       term *= mult;
        *T522 += 1./68640.*term;   term *= mult;
        *T522 += -3./12812800.*term; 
        *T522 /= 2;
        /* T512 terms */
        term = mult;
        *T512  = -1./3.*term;             term *= mult;
        *T512 += 1./20.*term;             term *= mult;
        *T512 += -1./336.*term;           term *= mult;
        *T512 += 1./10560.*term;          term *= mult;
        *T512 += -3/1601600.*term; term *= mult;
        *T512 += 1./39603200.*term;
        }
    log_exit("qfringe_T_matrix");
    }

/* This subroutine returns a stand-alone quadrpole fringe
 * matrix, in the event the user asks for a qfringe element
 * separate from a quadrupole
 */

VMATRIX *qfringe_matrix(
    double K1, double l, double tilt, long direction, long order, double fse
    )
{
    VMATRIX *M;

    log_entry("qfringe_matrix");

    M = quad_fringe(l, K1, order, (direction==-1?1:0), fse);

    M->C[4] = l;
    M->R[5][5] = M->R[4][4] = 1;

    tilt_matrices(M, tilt);

    log_exit("qfringe_matrix");
    return(M);
    }

VMATRIX *quse_matrix(double K1, double K2, double length, long maximum_order,
                           double tilt, double fse1, double fse2)
{
    VMATRIX *M;
    double *C, **R, ***T, ****U;
    double k, kl;

    K1 *= (1+fse1);
    
    if (K1==0 || length==0) {
      if (K2==0)
        return drift_matrix(length, maximum_order);
      return sextupole_matrix(K2, length, maximum_order, tilt, fse2, 0);
    }

    K2 *= (1+fse2);

    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order = MIN(3,maximum_order));
    C = M->C;
    R = M->R;

    if (K1>0) {
      double k2,k3,k4,k5,k6,k7,k8,l,l2;
      double cos_kl,cos_2kl,cos_3kl,sin_kl,sin_2kl,sin_3kl,cosh_kl,cosh_2kl,cosh_3kl,sinh_kl,sinh_2kl,sinh_3kl ;
      
      k = sqrt(K1);
      kl = k*length;
      l = length;
      
      k2 = ipow(k,2) ;
      k3 = ipow(k,3) ;
      k4 = ipow(k,4) ;
      k5 = ipow(k,5) ;
      k6 = ipow(k,6) ;
      k7 = ipow(k,7) ;
      k8 = ipow(k,8) ;
      l2 = ipow(l,2) ;
      cos_kl = cos(kl) ;
      cos_2kl = cos(2*kl) ;
      cos_3kl = cos(3*kl) ;
      sin_kl = sin(kl) ;
      sin_2kl = sin(2*kl) ;
      sin_3kl = sin(3*kl) ;
      cosh_kl = cosh(kl) ;
      cosh_2kl = cosh(2*kl) ;
      cosh_3kl = cosh(3*kl) ;
      sinh_kl = sinh(kl) ;
      sinh_2kl = sinh(2*kl) ;
      sinh_3kl = sinh(3*kl) ;

      C[4] = length;
      
      R[0][0] = cos_kl ;
      R[0][1] = sin_kl/k ;
      R[1][0] = -(k*sin_kl) ;
      R[1][1] = cos_kl ;
      R[2][2] = cosh_kl ;
      R[2][3] = sinh_kl/k ;
      R[3][2] = k*sinh_kl ;
      R[3][3] = cosh_kl ;
      R[4][4] = 1;
      R[5][5] = 1;
      if (M->order>1) {
        T = M->T;
        T[0][0][0] = (-2*K2*ipow(sin((kl)/2.),2))/(3.*k2) - (K2*cos_kl*ipow(sin((kl)/2.),2))/(3.*k2) ;
        T[0][1][0] = -(K2*sin_kl)/(6.*k3) + (K2*sin_2kl)/(12.*k3) ;
        T[0][1][1] = (-2*K2*ipow(sin((kl)/2.),4))/(3.*k4) ;
        T[0][2][2] = K2/(4.*k2) - (3*K2*cos_kl)/(10.*k2) + (K2*cosh_2kl)/(20.*k2) ;
        T[0][3][2] = -(K2*sin_kl)/(10.*k3) + (K2*sinh_2kl)/(20.*k3) ;
        T[0][3][3] = -K2/(4.*k4) + (K2*cos_kl)/(5.*k4) + (K2*cosh_2kl)/(20.*k4) ;
        T[0][5][0] = (kl*sin_kl)/2. ;
        T[0][5][1] = -(l*cos_kl)/2. + sin_kl/(2.*k) ;
        T[1][0][0] = (-2*K2*cos((kl)/2.)*sin((kl)/2.))/(3.*k) - (K2*cos((kl)/2.)*cos_kl*sin((kl)/2.))/(3.*k) + (K2*ipow(sin((kl)/2.),2)*sin_kl)/(3.*k) ;
        T[1][1][0] = -(K2*cos_kl)/(6.*k2) + (K2*cos_2kl)/(6.*k2) ;
        T[1][1][1] = (-4*K2*cos((kl)/2.)*ipow(sin((kl)/2.),3))/(3.*k3) ;
        T[1][2][2] = (3*K2*sin_kl)/(10.*k) + (K2*sinh_2kl)/(10.*k) ;
        T[1][3][2] = -(K2*cos_kl)/(10.*k2) + (K2*cosh_2kl)/(10.*k2) ;
        T[1][3][3] = -(K2*sin_kl)/(5.*k3) + (K2*sinh_2kl)/(10.*k3) ;
        T[1][5][0] = (k2*l*cos_kl)/2. + (k*sin_kl)/2. ;
        T[1][5][1] = (kl*sin_kl)/2. ;
        T[2][2][0] = (K2*cosh_kl)/(5.*k2) - (K2*cos_kl*cosh_kl)/(5.*k2) + (2*K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[2][2][1] = -(K2*cosh_kl*sin_kl)/(5.*k3) + (3*K2*sinh_kl)/(5.*k3) - (2*K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[2][3][0] = (2*K2*cosh_kl*sin_kl)/(5.*k3) - (K2*sinh_kl)/(5.*k3) - (K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[2][3][1] = (2*K2*cosh_kl)/(5.*k4) - (2*K2*cos_kl*cosh_kl)/(5.*k4) - (K2*sin_kl*sinh_kl)/(5.*k4) ;
        T[2][5][2] = -(kl*sinh_kl)/2. ;
        T[2][5][3] = -(l*cosh_kl)/2. + sinh_kl/(2.*k) ;
        T[3][2][0] = (3*K2*cosh_kl*sin_kl)/(5.*k) + (K2*sinh_kl)/(5.*k) + (K2*cos_kl*sinh_kl)/(5.*k) ;
        T[3][2][1] = (3*K2*cosh_kl)/(5.*k2) - (3*K2*cos_kl*cosh_kl)/(5.*k2) + (K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[3][3][0] = -(K2*cosh_kl)/(5.*k2) + (K2*cos_kl*cosh_kl)/(5.*k2) + (3*K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[3][3][1] = (K2*cosh_kl*sin_kl)/(5.*k3) + (2*K2*sinh_kl)/(5.*k3) - (3*K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[3][5][2] = -(k2*l*cosh_kl)/2. - (k*sinh_kl)/2. ;
        T[3][5][3] = -(kl*sinh_kl)/2. ;
        T[4][0][0] = (k2*l)/4. - (k*sin_2kl)/8. ;
        T[4][1][0] = -ipow(sin_kl,2)/2. ;
        T[4][1][1] = l/4. + sin_2kl/(8.*k) ;
        T[4][2][2] = -(k2*l)/4. + (k*sinh_2kl)/8. ;
        T[4][3][2] = ipow(sinh_kl,2)/2. ;
        T[4][3][3] = l/4. + sinh_2kl/(8.*k) ;
      }
      if (M->order>2) {
        U = M->Q;
        U[0][0][0][0] = -ipow(K2,2)/(12.*k4) + (3*k2*cos_kl)/64. + (29*ipow(K2,2)*cos_kl)/(576.*k4) + (ipow(K2,2)*cos_2kl)/(36.*k4) - (3*k2*cos_3kl)/64. + (ipow(K2,2)*cos_3kl)/(192.*k4) - (3*k3*l*sin_kl)/16. + (5*ipow(K2,2)*l*sin_kl)/(48.*k3) ;
        U[0][1][0][0] = (-3*k2*l*cos_kl)/16. - (ipow(K2,2)*l*cos_kl)/(16.*k4) + (21*k*sin_kl)/64. + (3*ipow(K2,2)*sin_kl)/(64.*k5) - (3*k*sin_3kl)/64. + (ipow(K2,2)*sin_3kl)/(192.*k5) ;
        U[0][1][1][0] = (-9*kl*sin_kl)/16. + (ipow(K2,2)*l*sin_kl)/(16.*k5) - (ipow(K2,2)*ipow(sin_kl,2))/(12.*k6) - (3*sin_kl*sin_2kl)/32. + (ipow(K2,2)*sin_kl*sin_2kl)/(96.*k6) ;
        U[0][1][1][1] = (3*l*cos_kl)/16. - (5*ipow(K2,2)*l*cos_kl)/(48.*k6) - (21*sin_kl)/(64.*k) + (5*ipow(K2,2)*sin_kl)/(576.*k7) + (ipow(K2,2)*sin_2kl)/(18.*k7) + (3*sin_3kl)/(64.*k) - (ipow(K2,2)*sin_3kl)/(192.*k7) ;
        U[0][2][2][0] = ipow(K2,2)/(8.*k4) - (k2*cos_kl)/32. - (141*ipow(K2,2)*cos_kl)/(1600.*k4) - (ipow(K2,2)*cos_2kl)/(40.*k4) + (ipow(K2,2)*cosh_2kl)/(100.*k4) + (k2*cos_kl*cosh_2kl)/32. - (7*ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k4) + (k3*l*sin_kl)/8. - (7*ipow(K2,2)*l*sin_kl)/(80.*k3) - (3*k2*sin_kl*sinh_2kl)/32. + (ipow(K2,2)*sin_kl*sinh_2kl)/(320.*k4) ;
        U[0][2][2][1] = -(k2*l*cos_kl)/8. + (7*ipow(K2,2)*l*cos_kl)/(80.*k4) - (3*k*sin_kl)/32. - (111*ipow(K2,2)*sin_kl)/(1600.*k5) + (k*cosh_2kl*sin_kl)/32. - (7*ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k5) - (ipow(K2,2)*sin_2kl)/(40.*k5) + (3*ipow(K2,2)*sinh_2kl)/(100.*k5) + (3*k*cos_kl*sinh_2kl)/32. - (ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k5) ;
        U[0][3][2][0] = (k2*l*cos_kl)/4. + (ipow(K2,2)*l*cos_kl)/(20.*k4) - (7*k*sin_kl)/32. - (61*ipow(K2,2)*sin_kl)/(4800.*k5) - (3*k*cosh_2kl*sin_kl)/32. + (ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k5) - (ipow(K2,2)*sin_2kl)/(120.*k5) + (ipow(K2,2)*sinh_2kl)/(100.*k5) + (k*cos_kl*sinh_2kl)/32. - (7*ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k5) ;
        U[0][3][2][1] = -ipow(K2,2)/(8.*k6) - (3*cos_kl)/32. + (431*ipow(K2,2)*cos_kl)/(4800.*k6) + (ipow(K2,2)*cos_2kl)/(120.*k6) + (3*ipow(K2,2)*cosh_2kl)/(100.*k6) + (3*cos_kl*cosh_2kl)/32. - (ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k6) + (kl*sin_kl)/4. + (ipow(K2,2)*l*sin_kl)/(20.*k5) + (sin_kl*sinh_2kl)/32. - (7*ipow(K2,2)*sin_kl*sinh_2kl)/(320.*k6) ;
        U[0][3][3][0] = -cos_kl/32. + (73*ipow(K2,2)*cos_kl)/(4800.*k6) + (ipow(K2,2)*cos_2kl)/(60.*k6) - (ipow(K2,2)*cosh_2kl)/(100.*k6) + (cos_kl*cosh_2kl)/32. - (7*ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k6) - (kl*sin_kl)/8. + (7*ipow(K2,2)*l*sin_kl)/(80.*k5) - (3*sin_kl*sinh_2kl)/32. + (ipow(K2,2)*sin_kl*sinh_2kl)/(320.*k6) ;
        U[0][3][3][1] = (l*cos_kl)/8. - (7*ipow(K2,2)*l*cos_kl)/(80.*k6) - (11*sin_kl)/(32.*k) + (203*ipow(K2,2)*sin_kl)/(4800.*k7) + (cosh_2kl*sin_kl)/(32.*k) - (7*ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k7) + (ipow(K2,2)*sin_2kl)/(60.*k7) + (ipow(K2,2)*sinh_2kl)/(50.*k7) + (3*cos_kl*sinh_2kl)/(32.*k) - (ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k7) ;
        U[0][5][0][0] = (K2*cos_kl)/(18.*k2) - (K2*cos_2kl)/(18.*k2) + (K2*l*sin_kl)/(12.*k) + (K2*l*sin_2kl)/(24.*k) ;
        U[0][5][1][0] = (K2*l)/(8.*k2) + (K2*l*cos_kl)/(12.*k2) - (K2*l*cos_2kl)/(24.*k2) - (5*K2*sin_kl)/(36.*k3) - (K2*sin_2kl)/(72.*k3) ;
        U[0][5][1][1] = -K2/(8.*k4) + (K2*cos_kl)/(9.*k4) + (K2*cos_2kl)/(72.*k4) + (K2*l*sin_kl)/(6.*k3) - (K2*l*sin_2kl)/(24.*k3) ;
        U[0][5][2][2] = (K2*cos_kl)/(50.*k2) - (K2*cosh_2kl)/(50.*k2) - (3*K2*l*sin_kl)/(20.*k) - (K2*l*sinh_2kl)/(40.*k) ;
        U[0][5][3][2] = -(K2*l)/(8.*k2) + (K2*l*cos_kl)/(20.*k2) - (K2*l*cosh_2kl)/(40.*k2) + (9*K2*sin_kl)/(100.*k3) + (K2*sinh_2kl)/(200.*k3) ;
        U[0][5][3][3] = -K2/(8.*k4) + (3*K2*cos_kl)/(25.*k4) + (K2*cosh_2kl)/(200.*k4) + (K2*l*sin_kl)/(10.*k3) - (K2*l*sinh_2kl)/(40.*k3) ;
        U[0][5][5][0] = -(k2*l2*cos_kl)/8. - (3*kl*sin_kl)/8. ;
        U[0][5][5][1] = (l*cos_kl)/8. - sin_kl/(8.*k) - (k*l2*sin_kl)/8. ;
        U[1][0][0][0] = (-3*k4*l*cos_kl)/16. + (5*ipow(K2,2)*l*cos_kl)/(48.*k2) - (15*k3*sin_kl)/64. + (31*ipow(K2,2)*sin_kl)/(576.*k3) - (ipow(K2,2)*sin_2kl)/(18.*k3) + (9*k3*sin_3kl)/64. - (ipow(K2,2)*sin_3kl)/(64.*k3) ;
        U[1][1][0][0] = (9*k2*cos_kl)/64. - (ipow(K2,2)*cos_kl)/(64.*k4) - (9*k2*cos_3kl)/64. + (ipow(K2,2)*cos_3kl)/(64.*k4) + (3*k3*l*sin_kl)/16. + (ipow(K2,2)*l*sin_kl)/(16.*k3) ;
        U[1][1][1][0] = (-9*k2*l*cos_kl)/16. + (ipow(K2,2)*l*cos_kl)/(16.*k4) - (9*k*sin_kl)/16. + (ipow(K2,2)*sin_kl)/(16.*k5) - (ipow(K2,2)*cos_kl*sin_kl)/(6.*k5) - (3*k*cos_2kl*sin_kl)/16. + (ipow(K2,2)*cos_2kl*sin_kl)/(48.*k5) - (3*k*cos_kl*sin_2kl)/32. + (ipow(K2,2)*cos_kl*sin_2kl)/(96.*k5) ;
        U[1][1][1][1] = (-9*cos_kl)/64. - (55*ipow(K2,2)*cos_kl)/(576.*k6) + (ipow(K2,2)*cos_2kl)/(9.*k6) + (9*cos_3kl)/64. - (ipow(K2,2)*cos_3kl)/(64.*k6) - (3*kl*sin_kl)/16. + (5*ipow(K2,2)*l*sin_kl)/(48.*k5) ;
        U[1][2][2][0] = (k4*l*cos_kl)/8. - (7*ipow(K2,2)*l*cos_kl)/(80.*k2) + (5*k3*sin_kl)/32. + (ipow(K2,2)*sin_kl)/(1600.*k3) - (7*k3*cosh_2kl*sin_kl)/32. + (9*ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k3) + (ipow(K2,2)*sin_2kl)/(20.*k3) + (ipow(K2,2)*sinh_2kl)/(50.*k3) - (k3*cos_kl*sinh_2kl)/32. - (13*ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k3) ;
        U[1][2][2][1] = (-7*k2*cos_kl)/32. + (29*ipow(K2,2)*cos_kl)/(1600.*k4) - (ipow(K2,2)*cos_2kl)/(20.*k4) + (3*ipow(K2,2)*cosh_2kl)/(50.*k4) + (7*k2*cos_kl*cosh_2kl)/32. - (9*ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k4) + (k3*l*sin_kl)/8. - (7*ipow(K2,2)*l*sin_kl)/(80.*k3) - (k2*sin_kl*sinh_2kl)/32. - (13*ipow(K2,2)*sin_kl*sinh_2kl)/(320.*k4) ;
        U[1][3][2][0] = (k2*cos_kl)/32. + (179*ipow(K2,2)*cos_kl)/(4800.*k4) - (ipow(K2,2)*cos_2kl)/(60.*k4) + (ipow(K2,2)*cosh_2kl)/(50.*k4) - (k2*cos_kl*cosh_2kl)/32. - (13*ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k4) - (k3*l*sin_kl)/4. - (ipow(K2,2)*l*sin_kl)/(20.*k3) - (7*k2*sin_kl*sinh_2kl)/32. + (9*ipow(K2,2)*sin_kl*sinh_2kl)/(320.*k4) ;
        U[1][3][2][1] = (k2*l*cos_kl)/4. + (ipow(K2,2)*l*cos_kl)/(20.*k4) + (11*k*sin_kl)/32. - (191*ipow(K2,2)*sin_kl)/(4800.*k5) - (k*cosh_2kl*sin_kl)/32. - (13*ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k5) - (ipow(K2,2)*sin_2kl)/(60.*k5) + (3*ipow(K2,2)*sinh_2kl)/(50.*k5) + (7*k*cos_kl*sinh_2kl)/32. - (9*ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k5) ;
        U[1][3][3][0] = -(k2*l*cos_kl)/8. + (7*ipow(K2,2)*l*cos_kl)/(80.*k4) - (3*k*sin_kl)/32. + (347*ipow(K2,2)*sin_kl)/(4800.*k5) - (7*k*cosh_2kl*sin_kl)/32. + (9*ipow(K2,2)*cosh_2kl*sin_kl)/(320.*k5) - (ipow(K2,2)*sin_2kl)/(30.*k5) - (ipow(K2,2)*sinh_2kl)/(50.*k5) - (k*cos_kl*sinh_2kl)/32. - (13*ipow(K2,2)*cos_kl*sinh_2kl)/(320.*k5) ;
        U[1][3][3][1] = (-7*cos_kl)/32. - (217*ipow(K2,2)*cos_kl)/(4800.*k6) + (ipow(K2,2)*cos_2kl)/(30.*k6) + (ipow(K2,2)*cosh_2kl)/(25.*k6) + (7*cos_kl*cosh_2kl)/32. - (9*ipow(K2,2)*cos_kl*cosh_2kl)/(320.*k6) - (kl*sin_kl)/8. + (7*ipow(K2,2)*l*sin_kl)/(80.*k5) - (sin_kl*sinh_2kl)/32. - (13*ipow(K2,2)*sin_kl*sinh_2kl)/ (320.*k6) ;
        U[1][5][0][0] = (K2*l*cos_kl)/12. + (K2*l*cos_2kl)/12. + (K2*sin_kl)/(36.*k) + (11*K2*sin_2kl)/(72.*k) ;
        U[1][5][1][0] = K2/(8.*k2) - (K2*cos_kl)/(18.*k2) - (5*K2*cos_2kl)/(72.*k2) - (K2*l*sin_kl)/(12.*k) + (K2*l*sin_2kl)/(12.*k) ;
        U[1][5][1][1] = (K2*l*cos_kl)/(6.*k2) - (K2*l*cos_2kl)/(12.*k2) + (K2*sin_kl)/(18.*k3) - (5*K2*sin_2kl)/(72.*k3) ;
        U[1][5][2][2] = (-3*K2*l*cos_kl)/20. - (K2*l*cosh_2kl)/20. - (17*K2*sin_kl)/(100.*k) - (13*K2*sinh_2kl)/(200.*k) ;
        U[1][5][3][2] = -K2/(8.*k2) + (7*K2*cos_kl)/(50.*k2) - (3*K2*cosh_2kl)/(200.*k2) - (K2*l*sin_kl)/(20.*k) - (K2*l*sinh_2kl)/(20.*k) ;
        U[1][5][3][3] = (K2*l*cos_kl)/(10.*k2) - (K2*l*cosh_2kl)/(20.*k2) - (K2*sin_kl)/(50.*k3) - (3*K2*sinh_2kl)/(200.*k3) ;
        U[1][5][5][0] = (-5*k2*l*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
        U[1][5][5][1] = -(k2*l2*cos_kl)/8. - (3*kl*sin_kl)/8. ;
        U[2][2][0][0] = (k2*cosh_kl*ipow(sin((kl)/2.),2))/8. + (143*ipow(K2,2)*cosh_kl*ipow(sin((kl)/2.),2))/(600.*k4) + (k2*cos_kl*cosh_kl*ipow(sin((kl)/2.),2))/8. + (11*ipow(K2,2)*cos_kl*cosh_kl*ipow(sin((kl)/2.),2))/ (120.*k4) + (k3*l*sinh_kl)/8. - (7*ipow(K2,2)*l*sinh_kl)/(40.*k3) + (11*ipow(K2,2)*sin_kl*sinh_kl)/(75.*k4) - (3*k2*sin_2kl*sinh_kl)/32. - (13*ipow(K2,2)*sin_2kl*sinh_kl)/(480.*k4) ;
        U[2][2][1][0] = -(k2*l*cosh_kl)/4. - (ipow(K2,2)*l*cosh_kl)/(10.*k4) + (41*ipow(K2,2)*cosh_kl*sin_kl)/(150.*k5) - (k*cosh_kl*sin_2kl)/32. - (11*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k5) + (7*k*sinh_kl)/32. - (81*ipow(K2,2)*sinh_kl)/(800.*k5) - (4*ipow(K2,2)*cos_kl*sinh_kl)/(75.*k5) + (3*k*cos_2kl*sinh_kl)/32. + (13*ipow(K2,2)*cos_2kl*sinh_kl)/(480.*k5) ;
        U[2][2][1][1] = -(cosh_kl*ipow(sin((kl)/2.),2))/8. + (313*ipow(K2,2)*cosh_kl*ipow(sin((kl)/2.),2))/(600.*k6) - (cos_kl*cosh_kl*ipow(sin((kl)/2.),2))/8. - (11*ipow(K2,2)*cos_kl*cosh_kl*ipow(sin((kl)/2.),2))/ (120.*k6) + (kl*sinh_kl)/8. - (7*ipow(K2,2)*l*sinh_kl)/(40.*k5) + (ipow(K2,2)*sin_kl*sinh_kl)/(75.*k6) + (3*sin_2kl*sinh_kl)/32. + (13*ipow(K2,2)*sin_2kl*sinh_kl)/(480.*k6) ;
        U[2][2][2][2] = (-3*k2*cosh_kl)/64. - (101*ipow(K2,2)*cosh_kl)/(1600.*k4) + (3*ipow(K2,2)*cos_kl*cosh_kl)/(50.*k4) + (3*k2*cosh_3kl)/64. + (ipow(K2,2)*cosh_3kl)/(320.*k4) - (3*k3*l*sinh_kl)/16. + (11*ipow(K2,2)*l*sinh_kl)/(80.*k3) - (3*ipow(K2,2)*sin_kl*sinh_kl)/(25.*k4) ;
        U[2][3][0][0] = (k2*l*cosh_kl)/8. - (7*ipow(K2,2)*l*cosh_kl)/(40.*k4) - (ipow(K2,2)*cosh_kl*sin_kl)/(75.*k5) - (3*k*cosh_kl*sin_2kl)/32. - (13*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k5) + (3*k*sinh_kl)/32. + (207*ipow(K2,2)*sinh_kl)/(800.*k5) + (ipow(K2,2)*cos_kl*sinh_kl)/(150.*k5) - (k*cos_2kl*sinh_kl)/32. - (11*ipow(K2,2)*cos_2kl*sinh_kl)/(480.*k5) ;
        U[2][3][1][0] = (ipow(K2,2)*cosh_kl)/(75.*k6) - (ipow(K2,2)*cos_kl*cosh_kl)/(75.*k6) - (3*cosh_kl*ipow(sin_kl,2))/16. - (13*ipow(K2,2)*cosh_kl*ipow(sin_kl,2))/(240.*k6) - (kl*sinh_kl)/4. - (ipow(K2,2)*l*sinh_kl)/(10.*k5) + (29*ipow(K2,2)*sin_kl*sinh_kl)/(150.*k6) - (sin_2kl*sinh_kl)/32. - (11*ipow(K2,2)*sin_2kl*sinh_kl)/ (480.*k6) ;
        U[2][3][1][1] = (l*cosh_kl)/8. - (7*ipow(K2,2)*l*cosh_kl)/(40.*k6) + (4*ipow(K2,2)*cosh_kl*sin_kl)/(75.*k7) + (3*cosh_kl*sin_2kl)/(32.*k) + (13*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k7) - (11*sinh_kl)/(32.*k) + (217*ipow(K2,2)*sinh_kl)/(800.*k7) - (17*ipow(K2,2)*cos_kl*sinh_kl)/(75.*k7) + (cos_2kl*sinh_kl)/(32.*k) + (11*ipow(K2,2)*cos_2kl*sinh_kl)/(480.*k7) ;
        U[2][3][2][2] = (3*k2*l*cosh_kl)/16. + (9*ipow(K2,2)*l*cosh_kl)/(80.*k4) - (3*ipow(K2,2)*cosh_kl*sin_kl)/(25.*k5) - (21*k*sinh_kl)/64. - (99*ipow(K2,2)*sinh_kl)/(1600.*k5) + (3*ipow(K2,2)*cos_kl*sinh_kl)/(50.*k5) + (3*k*sinh_3kl)/64. + (ipow(K2,2)*sinh_3kl)/(320.*k5) ;
        U[2][3][3][2] = (-2*ipow(K2,2)*cosh_kl*ipow(sin((kl)/2.),2))/(25.*k6) + (9*kl*sinh_kl)/16. - (ipow(K2,2)*l*sinh_kl)/(80.*k5) + (ipow(K2,2)*sin_kl*sinh_kl)/(50.*k6) + (3*sinh_kl*sinh_2kl)/32. + (ipow(K2,2)*sinh_kl*sinh_2kl)/(160.*k6) ;
        U[2][3][3][3] = (3*l*cosh_kl)/16. - (11*ipow(K2,2)*l*cosh_kl)/(80.*k6) + (2*ipow(K2,2)*cosh_kl*sin_kl)/(25.*k7) - (9*sinh_kl)/(32.*k) + (73*ipow(K2,2)*sinh_kl)/(800.*k7) - (ipow(K2,2)*cos_kl*sinh_kl)/(25.*k7) + (3*cosh_2kl*sinh_kl)/(32.*k) + (ipow(K2,2)*cosh_2kl*sinh_kl)/(160.*k7) ;
        U[2][5][2][0] = -(K2*cosh_kl)/(25.*k2) + (K2*cos_kl*cosh_kl)/(25.*k2) - (K2*l*cosh_kl*sin_kl)/(5.*k) - (K2*l*sinh_kl)/(10.*k) + (K2*l*cos_kl*sinh_kl)/(10.*k) - (7*K2*sin_kl*sinh_kl)/(25.*k2) ;
        U[2][5][2][1] = (-3*K2*l*cosh_kl)/(10.*k2) + (K2*l*cos_kl*cosh_kl)/(5.*k2) + (K2*cosh_kl*sin_kl)/(25.*k3) - (11*K2*sinh_kl)/(50.*k3) + (7*K2*cos_kl*sinh_kl)/(25.*k3) + (K2*l*sin_kl*sinh_kl)/(10.*k2) ;
        U[2][5][3][0] = (K2*l*cosh_kl)/(10.*k2) + (K2*l*cos_kl*cosh_kl)/(10.*k2) - (2*K2*cosh_kl*sin_kl)/(25.*k3) - (3*K2*sinh_kl)/(50.*k3) - (3*K2*cos_kl*sinh_kl)/(50.*k3) - (K2*l*sin_kl*sinh_kl)/(5.*k2) ;
        U[2][5][3][1] = (-2*K2*cosh_kl)/(25.*k4) + (2*K2*cos_kl*cosh_kl)/(25.*k4) + (K2*l*cosh_kl*sin_kl)/(10.*k3) - (K2*l*sinh_kl)/(5.*k3) + (K2*l*cos_kl*sinh_kl)/(5.*k3) - (3*K2*sin_kl*sinh_kl)/(50.*k4) ;
        U[2][5][5][2] = (k2*l2*cosh_kl)/8. + (3*kl*sinh_kl)/8. ;
        U[2][5][5][3] = (l*cosh_kl)/8. - sinh_kl/(8.*k) + (k*l2*sinh_kl)/8. ;
        U[3][2][0][0] = (k4*l*cosh_kl)/8. - (7*ipow(K2,2)*l*cosh_kl)/(40.*k2) + (k3*cos((kl)/2.)*cosh_kl*sin((kl)/2.))/8. + (143*ipow(K2,2)*cos((kl)/2.)*cosh_kl*sin((kl)/2.))/(600.*k3) + (k3*cos((kl)/2.)*cos_kl*cosh_kl*sin((kl)/2.))/8. + (11*ipow(K2,2)*cos((kl)/2.)*cos_kl*cosh_kl*sin((kl)/2.))/ (120.*k3) + (11*ipow(K2,2)*cosh_kl*sin_kl)/(75.*k3) - (k3*cosh_kl*ipow(sin((kl)/2.),2)*sin_kl)/8. - (11*ipow(K2,2)*cosh_kl*ipow(sin((kl)/2.),2)*sin_kl)/ (120.*k3) - (3*k3*cosh_kl*sin_2kl)/32. - (13*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k3) + (k3*sinh_kl)/8. - (7*ipow(K2,2)*sinh_kl)/(40.*k3) + (11*ipow(K2,2)*cos_kl*sinh_kl)/(75.*k3) - (3*k3*cos_2kl*sinh_kl)/16. - (13*ipow(K2,2)*cos_2kl*sinh_kl)/(240.*k3) + (k3*ipow(sin((kl)/2.),2)*sinh_kl)/8. + (143*ipow(K2,2)*ipow(sin((kl)/2.),2)*sinh_kl)/(600.*k3) + (k3*cos_kl*ipow(sin((kl)/2.),2)*sinh_kl)/8. + (11*ipow(K2,2)*cos_kl*ipow(sin((kl)/2.),2)*sinh_kl)/(120.*k3) ;
        U[3][2][1][0] = -(k2*cosh_kl)/32. - (161*ipow(K2,2)*cosh_kl)/(800.*k4) + (11*ipow(K2,2)*cos_kl*cosh_kl)/(50.*k4) + (k2*cos_2kl*cosh_kl)/32. - (3*ipow(K2,2)*cos_2kl*cosh_kl)/(160.*k4) - (k3*l*sinh_kl)/4. - (ipow(K2,2)*l*sinh_kl)/(10.*k3) + (49*ipow(K2,2)*sin_kl*sinh_kl)/(150.*k4) - (7*k2*sin_2kl*sinh_kl)/32. - (37*ipow(K2,2)*sin_2kl*sinh_kl)/(480.*k4) ;
        U[3][2][1][1] = (k2*l*cosh_kl)/8. - (7*ipow(K2,2)*l*cosh_kl)/(40.*k4) - (k*cos((kl)/2.)*cosh_kl*sin((kl)/2.))/8. + (313*ipow(K2,2)*cos((kl)/2.)*cosh_kl*sin((kl)/2.))/(600.*k5) - (k*cos((kl)/2.)*cos_kl*cosh_kl*sin((kl)/2.))/8. - (11*ipow(K2,2)*cos((kl)/2.)*cos_kl*cosh_kl*sin((kl)/2.))/ (120.*k5) + (ipow(K2,2)*cosh_kl*sin_kl)/(75.*k5) + (k*cosh_kl*ipow(sin((kl)/2.),2)*sin_kl)/8. + (11*ipow(K2,2)*cosh_kl*ipow(sin((kl)/2.),2)*sin_kl)/ (120.*k5) + (3*k*cosh_kl*sin_2kl)/32. + (13*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k5) + (k*sinh_kl)/8. - (7*ipow(K2,2)*sinh_kl)/(40.*k5) + (ipow(K2,2)*cos_kl*sinh_kl)/(75.*k5) + (3*k*cos_2kl*sinh_kl)/16. + (13*ipow(K2,2)*cos_2kl*sinh_kl)/(240.*k5) - (k*ipow(sin((kl)/2.),2)*sinh_kl)/8. + (313*ipow(K2,2)*ipow(sin((kl)/2.),2)*sinh_kl)/(600.*k5) - (k*cos_kl*ipow(sin((kl)/2.),2)*sinh_kl)/8. - (11*ipow(K2,2)*cos_kl*ipow(sin((kl)/2.),2)*sinh_kl)/(120.*k5) ;
        U[3][2][2][2] = (-3*k4*l*cosh_kl)/16. + (11*ipow(K2,2)*l*cosh_kl)/(80.*k2) - (9*ipow(K2,2)*cosh_kl*sin_kl)/(50.*k3) - (15*k3*sinh_kl)/64. + (119*ipow(K2,2)*sinh_kl)/(1600.*k3) - (3*ipow(K2,2)*cos_kl*sinh_kl)/(50.*k3) + (9*k3*sinh_3kl)/64. + (3*ipow(K2,2)*sinh_3kl)/(320.*k3) ;
        U[3][3][0][0] = (7*k2*cosh_kl)/32. + (67*ipow(K2,2)*cosh_kl)/(800.*k4) - (ipow(K2,2)*cos_kl*cosh_kl)/(150.*k4) - (7*k2*cos_2kl*cosh_kl)/32. - (37*ipow(K2,2)*cos_2kl*cosh_kl)/(480.*k4) + (k3*l*sinh_kl)/8. - (7*ipow(K2,2)*l*sinh_kl)/(40.*k3) - (ipow(K2,2)*sin_kl*sinh_kl)/(50.*k4) - (k2*sin_2kl*sinh_kl)/32. + (3*ipow(K2,2)*sin_2kl*sinh_kl)/(160.*k4) ;
        U[3][3][1][0] = -(k2*l*cosh_kl)/4. - (ipow(K2,2)*l*cosh_kl)/(10.*k4) + (31*ipow(K2,2)*cosh_kl*sin_kl)/(150.*k5) - (3*k*cos_kl*cosh_kl*sin_kl)/8. - (13*ipow(K2,2)*cos_kl*cosh_kl*sin_kl)/(120.*k5) - (k*cosh_kl*sin_2kl)/32. - (11*ipow(K2,2)*cosh_kl*sin_2kl)/(480.*k5) - (k*sinh_kl)/4. - (13*ipow(K2,2)*sinh_kl)/(150.*k5) + (9*ipow(K2,2)*cos_kl*sinh_kl)/(50.*k5) - (k*cos_2kl*sinh_kl)/16. - (11*ipow(K2,2)*cos_2kl*sinh_kl)/(240.*k5) - (3*k*ipow(sin_kl,2)*sinh_kl)/16. - (13*ipow(K2,2)*ipow(sin_kl,2)*sinh_kl)/(240.*k5) ;
        U[3][3][1][1] = (-7*cosh_kl)/32. + (77*ipow(K2,2)*cosh_kl)/(800.*k6) - (13*ipow(K2,2)*cos_kl*cosh_kl)/(75.*k6) + (7*cos_2kl*cosh_kl)/32. + (37*ipow(K2,2)*cos_2kl*cosh_kl)/(480.*k6) + (kl*sinh_kl)/8. - (7*ipow(K2,2)*l*sinh_kl)/(40.*k5) + (7*ipow(K2,2)*sin_kl*sinh_kl)/(25.*k6) + (sin_2kl*sinh_kl)/32. - (3*ipow(K2,2)*sin_2kl*sinh_kl)/ (160.*k6) ;
        U[3][3][2][2] = (-9*k2*cosh_kl)/64. + (81*ipow(K2,2)*cosh_kl)/(1600.*k4) - (3*ipow(K2,2)*cos_kl*cosh_kl)/(50.*k4) + (9*k2*cosh_3kl)/64. + (3*ipow(K2,2)*cosh_3kl)/(320.*k4) + (3*k3*l*sinh_kl)/16. + (9*ipow(K2,2)*l*sinh_kl)/(80.*k3) - (9*ipow(K2,2)*sin_kl*sinh_kl)/(50.*k4) ;
        U[3][3][3][2] = (9*k2*l*cosh_kl)/16. - (ipow(K2,2)*l*cosh_kl)/(80.*k4) - (2*ipow(K2,2)*cos((kl)/2.)*cosh_kl*sin((kl)/2.))/(25.*k5) + (ipow(K2,2)*cosh_kl*sin_kl)/(50.*k5) + (9*k*sinh_kl)/16. - (ipow(K2,2)*sinh_kl)/(80.*k5) + (ipow(K2,2)*cos_kl*sinh_kl)/(50.*k5) + (3*k*cosh_2kl*sinh_kl)/16. + (ipow(K2,2)*cosh_2kl*sinh_kl)/(80.*k5) - (2*ipow(K2,2)*ipow(sin((kl)/2.),2)*sinh_kl)/(25.*k5) + (3*k*cosh_kl*sinh_2kl)/32. + (ipow(K2,2)*cosh_kl*sinh_2kl)/(160.*k5) ;
        U[3][3][3][3] = (-3*cosh_kl)/32. - (37*ipow(K2,2)*cosh_kl)/(800.*k6) + (ipow(K2,2)*cos_kl*cosh_kl)/(25.*k6) + (3*cosh_kl*cosh_2kl)/32. + (ipow(K2,2)*cosh_kl*cosh_2kl)/(160.*k6) + (3*kl*sinh_kl)/16. - (11*ipow(K2,2)*l*sinh_kl)/(80.*k5) + (3*ipow(K2,2)*sin_kl*sinh_kl)/(25.*k6) + (3*sinh_kl*sinh_2kl)/16. + (ipow(K2,2)*sinh_kl*sinh_2kl)/(80.*k6) ;
        U[3][5][2][0] = -(K2*l*cosh_kl)/10. - (K2*l*cos_kl*cosh_kl)/10. - (13*K2*cosh_kl*sin_kl)/(25.*k) - (7*K2*sinh_kl)/(50.*k) - (7*K2*cos_kl*sinh_kl)/(50.*k) - (3*K2*l*sin_kl*sinh_kl)/10. ;
        U[3][5][2][1] = (-13*K2*cosh_kl)/(25.*k2) + (13*K2*cos_kl*cosh_kl)/(25.*k2) - (K2*l*cosh_kl*sin_kl)/(10.*k) - (3*K2*l*sinh_kl)/(10.*k) + (3*K2*l*cos_kl*sinh_kl)/(10.*k) - (7*K2*sin_kl*sinh_kl)/(50.*k2) ;
        U[3][5][3][0] = (K2*cosh_kl)/(25.*k2) - (K2*cos_kl*cosh_kl)/(25.*k2) - (3*K2*l*cosh_kl*sin_kl)/(10.*k) + (K2*l*sinh_kl)/(10.*k) - (K2*l*cos_kl*sinh_kl)/(10.*k) - (11*K2*sin_kl*sinh_kl)/(50.*k2) ;
        U[3][5][3][1] = -(K2*l*cosh_kl)/(5.*k2) + (3*K2*l*cos_kl*cosh_kl)/(10.*k2) - (K2*cosh_kl*sin_kl)/(25.*k3) - (7*K2*sinh_kl)/(25.*k3) + (11*K2*cos_kl*sinh_kl)/(50.*k3) - (K2*l*sin_kl*sinh_kl)/(10.*k2) ;
        U[3][5][5][2] = (5*k2*l*cosh_kl)/8. + (3*k*sinh_kl)/8. + (k3*l2*sinh_kl)/8. ;
        U[3][5][5][3] = (k2*l2*cosh_kl)/8. + (3*kl*sinh_kl)/8.  ;
      }
    } else {
      double k2,k3,k4,k5,k6,k7,k8,l,l2;
      double cos_kl,cos_2kl,cos_3kl,sin_kl,sin_2kl,sin_3kl,cosh_kl,cosh_2kl,cosh_3kl,sinh_kl,sinh_2kl,sinh_3kl ;

      k = sqrt(-K1);
      kl = k*length;
      l = length;
      
      C[4] = length;

      k2 = ipow(k,2) ;
      k3 = ipow(k,3) ;
      k4 = ipow(k,4) ;
      k5 = ipow(k,5) ;
      k6 = ipow(k,6) ;
      k7 = ipow(k,7) ;
      k8 = ipow(k,8) ;
      l2 = ipow(l,2) ;
      cos_kl = cos(kl) ;
      cos_2kl = cos(2*kl) ;
      cos_3kl = cos(3*kl) ;
      sin_kl = sin(kl) ;
      sin_2kl = sin(2*kl) ;
      sin_3kl = sin(3*kl) ;
      cosh_kl = cosh(kl) ;
      cosh_2kl = cosh(2*kl) ;
      cosh_3kl = cosh(3*kl) ;
      sinh_kl = sinh(kl) ;
      sinh_2kl = sinh(2*kl) ;
      sinh_3kl = sinh(3*kl) ;
      
      R[0][0] = cosh_kl ;
      R[0][1] = sinh_kl/k ;
      R[1][0] = k*sinh_kl ;
      R[1][1] = cosh_kl ;
      R[2][2] = cos_kl ;
      R[2][3] = sin_kl/k ;
      R[3][2] = -(k*sin_kl) ;
      R[3][3] = cos_kl ;
      R[4][4] = 1;
      R[5][5] = 1 ;
      
      if (M->order>1) {
        T = M->T;
        T[0][0][0] = (-2*K2*pow(sinh((kl)/2.),2))/(3.*k2) - (K2*cosh_kl*pow(sinh((kl)/2.),2))/(3.*k2) ;
        T[0][1][0] = (K2*sinh_kl)/(6.*k3) - (K2*sinh_2kl)/(12.*k3) ;
        T[0][1][1] = (-2*K2*pow(sinh((kl)/2.),4))/(3.*k4) ;
        T[0][2][2] = -K2/(4.*k2) - (K2*cos_2kl)/(20.*k2) + (3*K2*cosh_kl)/(10.*k2) ;
        T[0][3][2] = -(K2*sin_2kl)/(20.*k3) + (K2*sinh_kl)/(10.*k3) ;
        T[0][3][3] = -K2/(4.*k4) + (K2*cos_2kl)/(20.*k4) + (K2*cosh_kl)/(5.*k4) ;
        T[0][5][0] = -(kl*sinh_kl)/2. ;
        T[0][5][1] = -(l*cosh_kl)/2. + sinh_kl/(2.*k) ;
        T[1][0][0] = (-2*K2*cosh((kl)/2.)*sinh((kl)/2.))/(3.*k) - (K2*cosh((kl)/2.)*cosh_kl*sinh((kl)/2.))/(3.*k) - (K2*pow(sinh((kl)/2.),2)*sinh_kl)/(3.*k) ;
        T[1][1][0] = (K2*cosh_kl)/(6.*k2) - (K2*cosh_2kl)/(6.*k2) ;
        T[1][1][1] = (-4*K2*cosh((kl)/2.)*pow(sinh((kl)/2.),3))/(3.*k3) ;
        T[1][2][2] = (K2*sin_2kl)/(10.*k) + (3*K2*sinh_kl)/(10.*k) ;
        T[1][3][2] = -(K2*cos_2kl)/(10.*k2) + (K2*cosh_kl)/(10.*k2) ;
        T[1][3][3] = -(K2*sin_2kl)/(10.*k3) + (K2*sinh_kl)/(5.*k3) ;
        T[1][5][0] = -(k2*l*cosh_kl)/2. - (k*sinh_kl)/2. ;
        T[1][5][1] = -(kl*sinh_kl)/2. ;
        T[2][2][0] = -(K2*cos_kl)/(5.*k2) + (K2*cos_kl*cosh_kl)/(5.*k2) + (2*K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[2][2][1] = (-3*K2*sin_kl)/(5.*k3) + (2*K2*cosh_kl*sin_kl)/(5.*k3) + (K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[2][3][0] = (K2*sin_kl)/(5.*k3) + (K2*cosh_kl*sin_kl)/(5.*k3) - (2*K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[2][3][1] = (2*K2*cos_kl)/(5.*k4) - (2*K2*cos_kl*cosh_kl)/(5.*k4) + (K2*sin_kl*sinh_kl)/(5.*k4) ;
        T[2][5][2] = (kl*sin_kl)/2. ;
        T[2][5][3] = -(l*cos_kl)/2. + sin_kl/(2.*k) ;
        T[3][2][0] = (K2*sin_kl)/(5.*k) + (K2*cosh_kl*sin_kl)/(5.*k) + (3*K2*cos_kl*sinh_kl)/(5.*k) ;
        T[3][2][1] = (-3*K2*cos_kl)/(5.*k2) + (3*K2*cos_kl*cosh_kl)/(5.*k2) + (K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[3][3][0] = (K2*cos_kl)/(5.*k2) - (K2*cos_kl*cosh_kl)/(5.*k2) + (3*K2*sin_kl*sinh_kl)/(5.*k2) ;
        T[3][3][1] = (-2*K2*sin_kl)/(5.*k3) + (3*K2*cosh_kl*sin_kl)/(5.*k3) - (K2*cos_kl*sinh_kl)/(5.*k3) ;
        T[3][5][2] = (k2*l*cos_kl)/2. + (k*sin_kl)/2. ;
        T[3][5][3] = (kl*sin_kl)/2. ;
        T[4][0][0] = -(k2*l)/4. + (k*sinh_2kl)/8. ;
        T[4][1][0] = pow(sinh_kl,2)/2. ;
        T[4][1][1] = l/4. + sinh_2kl/(8.*k) ;
        T[4][2][2] = (k2*l)/4. - (k*sin_2kl)/8. ;
        T[4][3][2] = -pow(sin_kl,2)/2. ;
        T[4][3][3] = l/4. + sin_2kl/(8.*k) ;
      }
      
      if (M->order>2) {
        U = M->Q;
        U[0][0][0][0] = -pow(K2,2)/(12.*k4) - (3*k2*cosh_kl)/64. + (29*pow(K2,2)*cosh_kl)/(576.*k4) + (pow(K2,2)*cosh_2kl)/(36.*k4) + (3*k2*cosh_3kl)/64. + (pow(K2,2)*cosh_3kl)/(192.*k4) - (3*k3*l*sinh_kl)/16. - (5*pow(K2,2)*l*sinh_kl)/(48.*k3) ;
        U[0][1][0][0] = (3*k2*l*cosh_kl)/16. - (pow(K2,2)*l*cosh_kl)/(16.*k4) - (9*k*sinh_kl)/32. + (5*pow(K2,2)*sinh_kl)/(96.*k5) + (3*k*cosh_2kl*sinh_kl)/32. + (pow(K2,2)*cosh_2kl*sinh_kl)/(96.*k5) ;
        U[0][1][1][0] = (9*kl*sinh_kl)/16. + (pow(K2,2)*l*sinh_kl)/(16.*k5) - (pow(K2,2)*pow(sinh_kl,2))/(12.*k6) + (3*sinh_kl*sinh_2kl)/32. + (pow(K2,2)*sinh_kl*sinh_2kl)/(96.*k6) ;
        U[0][1][1][1] = (3*l*cosh_kl)/16. + (5*pow(K2,2)*l*cosh_kl)/(48.*k6) - (21*sinh_kl)/(64.*k) - (5*pow(K2,2)*sinh_kl)/(576.*k7) - (pow(K2,2)*sinh_2kl)/(18.*k7) + (3*sinh_3kl)/(64.*k) + (pow(K2,2)*sinh_3kl)/(192.*k7) ;
        U[0][2][2][0] = pow(K2,2)/(8.*k4) + (pow(K2,2)*cos_2kl)/(100.*k4) + (k2*cosh_kl)/32. - (141*pow(K2,2)*cosh_kl)/ (1600.*k4) - (k2*cos_2kl*cosh_kl)/32. - (7*pow(K2,2)*cos_2kl*cosh_kl)/(320.*k4) - (pow(K2,2)*cosh_2kl)/(40.*k4) + (k3*l*sinh_kl)/8. + (7*pow(K2,2)*l*sinh_kl)/(80.*k3) - (3*k2*sin_2kl*sinh_kl)/32. - (pow(K2,2)*sin_2kl*sinh_kl)/(320.*k4) ;
        U[0][2][2][1] = (k2*l*cosh_kl)/8. + (7*pow(K2,2)*l*cosh_kl)/(80.*k4) + (3*pow(K2,2)*sin_2kl)/(100.*k5) - (3*k*cosh_kl*sin_2kl)/32. - (pow(K2,2)*cosh_kl*sin_2kl)/(320.*k5) + (3*k*sinh_kl)/32. - (111*pow(K2,2)*sinh_kl)/(1600.*k5) - (k*cos_2kl*sinh_kl)/32. - (7*pow(K2,2)*cos_2kl*sinh_kl)/(320.*k5) - (pow(K2,2)*cosh_kl*sinh_kl)/(20.*k5) ;
        U[0][3][2][0] = -(k2*l*cosh_kl)/4. + (pow(K2,2)*l*cosh_kl)/(20.*k4) + (pow(K2,2)*sin_2kl)/(100.*k5) - (k*cosh_kl*sin_2kl)/32. - (7*pow(K2,2)*cosh_kl*sin_2kl)/(320.*k5) + (7*k*sinh_kl)/32. - (61*pow(K2,2)*sinh_kl)/(4800.*k5) + (3*k*cos_2kl*sinh_kl)/32. + (pow(K2,2)*cos_2kl*sinh_kl)/(320.*k5) - (pow(K2,2)*sinh_2kl)/(120.*k5) ;
        U[0][3][2][1] = pow(K2,2)/(8.*k6) - (3*pow(K2,2)*cos_2kl)/(100.*k6) - (3*cosh_kl)/32. - (431*pow(K2,2)*cosh_kl)/(4800.*k6) + (3*cos_2kl*cosh_kl)/32. + (pow(K2,2)*cos_2kl*cosh_kl)/(320.*k6) - (pow(K2,2)*cosh_2kl)/(120.*k6) - (kl*sinh_kl)/4. + (pow(K2,2)*l*sinh_kl)/(20.*k5) - (sin_2kl*sinh_kl)/32. - (7*pow(K2,2)*sin_2kl*sinh_kl)/(320.*k6) ;
        U[0][3][3][0] = (pow(K2,2)*cos_2kl)/(100.*k6) + (pow(K2,2)*cosh_kl)/(150.*k6) - (pow(K2,2)*cosh_2kl)/(60.*k6) - (cosh_kl*pow(sin_kl,2))/16. - (7*pow(K2,2)*cosh_kl*pow(sin_kl,2))/(160.*k6) + (kl*sinh_kl)/8. + (7*pow(K2,2)*l*sinh_kl)/(80.*k5) + (3*sin_2kl*sinh_kl)/32. + (pow(K2,2)*sin_2kl*sinh_kl)/(320.*k6) ;
        U[0][3][3][1] = (l*cosh_kl)/8. + (7*pow(K2,2)*l*cosh_kl)/(80.*k6) - (pow(K2,2)*sin_2kl)/(50.*k7) + (3*cosh_kl*sin_2kl)/(32.*k) + (pow(K2,2)*cosh_kl*sin_2kl)/(320.*k7) - (11*sinh_kl)/(32.*k) - (203*pow(K2,2)*sinh_kl)/(4800.*k7) + (cos_2kl*sinh_kl)/(32.*k) + (7*pow(K2,2)*cos_2kl*sinh_kl)/(320.*k7) - (pow(K2,2)*sinh_2kl)/(60.*k7) ;
        U[0][5][0][0] = -(K2*cosh_kl)/(18.*k2) + (K2*cosh_2kl)/(18.*k2) + (K2*l*sinh_kl)/(12.*k) + (K2*l*sinh_2kl)/(24.*k) ;
        U[0][5][1][0] = -(K2*l)/(8.*k2) - (K2*l*cosh_kl)/(12.*k2) + (K2*l*cosh_2kl)/(24.*k2) + (5*K2*sinh_kl)/(36.*k3) + (K2*sinh_2kl)/(72.*k3) ;
        U[0][5][1][1] = -K2/(8.*k4) + (K2*cosh_kl)/(9.*k4) + (K2*cosh_2kl)/(72.*k4) - (K2*l*sinh_kl)/(6.*k3) + (K2*l*sinh_2kl)/(24.*k3) ;
        U[0][5][2][2] = (K2*cos_2kl)/(50.*k2) - (K2*cosh_kl)/(50.*k2) - (K2*l*sin_2kl)/(40.*k) - (3*K2*l*sinh_kl)/(20.*k) ;
        U[0][5][3][2] = (K2*l)/(8.*k2) + (K2*l*cos_2kl)/(40.*k2) - (K2*l*cosh_kl)/(20.*k2) - (K2*sin_2kl)/(200.*k3) - (9*K2*sinh_kl)/(100.*k3) ;
        U[0][5][3][3] = -K2/(8.*k4) + (K2*cos_2kl)/(200.*k4) + (3*K2*cosh_kl)/(25.*k4) + (K2*l*sin_2kl)/(40.*k3) - (K2*l*sinh_kl)/(10.*k3) ;
        U[0][5][5][0] = (k2*l2*cosh_kl)/8. + (3*kl*sinh_kl)/8. ;
        U[0][5][5][1] = (l*cosh_kl)/8. - sinh_kl/(8.*k) + (k*l2*sinh_kl)/8. ;
        U[1][0][0][0] = (-3*k4*l*cosh_kl)/16. - (5*pow(K2,2)*l*cosh_kl)/(48.*k2) - (15*k3*sinh_kl)/64. - (31*pow(K2,2)*sinh_kl)/(576.*k3) + (pow(K2,2)*sinh_2kl)/(18.*k3) + (9*k3*sinh_3kl)/64. + (pow(K2,2)*sinh_3kl)/(64.*k3) ;
        U[1][1][0][0] = (-3*k2*cosh_kl)/32. - (pow(K2,2)*cosh_kl)/(96.*k4) + (3*k2*cosh_kl*cosh_2kl)/32. + (pow(K2,2)*cosh_kl*cosh_2kl)/(96.*k4) + (3*k3*l*sinh_kl)/16. - (pow(K2,2)*l*sinh_kl)/(16.*k3) + (3*k2*sinh_kl*sinh_2kl)/16. + (pow(K2,2)*sinh_kl*sinh_2kl)/(48.*k4) ;
        U[1][1][1][0] = (9*k2*l*cosh_kl)/16. + (pow(K2,2)*l*cosh_kl)/(16.*k4) + (9*k*sinh_kl)/16. + (pow(K2,2)*sinh_kl)/(16.*k5) - (pow(K2,2)*cosh_kl*sinh_kl)/(6.*k5) + (3*k*cosh_2kl*sinh_kl)/16. + (pow(K2,2)*cosh_2kl*sinh_kl)/(48.*k5) + (3*k*cosh_kl*sinh_2kl)/32. + (pow(K2,2)*cosh_kl*sinh_2kl)/(96.*k5) ;
        U[1][1][1][1] = (-9*cosh_kl)/64. + (55*pow(K2,2)*cosh_kl)/(576.*k6) - (pow(K2,2)*cosh_2kl)/(9.*k6) + (9*cosh_3kl)/64. + (pow(K2,2)*cosh_3kl)/(64.*k6) + (3*kl*sinh_kl)/16. + (5*pow(K2,2)*l*sinh_kl)/(48.*k5) ;
        U[1][2][2][0] = (k4*l*cosh_kl)/8. + (7*pow(K2,2)*l*cosh_kl)/(80.*k2) - (pow(K2,2)*sin_2kl)/(50.*k3) - (k3*cosh_kl*sin_2kl)/32. + (13*pow(K2,2)*cosh_kl*sin_2kl)/(320.*k3) + (5*k3*sinh_kl)/32. - (pow(K2,2)*sinh_kl)/(1600.*k3) - (7*k3*cos_2kl*sinh_kl)/32. - (9*pow(K2,2)*cos_2kl*sinh_kl)/(320.*k3) - (pow(K2,2)*sinh_2kl)/(20.*k3) ;
        U[1][2][2][1] = (3*pow(K2,2)*cos_2kl)/(50.*k4) + (7*k2*cosh_kl)/32. + (29*pow(K2,2)*cosh_kl)/(1600.*k4) - (7*k2*cos_2kl*cosh_kl)/32. - (9*pow(K2,2)*cos_2kl*cosh_kl)/(320.*k4) - (pow(K2,2)*pow(cosh_kl,2))/(20.*k4) + (k3*l*sinh_kl)/8. + (7*pow(K2,2)*l*sinh_kl)/(80.*k3) - (k2*sin_2kl*sinh_kl)/32. + (13*pow(K2,2)*sin_2kl*sinh_kl)/(320.*k4) - (pow(K2,2)*pow(sinh_kl,2))/(20.*k4) ;
        U[1][3][2][0] = (pow(K2,2)*cos_2kl)/(50.*k4) - (k2*cosh_kl)/32. + (179*pow(K2,2)*cosh_kl)/(4800.*k4) + (k2*cos_2kl*cosh_kl)/32. - (13*pow(K2,2)*cos_2kl*cosh_kl)/(320.*k4) - (pow(K2,2)*cosh_2kl)/(60.*k4) - (k3*l*sinh_kl)/4. + (pow(K2,2)*l*sinh_kl)/(20.*k3) - (7*k2*sin_2kl*sinh_kl)/32. - (9*pow(K2,2)*sin_2kl*sinh_kl)/(320.*k4) ;
        U[1][3][2][1] = -(k2*l*cosh_kl)/4. + (pow(K2,2)*l*cosh_kl)/(20.*k4) + (3*pow(K2,2)*sin_2kl)/(50.*k5) - (7*k*cosh_kl*sin_2kl)/32. - (9*pow(K2,2)*cosh_kl*sin_2kl)/(320.*k5) - (11*k*sinh_kl)/32. - (191*pow(K2,2)*sinh_kl)/(4800.*k5) + (k*cos_2kl*sinh_kl)/32. - (13*pow(K2,2)*cos_2kl*sinh_kl)/(320.*k5) - (pow(K2,2)*sinh_2kl)/(60.*k5) ;
        U[1][3][3][0] = (k2*l*cosh_kl)/8. + (7*pow(K2,2)*l*cosh_kl)/(80.*k4) - (k*cos_kl*cosh_kl*sin_kl)/8. - (7*pow(K2,2)*cos_kl*cosh_kl*sin_kl)/(80.*k5) - (pow(K2,2)*sin_2kl)/(50.*k5) + (3*k*cosh_kl*sin_2kl)/32. + (pow(K2,2)*cosh_kl*sin_2kl)/(320.*k5) + (k*sinh_kl)/8. + (113*pow(K2,2)*sinh_kl)/(1200.*k5) + (3*k*cos_2kl*sinh_kl)/16. + (pow(K2,2)*cos_2kl*sinh_kl)/(160.*k5) - (k*pow(sin_kl,2)*sinh_kl)/16. - (7*pow(K2,2)*pow(sin_kl,2)*sinh_kl)/(160.*k5) - (pow(K2,2)*sinh_2kl)/(30.*k5) ;
        U[1][3][3][1] = -(pow(K2,2)*cos_2kl)/(25.*k6) - (7*cosh_kl)/32. + (217*pow(K2,2)*cosh_kl)/(4800.*k6) + (7*cos_2kl*cosh_kl)/32. + (9*pow(K2,2)*cos_2kl*cosh_kl)/(320.*k6) - (pow(K2,2)*cosh_2kl)/(30.*k6) + (kl*sinh_kl)/8. + (7*pow(K2,2)*l*sinh_kl)/(80.*k5) + (sin_2kl*sinh_kl)/32. - (13*pow(K2,2)*sin_2kl*sinh_kl)/(320.*k6) ;
        U[1][5][0][0] = (K2*l*cosh_kl)/12. + (K2*l*cosh_2kl)/12. + (K2*sinh_kl)/(36.*k) + (11*K2*sinh_2kl)/(72.*k) ;
        U[1][5][1][0] = -K2/(8.*k2) + (K2*cosh_kl)/(18.*k2) + (5*K2*cosh_2kl)/(72.*k2) - (K2*l*sinh_kl)/(12.*k) + (K2*l*sinh_2kl)/(12.*k) ;
        U[1][5][1][1] = -(K2*l*cosh_kl)/(6.*k2) + (K2*l*cosh_2kl)/(12.*k2) - (K2*sinh_kl)/(18.*k3) + (5*K2*sinh_2kl)/(72.*k3) ;
        U[1][5][2][2] = -(K2*l*cos_2kl)/20. - (3*K2*l*cosh_kl)/20. - (13*K2*sin_2kl)/(200.*k) - (17*K2*sinh_kl)/(100.*k) ;
        U[1][5][3][2] = K2/(8.*k2) + (3*K2*cos_2kl)/(200.*k2) - (7*K2*cosh_kl)/(50.*k2) - (K2*l*sin_2kl)/(20.*k) - (K2*l*sinh_kl)/(20.*k) ;
        U[1][5][3][3] = (K2*l*cos_2kl)/(20.*k2) - (K2*l*cosh_kl)/(10.*k2) + (3*K2*sin_2kl)/(200.*k3) + (K2*sinh_kl)/(50.*k3) ;
        U[1][5][5][0] = (5*k2*l*cosh_kl)/8. + (3*k*sinh_kl)/8. + (k3*l2*sinh_kl)/8. ;
        U[1][5][5][1] = (k2*l2*cosh_kl)/8. + (3*kl*sinh_kl)/8. ;
        U[2][2][0][0] = (k3*l*sin_kl)/8. + (7*pow(K2,2)*l*sin_kl)/(40.*k3) + (k2*cos_kl*pow(sinh((kl)/2.),2))/8. - (143*pow(K2,2)*cos_kl*pow(sinh((kl)/2.),2))/(600.*k4) + (k2*cos_kl*cosh_kl*pow(sinh((kl)/2.),2))/8. - (11*pow(K2,2)*cos_kl*cosh_kl*pow(sinh((kl)/2.),2))/ (120.*k4) - (11*pow(K2,2)*sin_kl*sinh_kl)/(75.*k4) - (3*k2*sin_kl*sinh_2kl)/32. + (13*pow(K2,2)*sin_kl*sinh_2kl)/(480.*k4) ;
        U[2][2][1][0] = (k2*l*cos_kl)/4. - (pow(K2,2)*l*cos_kl)/(10.*k4) - (7*k*sin_kl)/32. - (81*pow(K2,2)*sin_kl)/(800.*k5) - (4*pow(K2,2)*cosh_kl*sin_kl)/(75.*k5) - (3*k*cosh_2kl*sin_kl)/32. + (13*pow(K2,2)*cosh_2kl*sin_kl)/(480.*k5) + (41*pow(K2,2)*cos_kl*sinh_kl)/(150.*k5) + (k*cos_kl*sinh_2kl)/32. - (11*pow(K2,2)*cos_kl*sinh_2kl)/(480.*k5) ;
        U[2][2][1][1] = -cos_kl/32. - (227*pow(K2,2)*cos_kl)/(800.*k6) + (23*pow(K2,2)*cos_kl*cosh_kl)/(75.*k6) + (cos_kl*cosh_2kl)/32. - (11*pow(K2,2)*cos_kl*cosh_2kl)/ (480.*k6) - (kl*sin_kl)/8. - (7*pow(K2,2)*l*sin_kl)/(40.*k5) + (pow(K2,2)*sin_kl*sinh_kl)/(75.*k6) - (3*sin_kl*sinh_2kl)/32. + (13*pow(K2,2)*sin_kl*sinh_2kl)/(480.*k6) ;
        U[2][2][2][2] = (3*k2*cos_kl)/64. - (101*pow(K2,2)*cos_kl)/(1600.*k4) - (3*k2*cos_3kl)/64. + (pow(K2,2)*cos_3kl)/(320.*k4) + (3*pow(K2,2)*cos_kl*cosh_kl)/(50.*k4) - (3*k3*l*sin_kl)/16. - (11*pow(K2,2)*l*sin_kl)/(80.*k3) + (3*pow(K2,2)*sin_kl*sinh_kl)/(25.*k4) ;
        U[2][3][0][0] = -(k2*l*cos_kl)/8. - (7*pow(K2,2)*l*cos_kl)/(40.*k4) - (3*k*sin_kl)/32. + (207*pow(K2,2)*sin_kl)/(800.*k5) + (pow(K2,2)*cosh_kl*sin_kl)/(150.*k5) + (k*cosh_2kl*sin_kl)/32. - (11*pow(K2,2)*cosh_2kl*sin_kl)/(480.*k5) - (pow(K2,2)*cos_kl*sinh_kl)/(75.*k5) + (3*k*cos_kl*sinh_2kl)/32. - (13*pow(K2,2)*cos_kl*sinh_2kl)/(480.*k5) ;
        U[2][3][1][0] = (kl*sin_kl)/4. - (pow(K2,2)*l*sin_kl)/(10.*k5) + (3*cos_kl*pow(sinh((kl)/2.),2))/8. - (49*pow(K2,2)*cos_kl*pow(sinh((kl)/2.),2))/(600.*k6) + (3*cos_kl*cosh_kl*pow(sinh((kl)/2.),2))/8. - (13*pow(K2,2)*cos_kl*cosh_kl*pow(sinh((kl)/2.),2))/ (120.*k6) + (29*pow(K2,2)*sin_kl*sinh_kl)/ (150.*k6) + (sin_kl*sinh_2kl)/32. - (11*pow(K2,2)*sin_kl*sinh_2kl)/(480.*k6) ;
        U[2][3][1][1] = (l*cos_kl)/8. + (7*pow(K2,2)*l*cos_kl)/(40.*k6) - (11*sin_kl)/(32.*k) - (217*pow(K2,2)*sin_kl)/(800.*k7) + (17*pow(K2,2)*cosh_kl*sin_kl)/(75.*k7) + (cosh_2kl*sin_kl)/(32.*k) - (11*pow(K2,2)*cosh_2kl*sin_kl)/(480.*k7) - (4*pow(K2,2)*cos_kl*sinh_kl)/(75.*k7) + (3*cos_kl*sinh_2kl)/(32.*k) - (13*pow(K2,2)*cos_kl*sinh_2kl)/(480.*k7) ;
        U[2][3][2][2] = (-3*k2*l*cos_kl)/16. + (9*pow(K2,2)*l*cos_kl)/(80.*k4) + (21*k*sin_kl)/64. - (99*pow(K2,2)*sin_kl)/(1600.*k5) + (3*pow(K2,2)*cosh_kl*sin_kl)/(50.*k5) - (3*k*sin_3kl)/64. + (pow(K2,2)*sin_3kl)/(320.*k5) - (3*pow(K2,2)*cos_kl*sinh_kl)/(25.*k5) ;
        U[2][3][3][2] = (pow(K2,2)*cos_kl)/(25.*k6) - (pow(K2,2)*cos_kl*cosh_kl)/(25.*k6) - (9*kl*sin_kl)/16. - (pow(K2,2)*l*sin_kl)/(80.*k5) - (3*sin_kl*sin_2kl)/32. + (pow(K2,2)*sin_kl*sin_2kl)/(160.*k6) + (pow(K2,2)*sin_kl*sinh_kl)/(50.*k6) ;
        U[2][3][3][3] = (3*l*cos_kl)/16. + (11*pow(K2,2)*l*cos_kl)/(80.*k6) - (21*sin_kl)/(64.*k) - (141*pow(K2,2)*sin_kl)/(1600.*k7) + (pow(K2,2)*cosh_kl*sin_kl)/(25.*k7) + (3*sin_3kl)/(64.*k) - (pow(K2,2)*sin_3kl)/(320.*k7) - (2*pow(K2,2)*cos_kl*sinh_kl)/(25.*k7) ;
        U[2][5][2][0] = (K2*cos_kl)/(25.*k2) - (K2*cos_kl*cosh_kl)/(25.*k2) - (K2*l*sin_kl)/(10.*k) + (K2*l*cosh_kl*sin_kl)/(10.*k) - (K2*l*cos_kl*sinh_kl)/(5.*k) - (7*K2*sin_kl*sinh_kl)/(25.*k2) ;
        U[2][5][2][1] = (3*K2*l*cos_kl)/(10.*k2) - (K2*l*cos_kl*cosh_kl)/(5.*k2) + (11*K2*sin_kl)/(50.*k3) - (7*K2*cosh_kl*sin_kl)/(25.*k3) - (K2*cos_kl*sinh_kl)/(25.*k3) + (K2*l*sin_kl*sinh_kl)/(10.*k2) ;
        U[2][5][3][0] = -(K2*l*cos_kl)/(10.*k2) - (K2*l*cos_kl*cosh_kl)/(10.*k2) + (3*K2*sin_kl)/(50.*k3) + (3*K2*cosh_kl*sin_kl)/(50.*k3) + (2*K2*cos_kl*sinh_kl)/(25.*k3) - (K2*l*sin_kl*sinh_kl)/(5.*k2) ;
        U[2][5][3][1] = (-2*K2*cos_kl)/(25.*k4) + (2*K2*cos_kl*cosh_kl)/(25.*k4) + (K2*l*sin_kl)/(5.*k3) - (K2*l*cosh_kl*sin_kl)/(5.*k3) - (K2*l*cos_kl*sinh_kl)/(10.*k3) + (3*K2*sin_kl*sinh_kl)/(50.*k4) ;
        U[2][5][5][2] = -(k2*l2*cos_kl)/8. - (3*kl*sin_kl)/8. ;
        U[2][5][5][3] = (l*cos_kl)/8. - sin_kl/(8.*k) - (k*l2*sin_kl)/8. ;
        U[3][2][0][0] = (k4*l*cos_kl)/8. + (7*pow(K2,2)*l*cos_kl)/(40.*k2) + (k3*sin_kl)/8. + (7*pow(K2,2)*sin_kl)/(40.*k3) - (11*pow(K2,2)*cosh_kl*sin_kl)/(75.*k3) - (3*k3*cosh_2kl*sin_kl)/16. + (13*pow(K2,2)*cosh_2kl*sin_kl)/(240.*k3) + (k3*cos_kl*cosh((kl)/2.)*sinh((kl)/2.))/8. - (143*pow(K2,2)*cos_kl*cosh((kl)/2.)*sinh((kl)/2.))/ (600.*k3) + (k3*cos_kl*cosh((kl)/2.)*cosh_kl* sinh((kl)/2.))/8. - (11*pow(K2,2)*cos_kl*cosh((kl)/2.)*cosh_kl* sinh((kl)/2.))/(120.*k3) - (k3*sin_kl*pow(sinh((kl)/2.),2))/8. + (143*pow(K2,2)*sin_kl*pow(sinh((kl)/2.),2))/(600.*k3) - (k3*cosh_kl*sin_kl*pow(sinh((kl)/2.),2))/8. + (11*pow(K2,2)*cosh_kl*sin_kl*pow(sinh((kl)/2.),2))/ (120.*k3) - (11*pow(K2,2)*cos_kl*sinh_kl)/(75.*k3) + (k3*cos_kl*pow(sinh((kl)/2.),2)*sinh_kl)/8. - (11*pow(K2,2)*cos_kl*pow(sinh((kl)/2.),2)*sinh_kl)/ (120.*k3) - (3*k3*cos_kl*sinh_2kl)/32. + (13*pow(K2,2)*cos_kl*sinh_2kl)/(480.*k3) ;
        U[3][2][1][0] = (k2*cos_kl)/32. - (161*pow(K2,2)*cos_kl)/(800.*k4) + (11*pow(K2,2)*cos_kl*cosh_kl)/(50.*k4) - (k2*cos_kl*cosh_2kl)/32. - (3*pow(K2,2)*cos_kl*cosh_2kl)/(160.*k4) - (k3*l*sin_kl)/4. + (pow(K2,2)*l*sin_kl)/(10.*k3) - (49*pow(K2,2)*sin_kl*sinh_kl)/(150.*k4) - (7*k2*sin_kl*sinh_2kl)/32. + (37*pow(K2,2)*sin_kl*sinh_2kl)/(480.*k4) ;
        U[3][2][1][1] = -(k2*l*cos_kl)/8. - (7*pow(K2,2)*l*cos_kl)/(40.*k4) - (3*k*sin_kl)/32. + (87*pow(K2,2)*sin_kl)/(800.*k5) - (22*pow(K2,2)*cosh_kl*sin_kl)/(75.*k5) - (7*k*cosh_2kl*sin_kl)/32. + (37*pow(K2,2)*cosh_2kl*sin_kl)/(480.*k5) + (8*pow(K2,2)*cos_kl*sinh_kl)/(25.*k5) - (k*cos_kl*sinh_2kl)/32. - (3*pow(K2,2)*cos_kl*sinh_2kl)/(160.*k5) ;
        U[3][2][2][2] = (-3*k4*l*cos_kl)/16. - (11*pow(K2,2)*l*cos_kl)/(80.*k2) - (15*k3*sin_kl)/64. - (119*pow(K2,2)*sin_kl)/(1600.*k3) + (3*pow(K2,2)*cosh_kl*sin_kl)/(50.*k3) + (9*k3*sin_3kl)/64. - (3*pow(K2,2)*sin_3kl)/(320.*k3) + (9*pow(K2,2)*cos_kl*sinh_kl)/(50.*k3) ;
        U[3][3][0][0] = (-7*k2*cos_kl)/32. + (67*pow(K2,2)*cos_kl)/(800.*k4) - (pow(K2,2)*cos_kl*cosh_kl)/(150.*k4) + (7*k2*cos_kl*cosh_2kl)/32. - (37*pow(K2,2)*cos_kl*cosh_2kl)/(480.*k4) + (k3*l*sin_kl)/8. + (7*pow(K2,2)*l*sin_kl)/(40.*k3) + (pow(K2,2)*sin_kl*sinh_kl)/(50.*k4) - (k2*sin_kl*sinh_2kl)/32. - (3*pow(K2,2)*sin_kl*sinh_2kl)/(160.*k4) ;
        U[3][3][1][0] = (k2*l*cos_kl)/4. - (pow(K2,2)*l*cos_kl)/(10.*k4) + (k*sin_kl)/4. - (pow(K2,2)*sin_kl)/(10.*k5) + (29*pow(K2,2)*cosh_kl*sin_kl)/(150.*k5) + (k*cosh_2kl*sin_kl)/16. - (11*pow(K2,2)*cosh_2kl*sin_kl)/(240.*k5) + (3*k*cos_kl*cosh((kl)/2.)*sinh((kl)/2.))/8. - (49*pow(K2,2)*cos_kl*cosh((kl)/2.)*sinh((kl)/2.))/(600.*k5) + (3*k*cos_kl*cosh((kl)/2.)*cosh_kl*sinh((kl)/2.))/8. - (13*pow(K2,2)*cos_kl*cosh((kl)/2.)*cosh_kl*sinh((kl)/2.))/ (120.*k5) - (3*k*sin_kl*pow(sinh((kl)/2.),2))/8. + (49*pow(K2,2)*sin_kl*pow(sinh((kl)/2.),2))/(600.*k5) - (3*k*cosh_kl*sin_kl*pow(sinh((kl)/2.),2))/8. + (13*pow(K2,2)*cosh_kl*sin_kl*pow(sinh((kl)/2.),2))/ (120.*k5) + (29*pow(K2,2)*cos_kl*sinh_kl)/ (150.*k5) + (3*k*cos_kl*pow(sinh((kl)/2.),2)*sinh_kl)/8. - (13*pow(K2,2)*cos_kl*pow(sinh((kl)/2.),2)*sinh_kl)/ (120.*k5) + (k*cos_kl*sinh_2kl)/32. - (11*pow(K2,2)*cos_kl*sinh_2kl)/(480.*k5) ;
        U[3][3][1][1] = (-7*cos_kl)/32. - (77*pow(K2,2)*cos_kl)/(800.*k6) + (13*pow(K2,2)*cos_kl*cosh_kl)/(75.*k6) + (7*cos_kl*cosh_2kl)/32. - (37*pow(K2,2)*cos_kl*cosh_2kl)/(480.*k6) - (kl*sin_kl)/8. - (7*pow(K2,2)*l*sin_kl)/(40.*k5) + (7*pow(K2,2)*sin_kl*sinh_kl)/(25.*k6) - (sin_kl*sinh_2kl)/32. - (3*pow(K2,2)*sin_kl*sinh_2kl)/ (160.*k6) ;
        U[3][3][2][2] = (9*k2*cos_kl)/64. + (81*pow(K2,2)*cos_kl)/(1600.*k4) - (9*k2*cos_3kl)/64. + (3*pow(K2,2)*cos_3kl)/(320.*k4) - (3*pow(K2,2)*cos_kl*cosh_kl)/(50.*k4) + (3*k3*l*sin_kl)/16. - (9*pow(K2,2)*l*sin_kl)/(80.*k3) + (9*pow(K2,2)*sin_kl*sinh_kl)/(50.*k4) ;
        U[3][3][3][2] = (-9*k2*l*cos_kl)/16. - (pow(K2,2)*l*cos_kl)/(80.*k4) - (9*k*sin_kl)/16. - (21*pow(K2,2)*sin_kl)/(400.*k5) - (3*k*cos_2kl*sin_kl)/16. + (pow(K2,2)*cos_2kl*sin_kl)/(80.*k5) + (3*pow(K2,2)*cosh_kl*sin_kl)/(50.*k5) - (3*k*cos_kl*sin_2kl)/32. + (pow(K2,2)*cos_kl*sin_2kl)/(160.*k5) - (pow(K2,2)*cos_kl*sinh_kl)/(50.*k5) ;
        U[3][3][3][3] = (-9*cos_kl)/64. + (79*pow(K2,2)*cos_kl)/(1600.*k6) + (9*cos_3kl)/64. - (3*pow(K2,2)*cos_3kl)/(320.*k6) - (pow(K2,2)*cos_kl*cosh_kl)/(25.*k6) - (3*kl*sin_kl)/16. - (11*pow(K2,2)*l*sin_kl)/(80.*k5) + (3*pow(K2,2)*sin_kl*sinh_kl)/(25.*k6) ;
        U[3][5][2][0] = -(K2*l*cos_kl)/10. - (K2*l*cos_kl*cosh_kl)/10. - (7*K2*sin_kl)/(50.*k) - (7*K2*cosh_kl*sin_kl)/(50.*k) - (13*K2*cos_kl*sinh_kl)/(25.*k) + (3*K2*l*sin_kl*sinh_kl)/10. ;
        U[3][5][2][1] = (13*K2*cos_kl)/(25.*k2) - (13*K2*cos_kl*cosh_kl)/(25.*k2) - (3*K2*l*sin_kl)/(10.*k) + (3*K2*l*cosh_kl*sin_kl)/(10.*k) - (K2*l*cos_kl*sinh_kl)/(10.*k) - (7*K2*sin_kl*sinh_kl)/(50.*k2) ;
        U[3][5][3][0] = -(K2*cos_kl)/(25.*k2) + (K2*cos_kl*cosh_kl)/(25.*k2) + (K2*l*sin_kl)/(10.*k) - (K2*l*cosh_kl*sin_kl)/(10.*k) - (3*K2*l*cos_kl*sinh_kl)/(10.*k) - (11*K2*sin_kl*sinh_kl)/(50.*k2) ;
        U[3][5][3][1] = (K2*l*cos_kl)/(5.*k2) - (3*K2*l*cos_kl*cosh_kl)/(10.*k2) + (7*K2*sin_kl)/(25.*k3) - (11*K2*cosh_kl*sin_kl)/(50.*k3) + (K2*cos_kl*sinh_kl)/(25.*k3) - (K2*l*sin_kl*sinh_kl)/(10.*k2) ;
        U[3][5][5][2] = (-5*k2*l*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
        U[3][5][5][3] = -(k2*l2*cos_kl)/8. - (3*kl*sin_kl)/8.  ;
      }
    }
    
    tilt_matrices(M, tilt);
    return(M);
  }
