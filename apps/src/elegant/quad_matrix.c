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

VMATRIX *quadrupole_matrix(double K1, double lHC, long maximum_order,
                           double tilt, double ffringe, double fse,
                           double xkick, double ykick,
                           char *fringeType)
{
    VMATRIX *M;
    VMATRIX *Mfringe, *Mtot, *Md, *tmp;
    double *C, **R, ***T, ****U;
    double kl, k, sin_kl, cos_kl, cosh_kl, sinh_kl;
    double lNominal, lEdge=0;
    static char *fringeTypeOpt[2] = {"inset", "fixed-strength"};
    long fixedStrengthFringe = 0;
    
    K1 *= (1+fse);
    
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
      if (ffringe) {
        /* If mode is fixedStrength, then the sloped area is symmetric about the nominal
         * entrance and exit points.  This means that the integrated strength is not
         * changed. 
         * If the mode is "inset", then the sloped areas end at the nominal ends of
         * the quad.  The integrated strength changes as the fringe fraction changes.
         */
        if (fringeType) {
          if ((fixedStrengthFringe = match_string(fringeType, fringeTypeOpt, 2, 0))<0)
            bombElegant("unrecognized fringe type for QUAD", NULL);
        }
        /* length of each edge */
        lEdge = lNominal*ffringe/2;
        if (fixedStrengthFringe)
          /* only half the total edge-field length is inside the nominal length */
          lHC = lNominal-lEdge;
        else 
          lHC = lNominal-2*lEdge;
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
            U[0][3][2][0] = -(k2*lHC*cos_kl)/4. + (9*k*sin_kl)/32. - (3*k*cosh_2kl*sin_kl)/32. + (k*cos_kl*sinh_2kl)/32. ;
            U[0][3][2][1] = -(k*lHC*sin_kl)/4. + (3*cos_kl*pow(sinh_kl,2))/16. + (sin_kl*sinh_2kl)/32. ;
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
            U[1][3][2][0] = (k2*cos_kl)/32. - (k2*cos_kl*cosh_2kl)/32. + (k3*lHC*sin_kl)/4. - (7*k2*sin_kl*sinh_2kl)/32. ;
            U[1][3][2][1] = -(k2*lHC*cos_kl)/4. - (k*sin_kl)/4. + (k*cosh_2kl*sin_kl)/16. + (3*k*cos_kl*cosh_kl*sinh_kl)/8. - (3*k*sin_kl*pow(sinh_kl,2))/16. + (k*cos_kl*sinh_2kl)/32. ;
            U[1][3][3][0] = -(k2*lHC*cos_kl)/8. - (k*sin_kl)/8. - (3*k*cosh_2kl*sin_kl)/16. + (k*cos_kl*cosh_kl*sinh_kl)/8. - (k*sin_kl*pow(sinh_kl,2))/16. - (3*k*cos_kl*sinh_2kl)/32. ;
            U[1][3][3][1] = (-7*cos_kl)/32. + (7*cos_kl*cosh_2kl)/32. - (k*lHC*sin_kl)/8. - (sin_kl*sinh_2kl)/32. ;
            U[1][5][5][0] = (-5*k2*lHC*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
            U[1][5][5][1] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8. ;
            U[2][2][2][2] = (-3*k3*lHC*sinh_kl)/16. + (3*k2*sinh_kl*sinh_2kl)/32. ;
            U[2][3][2][2] = (3*k2*lHC*cosh_kl)/16. - (21*k*sinh_kl)/64. + (3*k*sinh_3kl)/64. ;
            U[2][3][3][2] = (9*k*lHC*sinh_kl)/16. + (3*sinh_kl*sinh_2kl)/32. ;
            U[2][3][3][3] = (3*lHC*cosh_kl)/16. - (21*sinh_kl)/(64.*k) + (3*sinh_3kl)/(64.*k) ;
            U[2][5][5][2] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8. ;
            U[2][5][5][3] = (lHC*cosh_kl)/8. - sinh_kl/(8.*k) + (k*l2*sinh_kl)/8. ;
            U[3][2][2][2] = (-3*k4*lHC*cosh_kl)/16. - (3*k3*sinh_kl)/16. + (3*k3*cosh_2kl*sinh_kl)/16. + (3*k3*cosh_kl*sinh_2kl)/32. ;
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
            U[0][3][2][0] = (k2*lHC*cosh_kl)/4. - (k*cosh_kl*sin_2kl)/32. - (9*k*sinh_kl)/32. + (3*k*cos_2kl*sinh_kl)/32. ;
            U[0][3][2][1] = (-3*cosh_kl)/32. + (3*cos_2kl*cosh_kl)/32. + (k*lHC*sinh_kl)/4. - (sin_2kl*sinh_kl)/32. ;
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
            U[1][3][2][0] = -(k2*cosh_kl)/32. + (k2*cos_2kl*cosh_kl)/32. + (k3*lHC*sinh_kl)/4. - (7*k2*sin_2kl*sinh_kl)/32. ;
            U[1][3][2][1] = (k2*lHC*cosh_kl)/4. - (7*k*cosh_kl*sin_2kl)/32. + (5*k*sinh_kl)/32. + (k*cos_2kl*sinh_kl)/32. ;
            U[1][3][3][0] = (k2*lHC*cosh_kl)/8. - (k*cos_kl*cosh_kl*sin_kl)/8. + (3*k*cosh_kl*sin_2kl)/32. + (k*sinh_kl)/8. + (3*k*cos_2kl*sinh_kl)/16. - (k*pow(sin_kl,2)*sinh_kl)/16. ;
            U[1][3][3][1] = (-7*cosh_kl)/32. + (7*cos_2kl*cosh_kl)/32. + (k*lHC*sinh_kl)/8. + (sin_2kl*sinh_kl)/32. ;
            U[1][5][5][0] = (5*k2*lHC*cosh_kl)/8. + (3*k*sinh_kl)/8. + (k3*l2*sinh_kl)/8. ;
            U[1][5][5][1] = (k2*l2*cosh_kl)/8. + (3*k*lHC*sinh_kl)/8. ;
            U[2][2][2][2] = (-3*k3*lHC*sin_kl)/16. + (3*k2*sin_kl*sin_2kl)/32. ;
            U[2][3][2][2] = (-3*k2*lHC*cos_kl)/16. + (21*k*sin_kl)/64. - (3*k*sin_3kl)/64. ;
            U[2][3][3][2] = (-9*k*lHC*sin_kl)/16. - (3*sin_kl*sin_2kl)/32. ;
            U[2][3][3][3] = (3*lHC*cos_kl)/16. - (21*sin_kl)/(64.*k) + (3*sin_3kl)/(64.*k) ;
            U[2][5][5][2] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8. ;
            U[2][5][5][3] = (lHC*cos_kl)/8. - sin_kl/(8.*k) - (k*l2*sin_kl)/8. ;
            U[3][2][2][2] = (-3*k4*lHC*cos_kl)/16. - (3*k3*sin_kl)/16. + (3*k3*cos_2kl*sin_kl)/16. + (3*k3*cos_kl*sin_2kl)/32. ;
            U[3][3][2][2] = (9*k2*cos_kl)/64. - (9*k2*cos_3kl)/64. + (3*k3*lHC*sin_kl)/16. ;
            U[3][3][3][2] = (-9*k2*lHC*cos_kl)/16. - (9*k*sin_kl)/16. - (3*k*cos_2kl*sin_kl)/16. - (3*k*cos_kl*sin_2kl)/32. ;
            U[3][3][3][3] = (-9*cos_kl)/64. + (9*cos_3kl)/64. - (3*k*lHC*sin_kl)/16. ;
            U[3][5][5][2] = (-5*k2*lHC*cos_kl)/8. - (3*k*sin_kl)/8. + (k3*l2*sin_kl)/8. ;
            U[3][5][5][3] = -(k2*l2*cos_kl)/8. - (3*k*lHC*sin_kl)/8.  ;
          }
        }
      }

      if (lEdge && K1) {
        Md = NULL;
        Mtot = tmalloc(sizeof(*Mtot));
        initialize_matrices(Mtot, M->order);
        
        /* entrance fringe fields */
        Mfringe = quad_fringe(lEdge, K1, M->order, 0, 0.0);
        
        if (fixedStrengthFringe) {
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
        
        if (fixedStrengthFringe) {
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
    double *C, **R, ***T;
    double k, kl, K1L, cos_kl, sin_kl, cosh_kl, sinh_kl, l2, sin2_kl, sin_2kl, cos_2kl;
    double sinh_2kl, sinh2_kl, K1_1p5, K1_2;

    K1 *= (1+fse1);
    
    if (K1==0 || length==0) {
      if (K2==0)
        return drift_matrix(length, maximum_order);
      return sextupole_matrix(K2, length, maximum_order, tilt, fse2);
    }

    K2 *= (1+fse2);

    M = tmalloc(sizeof(*M));
    initialize_matrices(M, M->order = MIN(2,maximum_order));
    C = M->C;
    R = M->R;
    T = M->T;

    if (K1>0) {
      k = sqrt(K1);
      kl = k*length;
      K1L = K1*length;
      cos_kl = cos(kl);
      sin_kl = sin(kl);
      cosh_kl = cosh(kl);
      sinh_kl = sinh(kl);
      l2 = ipow(length, 2);
      sin2_kl = ipow(sin_kl, 2);
      sin_2kl = sin(2*kl);
      sinh_2kl = sinh(2*kl);
      sinh2_kl = ipow(sinh_kl, 2);
      
      K1_1p5 = pow(K1, 1.5);
      K1_2   = ipow(K1, 2);

      C[4] = length;
      
      R[0][0] = cos_kl;
      R[0][1] = sin_kl/k;
      R[1][0] = -(k*sin_kl);
      R[1][1] = cos_kl;
      R[2][2] = cosh_kl;
      R[2][3] = sinh_kl/k;
      R[3][2] = k*sinh_kl;
      R[3][3] = cosh_kl;
      R[4][4] = 1;
      R[5][5] = 1;
      if (M->order>1) {
        T[0][0][0] = -(K2*(K1*l2 + sin2_kl))/(8.*K1);
        T[0][1][0] = (K2*(-2*kl + sin_2kl))/(8.*K1_1p5);
        T[0][1][1] = (K2*(-(K1*l2) + sin2_kl))/(8.*K1_2);
        T[0][2][2] = (K2*(K1*l2 + sinh2_kl))/(8.*K1);
        T[0][3][2] = (K2*(-2*length + sinh_2kl/k))/(8.*K1);
        T[0][3][3] = (K2*(-(K1*l2) + sinh2_kl))/(8.*K1_2);
        T[0][5][0] = (kl*sin_kl)/2.;
        T[0][5][1] = (-(length*cos_kl) + sin_kl/k)/2.;
        T[1][0][0] = (K2*(-2*length - sin_2kl/k))/8.;
        T[1][1][0] = -(K2*sin2_kl)/(2.*K1);
        T[1][1][1] = (K2*(-2*kl + sin_2kl))/(8.*K1_1p5);
        T[1][2][2] = (K2*(2*length + sinh_2kl/k))/8.;
        T[1][3][2] = (K2*sinh2_kl)/(2.*K1);
        T[1][3][3] = (K2*(-2*kl + sinh_2kl))/(8.*K1_1p5);
        T[1][5][0] = (K1*length*cos_kl + k*sin_kl)/2.;
        T[1][5][1] = (kl*sin_kl)/2.;
        T[2][2][0] = (K2*sin_kl*sinh_kl)/(2.*K1);
        T[2][2][1] = (K2*(K1*length - k*cos_kl*sinh_kl))/(2.*K1_2);
        T[2][3][0] = (K2*(-(K1*length) + k*cosh_kl*sin_kl))/(2.*K1_2);
        if (kl<1e-3)
          T[2][3][1] = K2*(ipow(kl,4)/6)/(2.*K1_2);
        else
          T[2][3][1] = K2*(1 - cos_kl*cosh_kl)/(2.*K1_2);
        T[2][5][2] = -(kl*sinh_kl)/2.;
        T[2][5][3] = (-(length*cosh_kl) + sinh_kl/k)/2.;
        T[3][2][0] = (K2*(cosh_kl*sin_kl + cos_kl*sinh_kl))/(2.*k);
        T[3][2][1] = (K2 - K2*cos_kl*cosh_kl +  K2*sin_kl*sinh_kl)/(2.*K1);
        T[3][3][0] = (K2*(-1 + cos_kl*cosh_kl + sin_kl*sinh_kl))/(2.*K1);
        T[3][3][1] = (K2*(cosh_kl*sin_kl - cos_kl*sinh_kl))/(2.*K1_1p5);
        T[3][5][2] = (-(K1*length*cosh_kl) - k*sinh_kl)/2.;
        T[3][5][3] = -(kl*sinh_kl)/2.;
        T[4][0][0] = (2*K1*length - k*sin_2kl)/8.;
        T[4][1][0] = -sin2_kl/2.;
        T[4][1][1] = (2*length + sin_2kl/k)/8.;
        T[4][2][2] = (-2*K1*length + k*sinh_2kl)/8.;
        T[4][3][2] = sinh2_kl/2.;
        T[4][3][3] = (2*length + sinh_2kl/k)/8.;
      }
    } else {
      k = sqrt(-K1);
      kl = k*length;
      K1_2 = K1*K1;
      
      cos_kl = cos(kl);
      sin_kl = sin(kl);
      cosh_kl = cosh(kl);
      sinh_kl = sinh(kl);
      l2 = ipow(length, 2);
      sin2_kl = ipow(sin_kl, 2);
      sin_2kl = sin(2*kl);
      cos_2kl = cos(2*kl);
      sinh_2kl = sinh(2*kl);
      sinh2_kl = ipow(sinh_kl, 2);
      
      C[4] = length;

      R[0][0] = cosh_kl;
      R[0][1] = sinh_kl/k;
      R[1][0] = k*sinh_kl;
      R[1][1] = cosh_kl;
      R[2][2] = cos_kl;
      R[2][3] = sin_kl/k;
      R[3][2] = -(k*sin_kl);
      R[3][3] = cos_kl;
      R[4][4] = 1;
      R[5][5] = 1;
      if (M->order>1) {
        T[0][0][0] = (K2*(-(K1*l2) + sinh2_kl))/(8.*K1);
        T[0][1][0] = (K2*(-2*length + sinh_2kl/k))/(8.*K1);
        T[0][1][1] = -(K2*(K1*l2 + sinh2_kl))/(8.*K1_2);
        T[0][2][2] = (K2*(-1 + 2*K1*l2 + cos_2kl))/(16.*K1);
        T[0][3][2] = (K2*(-2*kl + sin_2kl))/(8.*K1*k);
        T[0][3][3] = (K2*(-(K1*l2) - sin2_kl))/(8.*K1_2);
        T[0][5][0] = -(kl*sinh_kl)/2.;
        T[0][5][1] = (-(length*cosh_kl) + sinh_kl/k)/2.;
        T[1][0][0] = (K2*(-2*K1*length + k*sinh_2kl))/(8.*K1);
        T[1][1][0] = (K2*sinh2_kl)/(2.*K1);
        T[1][1][1] = -(K2*(2*K1*length + k*sinh_2kl))/(8.*K1_2);
        T[1][2][2] = (K2*(2*K1*length - k*sin_2kl))/(8.*K1);
        T[1][3][2] = -(K2*sin2_kl)/(2.*K1);
        T[1][3][3] = (K2*(-2*kl + sin_2kl))/(8.*K1*k);
        T[1][5][0] = (K1*length*cosh_kl - k*sinh_kl)/2.;
        T[1][5][1] = -(kl*sinh_kl)/2.;
        T[2][2][0] = -(K2*sinh_kl*sin_kl)/(2.*K1);
        T[2][2][1] = (K2*(K1*length + k*cosh_kl*sin_kl))/(2.*K1_2);
        T[2][3][0] = (K2*(-(K1*length) - k*cos_kl*sinh_kl))/(2.*K1_2);
        T[2][3][1] = (K2 - K2*cosh_kl*cos_kl)/(2.*K1_2);
        T[2][5][2] = (kl*sin_kl)/2.;
        T[2][5][3] = (-(length*cos_kl) + sin_kl/k)/2.;
        T[3][2][0] = (K2*(cos_kl*sinh_kl + cosh_kl*sin_kl))/(2.*k);
        T[3][2][1] = (K2 - K2*cosh_kl*cos_kl - K2*sinh_kl*sin_kl)/(2.*K1);
        T[3][3][0] = (K2*(-1 + cosh_kl*cos_kl - sinh_kl*sin_kl))/(2.*K1);
        T[3][3][1] = (K2*(cos_kl*sinh_kl - cosh_kl*sin_kl))/(2.*K1*k);
        T[3][5][2] = (-(K1*length*cos_kl) + k*sin_kl)/2.;
        T[3][5][3] = (kl*sin_kl)/2.;
        T[4][0][0] = (2*K1*length + k*sinh_2kl)/8.;
        T[4][1][0] = sinh2_kl/2.;
        T[4][1][1] = (2*length + sinh_2kl/k)/8.;
        T[4][2][2] = (-2*K1*length - k*sin_2kl)/8.;
        T[4][3][2] = -sin2_kl/2.;
        T[4][3][3] = (2*length + sin_2kl/k)/8.;
      }
    }
    
    tilt_matrices(M, tilt);
    return(M);
    }
