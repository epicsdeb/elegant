/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: twiss.c
 * purpose: computation of Twiss parameters
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "matlib.h"
#include "chromDefs.h"
#include "fftpackC.h"
#include "twiss.h"
#include <stddef.h>

void computeDrivingTerms(DRIVING_TERMS *drivingTerms, ELEMENT_LIST *eptr, TWISS *twiss0, double *tune);
void copy_doubles(double *target, double *source, long n);
double find_acceptance(ELEMENT_LIST *elem, long plane, RUN *run, char **name, double *end_pos);
void modify_rfca_matrices(ELEMENT_LIST *eptr, long order);
void reset_rfca_matrices(ELEMENT_LIST *eptr, long order);
void incrementRadIntegrals(RADIATION_INTEGRALS *radIntegrals, double *dI, 
                           ELEMENT_LIST *elem, 
                           double beta0, double alpha0, double gamma0, double eta0, double etap0, 
                           double betay0, double alphay0, double gammay0, double etay0, double etapy0, 
                           double *coord);
void AddWigglerRadiationIntegrals(double length, long periods, double radius,
				  double eta, double etap, 
				  double beta, double alpha,
				  double *I1, double *I2, double *I3, double *I4, double *I5);
void LoadStartingTwissFromFile(double *betax, double *betay, double *alphax, double *alphay,
                               double *etax, double *etay, double *etaxp, double *etayp,
                               char *filename, char *elementName, long elementOccurrence);
void computeTuneShiftWithAmplitude(double dnuxdA[N_TSWA][N_TSWA], double dnuydA[N_TSWA][N_TSWA],
				   double *nuxExtrema, double *nuyExtrema,
                                   TWISS *twiss, double *tune, VMATRIX *M, LINE_LIST *beamline,
				   RUN *run, double *startingCoord);
void computeTuneShiftWithAmplitudeM(double dnuxdA[N_TSWA][N_TSWA], double dnuydA[N_TSWA][N_TSWA],
                                    TWISS *twiss, double *tune, VMATRIX *M);
void processTwissAnalysisRequests(ELEMENT_LIST *elem);
double adjustTuneHalfPlane(double frequency, double phase0, double phase1);

static long twissConcatOrder = 3;
static long doTuneShiftWithAmplitude = 0;
static SDDS_DATASET SDDSTswaTunes;
static long linearChromaticTrackingInitialized = 0;
void setLinearChromaticTrackingValues(LINE_LIST *beamline) ;

#define TWISS_ANALYSIS_QUANTITIES 8
static char *twissAnalysisQuantityName[TWISS_ANALYSIS_QUANTITIES] = {"betax", "betay", "etax", "etay", "alphax", "alphay", "etaxp", "etayp"};
static long twissAnalysisQuantityOffset[TWISS_ANALYSIS_QUANTITIES] = {
  offsetof(TWISS, betax), offsetof(TWISS, betay), offsetof(TWISS, etax), offsetof(TWISS, etay),
  offsetof(TWISS, alphax), offsetof(TWISS, alphay), 
  offsetof(TWISS, etapx), offsetof(TWISS, etapx), 
  };
#define TWISS_ANALYSIS_AVE 0
#define TWISS_ANALYSIS_MIN 1
#define TWISS_ANALYSIS_MAX 2
#define TWISS_ANALYSIS_STATS 3
static char *twissAnalysisStatName[TWISS_ANALYSIS_STATS] = {"ave", "min", "max" };
static long twissAnalysisStatCode[TWISS_ANALYSIS_STATS] = {
  TWISS_ANALYSIS_AVE, TWISS_ANALYSIS_MIN, TWISS_ANALYSIS_MAX };

typedef struct {
  char *startName, *endName, *matchName, *tag;
  double sStart, sEnd;
  long startOccurence, endOccurence;
  short initialized;
  long count;
  double twissMem[TWISS_ANALYSIS_STATS][TWISS_ANALYSIS_QUANTITIES];
} TWISS_ANALYSIS_REQUEST;

static long twissAnalysisRequests = 0;
static TWISS_ANALYSIS_REQUEST *twissAnalysisRequest = NULL;

static short mustResetRfcaMatrices = 0;

VMATRIX *compute_periodic_twiss(
                                double *betax, double *alphax, double *etax, double *etapx, double *NUx,
                                double *betay, double *alphay, double *etay, double *etapy, double *NUy,
                                ELEMENT_LIST *elem, double *clorb, RUN *run, unsigned long *unstable,
                                double *eta2, double *eta3)
{
  VMATRIX *M, *M1;
  double cos_phi, sin_phi, **R, beta[2], alpha[2], phi[2];
  double ***T, ****Q, eta2f[6];
  long i, j, k;
  MATRIX *dispR, *dispM, *dispMInv, *dispEta;
  
  log_entry("compute_periodic_twiss");

  *unstable = 0;

  if ((i = fill_in_matrices(elem, run)))
    fprintf(stdout, "%ld matrices recomputed for periodic Twiss parameter computation\n", i);
    fflush(stdout);

  if (cavities_are_drifts_if_matched) {
    if (run->always_change_p0)
      bomb("can't have run_setup/always_change_p0=1 and twiss_output/cavities_are_drifts_if_matched=1", NULL);
    modify_rfca_matrices(elem, run->default_order);  /* replace rf cavities with drifts */
  }
  
  if (clorb) {
    /* use the closed orbit to compute the on-orbit R matrix */
    M1 = tmalloc(sizeof(*M1));
    initialize_matrices(M1, 1);
    for (i=0; i<6; i++) {
      M1->C[i] = clorb[i];
      M1->R[i][i] = 1;
    }
    M = append_full_matrix(elem, run, M1, twissConcatOrder);
    free_matrices(M1);
/*
    fprintf(stdout, "Computed revolution matrix on closed orbit to %ld order\n",
            twissConcatOrder);
    fflush(stdout);
    fprintf(stdout, "matrix concatenation for periodic Twiss computation:\n");
    fflush(stdout);
    fprintf(stdout, "closed orbit at input:\n  ");
    fflush(stdout);
    for (i=0; i<6; i++)
      fprintf(stdout, "%14.6e ", clorb[i]);
      fflush(stdout);
    fprintf(stdout, "\nclosed orbit at output:\n  ");
    fflush(stdout);
    for (i=0; i<6; i++)
      fprintf(stdout, "%14.6e ", M->C[i]);
      fflush(stdout);
    fprintf(stdout, "\nR matrix:\n");
    fflush(stdout);
    for (i=0; i<6; i++) {
      fprintf(stdout, "  ");
      fflush(stdout);
      for (j=0; j<6; j++)
        fprintf(stdout, "%14.6e ", M->R[i][j]);
        fflush(stdout);
      fprintf(stdout, "\n");
      fflush(stdout);
    }
*/
  }
  else {
    M = full_matrix(elem, run, twissConcatOrder);
/*
    fprintf(stdout, "Computed revolution matrix to %ld order\n", twissConcatOrder);
    fflush(stdout);
*/
  }

  R = M->R;
  T = M->T;
  Q = M->Q;
  
  /* allocate matrices for computing dispersion, which I do
   * in 4-d using 
   * eta[i] = Inv(I - R)[i][j] R[j][5]
   * eta2[i] = Inv(I-R)[i][j] Sum[0<=k<=j<=5] T[i][j][k]*eta[i]*eta[k]
   *  with eta[4]=0 and eta[5]=1
   * The dispersion is defined by, for example, 
   *  x = x(delta=0) + delta*eta + delta^2*eta2 + delta^3*eta3 ... 
   */
  m_alloc(&dispM, 4, 4);
  m_alloc(&dispMInv, 4, 4);
  m_alloc(&dispR, 4, 1);
  m_alloc(&dispEta, 4, 1);
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      dispM->a[i][j] = (i==j?1:0) - R[i][j];
    }
    dispR->a[i][0] = R[i][5];
  }
  if (m_det(dispM)==0) {
    fprintf(stdout, "Unable to compute dispersion: unstable\n");
    fflush(stdout);
    *unstable = 3; /* both planes */
  } else {
    double *eta;
    m_invert(dispMInv, dispM);
    m_mult(dispEta, dispMInv, dispR);
    eta = tmalloc(sizeof(*eta)*6);
    eta[0] = *etax = dispEta->a[0][0];
    eta[1] = *etapx = dispEta->a[1][0];
    eta[2] = *etay = dispEta->a[2][0];    
    eta[3] = *etapy = dispEta->a[3][0];
    eta[4] = 0;
    eta[5] = 1;
    if (eta2) {
      /* now do the second-order dispersion */
      for (i=0; i<4; i++) {
        dispR->a[i][0] = 0;
        if (T) {
          for (j=0; j<6; j++) 
            for (k=0; k<=j; k++) 
              dispR->a[i][0] += T[i][j][k]*eta[j]*eta[k];
        }
      }
      m_mult(dispEta, dispMInv, dispR);
      for (i=0; i<4; i++)
        eta2f[i] = eta2[i] = dispEta->a[i][0];
      eta2f[4] = 0;
      eta2f[5] = 0;

      if (eta3) {
        long l;
        /* third-order dispersion */
        for (i=0; i<4; i++) {
          dispR->a[i][0] = 0;
          if (T) {
            for (j=0; j<6; j++) 
              for (k=0; k<=j; k++)  {
                dispR->a[i][0] += T[i][j][k]*(eta2f[j]*eta[k]+eta[j]*eta2f[k]);
                if (Q)
                  for (l=0; l<=k; l++) 
                    dispR->a[i][0] += Q[i][j][k][l]*eta[j]*eta[k]*eta[l];
              }
          }
        }
        m_mult(dispEta, dispMInv, dispR);
        for (i=0; i<4; i++)
          eta3[i] = dispEta->a[i][0];
      }
    }
    free(eta);
  }
  m_free(&dispM);
  m_free(&dispEta);
  m_free(&dispR);
  m_free(&dispMInv);
  
  for (i=0; i<4; i+=2 ) {
    if (fabs(cos_phi = (R[i][i] + R[i+1][i+1])/2)>1) {
      fprintf(stdout, "warning: beamline unstable for %c plane--can't match beta functions.\n", (char)('x'+i/2));
      fflush(stdout);
      *unstable |= (i==0?1:2);
      sin_phi = 1e-6;
    }
    beta[i/2] = fabs(R[i][i+1]/sin(acos(cos_phi)));
    sin_phi   = R[i][i+1]/beta[i/2];
    phi[i/2]  = atan2(sin_phi, cos_phi);
    if (phi[i/2]<0)
      phi[i/2] += PIx2;
    alpha[i/2] = (R[i][i]-R[i+1][i+1])/(2*sin_phi);
  }
  
  *betax = beta[0];
  *betay = beta[1];
  *alphax = alpha[0];
  *alphay = alpha[1];
  *NUx = phi[0]/PIx2;
  *NUy = phi[1]/PIx2;

  log_exit("compute_periodic_twiss");
  return(M);
}

void propagate_twiss_parameters(TWISS *twiss0, double *tune, long *waists,
                                RADIATION_INTEGRALS *radIntegrals, 
                                ELEMENT_LIST *elem,  RUN *run, double *traj,
				double *couplingFactor
				)
{
  double beta[2], alpha[2], phi[2], eta[2], etap[2], gamma[2], refAlpha[2];
  double *func, path[6], path0[6], detR[2], length, sTotal;
  double **R=NULL, C[2], S[2], Cp[2], Sp[2], D[2], Dp[2], sin_dphi, cos_dphi, dphi;
  long n_mat_computed, i, j, plane, otherPlane, hasMatrix;
  VMATRIX *M1, *M2;
  MATRIX *dispM, *dispOld, *dispNew;
  ELEMENT_LIST *elemOrig;
  static long asinWarning = 50;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  double complex kappa;  
#else
  doublecomplex_sdds kappa;  
#endif

  if (!twiss0)
    bomb("initial Twiss parameters not given (propagate_twiss_parameters())", NULL);
  elemOrig = elem;
  
  m_alloc(&dispM, 4, 4);
  m_alloc(&dispOld, 4, 1);
  m_alloc(&dispNew, 4, 1);
  
  for (plane=0; plane<2; plane++) {
    beta[plane]  = *(&twiss0->betax+plane*TWISS_Y_OFFSET);
    refAlpha[plane] = alpha[plane] = *(&twiss0->alphax+plane*TWISS_Y_OFFSET);
    phi[plane]   = *(&twiss0->phix+plane*TWISS_Y_OFFSET);
    eta[plane]   = *(&twiss0->etax+plane*TWISS_Y_OFFSET);
    etap[plane]  = *(&twiss0->etapx+plane*TWISS_Y_OFFSET);
  }

  M1 = tmalloc(sizeof(*M1));
  M2 = tmalloc(sizeof(*M2));
  initialize_matrices(M1, 1);
  initialize_matrices(M2, twissConcatOrder);
  if (traj) {
    for (i=0; i<6; i++) {
      path[i] = traj[i];
      M1->R[i][i] = 1;
    }
  }
  else {
    for (i=0; i<6; i++) {
      path[i] = 0;
      M1->R[i][i] = 1;
    }
  }
  
  n_mat_computed = 0;

/*
  fprintf(stdout, "Twiss parameter computation on path %e, %e, %e, %e, %e, %e\n",
          path[0], path[1], path[2], path[3], path[4], path[5]);
*/
  
  if (radIntegrals) {
    for (i=0; i<6; i++) 
      radIntegrals->RI[i] = 0;
  }
  waists[0] = waists[1] = 0;

  elem = elemOrig;
  sTotal = 0;
  while (elem) {
    for (plane=0; plane<2; plane++) 
      gamma[plane] = (1+sqr(alpha[plane]))/beta[plane];
    if (entity_description[elem->type].flags&HAS_MATRIX) {
      if (!elem->matrix || !(elem->matrix->R) || 
          (elem->pred && elem->pred->Pref_output!=elem->Pref_input) ||
          elem->type==T_TWISSELEMENT) {
        if (elem->type==T_TWISSELEMENT) {
          if (!((TWISSELEMENT*)elem->p_elem)->computeOnce || !((TWISSELEMENT*)elem->p_elem)->transformComputed) {
            TWISS twissInput;
            twissInput.betax = beta[0];
            twissInput.betay = beta[1];
            twissInput.alphax = alpha[0];
            twissInput.alphay = alpha[1];
            twissInput.etax   = eta[0];
            twissInput.etapx  = etap[0];
            twissInput.etay   = eta[1];
            twissInput.etapy  = etap[1];
            if (elem->matrix) {
              free_matrices(elem->matrix);
              free(elem->matrix);
            }
            if (((TWISSELEMENT*)elem->p_elem)->verbose) {
              printf("Computing twiss transformation matrix for %s at z=%e m from lattice twiss parameters\n", elem->name, elem->end_pos);
              printf("  * Initial twiss parameters:\n");
              printf("  betax = %le  alphax = %le  etax = %le, etaxp = %le\n",
                     twissInput.betax, twissInput.alphax, twissInput.etax, twissInput.etapx);
              printf("  betay = %le  alphay = %le  etay = %le, etayp = %le\n",
                     twissInput.betay, twissInput.alphay, twissInput.etay, twissInput.etapy);
              printf("  * Final twiss parameters:\n");
              printf("  betax = %le  alphax = %le  etax = %le, etaxp = %le\n",
                     ((TWISSELEMENT*)elem->p_elem)->twiss.betax,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.alphax,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.etax,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.etapx);
              printf("  betax = %le  alphax = %le  etax = %le, etaxp = %le\n",
                     ((TWISSELEMENT*)elem->p_elem)->twiss.betay,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.alphay,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.etay,
                     ((TWISSELEMENT*)elem->p_elem)->twiss.etapy);
            }
            elem->matrix = twissTransformMatrix((TWISSELEMENT*)elem->p_elem,
                                                &twissInput);
            if (((TWISSELEMENT*)elem->p_elem)->fromBeam) 
              /* If user wants the transform from the beam, we don't regard the transform as having
               * been computed in final form.  Hence, we set this flag to zero.  This allows twiss and moments
               * computations to be done with a possibly valid matrix prior to tracking
               */
              ((TWISSELEMENT*)elem->p_elem)->transformComputed = 0;
            else 
              ((TWISSELEMENT*)elem->p_elem)->transformComputed = 1;
          }
        } else {
          if (elem->matrix) {
            free_matrices(elem->matrix);
            free(elem->matrix);
          }
          elem->matrix = compute_matrix(elem, run, NULL);
        }
        n_mat_computed++;
      }
      hasMatrix = 1;
      /* use matrix concatenation to include effect of beam path */
      for (i=0; i<6; i++) 
        path0[i] = M1->C[i] = path[i];
      concat_matrices(M2, elem->matrix, M1,
                      entity_description[elem->type].flags&HAS_RF_MATRIX?
                      CONCAT_EXCLUDE_S0:0);
      R = M2->R;
      /* record new centroids for beam path */
      for (i=0; i<6; i++)
        path[i] = M2->C[i];
      for (plane=0; plane<2; plane++) {
        C[plane]  = R[0+2*plane][0+2*plane]; 
        S[plane]  = R[0+2*plane][1+2*plane];
        Cp[plane] = R[1+2*plane][0+2*plane]; 
        Sp[plane] = R[1+2*plane][1+2*plane];
        D[plane]  = R[0+2*plane][5]; 
        Dp[plane] = R[1+2*plane][5];
      }
    }
    else {
      hasMatrix = 0;
      if (elem->pred)
        elem->Pref_input = elem->Pref_output = elem->pred->Pref_output;
      for (plane=0; plane<2; plane++) {
        C[plane] = Sp[plane] = 1;
        Cp[plane] = D[plane] = Dp[plane] = 0;
      }
      if (entity_description[elem->type].flags&HAS_LENGTH) {
        S[0] = S[1] = *((double*)elem->p_elem);
        path[0] += path[1]*S[0];
        path[2] += path[3]*S[1];
      } else {
        S[0] = S[1] = 0;
      }
    }
    if (!elem->twiss)
      elem->twiss = tmalloc(sizeof(*elem->twiss));
    elem->twiss->periodic = matched;
    if (radIntegrals)
      incrementRadIntegrals(radIntegrals, elem->twiss->dI, 
                            elem, 
                            beta[0], alpha[0], gamma[0], eta[0], etap[0], 
                            beta[1], alpha[1], gamma[1], eta[1], etap[1], 
                            path0);
    if (elem->type==T_ROTATE) {
      if (fabs(((ROTATE*)elem->p_elem)->tilt-PI/2.0)<1e-6 ||
          fabs(((ROTATE*)elem->p_elem)->tilt-3*PI/2.0)<1e-6 ||
          fabs(((ROTATE*)elem->p_elem)->tilt+PI/2.0)<1e-6 ||
          fabs(((ROTATE*)elem->p_elem)->tilt+3*PI/2.0)<1e-6) {
        elem->twiss->betax = beta[1];
        elem->twiss->betay = beta[0];
        elem->twiss->alphax = alpha[1];
        elem->twiss->alphay = alpha[0];
        elem->twiss->phix = phi[1];
        elem->twiss->phiy = phi[0];
        elem->twiss->etax = eta[1];
        elem->twiss->etay = eta[0];
        elem->twiss->etapx = etap[1];
        elem->twiss->etapy = etap[0];
        SWAP_DOUBLE(beta[0], beta[1]);
        SWAP_DOUBLE(alpha[0], alpha[1]);
        SWAP_DOUBLE(eta[0], eta[1]);
        SWAP_DOUBLE(etap[0], etap[1]);
        SWAP_DOUBLE(phi[0], phi[1]);
        elem = elem->succ;
        continue;
      }
    }
    for (plane=0; plane<2; plane++) {
      otherPlane = plane?0:1;
      detR[plane] = C[plane]*Sp[plane] - Cp[plane]*S[plane];
      /* set up pointers to Twiss functions */
      func = ((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0);
      /* store centroid position */
      *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET) = path[plane?2:0];
      /* calculate new beta and alpha */
      if (elem->type==T_REFLECT) {
        func[0] = beta[plane];
        func[1] = -alpha[plane];
        dphi = 0;
      } else {
        func[0] = (sqr(C[plane])*beta[plane] - 2*C[plane]*S[plane]*alpha[plane] 
                   + sqr(S[plane])*gamma[plane])/detR[plane];
        func[1] = (-C[plane]*Cp[plane]*beta[plane] + 
                   (Sp[plane]*C[plane]+S[plane]*Cp[plane])*alpha[plane] - 
                   S[plane]*Sp[plane]*gamma[plane])/detR[plane];
        /* use R12=S to find sin(dphi) */
        if ((sin_dphi=S[plane]/sqrt(beta[plane]*func[0]))>1) {
          if (asinWarning>0) {
            asinWarning--;
            fprintf(stdout, "warning: argument of asin > 1 by %f (propagate_twiss)\n", sin_dphi-1);
            fprintf(stdout, "element is %s at z=%em\n", elem->name, elem->end_pos);
            fprintf(stdout, "%c-plane matrix:  C = %e,  S = %e,  ",
                    (plane==0?'x':'y'), C[plane], S[plane]);
            fprintf(stdout, "C' = %e,  S' = %e\n", Cp[plane], Sp[plane]);
            fprintf(stdout, "beta0 = %e, func[0] = %e\n", beta[plane], func[0]);
            fflush(stdout);
          }
          sin_dphi = 1;
          cos_dphi = 0;
        }
        else if (sin_dphi<-1) {
          if (asinWarning>0) {
            asinWarning--;
            fprintf(stdout, "warning: argument of asin < -1 by %f (propagate_twiss)\n", sin_dphi+1);
            fprintf(stdout, "element is %s at z=%em\n", elem->name, elem->end_pos);
            fflush(stdout);
          }
          sin_dphi = -1;
          cos_dphi = 0;
        }
        else {
          /* use R11=C to find cos(dphi) */
          cos_dphi = sqrt(beta[plane]/func[0])*C[plane] - alpha[plane]*sin_dphi;
        }
        if (entity_description[elem->type].flags&HAS_LENGTH) {
          if (*((double*)elem->p_elem)>=0) {
            if ((dphi = atan2(sin_dphi, cos_dphi))<0) 
              dphi += PIx2;
          } else {
            if ((dphi=atan2(sin_dphi, cos_dphi))>0)
              dphi -= PIx2;
          }
        } else {
          if ((dphi = atan2(sin_dphi, cos_dphi))<0) 
            dphi += PIx2;
        }
      }
      
      phi[plane]  = func[2] = phi[plane] + dphi;
      beta[plane]  = fabs(func[0]);
    
      alpha[plane] = func[1];
      if (SIGN(alpha[plane])!=SIGN(refAlpha[plane]) && refAlpha[plane]!=0 && alpha[plane]!=0) {
        /* if sign of alpha changes, it is called a "waist". */
        waists[plane] ++;
        refAlpha[plane] = alpha[plane];
      }
    }
    
    /* compute dispersion function and slope */
    if (hasMatrix) {
      for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
          dispM->a[i][j] = R[i][j];
        }
      }
    } else {
      for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
          dispM->a[i][j] = i==j?1:0;
        }
      }
      dispM->a[0][1] = S[0];
      dispM->a[2][3] = S[1];
    }
    dispOld->a[0][0] = eta[0];
    dispOld->a[1][0] = etap[0];
    dispOld->a[2][0] = eta[1];
    dispOld->a[3][0] = etap[1];
    
    m_mult(dispNew, dispM, dispOld);
    
    plane = 0;
    func = ((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0);
    eta[plane] = func[3] = dispNew->a[0][0] + (hasMatrix?R[0][5]:0);
    etap[plane] = func[4] = dispNew->a[1][0] + (hasMatrix?R[1][5]:0);
    plane = 1;
    func = ((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0);
    eta[plane] = func[3] = dispNew->a[2][0] + (hasMatrix?R[2][5]:0);
    etap[plane] = func[4] = dispNew->a[3][0] + (hasMatrix?R[3][5]:0);

    if (elem->type==T_MARK && ((MARK*)elem->p_elem)->fitpoint)
      store_fitpoint_twiss_parameters((MARK*)elem->p_elem, elem->name, elem->occurence,
                                      elem->twiss);

    sTotal = elem->end_pos;
    elem = elem->succ;
  }
  
  tune[0] = phi[0]/PIx2;
  tune[1] = phi[1]/PIx2;

  if (couplingFactor) {
    /* Compute linear coupling based on 187 of Handbook of Accelerator Physics and Engineering 
     * and P.J. Bryant, "A Simple Theory for Weak Betatron Coupling", CERN 94-01, Vol 1., 207-217.
     */
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
    double complex integrand, phaseFactor;
#else
    doublecomplex_sdds integrand, phaseFactor, tmp;
    double tmp2;
#endif
    double ks, phase, y, K1, K1r, tilt;
    long q;
#ifdef DEBUG_COUPLING
    FILE *fp;
    fp = fopen("kappa.sdds", "w");
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&column name=ElementName type=string &end\n");
    fprintf(fp, "&column name=Position, type=double &end\n");
    fprintf(fp, "&column name=Length , type=double &end\n");
    fprintf(fp, "&column name=K1 , type=double &end\n");
    fprintf(fp, "&column name=Tilt , type=double &end\n");
    fprintf(fp, "&column name=ks , type=double &end\n");
    fprintf(fp, "&column name=phase type=double &end\n");
    fprintf(fp, "&column name=exp0, type=double &end\n");
    fprintf(fp, "&column name=exp1 , type=double &end\n");
    fprintf(fp, "&column name=i0 , type=double &end\n");
    fprintf(fp, "&column name=i1 , type=double &end\n");
    fprintf(fp, "&column name=kappa0 , type=double &end\n");
    fprintf(fp, "&column name=kappa1 , type=double &end\n");
    fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
#endif
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
    kappa = 0;   /* coupling factor */
#else
    kappa.r = 0;   /* coupling factor */
    kappa.i = 0;
#endif
    elem = elemOrig;
    q = tune[0] - tune[1] + 0.5;
    couplingFactor[1] = (tune[0] - tune[1]) - q;
    while (elem) {
      if ((elem->type==T_QUAD || elem->type==T_KQUAD || elem->type==T_SEXT ||
	   elem->type==T_KSEXT || elem->type==T_SOLE) && (length = *((double*)elem->p_elem))) {
        if (elem->pred) {
          beta[0] = (elem->twiss->betax + elem->pred->twiss->betax)/2;
          beta[1] = (elem->twiss->betay + elem->pred->twiss->betay)/2;
          alpha[0] = (elem->twiss->alphax + elem->pred->twiss->alphax)/2;
          alpha[1] = (elem->twiss->alphay + elem->pred->twiss->alphay)/2;
          phi[0]  = (elem->twiss->phix + elem->pred->twiss->phix)/2;
          phi[1]  = (elem->twiss->phiy + elem->pred->twiss->phiy)/2;
          y = (elem->twiss->Cy + elem->pred->twiss->Cy)/2;
        }
        else {
          beta[0] = (elem->twiss->betax + twiss0->betax)/2;
          beta[1] = (elem->twiss->betay + twiss0->betay)/2;
          alpha[0] = (elem->twiss->alphax + twiss0->alphax)/2;
          alpha[1] = (elem->twiss->alphay + twiss0->alphay)/2;
          phi[0]  = (elem->twiss->phix + twiss0->phix)/2;
          phi[1]  = (elem->twiss->phiy + twiss0->phiy)/2;
          y = (elem->twiss->Cy + twiss0->Cy)/2;
        }
	ks = K1r = K1 = 0;
	switch (elem->type) {
	case T_QUAD:
	  K1r = (K1=((QUAD*)(elem->p_elem))->k1) * sin(2*(tilt=((QUAD*)(elem->p_elem))->tilt));
	  break;
	case T_KQUAD:
	  K1r = (K1=((KQUAD*)(elem->p_elem))->k1) * sin(2*(tilt=((KQUAD*)(elem->p_elem))->tilt));
	  break;
	case T_SEXT:
	  K1r = ((SEXT*)(elem->p_elem))->k2 * (y - ((SEXT*)(elem->p_elem))->dy);
	  break;
	case T_KSEXT:
	  K1r = ((KSEXT*)(elem->p_elem))->k2 * (y - ((SEXT*)(elem->p_elem))->dy);
	  break;
        case T_SOLE:
          ks = -((SOLE*)(elem->p_elem))->ks;
          break;
	}
        if (K1r!=0 || ks!=0) {
          phase = phi[0] - phi[1] 
            - (tune[0] - tune[1] - q)*2*PI*
              (elem->end_pos-length/2)/sTotal;
          phaseFactor = cexpi(phase);
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
          integrand = K1r + (alpha[0]/beta[0] - alpha[1]/beta[1])*ks/2 - I*(1/beta[0]+1/beta[1])*ks/2;
          kappa = kappa + integrand*phaseFactor*sqrt(beta[0]*beta[1])/PIx2*length;
#else
          integrand.r = K1r + (alpha[0]/beta[0] - alpha[1]/beta[1])*ks/2;
          integrand.i = -(1/beta[0]+1/beta[1])*ks/2;

	  tmp.r = (integrand.r * phaseFactor.r - integrand.i * phaseFactor.i);
          tmp.i = (integrand.i * phaseFactor.r + integrand.r * phaseFactor.i);
          tmp2 = sqrt(beta[0]*beta[1])/PIx2*length;
          tmp.r = tmp.r * tmp2;
          tmp.i = tmp.i * tmp2;
	  kappa.r = kappa.r + tmp.r;
          kappa.i = kappa.i + tmp.i;
#endif

#ifdef DEBUG_COUPLING
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
          fprintf(fp, "%s %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  elem->name, elem->end_pos-length/2, length, K1r, tilt, ks,
                  phase, creal(phaseFactor), cimag(phaseFactor), creal(integrand), cimag(integrand),
                  creal(kappa), cimag(kappa));
#else
          fprintf(fp, "%s %e %e %e %e %e %e %e %e %e %e %e %e\n",
                  elem->name, elem->end_pos-length/2, length, K1r, tilt, ks,
                  phase, phaseFactor.r, phaseFactor.i, integrand.r, integrand.i,
                  kappa.r, kappa.i);
#endif
#endif

        }
      }
      elem = elem->succ;
    }

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
    if ((couplingFactor[0] = cabs(kappa)))
#else
    if (couplingFactor[0] = sqrt(kappa.r * kappa.r + kappa.i * kappa.i))
#endif
      couplingFactor[2] = sqr(couplingFactor[0])/(sqr(couplingFactor[0]) + sqr(couplingFactor[1]));
    else
      couplingFactor[2] = 0;

#ifdef DEBUG_COUPLING
    fclose(fp);
#endif
  }

  if (radIntegrals)
    radIntegrals->computed = 1;
  
/*
    fprintf(stdout, "beta, eta, alpha: %e, %e; %e, %e; %e, %e\n",
          beta[0], beta[1],
          eta[0], eta[1], 
          alpha[0], alpha[1]);
  if (radIntegrals)
    fprintf(stdout, "Radiation integrals: %e, %e, %e, %e, %e\n",
            radIntegrals->RI[0], 
            radIntegrals->RI[1], 
            radIntegrals->RI[2], 
            radIntegrals->RI[3], 
            radIntegrals->RI[4]);
*/
  
  m_free(&dispNew);
  m_free(&dispM);
  m_free(&dispOld);
  
  free_matrices(M1); tfree(M1);
  free_matrices(M2); tfree(M2);

  processTwissAnalysisRequests(elemOrig);
}


static long twiss_initialized = 0;
static long SDDS_twiss_initialized = 0;
static SDDS_TABLE SDDS_twiss;
static long twiss_count = 0;

#define IC_S 0
#define IC_BETAX 1
#define IC_ALPHAX 2
#define IC_PHIX 3
#define IC_ETAX 4
#define IC_ETAPX 5
#define IC_APX 6
#define IC_BETAY 7
#define IC_ALPHAY 8
#define IC_PHIY 9
#define IC_ETAY 10
#define IC_ETAPY 11
#define IC_APY 12
#define IC_PCENTRAL 13
#define N_DOUBLE_COLUMNS 14
#define IC_ELEMENT 14
#define IC_OCCURENCE 15
#define IC_TYPE 16
#define N_COLUMNS 17
#define IC_I1 N_COLUMNS
#define IC_I2 (N_COLUMNS+1)
#define IC_I3 (N_COLUMNS+2)
#define IC_I4 (N_COLUMNS+3)
#define IC_I5 (N_COLUMNS+4)
#define N_COLUMNS_WRI (IC_I5+1)
static SDDS_DEFINITION column_definition[N_COLUMNS_WRI] = {
{"s", "&column name=s, type=double, units=m, description=Distance &end"},
{"betax", "&column name=betax, type=double, units=m, symbol=\"$gb$r$bx$n\", description=\"Horizontal beta-function\" &end"},
{"alphax", "&column name=alphax, type=double, symbol=\"$ga$r$bx$n\", description=\"Horizontal alpha-function\" &end"},
{"psix", "&column name=psix, type=double, units=rad, symbol=\"$gy$r$bx$n\", description=\"Horizontal phase advance\" &end"},
{"etax", "&column name=etax, type=double, units=m, symbol=\"$gc$r$bx$n\", description=\"Horizontal dispersion\" &end"},
{"etaxp", "&column name=etaxp, type=double, symbol=\"$gc$r$bx$n$a'$n\", description=\"Slope of horizontal dispersion\" &end"},
{"xAperture", "&column name=xAperture, type=double, units=m, symbol=\"a$bx,eff$n\", description=\"Effective horizontal aperture\" &end"},
{"betay", "&column name=betay, type=double, units=m, symbol=\"$gb$r$by$n\", description=\"Vertical beta-function\" &end"},
{"alphay", "&column name=alphay, type=double, symbol=\"$ga$r$by$n\", description=\"Vertical alpha-function\" &end"},
{"psiy", "&column name=psiy, type=double, units=rad, symbol=\"$gy$r$by$n\", description=\"Vertical phase advance\" &end"},
{"etay", "&column name=etay, type=double, units=m, symbol=\"$gc$r$by$n\", description=\"Vertical dispersion\" &end"},
{"etayp", "&column name=etayp, type=double, symbol=\"$gc$r$by$n$a'$n\", description=\"Slope of vertical dispersion\" &end"},
{"yAperture", "&column name=yAperture, type=double, units=m, symbol=\"a$by,eff$n\", description=\"Effective vertical aperture\" &end"},
{"pCentral0", "&column name=pCentral0, type=double, units=\"m$be$nc\", symbol=\"p$bcent$n\", description=\"Initial central momentum\" &end"},
{"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
{"ElementOccurence", "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
{"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
{"dI1", "&column name=dI1, type=double, description=\"Contribution to radiation integral 1\", units=m &end"} ,
{"dI2", "&column name=dI2, type=double, description=\"Contribution to radiation integral 2\", units=1/m &end"} ,
{"dI3", "&column name=dI3, type=double, description=\"Contribution to radiation integral 3\", units=1/m$a2$n &end"} ,
{"dI4", "&column name=dI4, type=double, description=\"Contribution to radiation integral 4\", units=1/m &end"} ,
{"dI5", "&column name=dI5, type=double, description=\"Contribution to radiation integral 5\", units=1/m &end"} ,
};

#define IP_STEP 0
#define IP_NUX 1
#define IP_DNUXDP 2
#define IP_DNUXDP2 3
#define IP_DNUXDP3 4
#define IP_AX 5
#define IP_NUY 6
#define IP_DNUYDP 7
#define IP_DNUYDP2 8
#define IP_DNUYDP3 9
#define IP_AY 10
#define IP_DPHRANGE 11
#define IP_NUXUPPER 12
#define IP_NUXLOWER 13
#define IP_NUYUPPER 14
#define IP_NUYLOWER 15
#define IP_STAGE 16
#define IP_PCENTRAL 17
#define IP_DBETAXDP 18
#define IP_DBETAYDP 19
#define IP_ETAX2    20
#define IP_ETAY2    21
#define IP_ETAX3    22
#define IP_ETAY3    23
#define IP_BETAXMIN 24
#define IP_BETAXAVE 25
#define IP_BETAXMAX 26
#define IP_BETAYMIN 27
#define IP_BETAYAVE 28
#define IP_BETAYMAX 29
#define IP_ETAXMAX 30
#define IP_ETAYMAX 31
#define IP_WAISTSX 32
#define IP_WAISTSY 33
#define IP_DNUXDAX 34
#define IP_DNUXDAY 35
#define IP_DNUYDAX 36
#define IP_DNUYDAY 37
#define IP_DNUXDAX2 38
#define IP_DNUXDAY2 39
#define IP_DNUXDAXAY 40
#define IP_DNUYDAX2 41
#define IP_DNUYDAY2 42
#define IP_DNUYDAXAY 43
#define IP_NUXTSWAMIN 44
#define IP_NUXTSWAMAX 45
#define IP_NUYTSWAMIN 46
#define IP_NUYTSWAMAX 47
#define IP_COUPLINGINTEGRAL 48
#define IP_COUPLINGOFFSET 49
#define IP_EMITRATIO 50
#define IP_H11001 51
#define IP_H00111 52
#define IP_H20001 53
#define IP_H00201 54
#define IP_H10002 55
#define IP_H21000 56
#define IP_H30000 57
#define IP_H10110 58
#define IP_H10020 59
#define IP_H10200 60
#define IP_H22000 61
#define IP_H11110 62
#define IP_H00220 63
#define IP_H31000 64
#define IP_H40000 65
#define IP_H20110 66
#define IP_H11200 67
#define IP_H20020 68
#define IP_H20200 69
#define IP_H00310 70
#define IP_H00400 71 
#define IP_DNUXDJX 72
#define IP_DNUXDJY 73
#define IP_DNUYDJY 74
#define IP_ALPHAC2 75
#define IP_ALPHAC  76
/* IP_ALPHAC must be the last item before the radiation-integral-related
 * items!
 */
#define IP_I1 IP_ALPHAC+1
#define IP_I2 IP_ALPHAC+2
#define IP_I3 IP_ALPHAC+3
#define IP_I4 IP_ALPHAC+4
#define IP_I5 IP_ALPHAC+5
#define IP_EX0 IP_ALPHAC+6
#define IP_ENX0 IP_ALPHAC+7
#define IP_TAUX IP_ALPHAC+8
#define IP_JX IP_ALPHAC+9
#define IP_TAUY IP_ALPHAC+10
#define IP_JY IP_ALPHAC+11
#define IP_SIGMADELTA IP_ALPHAC+12
#define IP_TAUDELTA IP_ALPHAC+13
#define IP_JDELTA IP_ALPHAC+14
#define IP_U0 IP_ALPHAC+15
#define N_PARAMETERS IP_U0+1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
{"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
{"nux", "&parameter name=nux, symbol=\"$gn$r$bx$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal tune\" &end"},
{"dnux/dp", "&parameter name=dnux/dp, symbol=\"$gx$r$bx$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal chromaticity\" &end"},
{"dnux/dp2", "&parameter name=dnux/dp2, symbol=\"$gx$r$bx2$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal 2nd-order chromaticity\" &end"},
{"dnux/dp3", "&parameter name=dnux/dp3, symbol=\"$gx$r$bx3$n\", type=double, units=\"1/(2$gp$r)\", description=\"Horizontal 3rd-order chromaticity\" &end"},
{"Ax", "&parameter name=Ax, symbol=\"A$bx$n\", type=double, units=\"m\", description=\"Horizontal acceptance\" &end"},
{"nuy", "&parameter name=nuy, symbol=\"$gn$r$by$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical tune\" &end"},
{"dnuy/dp", "&parameter name=dnuy/dp, symbol=\"$gx$r$by$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical chromaticity\" &end"},
{"dnuy/dp2", "&parameter name=dnuy/dp2, symbol=\"$gx$r$by2$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical 2nd-order chromaticity\" &end"},
{"dnuy/dp3", "&parameter name=dnuy/dp3, symbol=\"$gx$r$by3$n\", type=double, units=\"1/(2$gp$r)\", description=\"Vertical 3rd-order chromaticity\" &end"},
{"Ay", "&parameter name=Ay, symbol=\"A$by$n\", type=double, units=\"m\", description=\"Vertical acceptance\" &end"},
{"deltaHalfRange", "&parameter name=deltaHalfRange, symbol=\"$gDd$r/2\", type=double, description=\"Half range of momentum offset for chromatic tune spread evaluation\" &end"},
{"nuxChromUpper", "&parameter name=nuxChromUpper, symbol=\"$gx$r$bu$n\", type=double, description=\"Upper limit of x tune due to chromaticity and deltaRange\" &end"},
{"nuxChromLower", "&parameter name=nuxChromLower, symbol=\"$gx$r$bu$n\", type=double, description=\"Lower limit of x tune due to chromaticity and deltaRange\" &end"},
{"nuyChromUpper", "&parameter name=nuyChromUpper, symbol=\"$gy$r$bu$n\", type=double, description=\"Upper limit of y tune due to chromaticity and deltaRange\" &end"},
{"nuyChromLower", "&parameter name=nuyChromLower, symbol=\"$gy$r$bu$n\", type=double, description=\"Lower limit of y tune due to chromaticity and deltaRange\" &end"},
{"Stage", "&parameter name=Stage, type=string, description=\"Stage of computation\" &end"},
{"pCentral", "&parameter name=pCentral, type=double, units=\"m$be$nc\", description=\"Central momentum\" &end"},
{"dbetax/dp", "&parameter name=dbetax/dp, units=m, type=double, description=\"Derivative of betax with momentum offset\" &end"},
{"dbetay/dp", "&parameter name=dbetay/dp, units=m, type=double, description=\"Derivative of betay with momentum offset\" &end"},
{"etax2", "&parameter name=etax2, symbol=\"$gc$r$bx2$n\", units=m, type=double, description=\"Second-order dispersion (for matched or periodic case only)\" &end"},
{"etay2", "&parameter name=etay2, symbol=\"$gc$r$by2$n\", units=m, type=double, description=\"Second-order dispersion (for matched or periodic case only)\" &end"},
{"etax3", "&parameter name=etax3, symbol=\"$gc$r$bx3$n\", units=m, type=double, description=\"Third-order dispersion (for matched or periodic case only)\" &end"},
{"etay3", "&parameter name=etay3, symbol=\"$gc$r$by3$n\", units=m, type=double, description=\"Third-order dispersion (for matched or periodic case only)\" &end"},
{"betaxMin", "&parameter name=betaxMin, type=double, units=m, description=\"Minimum betax\" &end"},
{"betaxAve", "&parameter name=betaxAve, type=double, units=m, description=\"Average betax\" &end"},
{"betaxMax", "&parameter name=betaxMax, type=double, units=m, description=\"Maximum betax\" &end"},
{"betayMin", "&parameter name=betayMin, type=double, units=m, description=\"Minimum betay\" &end"},
{"betayAve", "&parameter name=betayAve, type=double, units=m, description=\"Average betay\" &end"},
{"betayMax", "&parameter name=betayMax, type=double, units=m, description=\"Maximum betay\" &end"},
{"etaxMax", "&parameter name=etaxMax, type=double, units=m, description=\"Maximum absolute value of etax\" &end"},
{"etayMax", "&parameter name=etayMax, type=double, units=m, description=\"Maximum absolute value of etay\" &end"},
{"waistsx", "&parameter name=waistsx, type=long, description=\"Number of changes in the sign of alphax\" &end"},
{"waistsy", "&parameter name=waistsy, type=long, description=\"Number of changes in the sign of alphay\" &end"},
{"dnux/dAx", "&parameter name=dnux/dAx, type=double, description=\"Horizontal tune shift with horizontal amplitude\", units=1/m &end"},
{"dnux/dAy", "&parameter name=dnux/dAy, type=double, description=\"Horizontal tune shift with vertical amplitude\", units=1/m &end"},
{"dnuy/dAx", "&parameter name=dnuy/dAx, type=double, description=\"Vertical tune shift with horizontal amplitude\", units=1/m &end"},
{"dnuy/dAy", "&parameter name=dnuy/dAy, type=double, description=\"Vertical tune shift with vertical amplitude\", units=1/m &end"},
{"dnux/dAx2", "&parameter name=dnux/dAx2, type=double, description=\"Horizontal tune shift with horizontal amplitude\", units=1/m$a2$n &end"},
{"dnux/dAy2", "&parameter name=dnux/dAy2, type=double, description=\"Horizontal tune shift with vertical amplitude\", units=1/m$a2$n &end"},
{"dnux/dAxAy", "&parameter name=dnux/dAxAy, type=double, description=\"Horizontal tune shift with horizontal and vertical amplitude\", units=1/m$a2$n &end"},
{"dnuy/dAx2", "&parameter name=dnuy/dAx2, type=double, description=\"Vertical tune shift with horizontal amplitude\", units=1/m$a2$n &end"},
{"dnuy/dAy2", "&parameter name=dnuy/dAy2, type=double, description=\"Vertical tune shift with vertical amplitude\", units=1/m$a2$n &end"},
{"dnuy/dAxAy", "&parameter name=dnuy/dAxAy, type=double, description=\"Vertical tune shift with horizontal and vertical amplitude\", units=1/m$a2$n &end"},
{"nuxTswaLower", "&parameter name=nuxTswaLower, type=double, description=\"Minimum horizontal tune from tune-shift-with-amplitude calculations\", &end"},
{"nuxTswaUpper", "&parameter name=nuxTswaUpper, type=double, description=\"Maximum horizontal tune from tune-shift-with-amplitude calculations\", &end"},
{"nuyTswaLower", "&parameter name=nuyTswaLower, type=double, description=\"Minimum vertical tune from tune-shift-with-amplitude calculations\", &end"},
{"nuyTswaUpper", "&parameter name=nuyTswaUpper, type=double, description=\"Maximum vertical tune from tune-shift-with-amplitude calculations\", &end"},
{"couplingIntegral", "&parameter name=couplingIntegral, type=double, description=\"Coupling integral for difference resonance\" &end"},
{"couplingDelta", "&parameter name=couplingDelta, type=double, description=\"Distance from difference resonance\" &end"},
{"emittanceRatio", "&parameter name=emittanceRatio, type=double, description=\"Emittance ratio from coupling integral\" &end"},
{"h11001", "&parameter name=h11001, type=double, description=\"Magnitude of chromatic driving term (x chromaticity)\", &end"},
{"h00111", "&parameter name=h00111, type=double, description=\"Magnitude of chromatic driving term (y chromaticity)\", &end"},
{"h20001", "&parameter name=h20001, type=double, description=\"Magnitude of chromatic driving term (synchro-betatron resonances)\", &end"},
{"h00201", "&parameter name=h00201, type=double, description=\"Magnitude of chromatic driving term (momentum-dependence of beta functions)\", &end"},
{"h10002", "&parameter name=h10002, type=double, description=\"Magnitude of chromatic driving term (second order dispersion)\", units=\"1/m$a1/2$n\" &end"},
{"h21000", "&parameter name=h21000, type=double, description=\"Magnitude of geometric driving term (nux)\", units=\"1/m$a1/2$n\" &end"},
{"h30000", "&parameter name=h30000, type=double, description=\"Magnitude of geometric driving term (3 nux)\", units=\"1/m$a1/2$n\" &end"},
{"h10110", "&parameter name=h10110, type=double, description=\"Magnitude of geometric driving term (nux)\", units=\"1/m$a1/2$n\" &end"},
{"h10020", "&parameter name=h10020, type=double, description=\"Magnitude of geometric driving term (nux - 2 nuy)\", units=\"1/m$a1/2$n\" &end"},
{"h10200", "&parameter name=h10200, type=double, description=\"Magnitude of geometric driving term (nux + 2 nuy)\", units=\"1/m\" &end"},
{"h22000", "&parameter name=h22000, type=double, description=\"Magnitude of geometric driving term (amplitude-dependent tune)\", units=\"1/m\" &end"},
{"h11110", "&parameter name=h11110, type=double, description=\"Magnitude of geometric driving term (amplitude-dependent tune)\", units=\"1/m\" &end"},
{"h00220", "&parameter name=h00220, type=double, description=\"Magnitude of geometric driving term (amplitude-dependent tune)\", units=\"1/m\" &end"},
{"h31000", "&parameter name=h31000, type=double, description=\"Magnitude of geometric driving term (2 nux)\", units=\"1/m\" &end"},
{"h40000", "&parameter name=h40000, type=double, description=\"Magnitude of geometric driving term (4 nux)\", units=\"1/m\" &end"},
{"h20110", "&parameter name=h20110, type=double, description=\"Magnitude of geometric driving term (2 nux)\", units=\"1/m\" &end"},
{"h11200", "&parameter name=h11200, type=double, description=\"Magnitude of geometric driving term (2 nuy)\", units=\"1/m\" &end"},
{"h20020", "&parameter name=h20020, type=double, description=\"Magnitude of geometric driving term (2 nux - 2 nuy)\", units=\"1/m\" &end"},
{"h20200", "&parameter name=h20200, type=double, description=\"Magnitude of geometric driving term (2 nux + 2 nuy)\", units=\"1/m\" &end"},
{"h00310", "&parameter name=h00310, type=double, description=\"Magnitude of geometric driving term (2 nuy)\", units=\"1/m\" &end"},
{"h00400", "&parameter name=h00400, type=double, description=\"Magnitude of geometric driving term (4 nuy)\", units=\"1/m\" &end"},
{"dnux/dJx", "&parameter name=dnux/dJx, type=double, description=\"Horizontal tune shift with horizontal invariant\", units=\"1/m\" &end"},
{"dnux/dJy", "&parameter name=dnux/dJy, type=double, description=\"Horizontal tune shift with vertical invariant\", units=\"1/m\" &end"},
{"dnuy/dJy", "&parameter name=dnuy/dJy, type=double, description=\"Vertical tune shift with vertical invariant\", units=\"1/m\" &end"},
{"alphac2", "&parameter name=alphac2, symbol=\"$ga$r$bc2$n\", type=double, description=\"2nd-order momentum compaction factor\" &end"},
{"alphac", "&parameter name=alphac, symbol=\"$ga$r$bc$n\", type=double, description=\"Momentum compaction factor\" &end"},
{"I1", "&parameter name=I1, type=double, description=\"Radiation integral 1\", units=m &end"} ,
{"I2", "&parameter name=I2, type=double, description=\"Radiation integral 2\", units=1/m &end"} ,
{"I3", "&parameter name=I3, type=double, description=\"Radiation integral 3\", units=1/m$a2$n &end"} ,
{"I4", "&parameter name=I4, type=double, description=\"Radiation integral 4\", units=1/m &end"} ,
{"I5", "&parameter name=I5, type=double, description=\"Radiation integral 5\", units=1/m &end"} ,
{"ex0", "&parameter name=ex0, type=double, description=\"Damped horizontal emittance\", units=$gp$rm &end"},
{"enx0", "&parameter name=enx0, type=double, units=\"m$be$nc $gp$rm\", description=\"Damped normalized horizontal emittance\""},
{"taux", "&parameter name=taux, type=double, description=\"Horizontal damping time\", units=s &end"},
{"Jx", "&parameter name=Jx, type=double, description=\"Horizontal damping partition number\" &end"},
{"tauy", "&parameter name=tauy, type=double, description=\"Vertical damping time\", units=s &end"},
{"Jy", "&parameter name=Jy, type=double, description=\"Vertical damping partition number\" &end"},
{"Sdelta0", "&parameter name=Sdelta0, type=double, description=\"RMS fractional energy spread\" &end"},
{"taudelta", "&parameter name=taudelta, type=double, description=\"Longitudinal damping time\", units=s &end"},
{"Jdelta", "&parameter name=Jdelta, type=double, description=\"Longitudinal damping partition number\" &end"},
{"U0", "&parameter name=U0, type=double, units=MeV, description=\"Energy loss per turn\" &end"},
} ;

void dump_twiss_parameters(
  LINE_LIST *beamline,
  long n_elem,
  double *tune,
  RADIATION_INTEGRALS *radIntegrals,                           
  double *chromaticity,
  double *dbeta,
  double *acceptance,
  double *alphac,
  long final_values_only,
  long tune_corrected,
  RUN *run
  )
{
  double *data;
  long i, j, row_count;
  char *stage;
  TWISS twiss_ave, twiss_min, twiss_max;

  TWISS *twiss0;
  ELEMENT_LIST *elem;

  log_entry("dump_twiss_parameters");

  twiss0 = beamline->twiss0;
  elem = beamline->elem_twiss;
  
  if (!twiss0)
    bomb("Twiss data not computed prior to dump_twiss_parameters() call (1)", NULL);

  data = tmalloc(sizeof(*data)*N_DOUBLE_COLUMNS);

  if (tune_corrected==1)
    stage = "tunes corrected";
  else
    stage = "tunes uncorrected";

  if (!SDDS_StartTable(&SDDS_twiss, final_values_only?1:n_elem+1)) {
    SDDS_SetError("Problem starting SDDS table (dump_twiss_parameters)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
  if (!SDDS_SetParameters(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                          IP_STEP, twiss_count, IP_STAGE, stage, 
                          IP_NUX, tune[0], IP_DNUXDP, chromaticity[0], IP_AX, acceptance[0],
                          IP_NUY, tune[1], IP_DNUYDP, chromaticity[1], IP_AY, acceptance[1], 
                          IP_DNUXDP2, beamline->chrom2[0],
                          IP_DNUYDP2, beamline->chrom2[1],
                          IP_DNUXDP3, beamline->chrom3[0],
                          IP_DNUYDP3, beamline->chrom3[1],
                          IP_DPHRANGE, beamline->chromDeltaHalfRange,
                          IP_NUXUPPER, beamline->tuneChromUpper[0],
                          IP_NUYUPPER, beamline->tuneChromUpper[1],
                          IP_NUXLOWER, beamline->tuneChromLower[0],
                          IP_NUYLOWER, beamline->tuneChromLower[1],
                          IP_COUPLINGINTEGRAL, beamline->couplingFactor[0],
			  IP_COUPLINGOFFSET, beamline->couplingFactor[1],
	                  IP_EMITRATIO, beamline->couplingFactor[2],
                          IP_ALPHAC, alphac[0], IP_ALPHAC2, alphac[1], 
                          IP_DBETAXDP, dbeta[0], IP_DBETAYDP, dbeta[1],
                          IP_BETAXMIN, twiss_min.betax, IP_BETAXAVE, twiss_ave.betax, IP_BETAXMAX, twiss_max.betax, 
                          IP_BETAYMIN, twiss_min.betay, IP_BETAYAVE, twiss_ave.betay, IP_BETAYMAX, twiss_max.betay, 
                          IP_ETAXMAX, MAX(fabs(twiss_min.etax), fabs(twiss_max.etax)),
                          IP_ETAYMAX, MAX(fabs(twiss_min.etay), fabs(twiss_max.etay)),
                          IP_WAISTSX, beamline->waists[0],
                          IP_WAISTSY, beamline->waists[1],
                          IP_DNUXDAX, beamline->dnux_dA[1][0],
                          IP_DNUXDAY, beamline->dnux_dA[0][1],
                          IP_DNUYDAX, beamline->dnuy_dA[1][0],
                          IP_DNUYDAY, beamline->dnuy_dA[0][1],
                          IP_DNUXDAX2, beamline->dnux_dA[2][0],
                          IP_DNUXDAY2, beamline->dnux_dA[0][2],
                          IP_DNUXDAXAY, beamline->dnux_dA[1][1],
                          IP_DNUYDAX2, beamline->dnuy_dA[2][0],
                          IP_DNUYDAY2, beamline->dnuy_dA[0][2],
                          IP_DNUYDAXAY, beamline->dnuy_dA[1][1],
			  IP_NUXTSWAMIN, beamline->nuxTswaExtrema[0],
			  IP_NUXTSWAMAX, beamline->nuxTswaExtrema[1],
			  IP_NUYTSWAMIN, beamline->nuyTswaExtrema[0],
			  IP_NUYTSWAMAX, beamline->nuyTswaExtrema[1],
                          IP_H11001, beamline->drivingTerms.h11001,
                          IP_H00111, beamline->drivingTerms.h00111,
                          IP_H20001, beamline->drivingTerms.h20001,
                          IP_H00201, beamline->drivingTerms.h00201,
                          IP_H10002, beamline->drivingTerms.h10002,
                          IP_H21000, beamline->drivingTerms.h21000,
                          IP_H30000, beamline->drivingTerms.h30000,
                          IP_H10110, beamline->drivingTerms.h10110,
                          IP_H10020, beamline->drivingTerms.h10020,
                          IP_H10200, beamline->drivingTerms.h10200,
                          IP_H22000, beamline->drivingTerms.h22000,
                          IP_H11110, beamline->drivingTerms.h11110,
                          IP_H00220, beamline->drivingTerms.h00220,
                          IP_H31000, beamline->drivingTerms.h31000,
                          IP_H40000, beamline->drivingTerms.h40000,
                          IP_H20110, beamline->drivingTerms.h20110,
                          IP_H11200, beamline->drivingTerms.h11200,
                          IP_H20020, beamline->drivingTerms.h20020,
                          IP_H20200, beamline->drivingTerms.h20200,
                          IP_H00310, beamline->drivingTerms.h00310,
                          IP_H00400, beamline->drivingTerms.h00400,
                          IP_DNUXDJX, beamline->drivingTerms.dnux_dJx,
                          IP_DNUXDJY, beamline->drivingTerms.dnux_dJy,
                          IP_DNUYDJY, beamline->drivingTerms.dnuy_dJy,
                          IP_ETAX2, beamline->eta2[0],
                          IP_ETAY2, beamline->eta2[2],
                          IP_ETAX3, beamline->eta3[0],
                          IP_ETAY3, beamline->eta3[2],
                          IP_PCENTRAL, run->p_central, -1)) {
    SDDS_SetError("Problem setting SDDS parameters (dump_twiss_parameters 1)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (radIntegrals) {
    if (!SDDS_SetParameters(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                            IP_I1, radIntegrals->RI[0],
                            IP_I2, radIntegrals->RI[1],
                            IP_I3, radIntegrals->RI[2],
                            IP_I4, radIntegrals->RI[3],
                            IP_I5, radIntegrals->RI[4],
                            IP_EX0, radIntegrals->ex0,
                            IP_TAUX, radIntegrals->taux,
                            IP_JX, radIntegrals->Jx,
                            IP_TAUY, radIntegrals->tauy,
                            IP_JY, radIntegrals->Jy,
                            IP_SIGMADELTA, radIntegrals->sigmadelta,
                            IP_TAUDELTA, radIntegrals->taudelta,
                            IP_JDELTA, radIntegrals->Jdelta,
                            IP_U0, radIntegrals->Uo, 
                            IP_ENX0, radIntegrals->ex0*sqrt(sqr(run->p_central)+1), 
                            -1)) {
      SDDS_SetError("Problem setting SDDS parameters (dump_twiss_parameters 2)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }

  if (!final_values_only) {
    row_count = 0;
    data[0] = 0;     /* position */
    copy_doubles(data+1, (double*)twiss0, N_DOUBLE_COLUMNS-2);
    data[N_DOUBLE_COLUMNS-1] = elem->Pref_input;
    for (j=0; j<N_DOUBLE_COLUMNS; j++)
      if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count++, 
                           IC_ELEMENT, "_BEG_", IC_OCCURENCE, (long)1, IC_TYPE, "MARK", -1)) {
      SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    i = 0;
    while (elem) {
      data[0] = elem->end_pos;     /* position */
      data[N_DOUBLE_COLUMNS-1] = elem->Pref_output;
      if (!elem->twiss)
        bomb("Twiss data not computed prior to dump_twiss_parameters() call (2)", NULL);
      copy_doubles(data+1, (double*)elem->twiss, N_DOUBLE_COLUMNS-2);
      for (j=0; j<N_DOUBLE_COLUMNS; j++)
        if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
          SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, 
                             IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                             IC_TYPE, entity_name[elem->type], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (radIntegrals) {
        if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count,
                               IC_I1, elem->twiss->dI[0],
                               IC_I2, elem->twiss->dI[1],
                               IC_I3, elem->twiss->dI[2],
                               IC_I4, elem->twiss->dI[3], 
                               IC_I5, elem->twiss->dI[4],
                               -1)) {
          SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
      i++;      
      row_count++;
      elem = elem->succ;
    }
    if (i!=n_elem)
      bomb("element count error in dump_twiss_parameters()", NULL);
  }
  else {
    /* find final element */
    i = 0;
    while (1) {
      if (!elem->twiss)
        bomb("Twiss data not computed prior to dump_twiss_parameters() call (2)", NULL);
      i++;
      if (!elem->succ)
        break;
      elem = elem->succ;
    }
    if (i!=n_elem)
      bomb("element count error in dump_twiss_parameters()", NULL);
    data[0] = elem->end_pos;     /* position */
    data[N_DOUBLE_COLUMNS-1] = elem->Pref_output;
    copy_doubles(data+1, (double*)elem->twiss, N_DOUBLE_COLUMNS-2);
    for (j=0; j<N_DOUBLE_COLUMNS; j++)
      if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, j, data[j], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, 
                           IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                           IC_TYPE, entity_name[elem->type], -1)) {
      SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (radIntegrals) {
      if (!SDDS_SetRowValues(&SDDS_twiss, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0,
                             IC_I1, elem->twiss->dI[0],
                             IC_I2, elem->twiss->dI[1],
                             IC_I3, elem->twiss->dI[2],
                             IC_I4, elem->twiss->dI[3], 
                             IC_I5, elem->twiss->dI[4],
                             -1)) {
        SDDS_SetError("Problem setting SDDS rows (dump_twiss_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }

  if (!SDDS_WriteTable(&SDDS_twiss)) {
    SDDS_SetError("Unable to write Twiss parameter data (dump_twiss_parameters)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDS_twiss);
  if (!SDDS_EraseData(&SDDS_twiss)) {
    SDDS_SetError("Unable to erase Twiss parameter data (dump_twiss_parameters)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  free(data);
  if (tune_corrected)
    twiss_count++;

  log_exit("dump_twiss_parameters");
}

long get_twiss_mode(long *mode, double *x_twiss, double *y_twiss)
{
  if (!twiss_initialized)
    return(0);
  *mode = matched;
  x_twiss[0] = beta_x;
  x_twiss[1] = alpha_x;
  x_twiss[2] = eta_x;
  x_twiss[3] = etap_x;
  x_twiss[4] = 0;
  y_twiss[0] = beta_y;
  y_twiss[1] = alpha_y;
  y_twiss[2] = eta_y;
  y_twiss[3] = etap_y;
  y_twiss[4] = 0;
  return(1);
}

void setup_twiss_output(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, long *do_twiss_output,
                        long default_order)
{

  log_entry("setup_twiss_output");

  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&twiss_output, nltext);
  if (echoNamelists) print_namelist(stdout, &twiss_output);
  
#if USE_MPI
    if (!writePermitted)
      filename = NULL;
#endif 

  if (filename)
    filename = compose_filename(filename, run->rootname);
  twissConcatOrder = concat_order;
  if (twissConcatOrder<default_order)
    twissConcatOrder = default_order;
  *do_twiss_output = output_at_each_step;
  if (higher_order_chromaticity && !matched) 
    bomb("higher order chromaticity calculations available only for a matched (periodic) beamline", NULL);

  if (reference_file && matched)
    bomb("reference_file and matched=1 are incompatible", NULL);
  if (!matched) {
    if (reference_file) {
      if (reference_element && reference_element_occurrence<0)
        bomb("invalid value of reference_element_occurrence---use 0 for last occurrence, >=1 for specific occurrence.", NULL);
      LoadStartingTwissFromFile(&beta_x, &beta_y, &alpha_x, &alpha_y, 
                                &eta_x, &etap_x, &eta_y, &etap_y,
                                reference_file, reference_element,
                                reference_element_occurrence);
      if (reflect_reference_values) {
        alpha_x = -alpha_x;
        alpha_y = -alpha_y;
        etap_x = -etap_x;
        etap_y = -etap_y;
      }
      fprintf(stdout, "Starting twiss parameters from reference file:\nbeta, alpha x: %le, %le\nbeta, alpha y: %le, %le\n",
              beta_x, alpha_x, beta_y, alpha_y);
      fprintf(stdout, "eta, eta' x: %le, %le\neta, eta' y: %le, %le\n",
              eta_x, etap_x, eta_y, etap_y);
      fflush(stdout);
    }
    if (beta_x<=0 || beta_y<=0)
      bomb("invalid initial beta-functions given in twiss_output namelist", NULL);
  }

  if (filename) {
#if SDDS_MPI_IO
    SDDS_twiss.parallel_io = 0;
#endif
    SDDS_ElegantOutputSetup(&SDDS_twiss, filename, SDDS_BINARY, 1, "Twiss parameters",
                            run->runfile, run->lattice, parameter_definition, 
                            (radiation_integrals?N_PARAMETERS:IP_ALPHAC+1),
                            column_definition, (radiation_integrals?N_COLUMNS_WRI:N_COLUMNS), "setup_twiss_output",
                            SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    SDDS_twiss_initialized = 1;
    twiss_count = 0;
  }
  else
    SDDS_twiss_initialized = 0;
  twiss_initialized = 1;

  beamline->flags |= BEAMLINE_TWISS_WANTED;
  if (radiation_integrals)
    beamline->flags |= BEAMLINE_RADINT_WANTED;
  
  beamline->chromDeltaHalfRange = chromatic_tune_spread_half_range;

  log_exit("setup_twiss_output");
}

void finish_twiss_output(void)
{
  log_entry("finish_twiss_output");
  if (SDDS_twiss_initialized && !SDDS_Terminate(&SDDS_twiss)) {
    SDDS_SetError("Problem terminating SDDS output (finish_twiss_output)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  SDDS_twiss_initialized = twiss_count = 0;
  if (doTuneShiftWithAmplitude && tune_shift_with_amplitude_struct.tune_output)
    SDDS_Terminate(&SDDSTswaTunes);
  doTuneShiftWithAmplitude = 0;
  log_exit("finish_twiss_output");
}

long run_twiss_output(RUN *run, LINE_LIST *beamline, double *starting_coord, long tune_corrected)
{
  ELEMENT_LIST *eptr, *elast;
  long n_elem, last_n_elem;
  unsigned long unstable;
  TWISS twiss_ave, twiss_min, twiss_max;

  /*
    if (beamline->flags&BEAMLINE_TWISS_CURRENT)
    return;
    */

  log_entry("run_twiss_output");

#ifdef DEBUG
  fprintf(stdout, "now in run_twiss_output\n");
  fflush(stdout);
#endif

  if (tune_corrected==0 && !output_before_tune_correction) {
    log_exit("run_twiss_output");
    return 1;
  }

  eptr = beamline->elem_twiss = &(beamline->elem);
  n_elem = last_n_elem = beamline->n_elems;
  while (eptr) {
    if (eptr->type==T_RECIRC && matched) {
      last_n_elem = n_elem;
      beamline->elem_twiss = beamline->elem_recirc = eptr;
    }
    eptr = eptr->succ;
    n_elem --;
  }
  n_elem = last_n_elem;
  
  compute_twiss_parameters(run, beamline, starting_coord, matched, radiation_integrals,
                           beta_x, alpha_x, eta_x, etap_x,
                           beta_y, alpha_y, eta_y, etap_y, &unstable);
  elast = beamline->elast;

  if (twissConcatOrder) {
#ifdef DEBUG
    fprintf(stdout, "computing chromaticities\n");
    fflush(stdout);
#endif
    fprintf(stdout, "%s Twiss parameters (chromaticity valid for fully second-order calculation only!):\n",
            matched?"periodic":"final");
    fflush(stdout);
    fprintf(stdout, "         beta          alpha           nu           eta          eta'       dnu/d(dp/p)   dbeta/(dp/p)     accept.\n");
    fflush(stdout);
    fprintf(stdout, "          m                          1/2pi           m                         1/2pi            m          mm-mrad\n");
    fflush(stdout);
    fprintf(stdout, "--------------------------------------------------------------------------------------------------------------------\n");
    fflush(stdout);
    fprintf(stdout, "  x: %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",  
            elast->twiss->betax, elast->twiss->alphax, beamline->tune[0], elast->twiss->etax, elast->twiss->etapx,
            beamline->chromaticity[0], beamline->dbeta_dPoP[0], 1e6*beamline->acceptance[0]);
    fflush(stdout);
    if (statistics) {
      compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
      fprintf(stdout, "ave: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_ave.betax, twiss_ave.alphax, "", twiss_ave.etax, twiss_ave.etapx);
      fflush(stdout);
      fprintf(stdout, "min: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_min.betax, twiss_min.alphax, "", twiss_min.etax, twiss_min.etapx);
      fflush(stdout);
      fprintf(stdout, "max: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_max.betax, twiss_max.alphax, "", twiss_max.etax, twiss_max.etapx);
      fflush(stdout);
    }
    fprintf(stdout, "  y: %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",  
            elast->twiss->betay, elast->twiss->alphay, beamline->tune[1], elast->twiss->etay, elast->twiss->etapy,
            beamline->chromaticity[1], beamline->dbeta_dPoP[1], 1e6*beamline->acceptance[1]);
    fflush(stdout);
    if (statistics) {
      fprintf(stdout, "ave: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_ave.betay, twiss_ave.alphay, "", twiss_ave.etay, twiss_ave.etapy);
      fflush(stdout);
      fprintf(stdout, "min: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_min.betay, twiss_min.alphay, "", twiss_min.etay, twiss_min.etapy);
      fflush(stdout);
      fprintf(stdout, "max: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_max.betay, twiss_max.alphay, "", twiss_max.etay, twiss_max.etapy);
      fflush(stdout);
    }
  }
  else {
    fprintf(stdout, "%s Twiss parameters:\n", matched?"periodic":"final");
    fflush(stdout);
    fprintf(stdout, "         beta          alpha           nu           eta          eta'        accept.\n");
    fflush(stdout);
    fprintf(stdout, "          m                          1/2pi           m                       mm-mrad\n");
    fflush(stdout);
    fprintf(stdout, "---------------------------------------------------------------------------------------\n");
    fflush(stdout);
    fprintf(stdout, "  x: %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",  
            elast->twiss->betax, elast->twiss->alphax, beamline->tune[0], elast->twiss->etax, elast->twiss->etapx,
            1e6*beamline->acceptance[0]);
    fflush(stdout);
    if (statistics) {
      compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
      fprintf(stdout, "ave: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_ave.betax, twiss_ave.alphax, "", twiss_ave.etax, twiss_ave.etapx);
      fflush(stdout);
      fprintf(stdout, "min: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_min.betax, twiss_min.alphax, "", twiss_min.etax, twiss_min.etapx);
      fflush(stdout);
      fprintf(stdout, "max: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_max.betax, twiss_max.alphax, "", twiss_max.etax, twiss_max.etapx);
      fflush(stdout);
    }
    fprintf(stdout, "  y: %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",  
            elast->twiss->betay, elast->twiss->alphay, beamline->tune[1], elast->twiss->etay, elast->twiss->etapy,
            1e6*beamline->acceptance[1]);
    fflush(stdout);
    if (statistics) {
      fprintf(stdout, "ave: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_ave.betay, twiss_ave.alphay, "", twiss_ave.etay, twiss_ave.etapy);
      fflush(stdout);
      fprintf(stdout, "min: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_min.betay, twiss_min.alphay, "", twiss_min.etay, twiss_min.etapy);
      fflush(stdout);
      fprintf(stdout, "max: %13.6e %13.6e %-13s %13.6e %13.6e\n",
              twiss_max.betay, twiss_max.alphay, "", twiss_max.etay, twiss_max.etapy);
      fflush(stdout);
    }
  }

  if (beamline->acc_limit_name[0])
    fprintf(stdout, "x acceptance limited by %s ending at %e m\n", beamline->acc_limit_name[0], beamline->acceptance[2]);
    fflush(stdout);
  if (beamline->acc_limit_name[1])
    fprintf(stdout, "y acceptance limited by %s ending at %e m\n", beamline->acc_limit_name[1], beamline->acceptance[3]);
    fflush(stdout);

  if (SDDS_twiss_initialized) {
    dump_twiss_parameters(beamline, n_elem,
                          beamline->tune, 
                          radiation_integrals?&(beamline->radIntegrals):NULL, 
                          beamline->chromaticity, beamline->dbeta_dPoP,
                          beamline->acceptance, 
                          beamline->alpha, final_values_only, tune_corrected, run);
  }

  if (isnan(beamline->tune[0]) || 
      isnan(elast->twiss->betax) ||
      isnan(elast->twiss->etax) ||
      isnan(beamline->tune[1]) || 
      isnan(elast->twiss->betay) ||
      isnan(elast->twiss->etay))
    return 0;

  log_exit("run_twiss_output");
  return 1;
}

void compute_twiss_parameters(RUN *run, LINE_LIST *beamline, double *starting_coord, 
			      long periodic,
                              long radiation_integrals,
                              double betax, double alphax, double etax, double etapx, 
                              double betay, double alphay, double etay, double etapy, 
                              unsigned long *unstable)
{
  VMATRIX *M;
  double chromx, chromy, dbetax, dbetay, alpha1, alpha2, dalphax, dalphay;
  double chromx2, chromy2;
  double x_acc_z, y_acc_z;
  ELEMENT_LIST *eptr, *elast;
  char *x_acc_name, *y_acc_name;
  long i, j;
  
  log_entry("compute_twiss_parameters");

  *unstable = 0;
  
  if (!beamline->twiss0)
    beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));

  eptr = beamline->elem_twiss = &(beamline->elem);

  elast = eptr;
  while (eptr) {
    if (eptr->type==T_RECIRC && matched)
      beamline->elem_twiss = beamline->elem_recirc = eptr;
    elast = eptr;
    eptr = eptr->succ;
  }
  beamline->elast = elast;

  if (periodic) {
    if (beamline->matrix) {
      free_matrices(beamline->matrix);
      free(beamline->matrix);
    }
    beamline->matrix = compute_periodic_twiss(&betax, &alphax, &etax, &etapx, beamline->tune,
                                              &betay, &alphay, &etay, &etapy, beamline->tune+1, 
                                              beamline->elem_twiss, starting_coord, run,
                                              unstable, 
                                              beamline->eta2, beamline->eta3);
#ifdef DEBUG
    fprintf(stdout, "matched parameters computed--returned to compute_twiss_parameters\n");
    fflush(stdout);
    fprintf(stdout, "beamline matrix has order %ld\n", beamline->matrix->order);
    fflush(stdout);
    fprintf(stdout, "beamline matrix pointers:  %x, %x, %x, %x\n",
            beamline->matrix, beamline->matrix->C, beamline->matrix->R, beamline->matrix->T);
    fflush(stdout);
#endif
    if (twissConcatOrder>=2 && !(beamline->matrix->T))
      bomb("logic error: T matrix is NULL on return from compute_periodic_twiss", NULL);
  }
  else {
    VMATRIX *M1;
    if (twissConcatOrder>1 && starting_coord) {
      M1 = tmalloc(sizeof(*M1));
      initialize_matrices(M1, 1);
      M1->C[0] = starting_coord[0];
      M1->C[1] = starting_coord[0];
      M1->C[2] = starting_coord[0];
      M1->C[3] = starting_coord[0];
      M1->C[4] = starting_coord[0];
      M1->C[5] = starting_coord[0];
      for (i=0; i<6; i++)
        M1->R[i][i] = 1;
      fill_in_matrices(beamline->elem_twiss, run);
      beamline->matrix = append_full_matrix(beamline->elem_twiss, run, M1, twissConcatOrder);
    }
    else
      beamline->matrix = full_matrix(beamline->elem_twiss, run, twissConcatOrder);
  }

  beamline->twiss0->betax  = betax;
  beamline->twiss0->alphax = alphax;
  beamline->twiss0->phix   = 0;
  beamline->twiss0->etax   = etax;
  beamline->twiss0->etapx  = etapx;
  beamline->twiss0->betay  = betay;
  beamline->twiss0->alphay = alphay;
  beamline->twiss0->phiy   = 0;
  beamline->twiss0->etay   = etay;
  beamline->twiss0->etapy  = etapy;
  if (starting_coord) {
    beamline->twiss0->Cx  = starting_coord[0];
    beamline->twiss0->Cy  = starting_coord[2];
  } else
    beamline->twiss0->Cx  = beamline->twiss0->Cy = 0;
  
  if (twissConcatOrder>=2 && !(beamline->matrix->T))
    bomb("logic error: beamline T matrix is NULL in compute_twiss_parameters", NULL);

#ifdef DEBUG
  fprintf(stdout, "propagating parameters\n");
  fflush(stdout);
  fprintf(stdout, "beamline matrix pointers:  %x, %x, %x, %x\n",
          beamline->matrix, beamline->matrix->C, beamline->matrix->R, beamline->matrix->T);
  fflush(stdout);
#endif

  beamline->radIntegrals.computed = 0;
  
  propagate_twiss_parameters(beamline->twiss0, beamline->tune, beamline->waists,
                             (radiation_integrals?&(beamline->radIntegrals):NULL),
                             beamline->elem_twiss, run, starting_coord,
			     beamline->couplingFactor);
  
  if (radiation_integrals)
    completeRadiationIntegralComputation(&(beamline->radIntegrals), run->p_central,
                              beamline->revolution_length);
  
#ifdef DEBUG
  fprintf(stdout, "finding acceptance\n");
  fflush(stdout);
#endif
  beamline->acceptance[0] = find_acceptance(beamline->elem_twiss, 0, run, &x_acc_name, &x_acc_z);
  beamline->acceptance[1] = find_acceptance(beamline->elem_twiss, 1, run, &y_acc_name, &y_acc_z);
  beamline->acceptance[2] = x_acc_z;
  beamline->acceptance[3] = y_acc_z;
  if (x_acc_name)
    cp_str(&beamline->acc_limit_name[0], x_acc_name);
  else
    beamline->acc_limit_name[0] = NULL;
  if (y_acc_name)
    cp_str(&beamline->acc_limit_name[1], y_acc_name);
  else
    beamline->acc_limit_name[1] = NULL;

  beamline->twiss0->apx = beamline->elem_twiss->twiss->apx;
  beamline->twiss0->apy = beamline->elem_twiss->twiss->apy;

#ifdef DEBUG
  fprintf(stdout, "computing chromaticities\n");
  fflush(stdout);
#endif

  chromx = chromy = 0;
  dbetax = dbetay = 0;
  chromx2 = chromy2 = 0;

  for (i=0; i<N_TSWA; i++)
    for (j=0; j<N_TSWA; j++)
      beamline->dnux_dA[i][j] = beamline->dnuy_dA[i][j] = 0;

  beamline->drivingTerms.h21000 = beamline->drivingTerms.h30000 = beamline->drivingTerms.h10110 =
    beamline->drivingTerms.h10020 = beamline->drivingTerms.h10200 = 
      beamline->drivingTerms.h22000 = beamline->drivingTerms.h11110 = beamline->drivingTerms.h00220 = beamline->drivingTerms.h31000 = 
        beamline->drivingTerms.h40000 = beamline->drivingTerms.h20110 = beamline->drivingTerms.h11200 = beamline->drivingTerms.h20020 = 
          beamline->drivingTerms.h20200 = beamline->drivingTerms.h00310 = beamline->drivingTerms.h00400 = 0;
  
  if (periodic) {
    if (!(M = beamline->matrix))
      bomb("logic error: revolution matrix is NULL in compute_twiss_parameters", NULL);

    if (twissConcatOrder>1) {
#ifdef DEBUG
      fprintf(stdout, "computing chromaticities\n");
      fflush(stdout);
#endif
      if (!(M->T))
        bomb("logic error: T matrix is NULL in compute_twiss_parameters", NULL);
      computeChromaticities(&chromx, &chromy, 
                            &dbetax, &dbetay, &dalphax, &dalphay, beamline->twiss0, 
			    beamline->elast->twiss, M);
      beamline->chromaticity[0] = chromx;
      beamline->chromaticity[1] = chromy;
      if (twissConcatOrder>1 && higher_order_chromaticity)
	computeHigherOrderChromaticities(beamline, starting_coord, run, twissConcatOrder,
                                         higher_order_chromaticity_range/
                                         (higher_order_chromaticity_points-1),
                                         higher_order_chromaticity_points, quick_higher_order_chromaticity);
      computeChromaticTuneLimits(beamline);
      if (doTuneShiftWithAmplitude)
        computeTuneShiftWithAmplitude(beamline->dnux_dA, beamline->dnuy_dA,
				      beamline->nuxTswaExtrema, beamline->nuyTswaExtrema,
                                      beamline->twiss0, beamline->tune, M, beamline, run,
                                      starting_coord); 
      if (compute_driving_terms)
        computeDrivingTerms(&(beamline->drivingTerms), beamline->elem_twiss, beamline->twiss0, beamline->tune);
    }
  }
  else {
    M = beamline->matrix;
    elast = beamline->elast;
    
    if (!elast)
      bomb("logic error in compute_twiss_parameters--elast pointer is NULL", NULL);
    if (!elast->twiss)
      bomb("logic error in compute_twiss_parameters--elast->twiss pointer is NULL", NULL);

    betax = elast->twiss->betax;
    alphax = elast->twiss->alphax;
    etax = elast->twiss->etax;
    etapx = elast->twiss->etapx;
    betay = elast->twiss->betay;
    alphay = elast->twiss->alphay;
    etay = elast->twiss->etay;
    etapy = elast->twiss->etapy;

    if (twissConcatOrder>=2) {
      if (!(M->T))
        bomb("logic error: T matrix is NULL in compute_twiss_parameters", NULL);
      computeChromaticities(&chromx, &chromy, 
                            &dbetax, &dbetay, &dalphax, &dalphay, beamline->twiss0, 
			    beamline->elast->twiss, M);
      beamline->chromaticity[0] = chromx;
      beamline->chromaticity[1] = chromy;
      if (twissConcatOrder>1 && higher_order_chromaticity)
	computeHigherOrderChromaticities(beamline, starting_coord, run, twissConcatOrder,
                                         higher_order_chromaticity_range/
                                         (higher_order_chromaticity_points-1),
                                         higher_order_chromaticity_points, quick_higher_order_chromaticity);
      computeChromaticTuneLimits(beamline);
      if (doTuneShiftWithAmplitude)
        computeTuneShiftWithAmplitude(beamline->dnux_dA, beamline->dnuy_dA,
				      beamline->nuxTswaExtrema, beamline->nuyTswaExtrema,
                                      beamline->twiss0, beamline->tune, M, beamline, run,
                                      starting_coord); 
#ifdef DEBUG
      fprintf(stdout, "chomaticities: %e, %e\n", chromx, chromy);
      fflush(stdout);
#endif
    }
  }
  beamline->dbeta_dPoP[0] = dbetax;
  beamline->dbeta_dPoP[1] = dbetay;
  beamline->dalpha_dPoP[0] = dalphax;
  beamline->dalpha_dPoP[1] = dalphay;
  
  alpha1 = alpha2 = 0;
  if (beamline->matrix->C[4]!=0) {
    alpha1 = (beamline->matrix->R[4][5] +
              beamline->matrix->R[4][0]*elast->twiss->etax +
              beamline->matrix->R[4][1]*elast->twiss->etapx +
              beamline->matrix->R[4][2]*elast->twiss->etay +
              beamline->matrix->R[4][3]*elast->twiss->etapy)/beamline->matrix->C[4];
    if (beamline->matrix->T) {
      double eta[6];
      long j, k;
      eta[0] = elast->twiss->etax;
      eta[1] = elast->twiss->etapx;
      eta[2] = elast->twiss->etay;
      eta[3] = elast->twiss->etapy;
      eta[4] = 0;
      eta[5] = 1;
      for (j=0; j<4; j++)
        alpha2 += beamline->matrix->R[4][j]*beamline->eta2[j];
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++) 
          alpha2 += beamline->matrix->T[4][j][k]*eta[j]*eta[k];
      alpha2 /= beamline->matrix->C[4];
    }
  }
  beamline->alpha[0] = alpha1;
  beamline->alpha[1] = alpha2;

  beamline->flags |= BEAMLINE_TWISS_DONE+BEAMLINE_TWISS_CURRENT;
  if (radiation_integrals) 
    beamline->flags |= BEAMLINE_RADINT_DONE+BEAMLINE_RADINT_CURRENT;
  
  if (cavities_are_drifts_if_matched && mustResetRfcaMatrices)
    reset_rfca_matrices(&(beamline->elem), run->default_order);
  
  log_exit("compute_twiss_parameters");
}

void update_twiss_parameters(RUN *run, LINE_LIST *beamline, unsigned long *unstable)
{
  unsigned long unstable0;
  if (linearChromaticTrackingInitialized) {
    /* This is a kludge to fool SREffects into thinking we've computed twiss parameters.
     * First, propagate the users parameters around the beamline.
     * Second, restore the global parameters in case the beamline isn't right (which is likely).
     */
    compute_twiss_parameters(run, beamline, 
                             beamline->closed_orbit?beamline->closed_orbit->centroid:NULL, 0, 0,
                             beamline->twiss0->betax, beamline->twiss0->alphax, 
                             beamline->twiss0->etax, beamline->twiss0->etapx, 
                             beamline->twiss0->betay, beamline->twiss0->alphay, 
                             beamline->twiss0->etay, beamline->twiss0->etapy,
                             &unstable0);
    setLinearChromaticTrackingValues(beamline);
    
  }
  else {
    compute_twiss_parameters(run, beamline, 
                             beamline->closed_orbit?beamline->closed_orbit->centroid:NULL, matched, 
                             radiation_integrals,
                             beta_x, alpha_x, eta_x, etap_x, beta_y, alpha_y, eta_y, etap_y,
                             &unstable0);
  }
  if (unstable)
    *unstable = unstable0;
}

void copy_doubles(double *target, double *source, long n)
{
  log_entry("copy_doubles");
  while (n--)
    *target++ = *source++;
  log_exit("copy_doubles");
}


long has_aperture(ELEMENT_LIST *elem);

double find_acceptance(
                       ELEMENT_LIST *elem, long plane, RUN *run, char **name, double *z
                       )
{
  double beta, acceptance, tmp;
  double tube_aperture, aperture, centroid, offset;
  double other_centroid, a, b, a_tube, b_tube, aperture1;
  long aperture_set, elliptical_tube, tube_set;
  SCRAPER *scraper;
  ELEMENT_LIST *ap_elem;

  log_entry("find_acceptance");

  acceptance = tube_aperture = a_tube = b_tube = 0;
  elliptical_tube = tube_set = 0;
  *name = NULL;
  *z = -DBL_MAX;
  while (elem) {
    beta = *(((double*)elem->twiss) + (plane?TWISS_Y_OFFSET:0));
    centroid = *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET);
    other_centroid = *(((double*)elem->twiss) + (plane?0:1) + TWISS_CENT_OFFSET);
    aperture = 0;
    aperture_set = 0;
    if (!has_aperture(elem) && has_aperture(elem->succ))
      ap_elem = elem->succ;
    else
      ap_elem = elem;
    switch (ap_elem->type) {
    case T_MAXAMP:
      if (((MAXAMP*)ap_elem->p_elem)->elliptical==0) {
        elliptical_tube = tube_set = aperture_set = 0;
        if (plane)
          tube_aperture = ((MAXAMP*)ap_elem->p_elem)->y_max;
        else
          tube_aperture = ((MAXAMP*)ap_elem->p_elem)->x_max;
        if (tube_aperture>0) {
          aperture = tube_aperture - fabs(centroid);
          aperture_set = 1;
          tube_set = 1;
        }
      } else {
        elliptical_tube = 1;
        tube_set = aperture_set = 0;
        if (plane) {
          a_tube = ((MAXAMP*)ap_elem->p_elem)->y_max;
          b_tube = ((MAXAMP*)ap_elem->p_elem)->x_max;
        }
        else {
          a_tube = ((MAXAMP*)ap_elem->p_elem)->x_max;
          b_tube = ((MAXAMP*)ap_elem->p_elem)->y_max;
        }
        if (a_tube>0) {
          if (b_tube>0) {
            if ((aperture = sqr(a_tube)*(1 - sqr(other_centroid/b_tube)))<0)
              aperture = 0;
          } else
            aperture = sqr(a_tube);
          aperture = sqrt(aperture) - fabs(centroid);
          aperture_set = 1;
          tube_set = 1;
        }
      }
      break;
    case T_RCOL:
      if (plane) {
        if (((RCOL*)ap_elem->p_elem)->y_max>0) {
          aperture = ((RCOL*)ap_elem->p_elem)->y_max;
          offset = ((RCOL*)ap_elem->p_elem)->dy;
          aperture_set = 1;
        }
      } else {
        if (((RCOL*)ap_elem->p_elem)->x_max>0) {
          aperture = ((RCOL*)ap_elem->p_elem)->x_max;
          offset = ((RCOL*)ap_elem->p_elem)->dx;
          aperture_set = 1;
        }
      }
      if (aperture_set)
        aperture = aperture - fabs(centroid-offset);
      break;
    case T_ECOL:
      if (plane) {
        a = ((ECOL*)ap_elem->p_elem)->y_max;
        b = ((ECOL*)ap_elem->p_elem)->x_max;
        if (a>0) {
          centroid -= ((ECOL*)ap_elem->p_elem)->dy;
          other_centroid -= ((ECOL*)ap_elem->p_elem)->dx;
          aperture_set = 1;
        }
      } else {
        a = ((ECOL*)ap_elem->p_elem)->x_max;
        b = ((ECOL*)ap_elem->p_elem)->y_max;
        if (a>0) {
          centroid -= ((ECOL*)ap_elem->p_elem)->dx;
          other_centroid -= ((ECOL*)ap_elem->p_elem)->dy;
          aperture_set = 1;
        }
      }
      if (aperture_set) {
        if (b>0) {
          if ((aperture = sqr(a)*(1 - sqr(other_centroid/b)))<0)
            aperture = 0;
        } else
          aperture = sqr(a);
        aperture = sqrt(aperture) - fabs(centroid);
      }
      break;
    case T_SCRAPER:
      scraper = (SCRAPER*)ap_elem->p_elem;
      if (plane) {
        if (scraper->direction==1 || scraper->direction==3) {
          offset = scraper->dy;
          aperture = scraper->position + offset - centroid;
          if (scraper->direction==3)
            aperture = -aperture;
          aperture_set = 1;
        }
      } else {
        if (scraper->direction==0 || scraper->direction==2) {
          offset = scraper->dx;
          aperture = scraper->position + offset - centroid;
          if (scraper->direction==2)
            aperture = -aperture;
          aperture_set = 1;
        }
      }
      break;
    default:
      break;
    }
    if (run->apertureData.initialized) {
      /* Find aperture from the aperture data file */
      double dx, dy, xsize, ysize;
      if (interpolateApertureData(elem->end_pos, &(run->apertureData), &dx, &dy, &xsize, &ysize)) {
        centroid = *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET);
        if (plane) {
          aperture1 = ysize;
          offset = dy;
        } else {
          aperture1 = xsize;
          offset = dx;
        }
        aperture1 = aperture1 - fabs(centroid-offset);
        if (!aperture_set || aperture1<aperture) {
          aperture = aperture1;
          aperture_set = 1;
        }
      }
    }
    if (tube_set) {
      if (elliptical_tube) {
        centroid = *(((double*)elem->twiss) + (plane?1:0) + TWISS_CENT_OFFSET);
        other_centroid = *(((double*)elem->twiss) + (plane?0:1) + TWISS_CENT_OFFSET);
        if (b_tube>0) {
          if ((aperture1 = sqr(a_tube)*(1 - sqr(other_centroid/b_tube)))<0)
            aperture1 = 0;
        } else 
          aperture1 = sqr(a_tube);
        aperture1 = sqrt(aperture1) - fabs(centroid);
      } else
        aperture1 = tube_aperture - fabs(centroid);
      if (!aperture_set || aperture1<aperture) {
        aperture = aperture1;
        aperture_set = 1;
      }
    }
    
    if (aperture_set) {
      if (((tmp=sqr(aperture)/beta)<acceptance || !acceptance)) {
        *name = elem->name;
        *z = elem->end_pos;
        acceptance = tmp;
      }
    } else
      aperture = COORD_LIMIT;

    if (plane)
      elem->twiss->apy = aperture;
    else
      elem->twiss->apx = aperture;

    elem = elem->succ;
  }
  log_exit("find_acceptance");
  return(acceptance);
}

long has_aperture(ELEMENT_LIST *elem)
{
  if (!elem)
    return(0);
  switch (elem->type) {
  case T_MAXAMP: case T_RCOL: case T_ECOL: case T_SCRAPER:
    return(1);
  default:
    return(0);
  }
}

void modify_rfca_matrices(ELEMENT_LIST *eptr, long order)
/* Replace the matrices for rf cavities with drift matrices. 
   Used prior to periodic twiss parameter computations. */
{
  while (eptr) {
    if (eptr->type==T_RFCA || eptr->type==T_MODRF || eptr->type==T_RFCW) {
      if (eptr->matrix) {
        free_matrices(eptr->matrix);
        tfree(eptr->matrix);
      }
      switch (eptr->type) {
      case T_RFCA:
        if (((RFCA*)eptr->p_elem)->change_p0)
          bomb("can't treat cavities like drifts when CHANGE_P0=1", NULL);
        eptr->matrix = drift_matrix(((RFCA*)eptr->p_elem)->length, order);
        break;
      case T_RFCW:
        if (((RFCW*)eptr->p_elem)->change_p0)
          bomb("can't treat cavities like drifts when CHANGE_P0=1", NULL);
        eptr->matrix = drift_matrix(((RFCW*)eptr->p_elem)->length, order);
        break;
      case T_MODRF:
        eptr->matrix = drift_matrix(((MODRF*)eptr->p_elem)->length, order);
        break;
      }
    }
    eptr = eptr->succ;
  }
  mustResetRfcaMatrices = 1;
}

void reset_rfca_matrices(ELEMENT_LIST *eptr, long order)
/* Delete the rf cavity matrices so they'll get updated */
{
  while (eptr) {
    if (eptr->type==T_RFCA || eptr->type==T_MODRF || eptr->type==T_RFCW) {
      if (eptr->matrix) {
        free_matrices(eptr->matrix);
        tfree(eptr->matrix);
        eptr->matrix = NULL;
      }
    }
    eptr = eptr->succ;
  }
}

#define ASSIGN_MINMAX(min, max, value) ((min>value?min=value:1),(max<value?max=value:1))

void compute_twiss_statistics(LINE_LIST *beamline, TWISS *twiss_ave, TWISS *twiss_min, TWISS *twiss_max)
{
  ELEMENT_LIST *eptr;
  double dz, end_pos=0.0;
  long nElems;
  
  if (!twiss_ave) {
    fprintf(stdout, "error: NULL twiss_ave pointer in compute_twiss_statistics\n");
    fflush(stdout);
    abort();
  }
  if (!twiss_min) {
    fprintf(stdout, "error: NULL twiss_min pointer in compute_twiss_statistics\n");
    fflush(stdout);
    abort();
  }
  if (!twiss_max) {
    fprintf(stdout, "error: NULL twiss_max pointer in compute_twiss_statistics\n");
    fflush(stdout);
    abort();
  }
  
  zero_memory(twiss_ave, sizeof(*twiss_ave));
  twiss_min->betax = twiss_min->alphax = twiss_min->etax = twiss_min->etapx = DBL_MAX;
  twiss_min->betay = twiss_min->alphay = twiss_min->etay = twiss_min->etapy = DBL_MAX;
  twiss_max->betax = twiss_max->alphax = twiss_max->etax = twiss_max->etapx = -DBL_MAX;
  twiss_max->betay = twiss_max->alphay = twiss_max->etay = twiss_max->etapy = -DBL_MAX;

  eptr = beamline->elem_twiss;
  nElems = 0;
  while (eptr) {
    ASSIGN_MINMAX(twiss_min->betax, twiss_max->betax, eptr->twiss->betax);
    ASSIGN_MINMAX(twiss_min->alphax, twiss_max->alphax, eptr->twiss->alphax);
    ASSIGN_MINMAX(twiss_min->etax, twiss_max->etax, eptr->twiss->etax);
    ASSIGN_MINMAX(twiss_min->etapx, twiss_max->etapx, eptr->twiss->etapx);
    ASSIGN_MINMAX(twiss_min->betay, twiss_max->betay, eptr->twiss->betay);
    ASSIGN_MINMAX(twiss_min->alphay, twiss_max->alphay, eptr->twiss->alphay);
    ASSIGN_MINMAX(twiss_min->etay, twiss_max->etay, eptr->twiss->etay);
    ASSIGN_MINMAX(twiss_min->etapy, twiss_max->etapy, eptr->twiss->etapy);
    if (eptr->pred && eptr->pred->twiss && eptr->twiss) {
      dz = (eptr->end_pos - eptr->pred->end_pos)/2;
      twiss_ave->betax += (eptr->pred->twiss->betax + eptr->twiss->betax)*dz;
      twiss_ave->alphax += (eptr->pred->twiss->alphax + eptr->twiss->alphax)*dz;
      twiss_ave->etax += (eptr->pred->twiss->etax + eptr->twiss->etax)*dz;
      twiss_ave->etapx += (eptr->pred->twiss->etapx + eptr->twiss->etapx)*dz;
      twiss_ave->betay += (eptr->pred->twiss->betay + eptr->twiss->betay)*dz;
      twiss_ave->alphay += (eptr->pred->twiss->alphay + eptr->twiss->alphay)*dz;
      twiss_ave->etay += (eptr->pred->twiss->etay + eptr->twiss->etay)*dz;
      twiss_ave->etapy += (eptr->pred->twiss->etapy + eptr->twiss->etapy)*dz;
      nElems++;
    }
    end_pos = eptr->end_pos;
    eptr = eptr->succ;
  }
  if (nElems==0) {
    twiss_ave->betax  = (twiss_min->betax  + twiss_max->betax )/2;
    twiss_ave->alphax = (twiss_min->alphax + twiss_max->alphax)/2;
    twiss_ave->etax   = (twiss_min->etax   + twiss_max->etax  )/2;
    twiss_ave->etapx  = (twiss_min->etapx  + twiss_max->etapx )/2;
    twiss_ave->betay  = (twiss_min->betay  + twiss_max->betay )/2;
    twiss_ave->alphay = (twiss_min->alphay + twiss_max->alphay)/2;
    twiss_ave->etay   = (twiss_min->etay   + twiss_max->etay  )/2;
    twiss_ave->etapy  = (twiss_min->etapy  + twiss_max->etapy )/2;
  }
  else if (end_pos) {
    twiss_ave->betax  /= end_pos;
    twiss_ave->alphax /= end_pos;
    twiss_ave->etax   /= end_pos;
    twiss_ave->etapx  /= end_pos;
    twiss_ave->betay  /= end_pos;
    twiss_ave->alphay /= end_pos;
    twiss_ave->etay   /= end_pos;
    twiss_ave->etapy  /= end_pos;
  }
}

void incrementRadIntegrals(RADIATION_INTEGRALS *radIntegrals, double *dI, 
                           ELEMENT_LIST *elem, 
                           double beta0, double alpha0, double gamma0, double eta0, double etap0, 
                           double betay0, double alphay0, double gammay0, double etay0, double etapy0, 
                           double *coord)
{
  /* compute contribution to radiation integrals */
  long isBend;
  BEND *bptr;
  KSBEND *kbptr;
  CSBEND *cbptr;
  CSRCSBEND *csrbptr;
  QUAD *qptr;
  KQUAD *qptrk;
  SEXT *sptr;
  KSEXT *sptrk;
  double length=0.0, angle=0.0, E1=0.0, E2=0.0, K1=0.0;
  double k2, rho, k, kl, kick=0.0;
  double I1, I2, I3, I4, I5;
  double alpha1, gamma1, etap1, eta2, sin_kl, cos_kl;
  double etaAve=0, etaK1_rhoAve=0, HAve=0, h, K2=0.0, dx=0.0;

  I1 = I2 = I3 = I4 = I5 = 0;

#ifdef DEBUG
  fprintf(stderr, "incrementRadIntegrals: %s at %e: beta=%e, alpha=%e, eta=%e, etap=%e\n", 
	  elem->name, elem->end_pos, beta0, alpha0, eta0, etap0);
#endif
  if (elem->type==T_WIGGLER) {
    WIGGLER *wiggler;
    wiggler = (WIGGLER*)(elem->p_elem);
    AddWigglerRadiationIntegrals(wiggler->length, wiggler->poles, wiggler->radiusInternal,
				 eta0, etap0, 
				 beta0, alpha0,
				 &I1, &I2, &I3, &I4, &I5);
    radIntegrals->RI[0] += I1;
    radIntegrals->RI[1] += I2;
    radIntegrals->RI[2] += I3;
    radIntegrals->RI[3] += I4;
    radIntegrals->RI[4] += I5;
  } else if (elem->type==T_CWIGGLER) {
    CWIGGLER *wiggler;
    wiggler = (CWIGGLER*)(elem->p_elem);
    if (wiggler->BPeak[1]) {
      /* By */
      AddWigglerRadiationIntegrals(wiggler->length, wiggler->periods*2, 
                                   wiggler->radiusInternal[1],
                                   eta0, etap0, beta0, alpha0,
                                   &I1, &I2, &I3, &I4, &I5);
      radIntegrals->RI[0] += I1;
      radIntegrals->RI[1] += I2;
      radIntegrals->RI[2] += I3;
      radIntegrals->RI[3] += I4;
      radIntegrals->RI[4] += I5;
    }
    if (wiggler->BPeak[0]) {
      /* Bx */
      I1 = I2 = I3 = I4 = I5 = 0;
      AddWigglerRadiationIntegrals(wiggler->length, wiggler->periods*2, 
                                   wiggler->radiusInternal[0],
                                   etay0, etapy0, betay0, alphay0,
                                   &I1, &I2, &I3, &I4, &I5);
      radIntegrals->RI[1] += I2;
      radIntegrals->RI[2] += I3;
    }
  } else if (elem->type==T_UKICKMAP) {
    UKICKMAP *ukmap;
    ukmap = (UKICKMAP*)(elem->p_elem);
    if (ukmap->radiusInternal>0) {
      AddWigglerRadiationIntegrals(ukmap->length, 2*ukmap->periods, ukmap->radiusInternal,
                                   eta0, etap0,
                                   beta0, alpha0,
                                   &I1, &I2, &I3, &I4, &I5);
      radIntegrals->RI[0] += I1;
      radIntegrals->RI[1] += I2;
      radIntegrals->RI[2] += I3;
      radIntegrals->RI[3] += I4;
      radIntegrals->RI[4] += I5;
    }
  } else {
    isBend = 1;
    switch (elem->type) {
    case T_QUAD:
    case T_KQUAD:
      if (!coord && !coord[0]) {
	isBend = 0;
	break;
      }
      switch (elem->type) {
      case T_QUAD:
	qptr = (QUAD*)(elem->p_elem);
	length = qptr->length;
	kick = qptr->xkick*qptr->xKickCalibration;
	K1 = qptr->k1;
	dx = qptr->dx;
	break;
      case T_KQUAD:
	qptrk = (KQUAD*)(elem->p_elem);
	length = qptrk->length;
	kick = qptrk->xkick*qptrk->xKickCalibration;
	K1 = qptrk->k1;
	dx = qptrk->dx;
	break;
      }
      if (!(h = K1*(coord[0]-dx))) {
	isBend = 0;
	break;
      }
      /* The kick is subtracted because a positive angle bends toward the
       * right (negative x) 
       */
      angle = length*h - kick;
      E1 = E2 = 0;
      isBend = 1;
      break;
    case T_SEXT:
    case T_KSEXT:
      if (!coord && !coord[0]) {
	isBend = 0;
	break;
      }
      switch (elem->type) {
      case T_SEXT:
	sptr = (SEXT*)(elem->p_elem);
	length = sptr->length;
	K2 = sptr->k2;
	dx = sptr->dx;
	break;
      case T_KSEXT:
	sptrk = (KSEXT*)(elem->p_elem);
	length = sptrk->length;
	K2 = sptrk->k2;
	dx = sptrk->dx;
	break;
      }
      if (!(h = K2*sqr(coord[0]-dx)/2)) {
	isBend = 0;
	break;
      }
      K1 = K2*(coord[0]-dx);
      angle = length*h;
      E1 = E2 = 0;
      isBend = 1;
      break;
    case T_SBEN:
    case T_RBEN:
      bptr = (BEND*)(elem->p_elem);
      length = bptr->length;
      angle = bptr->angle;
      E1 = bptr->e1*(bptr->edgeFlags&BEND_EDGE1_EFFECTS?1:0);
      E2 = bptr->e2*(bptr->edgeFlags&BEND_EDGE2_EFFECTS?1:0);
      K1 = bptr->k1;
      break;
    case T_KSBEND:
      kbptr = (KSBEND*)(elem->p_elem);
      length = kbptr->length;
      angle = kbptr->angle;
      E1 = kbptr->e1*(kbptr->flags&BEND_EDGE1_EFFECTS?1:0);
      E2 = kbptr->e2*(kbptr->flags&BEND_EDGE2_EFFECTS?1:0);
      K1 = kbptr->k1;
      break;
    case T_CSBEND:
      cbptr = (CSBEND*)(elem->p_elem);
      length = cbptr->length;
      angle = cbptr->angle;
      E1 = cbptr->e1*(cbptr->edgeFlags&BEND_EDGE1_EFFECTS?1:0);
      E2 = cbptr->e2*(cbptr->edgeFlags&BEND_EDGE2_EFFECTS?1:0);
      K1 = cbptr->k1;
      break;
    case T_CSRCSBEND:
      csrbptr = (CSRCSBEND*)(elem->p_elem);
      length = csrbptr->length;
      angle = csrbptr->angle;
      E1 = csrbptr->e1*(csrbptr->edgeFlags&BEND_EDGE1_EFFECTS?1:0);
      E2 = csrbptr->e2*(csrbptr->edgeFlags&BEND_EDGE2_EFFECTS?1:0);
      K1 = csrbptr->k1;
      break;
    default:
      isBend = 0;
      break;
    }
    if (isBend && angle!=0) {
      /* This code seems like it should be here, but including it makes the
	 deviations from tracking results larger.
	 if (coord) {
	 K1 /= 1+coord[5];
	 angle /= 1+coord[5];
	 }
      */
      rho = length/angle;
      k2 = K1+1./(rho*rho);
      /* equations are from SLAC 1193 */
      k = sqrt(fabs(k2));
      kl = k*length;
      if (k2<0) {
	cos_kl = cosh(kl);
	sin_kl = sinh(kl);
      } else {
	sin_kl = sin(kl);
	cos_kl = cos(kl);
      }
      etap1 = etap0 + eta0/rho*tan(E1);
      eta2  = eta0*cos_kl + etap1*sin_kl/k + (1-cos_kl)/(rho*k2);
      alpha1 = alpha0 - beta0/rho*tan(E1);
      gamma1 = (1+sqr(alpha1))/beta0;
      if (kl>1e-2) {
	etaAve = eta0*sin_kl/kl + etap1*(1-cos_kl)/(k2*length) +
	  (kl-sin_kl)/(k2*kl*rho);
	etaK1_rhoAve =  -etaAve*K1/rho + (eta0*tan(E1)+eta2*tan(E2))/(2*length*sqr(rho));
	HAve = gamma1*sqr(eta0) + 2*alpha1*eta0*etap1 + beta0*sqr(etap1) 
	  + 2*angle*( -(gamma1*eta0+alpha1*etap1)*(kl-sin_kl)/(k*k2*sqr(length)) +
		      (alpha1*eta0+beta0*etap1)*(1-cos_kl)/(k2*sqr(length))
		      )
	  + sqr(angle)*(gamma1*(3*kl-4*sin_kl+sin_kl*cos_kl)/(2*k*k2*k2*ipow(length,3)) 
			- alpha1*sqr(1-cos_kl)/(k2*k2*ipow(length,3))
			+ beta0*(kl-cos_kl*sin_kl)/(2*kl*k2*sqr(length)));
      } 
      if (kl<=1e-2 || HAve<0) {
	double kl2, kl4, kl6, kl8, l2, l3, l4, h2;
	double T0, T2, T4, T6, T8;
	h = 1./rho;
	h2 = h*h;
	kl2 = kl*kl;
	kl4 = kl2*kl2;
	kl6 = kl4*kl2;
	kl8 = kl6*kl2;
	l2 = length*length;
	l3 = l2*length;
	l4 = l3*length;

	T0 = eta0 + (length*(3*etap1 + h*length))/6.;
	T2 = (-20*eta0 - length*(5*etap1 + h*length))/120.;
	T4 = (42*eta0 + length*(7*etap1 + h*length))/5040.;
	T6 =  (-72*eta0 - length*(9*etap1 + h*length))/362880.;
	T8 = 2.505210838544172e-8*(110.*eta0 + length*(11.*etap1 + h*length));
	etaAve = T0 + T2*kl2 + T4*kl4 + T6*kl6 + T8*kl8;
	etaK1_rhoAve =  -etaAve*K1/rho + (eta0*tan(E1)+eta2*tan(E2))/(2*length*sqr(rho));

	T0 = (ipow(eta0,2)*gamma1 - 
	      (eta0*gamma1*h*l2)/3. + 
	      (gamma1*h2*l4)/20. + 
	      beta0*(ipow(etap1,2) + etap1*h*length + 
		     (h2*l2)/3.) + 
	      alpha1*(eta0*(2*etap1 + h*length) - 
		      (h*l2*(4*etap1 + 3*h*length))/12.));
	T2 = -(h*l2*
	       (70*beta0*etap1 - 14*eta0*gamma1*length + 56*beta0*h*length + 
		5*gamma1*h*l3 + 
		7*alpha1*(10*eta0 - length*(2*etap1 + 5*h*length))))/(840.*length);
	T4 = (h*l2*
	      (-8*eta0*gamma1*length + 7*gamma1*h*l3 + 
	       8*beta0*(7*etap1 + 16*h*length) + 
	       alpha1*(56*eta0 - length*(8*etap1 + 63*h*length))))/(20160.*length);
	T6 =  (-2.505210838544172e-7*h*l2*
	       (-22.*eta0*gamma1*length + 51.*gamma1*h*l3 + 
		22.*beta0*(9.*etap1 + 64.*h*length) + 
		11.*alpha1*(18.*eta0 - 1.*length*(2.*etap1 + 51.*h*length))))/length;
	T8 = (9.635426302092969e-10*h*l2*
	      (-52.*eta0*gamma1*length + 341.*gamma1*h*l3 + 
	       52.*beta0*(11.*etap1 + 256.*h*length) + 
	       13.*alpha1*(44.*eta0 - 1.*length*(4.*etap1 + 341.*h*length))))/length;
	HAve = T0 + T2*kl2 + T4*kl4 + T6*kl6 + T8*kl8;
      }
      I1 = etaAve*length/rho;
      I2 = length/sqr(rho);
      I3 = I2/fabs(rho);
      I4 = I2/rho*etaAve - 2*length*etaK1_rhoAve;
      I5 = HAve*I3;
      radIntegrals->RI[0] += I1;
      radIntegrals->RI[1] += I2;
      radIntegrals->RI[2] += I3;
      radIntegrals->RI[3] += I4;
      radIntegrals->RI[4] += I5;
    }
  }
  if (dI) {
    dI[0] = I1;
    dI[1] = I2;
    dI[2] = I3;
    dI[3] = I4;
    dI[4] = I5;
  }
}


void completeRadiationIntegralComputation(RADIATION_INTEGRALS *RI, double Po, double revolutionLength)
{    
    double Rce, gamma;
    gamma = sqrt(sqr(Po)+1);
    Rce = sqr(particleCharge)/(1e7*particleMass);
    RI->Uo = particleMassMV*Rce*RI->RI[1]*2./3.*ipow(gamma,4);
    RI->Jx = 1 - RI->RI[3]/RI->RI[1];
    RI->Jdelta = 3 - RI->Jx;
    RI->Jy = 1;
    RI->tauy = 1./(Rce/3*ipow(gamma,3)*c_mks/revolutionLength*RI->RI[1]);
    RI->taux = RI->tauy*RI->Jy/RI->Jx;
    RI->taudelta = RI->tauy*RI->Jy/RI->Jdelta;
    RI->sigmadelta = gamma*sqrt(55./32./sqrt(3.)*hbar_mks/(particleMass*c_mks)*RI->RI[2]/(2*RI->RI[1]+RI->RI[3]));
    RI->ex0 = sqr(gamma)*55./32./sqrt(3.)*hbar_mks/(particleMass*c_mks)*RI->RI[4]/(RI->RI[1]-RI->RI[3]);
    RI->computed = 1;
  }

void LoadStartingTwissFromFile(double *betax, double *betay, double *alphax, double *alphay,
                               double *etax, double *etaxp, double *etay, double *etayp,
                               char *filename, char *elementName, long elementOccurrence)
{
  SDDS_DATASET SDDSin;
  long rows=0, rowOfInterest;
  double *betaxData, *betayData=NULL, *alphaxData=NULL, *alphayData=NULL;
  double *etaxData=NULL, *etayData=NULL, *etaxpData=NULL, *etaypData=NULL;
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename) || 
      SDDS_ReadPage(&SDDSin)!=1)
    SDDS_Bomb("problem reading Twiss reference file");
  if (SDDS_CheckColumn(&SDDSin, "betax", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "betay", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "alphax", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "alphay", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etax", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etay", "m", SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etaxp", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "etayp", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "ElementName", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK)
    SDDS_Bomb("invalid/missing columns in Twiss reference file");
  if (elementName) {
    if (!SDDS_SetRowFlags(&SDDSin, 1) ||
        (rows=SDDS_MatchRowsOfInterest(&SDDSin, "ElementName", elementName, SDDS_AND))<=0)
      SDDS_Bomb("Problem finding data for beta function reference.  Check for existence of element.");
    if (elementOccurrence>0 && elementOccurrence>rows)
      SDDS_Bomb("Too few occurrences of reference element in beta function reference file.");
  } 
  if ((rows=SDDS_CountRowsOfInterest(&SDDSin))<1)
    SDDS_Bomb("No data in beta function reference file.");
    
  if (!(betaxData=SDDS_GetColumnInDoubles(&SDDSin, "betax")) ||
      !(betayData=SDDS_GetColumnInDoubles(&SDDSin, "betay")) ||
      !(alphaxData=SDDS_GetColumnInDoubles(&SDDSin, "alphax")) ||
      !(alphayData=SDDS_GetColumnInDoubles(&SDDSin, "alphay")) || 
      !(etaxData=SDDS_GetColumnInDoubles(&SDDSin, "etax")) ||
      !(etayData=SDDS_GetColumnInDoubles(&SDDSin, "etay")) ||
      !(etaxpData=SDDS_GetColumnInDoubles(&SDDSin, "etaxp")) ||
      !(etaypData=SDDS_GetColumnInDoubles(&SDDSin, "etayp")) )
    SDDS_Bomb("Problem getting data for beta function reference.");
  if (elementName && elementOccurrence>0)
    rowOfInterest = elementOccurrence-1;
  else
    rowOfInterest = rows-1;
  *betax = betaxData[rowOfInterest];
  *betay = betayData[rowOfInterest];
  *alphax = alphaxData[rowOfInterest];
  *alphay = alphayData[rowOfInterest];
  *etax = etaxData[rowOfInterest];
  *etay = etayData[rowOfInterest];
  *etaxp = etaxpData[rowOfInterest];
  *etayp = etaypData[rowOfInterest];
  free(betaxData);
  free(betayData);
  free(alphaxData);
  free(alphayData);
  free(etaxData);
  free(etayData);
  free(etaxpData);
  free(etaypData);
}

void setupTuneShiftWithAmplitude(NAMELIST_TEXT *nltext, RUN *run)
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&tune_shift_with_amplitude, nltext);
  if (echoNamelists) print_namelist(stdout, &tune_shift_with_amplitude);

  if (tune_shift_with_amplitude_struct.turns<100 && tune_shift_with_amplitude_struct.turns!=0)
    bomb("too few turns requested (tune_shift_with_amplitude)", NULL);
  if (!tune_shift_with_amplitude_struct.turns && tune_shift_with_amplitude_struct.spread_only)
    bomb("spread_only requires turns!=0", NULL);
  if (tune_shift_with_amplitude_struct.turns) {
    if (tune_shift_with_amplitude_struct.x0<=0 ||
        tune_shift_with_amplitude_struct.y0<=0)
      bomb("x0 or y0 is zero or negative (tune_shift_with_amplitude)", NULL);
    if (!tune_shift_with_amplitude_struct.spread_only &&
	(tune_shift_with_amplitude_struct.x1<tune_shift_with_amplitude_struct.x0*10 ||
	 tune_shift_with_amplitude_struct.y1<tune_shift_with_amplitude_struct.y0*10))
      bomb("x1 or y1 is too small (tune_shift_with_amplitude)", NULL);
  }
  doTuneShiftWithAmplitude = 1;
  if (tune_shift_with_amplitude_struct.tune_output) {
    char *filename;
    filename = compose_filename(tune_shift_with_amplitude_struct.tune_output, run->rootname);
    if (!SDDS_InitializeOutput(&SDDSTswaTunes, SDDS_BINARY, 0, NULL, NULL, 
			       filename) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "x", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "y", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "Ax", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "Ay", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "nux", NULL, SDDS_DOUBLE) ||
	!SDDS_DefineSimpleColumn(&SDDSTswaTunes, "nuy", NULL, SDDS_DOUBLE) ||
	!SDDS_WriteLayout(&SDDSTswaTunes)) {
      fprintf(stdout, "Unable set up file %s\n", 
	      tune_shift_with_amplitude_struct.tune_output);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (filename!=tune_shift_with_amplitude_struct.tune_output)
      free(filename);
  }
}

double QElement(double ****Q, long i1, long i2, long i3, long i4)
{
  if (i3<i4)
    SWAP_LONG(i3, i4);
  if (i2<i3)
    SWAP_LONG(i2, i3);
  if (i3<i4)
    SWAP_LONG(i3, i4);
  if (i2<i3 || i3<i4)
    bomb("it didn't work",NULL);
  return Q[i1][i2][i3][i4];
}

void computeTuneShiftWithAmplitude(double dnux_dA[N_TSWA][N_TSWA], double dnuy_dA[N_TSWA][N_TSWA],
				   double *xTuneExtrema, double *yTuneExtrema, 
                                   TWISS *twiss, double *tune, VMATRIX *M, LINE_LIST *beamline, 
                                   RUN *run, double *startingCoord)
{
#define TSWA_TRACKING_EXTRA_PTS 15
#define TSWA_TRACKING_PTS (N_TSWA+TSWA_TRACKING_EXTRA_PTS)
  double tune1[2];
  double Ax[TSWA_TRACKING_PTS], Ay[TSWA_TRACKING_PTS];
  double x0[TSWA_TRACKING_PTS], y0[TSWA_TRACKING_PTS];
  double xTune[TSWA_TRACKING_PTS][TSWA_TRACKING_PTS], yTune[TSWA_TRACKING_PTS][TSWA_TRACKING_PTS];
  short lost[TSWA_TRACKING_PTS][TSWA_TRACKING_PTS];
  double tuneLowerLimit[2], tuneUpperLimit[2];
  double result, maxResult;
  double x, y;
  long ix, iy, tries, gridSize, nLost=0;
  double upperLimit[2], lowerLimit[2];
  MATRIX *AxAy, *Coef, *Nu, *AxAyTr, *Mf, *MfInv, *AxAyTrNu;
  long i, j, m, n, ix1, iy1, nTSWA, order;

  if (tune_shift_with_amplitude_struct.turns==0) {
    /* use the matrix only without tracking */
    xTuneExtrema[0] = xTuneExtrema[1] = 0;
    yTuneExtrema[0] = yTuneExtrema[1] = 0;
    computeTuneShiftWithAmplitudeM(dnux_dA, dnuy_dA, twiss, tune, M);
    return;
  }

  if ((gridSize=tune_shift_with_amplitude_struct.grid_size)>TSWA_TRACKING_PTS) {
    gridSize = tune_shift_with_amplitude_struct.grid_size = TSWA_TRACKING_PTS;
    fprintf(stdout, "Error: grid_size for TSWA is limited to %ld\n", gridSize);
    exit(1);
  }
  
  if (tune_shift_with_amplitude_struct.tune_output &&
      !SDDS_StartPage(&SDDSTswaTunes, gridSize*gridSize)) {
    fprintf(stdout, "Problem starting SDDS page for TSWA tune output\n");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  /* use tracking and NAFF */
  tries = tune_shift_with_amplitude_struct.scaling_iterations;
  upperLimit[0] = upperLimit[1] = 1;
  lowerLimit[0] = lowerLimit[1] = 0;
  while (tries--) {
    m = 0;  /* number of tune points */
    tuneLowerLimit[0] = tuneLowerLimit[1] = 0;
    tuneUpperLimit[0] = tuneUpperLimit[1] = 0;
    nLost = 0;
    for (ix=0; ix<gridSize; ix++) {
      x0[ix] = 
	x = sqrt(ix*
		 sqr((tune_shift_with_amplitude_struct.x1-tune_shift_with_amplitude_struct.x0))/(gridSize-1)
		 + sqr(tune_shift_with_amplitude_struct.x0));
      Ax[ix] = sqr(x)*(1+sqr(twiss->alphax))/twiss->betax;
      for (iy=0; iy<gridSize; iy++) {
	lost[ix][iy] = 0;
        xTune[ix][iy] = -1;
        yTune[ix][iy] = -1;
        if ((tune_shift_with_amplitude_struct.sparse_grid &&
             !(ix==0 || iy==0 || ix==iy)) ||
               (tune_shift_with_amplitude_struct.lines_only &&
                !(ix==0 || iy==0)))
          continue;
        m++;
        y0[iy] = 
	  y = sqrt(iy*
		   sqr((tune_shift_with_amplitude_struct.y1-tune_shift_with_amplitude_struct.y0))/(gridSize-1)
		   + sqr(tune_shift_with_amplitude_struct.y0));
        Ay[iy] = sqr(y)*(1+sqr(twiss->alphay))/twiss->betay;
	if (ix==0 && iy==0) {
	  tuneLowerLimit[0] = tuneLowerLimit[1] = 0;
	  tuneUpperLimit[0] = tuneUpperLimit[1] = 0;
	} else {
	  if ((ix1 = ix-1)<0)
	    ix1 = 0;
	  if ((iy1 = iy-1)<0)
	    iy1 = 0;
	  tuneLowerLimit[0] = xTune[ix1][iy1] - tune_shift_with_amplitude_struct.nux_roi_width/2;
	  tuneUpperLimit[0] = xTune[ix1][iy1] + tune_shift_with_amplitude_struct.nux_roi_width/2;
	  tuneLowerLimit[1] = yTune[ix1][iy1] - tune_shift_with_amplitude_struct.nux_roi_width/2;
	  tuneUpperLimit[1] = yTune[ix1][iy1] + tune_shift_with_amplitude_struct.nux_roi_width/2;
	}
        if (!computeTunesFromTracking(tune1, NULL, M, beamline, run, startingCoord,
                                      x, y,
                                      tune_shift_with_amplitude_struct.turns, 0, NULL,
				      tuneLowerLimit, tuneUpperLimit, 0)) {
	  lost[ix][iy] = 1;
	  nLost ++;
        }
        xTune[ix][iy] = tune1[0];
        yTune[ix][iy] = tune1[1];
        if (tune_shift_with_amplitude_struct.verbose)
          fprintf(stdout, "Tunes for TSWA: x=%e, y=%e, nux=%.15e, nuy=%.15e\n",
                  x, y, tune1[0], tune1[1]);
      }
    }
    if (tune_shift_with_amplitude_struct.verbose)
      fprintf(stdout, "All tunes computed for TSWA (%ld particles lost)\n", nLost);

    maxResult = -DBL_MAX;
    if (!tune_shift_with_amplitude_struct.spread_only) {
      if (nLost) {
	fprintf(stdout, "Amplitude too large for TSWA computation---particles lost\n");
      } else {
	for (ix=0; ix<gridSize; ix++) {
	  for (iy=0; iy<gridSize; iy++) {
	    if (ix==0 && iy==0)
	      continue;
            if ((tune_shift_with_amplitude_struct.sparse_grid &&
                 !(ix==0 || iy==0 || ix==iy)) ||
                (tune_shift_with_amplitude_struct.lines_only &&
                 !(ix==0 || iy==0)))
              continue;
	    result = fabs(xTune[ix][iy] - xTune[0][0]);
	    if (result>maxResult)
	      maxResult = result;
	    result = fabs(yTune[ix][iy] - yTune[0][0]);
	    if (result>maxResult)
	      maxResult = result;
	  }
	}
	if (tune_shift_with_amplitude_struct.verbose) 
	  fprintf(stdout, "maximum tune change: %e\n", maxResult);
      }
      if (nLost || maxResult>tune_shift_with_amplitude_struct.scale_down_limit) {
	if (upperLimit[0]>tune_shift_with_amplitude_struct.x1)
	  upperLimit[0] = tune_shift_with_amplitude_struct.x1;
	if (upperLimit[1]>tune_shift_with_amplitude_struct.y1)
	  upperLimit[1] = tune_shift_with_amplitude_struct.y1;
	tune_shift_with_amplitude_struct.x1 = (upperLimit[0] + 2*lowerLimit[0])/3;
	tune_shift_with_amplitude_struct.y1 = (upperLimit[1] + 2*lowerLimit[1])/3;
	fprintf(stdout, "Warning: the amplitude you specified for tune shift with amplitude is too large.\n");
	fprintf(stdout, "Reducing tune_shift_with_amplitude_struct.x1=%le and tune_shift_with_amplitude_struct.y1=%le\n",
		tune_shift_with_amplitude_struct.x1, tune_shift_with_amplitude_struct.y1);
	if (tries==0) 
	  tries = 1;   /* ensures we don't exit on amplitude too large */
	continue;
      }
    }
    break;
  }

  if (nLost)
    for (ix=0; ix<N_TSWA; ix++)
      for (iy=0; iy<N_TSWA; iy++)
        dnux_dA[ix][iy] = dnuy_dA[ix][iy] = sqrt(DBL_MAX);
  else
    for (ix=0; ix<N_TSWA; ix++)
      for (iy=0; iy<N_TSWA; iy++)
	dnux_dA[ix][iy] = dnuy_dA[ix][iy] = 0;

  if (tune_shift_with_amplitude_struct.tune_output) {
    for (ix=j=0; ix<gridSize; ix++) {
      for (iy=0; iy<gridSize; iy++) {
        if ((tune_shift_with_amplitude_struct.sparse_grid &&
             !(ix==0 || iy==0 || ix==iy)) ||
               (tune_shift_with_amplitude_struct.lines_only &&
                !(ix==0 || iy==0)))
          continue;
	if (!SDDS_SetRowValues(&SDDSTswaTunes, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, j,
			       "x", x0[ix], "Ax", Ax[ix],
			       "nux", xTune[ix][iy],
			       "y", y0[iy], "Ay", Ay[iy],
			       "nuy", yTune[ix][iy],
			       NULL)) {
	  fprintf(stdout, "Problem filling SDDS page for TSWA tune output\n");
	  SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
	  exit(1);
	}
        j++;
      }
    }
    if (!SDDS_WritePage(&SDDSTswaTunes)) {
      fprintf(stdout, "Problem writing SDDS page for TSWA tune output\n");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }

  xTuneExtrema[0] = -(xTuneExtrema[1] = -DBL_MAX);
  yTuneExtrema[0] = -(yTuneExtrema[1] = -DBL_MAX);
  for (ix=0; ix<gridSize; ix++) {
    for (iy=0; iy<gridSize; iy++) {
      if ((tune_shift_with_amplitude_struct.sparse_grid &&
           !(ix==0 || iy==0 || ix==iy)) ||
          (tune_shift_with_amplitude_struct.lines_only &&
           !(ix==0 || iy==0)))
        continue;
      if (lost[ix][iy] ||
	  xTune[ix][iy]<=0.0 || xTune[ix][iy]>=1.0 ||
	  yTune[ix][iy]<=0.0 || yTune[ix][iy]>=1.0)
	continue;
      if (xTuneExtrema[0]>xTune[ix][iy])
	xTuneExtrema[0] = xTune[ix][iy];
      if (xTuneExtrema[1]<xTune[ix][iy])
	xTuneExtrema[1] = xTune[ix][iy];
      if (yTuneExtrema[0]>yTune[ix][iy])
	yTuneExtrema[0] = yTune[ix][iy];
      if (yTuneExtrema[1]<yTune[ix][iy])
	yTuneExtrema[1] = yTune[ix][iy];
    }
  }
  if (nLost && !tune_shift_with_amplitude_struct.exclude_lost_particles) {
    fprintf(stdout, "%ld particles in tune tracking: setting spreads to 1\n", nLost);
    xTuneExtrema[0] = yTuneExtrema[0] = 0;
    xTuneExtrema[1] = yTuneExtrema[1] = 1.0;
  }

  if (tune_shift_with_amplitude_struct.verbose) {
      fprintf(stdout, "xTune extrema: %21.15e, %21.15e, delta = %21.15e\n",
	      xTuneExtrema[0], xTuneExtrema[1], 
	      fabs(xTuneExtrema[0]-xTuneExtrema[1]));
      fprintf(stdout, "yTune extrema: %21.15e, %21.15e, delta = %21.15e\n",
	      yTuneExtrema[0], yTuneExtrema[1], 
	      fabs(yTuneExtrema[0]-yTuneExtrema[1]));
  }

  if (!nLost && !tune_shift_with_amplitude_struct.spread_only) {
    m = gridSize*gridSize;
    order = tune_shift_with_amplitude_struct.order<=1?1:2;
    if (tune_shift_with_amplitude_struct.lines_only) {
      /* Perform 1-D fits vs Ax and Ay. 
       * No cross terms get included. 
       */
      double *tuneData, *coefData;
      tuneData = tmalloc(sizeof(*tuneData)*gridSize);
      coefData = tmalloc(sizeof(*coefData)*3);

      /* x tune vs x amplitude */
      for (ix=iy=0; ix<gridSize; ix++)
        tuneData[ix] = xTune[ix][iy];
      lsfn(Ax, tuneData, NULL, gridSize, order, coefData, NULL, NULL, NULL);
      for (ix=0; ix<(order+1); ix++) {
        dnux_dA[ix][0] = coefData[ix]*factorial(ix);
      }
        
      /* x tune vs y amplitude */
      for (ix=iy=0; iy<gridSize; iy++)
        tuneData[iy] = xTune[ix][iy];
      lsfn(Ay, tuneData, NULL, gridSize, order, coefData, NULL, NULL, NULL);
      for (iy=0; iy<(order+1); iy++)
        dnux_dA[0][iy] = coefData[iy]*factorial(iy);

      /* y tune vs x amplitude */
      for (ix=iy=0; ix<gridSize; ix++) 
        tuneData[ix] = yTune[ix][iy];
      lsfn(Ax, tuneData, NULL, gridSize, order, coefData, NULL, NULL, NULL);
      for (ix=0; ix<(order+1); ix++)
        dnuy_dA[ix][0] = coefData[ix]*factorial(ix);
            
      /* y tune vs y amplitude */
      for (ix=iy=0; iy<gridSize; iy++)
        tuneData[iy] = yTune[ix][iy];
      lsfn(Ay, tuneData, NULL, gridSize, order, coefData, NULL, NULL, NULL);
      for (iy=0; iy<(order+1); iy++) 
        dnuy_dA[0][iy] = coefData[iy]*factorial(iy);

      if (tune_shift_with_amplitude_struct.verbose) {
        fprintf(stdout, "dnux/(dAx^i):  ");
        for (ix=0; ix<(order+1); ix++)
          fprintf(stdout, "%10.3g ", dnux_dA[ix][0]);
        fputc('\n', stdout);
        fprintf(stdout, "dnux/(dAy^i):  ");
        for (iy=0; iy<(order+1); iy++)
          fprintf(stdout, "%10.3g ", dnux_dA[0][iy]);
        fputc('\n', stdout);
        fprintf(stdout, "dnuy/(dAx^i):  ");
        for (ix=0; ix<(order+1); ix++)
          fprintf(stdout, "%10.3g ", dnuy_dA[ix][0]);
        fputc('\n', stdout);
        fprintf(stdout, "dnuy/(dAy^i):  ");
        for (iy=0; iy<(order+1); iy++)
          fprintf(stdout, "%10.3g ", dnuy_dA[0][iy]);
        fputc('\n', stdout);
      }
      free(tuneData);
      free(coefData);
    } else {
      /* Perform 2d fits vs Ax and Ay. 
       * May have trouble getting dnux/Ay=dnuy/dAx 
       */
      if (order==1) {
        /* the expansion is
           nu = Ax^0 (TS00 + Ay*TS01) +
           Ax^1 (TS10 ) 
           where TSij = 1/(i!j!) dnu/(dAx^i dAy^j)
           */
        nTSWA = 2;
        n = 4;
      }
      else {
        /* the expansion is
           nu = Ax^0 (TS00 + Ay*TS01 + Ay^2*TS02) +
           Ax^1 (TS10 + Ay*TS11 ) +
           Ax^2 (TS20 )
           where TSij = 1/(i!j!) dnu/(dAx^i dAy^j)
           */
        nTSWA = N_TSWA;
        n = (N_TSWA*N_TSWA + N_TSWA)/2;
      }
      /* Nu = AxAy*Coef 
       * Nu is (m-1)x1, AxAy is (m-1)xn, Coef is nx1 */
      m --;  /* won't use [0][0] from tune arrays due to NAFF issues */
      m_alloc(&AxAy, m, n);
      m_alloc(&Coef, n, 1);
      m_alloc(&Nu, m, 1);
      m_alloc(&AxAyTr, n, m);
      m_alloc(&Mf, n, n);
      m_alloc(&MfInv, n, n);
      m_alloc(&AxAyTrNu, n, 1);
      for (ix=i=0; ix<gridSize; ix++) {
        for (iy=0; iy<gridSize; iy++) {
          if ((tune_shift_with_amplitude_struct.sparse_grid &&
               !(ix==0 || iy==0 || ix==iy)) ||
              (tune_shift_with_amplitude_struct.lines_only &&
               !(ix==0 || iy==0)))
            continue;
          if (ix==0 && iy==0)
            continue;
          for (ix1=j=0; ix1<nTSWA; ix1++) {
            for (iy1=0; (ix1+iy1)<=order; iy1++, j++) {
              if (j>=n)
                bomb("indexing problem for TSWA calculation", NULL);
              AxAy->a[i][j] = ipow(Ax[ix], ix1)*ipow(Ay[iy], iy1);
            }
          }
          i++;
        }
      }
      m_trans(AxAyTr, AxAy);
      m_mult(Mf, AxAyTr, AxAy);
      m_invert(MfInv, Mf);
      for (ix=i=0; ix<gridSize; ix++)
        for (iy=0; iy<gridSize; iy++) {
          if ((tune_shift_with_amplitude_struct.sparse_grid &&
               !(ix==0 || iy==0 || ix==iy)) ||
              (tune_shift_with_amplitude_struct.lines_only &&
               !(ix==0 || iy==0)))
            continue;
          if (ix==0 && iy==0)
            continue;
          Nu->a[i][0] = xTune[ix][iy];
          i++;
        }
      m_mult(AxAyTrNu, AxAyTr, Nu);
      m_mult(Coef, MfInv, AxAyTrNu);
      for (ix=i=0; ix<nTSWA; ix++) 
        for (iy=0; (ix+iy)<=order; iy++, i++)
          dnux_dA[ix][iy] = Coef->a[i][0]*factorial(ix)*factorial(iy);
      if (tune_shift_with_amplitude_struct.verbose) {
        fprintf(stdout, "dnux/(dAx^i dAy^j):\n");
        for (ix=0; ix<nTSWA; ix++) {
          for (iy=0; iy<nTSWA; iy++) {
            fprintf(stdout, "%10.3g%c", dnux_dA[ix][iy], iy==(nTSWA-1)?'\n':' ');
          }
        }
      }
      
      for (ix=i=0; ix<gridSize; ix++)
        for (iy=0; iy<gridSize; iy++) {
          if ((tune_shift_with_amplitude_struct.sparse_grid &&
               !(ix==0 || iy==0 || ix==iy)) ||
              (tune_shift_with_amplitude_struct.lines_only &&
               !(ix==0 || iy==0)))
            continue;
          if (ix==0 && iy==0)
            continue;
          Nu->a[i][0] = yTune[ix][iy];
          i++;
        }
      m_mult(AxAyTrNu, AxAyTr, Nu);
      m_mult(Coef, MfInv, AxAyTrNu);
      for (ix=i=0; ix<nTSWA; ix++) 
        for (iy=0; (ix+iy)<=order; iy++, i++)
          dnuy_dA[ix][iy] = Coef->a[i][0]*factorial(ix)*factorial(iy);
      if (tune_shift_with_amplitude_struct.verbose) {
        fprintf(stdout, "dnuy/(dAx^i dAy^j):\n");
        for (ix=0; ix<nTSWA; ix++) {
          for (iy=0; iy<nTSWA; iy++) {
            fprintf(stdout, "%10.3g%c", dnuy_dA[ix][iy], iy==(nTSWA-1)?'\n':' ');
          }
        }
      }
      
      m_free(&AxAy);
      m_free(&Coef);
      m_free(&Nu);
      m_free(&AxAyTr);
      m_free(&Mf);
      m_free(&MfInv);
      m_free(&AxAyTrNu);
    }
  }  
}

long computeTunesFromTracking(double *tune, double *amp, VMATRIX *M, LINE_LIST *beamline, RUN *run,
			      double *startingCoord, 
			      double xAmplitude, double yAmplitude, long turns,
                              long useMatrix,
                              double *endingCoord, double *tuneLowerLimit, double *tuneUpperLimit,
			      long allowLosses)
{
  double **oneParticle, dummy;
  double frequency[4], amplitude[4], phase[4];
  double *x, *y, *xp, *yp, p;
  long i, one=1;
  double CSave[6];
  
#ifdef DEBUG1
  static FILE *fpdeb = NULL;
  if (!fpdeb) {
    fpdeb = fopen("tuneTracking.sdds", "w");
    fprintf(fpdeb, "SDDS1\n");
    fprintf(fpdeb, "&column name=Pass type=long &end\n");
    fprintf(fpdeb, "&column name=x type=double &end\n");
    fprintf(fpdeb, "&column name=xp type=double &end\n");
    fprintf(fpdeb, "&column name=y type=double &end\n");
    fprintf(fpdeb, "&column name=yp type=double &end\n");
    fprintf(fpdeb, "&column name=t type=double &end\n");
    fprintf(fpdeb, "&column name=p type=double &end\n");
    fprintf(fpdeb, "&data mode=ascii &end\n");
  }
  fprintf(fpdeb, "%ld\n", turns);

#endif

#ifdef DEBUG
  fprintf(stdout, "In computeTunesFromTracking: turns=%ld, xAmp=%le, yAmp=%le, useMatrix=%ld\n",
          turns, xAmplitude, yAmplitude, useMatrix);
  fflush(stdout);
#endif
  x = xp = y = yp = NULL;
  if (useMatrix) {
    /* this is necessary because the concatenated matrix includes the closed orbit in 
     * C.  We don't want to put this in at each turn.
     */
    for (i=0; i<6; i++) {
      CSave[i] = M->C[i];
      M->C[i] = 0;
    }
  }
  oneParticle = (double**)czarray_2d(sizeof(**oneParticle), 1, 7);
  if (!startingCoord)
    fill_double_array(oneParticle[0], 7, 0.0);
  else {
    memcpy(oneParticle[0], startingCoord, 6*sizeof(**oneParticle));
    oneParticle[0][6] = 0;
  }
  oneParticle[0][0] += xAmplitude;
  oneParticle[0][2] += yAmplitude;

#ifdef DEBUG
  fprintf(stdout, "Starting coordinates: %21.15le, %21.15le, %21.15le, %21.15le, %21.15le, %21.15le\n",
          oneParticle[0][0], oneParticle[0][1],
          oneParticle[0][2], oneParticle[0][3],
          oneParticle[0][4], oneParticle[0][5]);
  fflush(stdout);
#endif

#ifdef DEBUG
  fprintf(stdout, "Doing malloc: turns=%ld\n", turns);
  fflush(stdout);
#endif


  if (!(x = malloc(sizeof(*x)*turns)) || !(y = malloc(sizeof(*y)*turns)) ||
      !(xp = malloc(sizeof(*xp)*turns)) || !(yp = malloc(sizeof(*yp)*turns)))
    bomb("memory allocation failure (computeTunesFromTracking)", NULL);

#ifdef DEBUG
  fprintf(stdout, "Did malloc\n");
  fflush(stdout);
#endif

  x[0] = oneParticle[0][0];
  xp[0] = oneParticle[0][1];
  y[0] = oneParticle[0][2];
  yp[0] = oneParticle[0][3];
  p = run->p_central;


#ifdef DEBUG
  fprintf(stdout, "Starting to track\n");
  fflush(stdout);
#endif
#ifdef DEBUG1
  fprintf(fpdeb, "0 %21.15e %21.15e %21.15e %21.15e %21.15e %21.15e\n",
            oneParticle[0][0], oneParticle[0][1], oneParticle[0][2], oneParticle[0][3], oneParticle[0][4], oneParticle[0][5]);
#endif

  for (i=1; i<turns; i++) {
    if (useMatrix)
      track_particles(oneParticle, M, oneParticle, one);
    else {
      if (!do_tracking(NULL, oneParticle, 1, NULL, beamline, &p,  (double**)NULL, 
		       (BEAM_SUMS**)NULL, (long*)NULL,
                       (TRAJECTORY*)NULL, run, 0, 
		       TEST_PARTICLES+(allowLosses?TEST_PARTICLE_LOSSES:0)+TIME_DEPENDENCE_OFF, 
                       1, i-1, NULL, NULL, NULL, NULL, NULL)) {
	if (!allowLosses) {
	  fprintf(stdout, "warning: test particle lost on turn %ld (computeTunesFromTracking)\n", i);
	  fflush(stdout);
	}
        return 0;
      }
    }
    if (isnan(oneParticle[0][0]) || isnan(oneParticle[0][1]) ||
        isnan(oneParticle[0][2]) || isnan(oneParticle[0][3]) ||
        isnan(oneParticle[0][4]) || isnan(oneParticle[0][5])) {
      if (!allowLosses) {
	fprintf(stdout, "warning: test particle lost on turn %ld (computeTunesFromTracking)\n", i);
	fflush(stdout);
      }
      return 0;
    }
    x[i] = oneParticle[0][0];
    xp[i] = oneParticle[0][1];
    y[i] = oneParticle[0][2];
    yp[i] = oneParticle[0][3];
#ifdef DEBUG1
    fprintf(fpdeb, "%ld %21.15e %21.15e %21.15e %21.15e %21.15e %21.15e\n", i, 
            oneParticle[0][0], oneParticle[0][1], oneParticle[0][2], oneParticle[0][3], oneParticle[0][4], oneParticle[0][5]);
#endif
  }
  if (endingCoord) {
    for (i=0; i<6; i++)
      endingCoord[i] = oneParticle[0][i];
  }
#ifdef DEBUG
  fprintf(stdout, "Ending coordinates: %21.15le, %21.15le, %21.15le, %21.15le, %21.15le, %21.15le\n",
          oneParticle[0][0], oneParticle[0][1],
          oneParticle[0][2], oneParticle[0][3],
          oneParticle[0][4], oneParticle[0][5]);
  fflush(stdout);
#endif


#ifdef DEBUG
  fprintf(stdout, "Performing NAFF (1)\n");
  fflush(stdout);
#endif
  if (PerformNAFF(&frequency[0], &amplitude[0], &phase[0], 
		  &dummy, 0.0, 1.0, x, turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,
		  tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
		  tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (x).\n");
    fprintf(stdout, "Limits: %e, %e\n",
	    tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
	    tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0);
    return 0;
  }
#ifdef DEBUG
  fprintf(stdout, "Performing NAFF (2)\n");
  fflush(stdout);
#endif
  if (PerformNAFF(&frequency[1], &amplitude[1], &phase[1], 
		  &dummy, 0.0, 1.0, xp, turns, 
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,
		  tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
		  tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (xp).\n");
    fprintf(stdout, "Limits: %e, %e\n",
	    tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
	    tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0);
    return 0;
  }
#ifdef DEBUG
  fprintf(stdout, "Performing NAFF (3)\n");
  fflush(stdout);
#endif
  if (PerformNAFF(&frequency[2], &amplitude[2], &phase[2], 
		  &dummy, 0.0, 1.0, y, turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,
		  tuneLowerLimit?(tuneLowerLimit[1]>0.5 ? 1-tuneLowerLimit[1] : tuneLowerLimit[1]):0,
		  tuneUpperLimit?(tuneUpperLimit[1]>0.5 ? 1-tuneUpperLimit[1] : tuneUpperLimit[1]):0)!=1) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (y).\n");
    fprintf(stdout, "Limits: %e, %e\n",
	    tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
	    tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0);
    return 0;
  }
#ifdef DEBUG
  fprintf(stdout, "Performing NAFF (4)\n");
  fflush(stdout);
#endif
  if (PerformNAFF(&frequency[3], &amplitude[3], &phase[3], 
		  &dummy, 0.0, 1.0, yp, turns,
		  NAFF_MAX_FREQUENCIES|NAFF_FREQ_CYCLE_LIMIT|NAFF_FREQ_ACCURACY_LIMIT,
		  0.0, 1, 200, 1e-12,
		  tuneLowerLimit?(tuneLowerLimit[1]>0.5 ? 1-tuneLowerLimit[1] : tuneLowerLimit[1]):0,
		  tuneUpperLimit?(tuneUpperLimit[1]>0.5 ? 1-tuneUpperLimit[1] : tuneUpperLimit[1]):0)!=1 ) {
    fprintf(stdout, "Warning: NAFF failed for tune analysis from tracking (yp).\n");
    fprintf(stdout, "Limits: %e, %e\n",
	    tuneLowerLimit?(tuneLowerLimit[0]>0.5 ? 1-tuneLowerLimit[0] : tuneLowerLimit[0]):0,
	    tuneUpperLimit?(tuneUpperLimit[0]>0.5 ? 1-tuneUpperLimit[0] : tuneUpperLimit[0]):0);
    return 0;
  }

#ifdef DEBUG
  fprintf(stdout, "NAFF done\n");
  for (i=0; i<4; i++)
    fprintf(stdout, "%ld: freq=%e, phase=%e, amp=%e\n",
	    i, frequency[i], phase[i], amplitude[i]);
  fflush(stdout);
#endif

  tune[0] = adjustTuneHalfPlane(frequency[0], phase[0], phase[1]);
  tune[1] = adjustTuneHalfPlane(frequency[2], phase[2], phase[3]);
  if (amp) {
    amp[0] = amplitude[0];
    amp[1] = amplitude[2];
  }
#ifdef DEBUG
  fprintf(stdout, "xtune = %e, ytune = %e\n", tune[0], tune[1]);
#endif

  free(x);
  free(y);
  free(xp);
  free(yp);
  free_czarray_2d((void*)oneParticle, 1, 7);
  if (useMatrix) {
    M->C[i] = CSave[i];
  }
  return 1;
}

double adjustTuneHalfPlane(double frequency, double phase0, double phase1)
{
  if (fabs(phase0-phase1)>PI) {
    if (phase0<phase1)
      phase0 += PIx2;
    else 
      phase1 += PIx2;
  }
  if (phase0<phase1)
    return frequency;
  return 1-frequency;
}

void computeTuneShiftWithAmplitudeM(double dnux_dA[N_TSWA][N_TSWA], double dnuy_dA[N_TSWA][N_TSWA],
                                    TWISS *twiss, double *tune, VMATRIX *M)
{
  VMATRIX M1, M2;
  long tplane, splane;
  long it, jt, is, js, ia;
  double alpha[2], beta[2], shift[2];
  double turns;  
  double C, S, theta;

  initialize_matrices(&M1, 3);
  initialize_matrices(&M2, 3);
  copy_matrices(&M2, M);
  /* the matrix M is already formed for any closed orbit, so we don't want 
   * to concatenate the closed orbit
   */
  for (it=0; it<6; it++)
    M2.C[it] = 0;
  concat_matrices(&M1, &M2, &M2, 0);
  turns = 2;

  while (turns<32768) {
    copy_matrices(&M2, &M1);
    concat_matrices(&M1, &M2, &M2, 0);
    turns *= 2;
  }

  free_matrices(&M2);

/*
  fprintf(stderr, "Tune check: \n%le %le\n%le %le\n",
          cos(PIx2*turns*tune[0]), (M1.R[0][0]+M1.R[1][1])/2,
          cos(PIx2*turns*tune[1]), (M1.R[2][2]+M1.R[3][3])/2);
*/
  
  alpha[0] = twiss->alphax;
  alpha[1] = twiss->alphay;
  beta[0] = twiss->betax;
  beta[1] = twiss->betay;

  for (it=0; it<N_TSWA; it++)
    for (is=0; is<N_TSWA; is++)
      dnux_dA[it][is] = dnuy_dA[it][is] = 0;
  dnux_dA[0][0] = tune[0];
  dnuy_dA[0][0] = tune[1];
  
  for (tplane=0; tplane<2; tplane++) {
    it = 2*tplane;
    jt = it+1;
    for (splane=0; splane<2; splane++) {
      is = 2*splane;
      js = is+1;
      shift[splane] = 0;
      /*
      fprintf(stderr, "tplane=%ld splane=%ld\n", tplane, splane);
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              it, it, is, is, QElement(M1.Q, it, it, is, is));
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              jt, jt, is, is, QElement(M1.Q, jt, jt, is, is));
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              it, it, is, js, QElement(M1.Q, it, it, is, js));
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              jt, jt, is, js, QElement(M1.Q, jt, jt, is, js));
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              it, it, js, js, QElement(M1.Q, it, it, js, js));
      fprintf(stderr, "Q%ld%ld%ld%ld = %e\n",
              jt, jt, js, js, QElement(M1.Q, jt, jt, js, js));
      fprintf(stderr, "\n");
      */    
      for (ia=0, theta=0; ia<72; ia++, theta+=PIx2/72) {
        C = cos(theta);
        S = sin(theta)*alpha[splane] + C;
        shift[splane] += 
          (QElement(M1.Q, it, it, is, is) + QElement(M1.Q, jt, jt, is, is))*beta[splane]*sqr(C) 
            - (QElement(M1.Q, it, it, is, js) + QElement(M1.Q, jt, jt, is, js))*C*S 
             + (QElement(M1.Q, it, it, js, js) + QElement(M1.Q, jt, jt, js, js))*sqr(S)/beta[splane];
      }
      shift[splane] /= 72;
    }
    if (tplane==0) {
      /* nux */
      dnux_dA[1][0] = shift[0]/(-2*PIx2*sin(PIx2*tune[0]*turns)*turns);
      dnux_dA[0][1] = shift[1]/(-2*PIx2*sin(PIx2*tune[0]*turns)*turns);
    }
    else {
      /* nuy */
      dnuy_dA[1][0] = shift[0]/(-2*PIx2*sin(PIx2*tune[1]*turns)*turns);
      dnuy_dA[0][1] = shift[1]/(-2*PIx2*sin(PIx2*tune[1]*turns)*turns);
    }
  }

  free_matrices(&M1);
}

void store_fitpoint_twiss_parameters(MARK *fpt, char *name, long occurence,TWISS *twiss)
{
  long i;
  static char *twiss_name_suffix[12] = {
    "betax", "alphax", "nux", "etax", "etapx", "etaxp",
    "betay", "alphay", "nuy", "etay", "etapy", "etayp",
    } ;
  static char s[200];
  if (!(fpt->init_flags&1)) {
    fpt->twiss_mem = tmalloc(sizeof(*(fpt->twiss_mem))*12);
    fpt->init_flags |= 1;
    for (i=0; i<12; i++) {
      sprintf(s, "%s#%ld.%s", name, occurence, twiss_name_suffix[i]);
      fpt->twiss_mem[i] = rpn_create_mem(s, 0);
    }
  }
  if (!twiss) {
    fprintf(stdout, "twiss parameter pointer unexpectedly NULL\n");
    fflush(stdout);
    abort();
  }
  for (i=0; i<5; i++) {
    rpn_store(*((&twiss->betax)+i)/(i==2?PIx2:1), NULL, fpt->twiss_mem[i]);
    rpn_store(*((&twiss->betay)+i)/(i==2?PIx2:1), NULL, fpt->twiss_mem[i+6]);
  }
  /* store etaxp and etayp in under two names each: etapx and etaxp */
  i = 4;
  rpn_store(*((&twiss->betax)+i), NULL, fpt->twiss_mem[i+1]);
  rpn_store(*((&twiss->betay)+i), NULL, fpt->twiss_mem[i+7]);
}


void clearTwissAnalysisRequests() 
{
  long i;
  for (i=0; i<twissAnalysisRequests; i++) {
    if (twissAnalysisRequest[i].startName)
      free(twissAnalysisRequest[i].startName);
    if (twissAnalysisRequest[i].endName)
      free(twissAnalysisRequest[i].endName);
    if (twissAnalysisRequest[i].matchName)
      free(twissAnalysisRequest[i].matchName);
    if (twissAnalysisRequest[i].tag)
      free(twissAnalysisRequest[i].tag);
  }
  free(twissAnalysisRequest);
  twissAnalysisRequests = 0;
}

void addTwissAnalysisRequest(char *tag, char *startName, char *endName, 
                             char *matchName,
                             long startOccurence, long endOccurence,
                             double sStart, double sEnd)
{
  long i;
  if (!tag || !strlen(tag))
    bomb("NULL or blank tag passed to addTwissAnalysisRequest", NULL);
  if (!(startName && strlen(startName) && endName && strlen(endName)) 
      && !(matchName && strlen(matchName)) && sStart==sEnd)
    bomb("must have both startName and endName, or matchName, or sStart!=sEnd (addTwissAnalysisRequest)", NULL);
  if ((startName || endName) && matchName)
    bomb("can't have startName or endName with matchName (addTwissAnalysisRequest)", NULL);
  for (i=0; i<twissAnalysisRequests; i++)
    if (strcmp(twissAnalysisRequest[i].tag, tag)==0)
      bomb("duplicate tag names seen (addTwissAnalysisRequest)", NULL);
  if (!(twissAnalysisRequest = 
        SDDS_Realloc(twissAnalysisRequest, sizeof(*twissAnalysisRequest)*(twissAnalysisRequests+1))) ||
      !SDDS_CopyString(&twissAnalysisRequest[twissAnalysisRequests].tag, tag))
    bomb("memory allocation failure (addTwissAnalysisRequest)", NULL);
  twissAnalysisRequest[twissAnalysisRequests].startName = 
    twissAnalysisRequest[twissAnalysisRequests].endName = 
      twissAnalysisRequest[twissAnalysisRequests].matchName = NULL;
  if ((startName &&
       !SDDS_CopyString(&twissAnalysisRequest[twissAnalysisRequests].startName, startName)) ||
      (endName &&
       !SDDS_CopyString(&twissAnalysisRequest[twissAnalysisRequests].endName, endName)) ||
      (matchName &&
       !SDDS_CopyString(&twissAnalysisRequest[twissAnalysisRequests].matchName, matchName)))
    bomb("memory allocation failure (addTwissAnalysisRequest)", NULL);
  if (matchName && 
      has_wildcards(twissAnalysisRequest[twissAnalysisRequests].matchName) &&
      strchr(twissAnalysisRequest[twissAnalysisRequests].matchName, '-'))
    twissAnalysisRequest[twissAnalysisRequests].matchName = 
      expand_ranges(twissAnalysisRequest[twissAnalysisRequests].matchName);
  twissAnalysisRequest[twissAnalysisRequests].sStart = sStart;
  twissAnalysisRequest[twissAnalysisRequests].sEnd = sEnd;
  twissAnalysisRequest[twissAnalysisRequests].startOccurence = startOccurence;
  twissAnalysisRequest[twissAnalysisRequests].endOccurence = endOccurence;
  twissAnalysisRequest[twissAnalysisRequests].initialized = 0;
  twissAnalysisRequests++;
}

void processTwissAnalysisRequests(ELEMENT_LIST *elem)
{
  long i, is, iq, count, firstTime;
  char buffer[1024];
  ELEMENT_LIST *elemOrig;
  double value, end_pos, start_pos;
  double twissData[TWISS_ANALYSIS_STATS][TWISS_ANALYSIS_QUANTITIES];

  elemOrig = elem;
  start_pos = 0;

  for (i=0; i<twissAnalysisRequests; i++) {
    /* initialize statistics buffers and rpn memories */
    firstTime = !twissAnalysisRequest[i].initialized;
    for (iq=0; iq<TWISS_ANALYSIS_QUANTITIES; iq++)  {
      for (is=0; is<TWISS_ANALYSIS_STATS; is++)
        if (!twissAnalysisRequest[i].initialized) {
          sprintf(buffer, "%s.%s.%s", twissAnalysisRequest[i].tag,
                  twissAnalysisStatName[is], twissAnalysisQuantityName[iq]);
          twissAnalysisRequest[i].twissMem[is][iq] = rpn_create_mem(buffer, 0);
        }
      twissData[TWISS_ANALYSIS_AVE][iq] = 0;
      twissData[TWISS_ANALYSIS_MIN][iq] = DBL_MAX;
      twissData[TWISS_ANALYSIS_MAX][iq] = -DBL_MAX;
    }
    twissAnalysisRequest[i].initialized = 1;
    
    count = end_pos = 0;
    while (elem) {
      if (twissAnalysisRequest[i].sStart<twissAnalysisRequest[i].sEnd &&
              elem->end_pos>twissAnalysisRequest[i].sEnd) 
        break;
      if (twissAnalysisRequest[i].sStart<twissAnalysisRequest[i].sEnd &&
          elem->end_pos<twissAnalysisRequest[i].sStart) {
        elem = elem->succ;
        continue;
      }
      if (twissAnalysisRequest[i].matchName || twissAnalysisRequest[i].startName) {
        if (twissAnalysisRequest[i].matchName) {
          if (!wild_match(elem->name, twissAnalysisRequest[i].matchName)) {
            elem = elem->succ;
            continue;
          } 
        } else {
          if (!count &&
              !(twissAnalysisRequest[i].startName && 
                strcmp(twissAnalysisRequest[i].startName, elem->name)==0 &&
                (twissAnalysisRequest[i].startOccurence<=1 ||
                 twissAnalysisRequest[i].startOccurence==elem->occurence))) {
            elem = elem->succ;
            continue;
          }
        }
      }
      count++;
      if (count==1) {
        if (elem->pred)
          start_pos = elem->pred->end_pos;
        else
          start_pos = 0;
      }
      if (twiss_analysis_struct.verbosity>1 && firstTime) {
        fprintf(stdout, "twiss analysis %s will include %s#%ld\n",
                twissAnalysisRequest[i].tag,
                elem->name, elem->occurence);
      }
      for (iq=0; iq<TWISS_ANALYSIS_QUANTITIES; iq++)  {
        value = *(double*)((char*)elem->twiss+twissAnalysisQuantityOffset[iq]);
        for (is=0; is<TWISS_ANALYSIS_STATS; is++) {
          switch (twissAnalysisStatCode[is]) {
          case TWISS_ANALYSIS_AVE:
            twissData[is][iq] += value;
            break;
          case TWISS_ANALYSIS_MIN:
            if (twissData[is][iq]>value)
              twissData[is][iq] = value;
            break;
          case TWISS_ANALYSIS_MAX:
            if (twissData[is][iq]<value)
              twissData[is][iq] = value;
            break;
          }
        }
      }
      if (twissAnalysisRequest[i].endName && 
           strcmp(twissAnalysisRequest[i].endName, elem->name)==0 &&
          (twissAnalysisRequest[i].endOccurence<=1 ||
           twissAnalysisRequest[i].endOccurence==elem->occurence))
        break;
      elem = elem->succ;
    }
    if (!count) {
      fprintf(stdout, "error: twiss analysis conditions never satisfied for request with tag %s\n",
              twissAnalysisRequest[i].tag);
      exit(1);
    }
    for (iq=0; iq<TWISS_ANALYSIS_QUANTITIES; iq++)  {
      twissData[TWISS_ANALYSIS_AVE][iq] /= count;
      for (is=0; is<TWISS_ANALYSIS_STATS; is++) {
        rpn_store(twissData[is][iq], NULL, twissAnalysisRequest[i].twissMem[is][iq]);
      }
    }
    if (twiss_analysis_struct.verbosity && firstTime) {
      fprintf(stdout, "%ld elements included in computations for twiss analysis request with tag %s\n",
              count, 
              twissAnalysisRequest[i].tag);
    }
    
#if DEBUG
    if (twissAnalysisRequest[i].matchName) {
      fprintf(stdout, "%ld matches for %s\n",
              count, twissAnalysisRequest[i].tag
              );
      for (iq=0; iq<TWISS_ANALYSIS_QUANTITIES; iq++) {
        fprintf(stdout, "%s: ", twissAnalysisQuantityName[iq]);
        for (is=0; is<TWISS_ANALYSIS_STATS; is++) 
          fprintf(stdout, "%s=%e%s",
                  twissAnalysisStatName[is],
                  twissData[is][iq],
                  (is==TWISS_ANALYSIS_STATS-1?"\n":", "));
      }
      fflush(stdout);
    }
#endif
    elem = elemOrig;
  }
}

void setupTwissAnalysisRequest(NAMELIST_TEXT *nltext, RUN *run, 
                               LINE_LIST *beamline)
{
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&twiss_analysis, nltext);
  if (echoNamelists) print_namelist(stdout, &twiss_analysis);

  if (twiss_analysis_struct.clear) {
    clearTwissAnalysisRequests();
    if (!(twiss_analysis_struct.start_name && twiss_analysis_struct.end_name) &&
        !twiss_analysis_struct.match_name &&
        twiss_analysis_struct.s_start==twiss_analysis_struct.s_end)
      return;
  }
  
  if (twiss_analysis_struct.start_name &&
      !strlen(trim_spaces(str_toupper(twiss_analysis_struct.start_name))))
    bomb("start_name is blank", NULL);
  if (twiss_analysis_struct.end_name &&
      !strlen(trim_spaces(str_toupper(twiss_analysis_struct.end_name))))
    bomb("end_name is blank", NULL);
  if (twiss_analysis_struct.match_name &&
      !strlen(trim_spaces(str_toupper(twiss_analysis_struct.match_name))))
    bomb("match_name is blank", NULL);
  if ((twiss_analysis_struct.tag &&
       !strlen(trim_spaces(twiss_analysis_struct.tag))) ||
      !twiss_analysis_struct.tag)
    bomb("tag is blank", NULL);
  
  if (!(twiss_analysis_struct.start_name && twiss_analysis_struct.end_name)
        && !twiss_analysis_struct.match_name &&
      twiss_analysis_struct.s_start==twiss_analysis_struct.s_end)
    bomb("you must give start_name and end_name, or match_name, or s_start different from s_end", NULL);
  if (twiss_analysis_struct.s_start>twiss_analysis_struct.s_end)
    bomb("s_start>s_end", NULL);
  addTwissAnalysisRequest(twiss_analysis_struct.tag,
                          twiss_analysis_struct.start_name,
                          twiss_analysis_struct.end_name,
                          twiss_analysis_struct.match_name,
                          twiss_analysis_struct.start_occurence,
                          twiss_analysis_struct.end_occurence,
                          twiss_analysis_struct.s_start,
                          twiss_analysis_struct.s_end);
}

void AddWigglerRadiationIntegrals(double length, long poles, double radius,
                                   double eta, double etap, 
                                   double beta, double alpha,
                                   double *I1, double *I2, double *I3, double *I4, double *I5)
{
  double h0, gamma;
  double Lp;
  long pole, fieldSign;
#ifdef DEBUG
  FILE *fpd = NULL;
#endif

  if (poles<2)
    bomb("wiggler must have at least 2 poles", NULL);
  if (radius<=0)
    bomb("wiggler must have positive, nonzero radius", NULL);

  /* length of each pole */
  Lp = length/poles;

#ifdef DEBUG
  fpd = fopen_e("wiggler.sdds", "w", 0);
  fprintf(fpd, "SDDS1\n&column name=Pole type=long &end\n");
  fprintf(fpd, "&column name=beta type=double units=m &end\n");
  fprintf(fpd, "&column name=alpha type=double &end\n");
  fprintf(fpd, "&column name=eta type=double units=m &end\n");
  fprintf(fpd, "&column name=etap type=double units=m &end\n");
  fprintf(fpd, "&column name=h0 type=double units=1/m &end\n");
  fprintf(fpd, "&column name=I1 type=double units=m &end\n");
  fprintf(fpd, "&column name=I2 type=double units=1/m &end\n");
  fprintf(fpd, "&column name=I3 type=double units=1/m$a2$n &end\n");
  fprintf(fpd, "&column name=I4 type=double units=1/m &end\n");
  fprintf(fpd, "&column name=I5 type=double units=1/m &end\n");
  fprintf(fpd, "&data mode=ascii no_row_counts=1 &end\n");
  fprintf(fpd, "0 %e %e %e %e 0 0 0 0 0 0\n", beta, alpha, eta, etap);
#endif
  
  gamma = (1+alpha*alpha)/beta;
  if (poles%2) {
    /* Odd number of poles: use half-strength end-poles to match */
    /* Integrate a half period at a time */
    fieldSign = 1;
    for (pole=0; pole<poles; pole++) {
      fieldSign *= -1;
      if (pole==0 || pole==poles-1) {
	h0 = fieldSign*0.5/radius;
      } else
	h0 = fieldSign/radius;
      
      *I1 += (h0*Lp*(h0*ipow(Lp,2) + 4*eta*PI + 2*etap*Lp*PI))/(2.*ipow(PI,2));
      *I2 += (ipow(h0,2)*Lp)/2.;
      *I3 += SIGN(h0)*(4*ipow(h0,3)*Lp)/(3.*PI);
      
      *I5 += SIGN(h0)*
	(ipow(h0,3)*Lp*(gamma*
			(-50625*eta*h0*ipow(Lp,2)*ipow(PI,3) +
			 72000*ipow(eta,2)*ipow(PI,4) +
			 32*ipow(h0,2)*ipow(Lp,4)*(1664 + 225*ipow(PI,2))) +
			225*ipow(PI,2)*
			(alpha*(-289*ipow(h0,2)*ipow(Lp,3) +
				5*h0*Lp*(128*eta - 45*etap*Lp)*PI + 640*eta*etap*ipow(PI,2)
				) + 64*beta*(6*ipow(h0,2)*ipow(Lp,2) + 10*etap*h0*Lp*PI +
					     5*ipow(etap,2)*ipow(PI,2)))))/(54000.*ipow(PI,5));
      
      beta  = beta - 2*Lp*alpha + sqr(Lp)*gamma;
      alpha = alpha - Lp*gamma;
      gamma = (1+alpha*alpha)/beta;
      eta   = eta + (etap + Lp/PI*h0)*Lp ;
      etap  = etap + 2*Lp/PI*h0;
      
#ifdef DEBUG
      fprintf(fpd, "%ld %e %e %e %e %e %e %e %e %e %e\n", pole+1, beta, alpha, eta, etap,
	      h0, *I1, *I2, *I3, *I4, *I5);
#endif
    }
    
  } else {
    /* Even number of poles: use half-length end-poles to match 
     * Integrate a full period at a time (starts and ends in the
     * middle of a pole).
     */
    double L;
    L = 2*Lp;
    h0 = 1./radius;
    for (pole=0; pole<poles; pole+=2) {
      *I5 += (ipow(h0,3)*Lp*(gamma*
                              (9000*ipow(eta,2)*ipow(PI,4) +
                               1125*eta*h0*ipow(Lp,2)*ipow(PI,2)*(16 + 3*PI) +
                               ipow(h0,2)*ipow(Lp,4)*(15656 + 2235*PI + 2250*ipow(PI,2))) +
                              225*ipow(PI,2)*(8*beta*
                                               (ipow(h0,2)*ipow(Lp,2) + 5*ipow(etap,2)*ipow(PI,2)) +
                                               alpha*(-16*ipow(h0,2)*ipow(Lp,3) +
                                                       5*etap*(16*eta*ipow(PI,2) + h0*ipow(Lp,2)*(16 + 3*PI))))))/
                                                         (3375.*ipow(PI,5));
      beta  = beta - 2*L*alpha + sqr(L)*gamma;
      alpha = alpha - L*gamma;
      gamma = (1+alpha*alpha)/beta;
      eta   = eta + 2*Lp*etap;
    }
    *I1 += -(poles/2)*((ipow(h0,2)*ipow(Lp,3))/ipow(PI,2));
    *I2 += (poles/2)*ipow(h0,2)*Lp;
    *I3 += (poles/2)*(8*ipow(h0,3)*Lp)/(3.*PI);
  }

#ifdef DEBUG
    fclose(fpd);
#endif
}


void setupLinearChromaticTracking(NAMELIST_TEXT *nltext, LINE_LIST *beamline)
{
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&setup_linear_chromatic_tracking, nltext);
  if (echoNamelists) print_namelist(stdout, &setup_linear_chromatic_tracking);

  if (setup_linear_chromatic_tracking_struct.nux[0]<0)
    bomb("nux < 0", NULL);
  if (setup_linear_chromatic_tracking_struct.nuy[0]<0)
    bomb("nuy < 0", NULL);
  
  beamline->twiss0 = malloc(sizeof(TWISS));

  setLinearChromaticTrackingValues(beamline);
  linearChromaticTrackingInitialized = 1;
}


void setLinearChromaticTrackingValues(LINE_LIST *beamline) 
{
  beamline->tune[0] = setup_linear_chromatic_tracking_struct.nux[0];
  beamline->chromaticity[0] = setup_linear_chromatic_tracking_struct.nux[1];
  beamline->chrom2[0] = setup_linear_chromatic_tracking_struct.nux[2];
  beamline->chrom3[0] = setup_linear_chromatic_tracking_struct.nux[3];

  beamline->twiss0->betax = setup_linear_chromatic_tracking_struct.betax[0];
  beamline->dbeta_dPoP[0] = setup_linear_chromatic_tracking_struct.betax[1];

  beamline->twiss0->alphax = setup_linear_chromatic_tracking_struct.alphax[0];
  beamline->dalpha_dPoP[0] = setup_linear_chromatic_tracking_struct.alphax[1];
  
  beamline->twiss0->etax = setup_linear_chromatic_tracking_struct.etax[0];
  beamline->eta2[0] = setup_linear_chromatic_tracking_struct.etax[1];
  beamline->eta3[0] = 0.0;
  beamline->twiss0->etapx = setup_linear_chromatic_tracking_struct.etapx[0];
  beamline->eta2[1] = setup_linear_chromatic_tracking_struct.etapx[1];
  beamline->eta3[1] = 0.0;

  beamline->tune[1] = setup_linear_chromatic_tracking_struct.nuy[0];
  beamline->chromaticity[1] = setup_linear_chromatic_tracking_struct.nuy[1];
  beamline->chrom2[1] = setup_linear_chromatic_tracking_struct.nuy[2];
  beamline->chrom3[1] = setup_linear_chromatic_tracking_struct.nuy[3];
  
  beamline->twiss0->betay = setup_linear_chromatic_tracking_struct.betay[0];
  beamline->dbeta_dPoP[1] = setup_linear_chromatic_tracking_struct.betay[1];
  
  beamline->twiss0->alphay = setup_linear_chromatic_tracking_struct.alphay[0];
  beamline->dalpha_dPoP[1] = setup_linear_chromatic_tracking_struct.alphay[1];
  
  beamline->twiss0->etay = setup_linear_chromatic_tracking_struct.etay[0];
  beamline->eta2[2] = setup_linear_chromatic_tracking_struct.etay[1];
  beamline->eta3[2] = 0.0;
  beamline->twiss0->etapy = setup_linear_chromatic_tracking_struct.etapy[0];
  beamline->eta2[3] = setup_linear_chromatic_tracking_struct.etapy[1];
  beamline->eta3[3] = 0.0;

  beamline->alpha[0] = setup_linear_chromatic_tracking_struct.alphac[0];
  beamline->alpha[1] = setup_linear_chromatic_tracking_struct.alphac[1];
}


void computeDrivingTerms(DRIVING_TERMS *d, ELEMENT_LIST *elem, TWISS *twiss0, double *tune)
/* Based on J. Bengtsson, SLS Note 9/97, March 7, 1997, with corrections per W. Guo (NSLS) */
{
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  double complex h11001, h00111, h20001, h00201, h10002;
  double complex h21000, h30000, h10110, h10020, h10200;
  double complex h22000, h11110, h00220, h31000, h40000, h12000;
  double complex h20110, h11200, h20020, h20200, h00310, h00400;
  double complex t1, t2, t3, t4;
#else
  doublecomplex_sdds h11001, h00111, h20001, h00201, h10002;
  doublecomplex_sdds h21000, h30000, h10110, h10020, h10200;
  doublecomplex_sdds h22000, h11110, h00220, h31000, h40000, h12000;
  doublecomplex_sdds h20110, h11200, h20020, h20200, h00310, h00400;
  doublecomplex_sdds t1, t2, t3, t4;
#endif
  double betax1, betay1, phix1, phiy1, etax1;
  double betax2, betay2, phix2, phiy2;
  double coef, b2L, b3L1, b3L2, sqrt_betax, sqrt3_betax, nux, nuy;
  double dnudA[2][2];
  ELEMENT_LIST *eptr1, *eptr2;
  
  /* accumulate real and imaginary parts */
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h11001 = h00111 = h20001 = h00201 = h10002 = 0;
  h21000 = h30000 = h10110 = h10020 = h10200 = 0;
#else
  h11001.r = h00111.r = h20001.r = h00201.r = h10002.r = 0;
  h21000.r = h30000.r = h10110.r = h10020.r = h10200.r = 0;
  h11001.i = h00111.i = h20001.i = h00201.i = h10002.i = 0;
  h21000.i = h30000.i = h10110.i = h10020.i = h10200.i = 0;
#endif

  eptr1 = elem;
  while (eptr1) {
    b2L = b3L1 = 0;
    switch (eptr1->type) {
    case T_SEXT:
      b3L1 = ((SEXT*)eptr1->p_elem)->k2 * ((SEXT*)eptr1->p_elem)->length/2;
      break;
    case T_KSEXT:
      b3L1 = ((KSEXT*)eptr1->p_elem)->k2 * ((KSEXT*)eptr1->p_elem)->length/2;
      break;
    case T_KQUSE:
      b2L = ((KQUSE*)eptr1->p_elem)->k1 * ((KQUSE*)eptr1->p_elem)->length;
      b3L1 = ((KQUSE*)eptr1->p_elem)->k2 * ((KQUSE*)eptr1->p_elem)->length/2;
      break;
    case T_SBEN:
    case T_RBEN:
      b2L = ((BEND*)eptr1->p_elem)->k1 * ((BEND*)eptr1->p_elem)->length;
      b3L1 = ((BEND*)eptr1->p_elem)->k2 * ((BEND*)eptr1->p_elem)->length/2;
      break;
    case T_CSBEND:
      b2L = ((CSBEND*)eptr1->p_elem)->k1 * ((CSBEND*)eptr1->p_elem)->length;
      b3L1 = ((CSBEND*)eptr1->p_elem)->k2 * ((CSBEND*)eptr1->p_elem)->length/2;
      break;
    case T_CSRCSBEND:
      b2L = ((CSRCSBEND*)eptr1->p_elem)->k1 * ((CSRCSBEND*)eptr1->p_elem)->length;
      b3L1 = ((CSRCSBEND*)eptr1->p_elem)->k2 * ((CSRCSBEND*)eptr1->p_elem)->length/2;
      break;
    case T_QUAD:
      b2L = ((QUAD*)eptr1->p_elem)->k1 * ((QUAD*)eptr1->p_elem)->length;
      break;
    case T_KQUAD:
      b2L = ((KQUAD*)eptr1->p_elem)->k1 * ((KQUAD*)eptr1->p_elem)->length;
      break;
    default:
      break;
    }      
    if (b2L || b3L1) {
      if (eptr1->pred) {
        betax1 = (eptr1->twiss->betax + eptr1->pred->twiss->betax)/2;
        etax1  = (eptr1->twiss->etax + eptr1->pred->twiss->etax)/2;
        phix1  = (eptr1->twiss->phix + eptr1->pred->twiss->phix)/2;
        betay1 = (eptr1->twiss->betay + eptr1->pred->twiss->betay)/2;
        phiy1  = (eptr1->twiss->phiy + eptr1->pred->twiss->phiy)/2;
      } else {
        betax1 = (eptr1->twiss->betax + twiss0->betax)/2;
        etax1  = (eptr1->twiss->etax + twiss0->etax)/2;
        phix1  = (eptr1->twiss->phix + twiss0->phix)/2;
        betay1 = (eptr1->twiss->betay + twiss0->betay)/2;
        phiy1  = (eptr1->twiss->phiy + twiss0->phiy)/2;
      }
      sqrt_betax = sqrt(betax1);
      sqrt3_betax = ipow(sqrt_betax, 3);

      /* first-order chromatic terms */
      /* h11001 and h00111 */
      coef = b2L-2*b3L1*etax1;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
      h11001 += coef*betax1/4;
      h00111 += -coef*betay1/4;
#else
      h11001.r += coef*betax1/4;
      h00111.r += -coef*betay1/4;
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
      /* h20001, h00201 */
      h20001 += coef/8*betax1*cos(2*phix1) + I*coef/8*betax1*sin(2*phix1);
      h00201 += -coef/8*betay1*cos(2*phiy1) + I*(-coef/8*betay1*sin(2*phiy1));
#else
      /* h20001, h00201 */
      h20001.r += coef/8*betax1*cos(2*phix1);
      h20001.i += coef/8*betax1*sin(2*phix1);
      h00201.r += -coef/8*betay1*cos(2*phiy1);
      h00201.i += -coef/8*betay1*sin(2*phiy1);
#endif

      /* h10002 */
      coef = b2L-b3L1*etax1;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
      h10002 += coef/2*etax1*sqrt_betax*cos(phix1) + I*coef/2*etax1*sqrt_betax*sin(phix1);
#else
      h10002.r += coef/2*etax1*sqrt_betax*cos(phix1);
      h10002.i += coef/2*etax1*sqrt_betax*sin(phix1);
#endif

      if (b3L1) {
        /* first-order geometric terms */
        /* h21000 */
        coef = -b3L1/8*sqrt3_betax;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
        h21000 += coef*cos(phix1) + I*coef*sin(phix1);
#else
        h21000.r += coef*cos(phix1);
        h21000.i += coef*sin(phix1);
#endif
        
        /* h30000 */
        coef = coef/3;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
        h30000 += coef*cos(3*phix1) + I*coef*sin(3*phix1);
#else
        h30000.r += coef*cos(3*phix1);
        h30000.i += coef*sin(3*phix1);
#endif
        
        /* h10110 */
        coef = b3L1/4*sqrt_betax*betay1;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
        h10110 += coef*cos(phix1) + I*coef*sin(phix1);
#else
        h10110.r += coef*cos(phix1);
        h10110.i += coef*sin(phix1);
#endif
        
        /* h10020 and h10200 */
        coef = coef/2;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
        h10020 += coef*cos(phix1-2*phiy1) + I*coef*sin(phix1-2*phiy1);
        h10200 += coef*cos(phix1+2*phiy1) + I*coef*sin(phix1+2*phiy1);
#else
        h10020.r += coef*cos(phix1-2*phiy1);
        h10020.i += coef*sin(phix1-2*phiy1);
        h10200.r += coef*cos(phix1+2*phiy1);
        h10200.i += coef*sin(phix1+2*phiy1);
#endif
      }
    }
    eptr1 = eptr1->succ;
  }

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  d->h11001 = cabs(h11001);
  d->h00111 = cabs(h00111);
  d->h20001 = cabs(h20001);
  d->h00201 = cabs(h00201);
  d->h10002 = cabs(h10002);

  d->h21000 = cabs(h21000);
  d->h30000 = cabs(h30000);
  d->h10110 = cabs(h10110);
  d->h10020 = cabs(h10020);
  d->h10200 = cabs(h10200);
#else
  d->h11001 = sqrt(h11001.r * h11001.r + h11001.i * h11001.i);
  d->h00111 = sqrt(h00111.r * h00111.r + h00111.i * h00111.i);
  d->h20001 = sqrt(h20001.r * h20001.r + h20001.i * h20001.i);
  d->h00201 = sqrt(h00201.r * h00201.r + h00201.i * h00201.i);
  d->h10002 = sqrt(h10002.r * h10002.r + h10002.i * h10002.i);

  d->h21000 = sqrt(h21000.r * h21000.r + h21000.i * h21000.i);
  d->h30000 = sqrt(h30000.r * h30000.r + h30000.i * h30000.i);
  d->h10110 = sqrt(h10110.r * h10110.r + h10110.i * h10110.i);
  d->h10020 = sqrt(h10020.r * h10020.r + h10020.i * h10020.i);
  d->h10200 = sqrt(h10200.r * h10200.r + h10200.i * h10200.i);
#endif

  /* compute second-order geomeric terms */
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h12000 = conj(h21000);
  h22000 = (3*sqr(cabs(h21000)) + sqr(cabs(h30000)))/64;
  d->h22000 = cabs(h22000);
#else
  h12000.r = h21000.r;
  h12000.i = -h21000.i;
  h22000.r = (3*sqr(sqrt(h21000.r * h21000.r + h21000.i * h21000.i)) + sqr(sqrt(h30000.r * h30000.r + h30000.i * h30000.i)))/64;
  h22000.i =0;
  d->h22000 = abs(h22000.r);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = 2*h21000*conj(h10110);
  t2 = h10020*conj(h10020);
  t3 = h10200*conj(h10200);
  h11110 = (t1+t2+t3)/16;
  d->h11110 = cabs(h11110);
#else
  t1.r = 2*(h21000.r*h10110.r + h21000.i*h10110.i);
  t1.i = 2*(h21000.i*h10110.r - h21000.r*h10110.i);
  t2.r = h10020.r*h10020.r + h10020.i*h10020.i;
  t2.i = h10020.i*h10020.r - h10020.r*h10020.i;
  t3.r = h10200.r*h10200.r + h10200.i*h10200.i;
  t3.i = h10200.i*h10200.r - h10200.r*h10200.i;
  h11110.r = (t1.r+t2.r+t3.r)/16;
  h11110.i = (t1.i+t2.i+t3.i)/16;
  d->h11110 = sqrt(h11110.r * h11110.r + h11110.i * h11110.i);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h00220 = (4*sqr(cabs(h10110)) + sqr(cabs(h10020)) + sqr(cabs(h10200)))/64;
  d->h00220 = cabs(h00220);
#else
  h00220.r = (4*sqr(sqrt(h10110.r * h10110.r + h10110.i * h10110.i)) + sqr(sqrt(h10020.r * h10020.r + h10020.i * h10020.i)) + sqr(sqrt(h10200.r * h10200.r + h10200.i * h10200.i)))/64;
  h00220.i = 0;
  d->h00220 = abs(h00220.r);
#endif
 
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h31000 = h30000*conj(h21000)/32;
  d->h31000 = cabs(h31000);
#else
  h31000.r = (h30000.r * h21000.r + h30000.i * h21000.i)/32;
  h31000.i = (h30000.i * h21000.r - h30000.r * h21000.i)/32;
  d->h31000 = sqrt(h31000.r * h31000.r + h31000.i * h31000.i);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h40000 = h30000*h21000/64;
  d->h40000 = cabs(h40000);
#else
  h40000.r = (h30000.r * h21000.r - h30000.i * h21000.i)/64;
  h40000.i = (h30000.i * h21000.r + h30000.r * h21000.i)/64;
  d->h40000 = sqrt(h40000.r * h40000.r + h40000.i * h40000.i);
#endif
  
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = 2*h30000*conj(h10110);
  t2 = 2*h21000*h10110;
  t3 = 4*h10200*h10020;
  h20110 = (t1+t2+t3)/64;
  d->h20110 = cabs(h20110);
#else
  t1.r = 2*(h30000.r*h10110.r + h30000.i*h10110.i);
  t1.i = 2*(h30000.i*h10110.r - h30000.r*h10110.i);
  t2.r = 2*(h21000.r*h10110.r - h21000.i*h10110.i);
  t2.i = 2*(h21000.i*h10110.r + h21000.r*h10110.i);
  t3.r = 4*(h10200.r*h10020.r - h10200.i*h10020.i);
  t3.i = 4*(h10200.i*h10020.r + h10200.r*h10020.i);
  h20110.r = (t1.r+t2.r+t3.r)/64;
  h20110.i = (t1.i+t2.i+t3.i)/64;
  d->h20110 = sqrt(h20110.r * h20110.r + h20110.i * h20110.i);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = 2*h10200*h12000;
  t2 = 2*h21000*conj(h10020);
  t3 = 4*h10200*conj(h10110);
  t4 = 4*h10110*conj(h10020);
  h11200 = (t1+t2+t3+t4)/64;
  d->h11200 = cabs(h11200);
#else
  t1.r = 2*(h10200.r*h12000.r - h10200.i*h12000.i);
  t1.i = 2*(h10200.i*h12000.r + h10200.r*h12000.i);
  t2.r = 2*(h21000.r*h10020.r + h21000.i*h10020.i);
  t2.i = 2*(h21000.i*h10020.r - h21000.r*h10020.i);
  t3.r = 4*(h10200.r*h10110.r + h10200.i*h10110.i);
  t3.i = 4*(h10200.i*h10110.r - h10200.r*h10110.i);
  t4.r = 4*(h10110.r*h10020.r + h10110.i*h10020.i);
  t4.i = 4*(h10110.i*h10020.r - h10110.r*h10020.i);
  h11200.r = (t1.r+t2.r+t3.r+t4.r)/64;
  h11200.i = (t1.i+t2.i+t3.i+t4.i)/64;
  d->h11200 = sqrt(h11200.r * h11200.r + h11200.i * h11200.i);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = h21000*h10020;
  t2 = h30000*conj(h10200);
  t3 = 4*h10110*h10020;
  h20020 = (t1+t2+t3)/64;
  d->h20020 = cabs(h20020);
#else
  t1.r = h21000.r*h10020.r - h21000.i*h10020.i;
  t1.i = h21000.i*h10020.r + h21000.r*h10020.i;
  t2.r = h30000.r*h10200.r + h30000.i*h10200.i;
  t2.i = h30000.i*h10200.r - h30000.r*h10200.i;
  t3.r = 4*(h10110.r*h10020.r - h10110.i*h10020.i);
  t3.i = 4*(h10110.i*h10020.r + h10110.r*h10020.i);
  h20020.r = (t1.r+t2.r+t3.r)/64;
  h20020.i = (t1.i+t2.i+t3.i)/64;
  d->h20020 = sqrt(h20020.r * h20020.r + h20020.i * h20020.i);
#endif

#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = h30000*conj(h10020)/64;
  t2 = h10200*h21000/64;
  t3 = h10110*h10200/16;
  h20200 = t1+t2+t3;
  d->h20200 = cabs(h20200);
#else
  t1.r = (h30000.r*h10020.r + h30000.i*h10020.i)/64;
  t1.i = (h30000.i*h10020.r - h30000.r*h10020.i)/64;
  t2.r = (h10200.r*h21000.r - h10200.i*h21000.i)/64;
  t2.i = (h10200.i*h21000.r + h10200.r*h21000.i)/64;
  t3.r = (h10110.r*h10200.r - h10110.i*h10200.i)/16;
  t3.i = (h10110.i*h10200.r + h10110.r*h10200.i)/16;
  h20200.r = t1.r+t2.r+t3.r;
  h20200.i = t1.i+t2.i+t3.i;
  d->h20200 = sqrt(h20200.r * h20200.r + h20200.i * h20200.i);
#endif
  
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  t1 = h10200*conj(h10110);
  t2 = h10110*conj(h10020);
  h00310 = (t1+t2)/32;
  d->h00310 = cabs(h00310);
#else
  t1.r = h10200.r*h10110.r + h10200.i*h10110.i;
  t1.i = h10200.i*h10110.r - h10200.r*h10110.i;
  t2.r = h10110.r*h10020.r + h10110.i*h10020.i;
  t2.i = h10110.i*h10020.r - h10110.r*h10020.i;
  h00310.r = (t1.r+t2.r)/32;
  h00310.i = (t1.i+t2.i)/32;
  d->h00310 = sqrt(h00310.r * h00310.r + h00310.i * h00310.i);
#endif
  
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  h00400 = h10200*conj(h10020)/64;
  d->h00400 = cabs(h00400);
#else
  h00400.r = (h10200.r * h10020.r + h10200.i * h10020.i)/64;
  h00400.i = (h10200.i * h10020.r - h10200.r * h10020.i)/64;
  d->h00400 = sqrt(h00400.r * h00400.r + h00400.i * h00400.i);
#endif
  
  nux = tune[0];
  nuy = tune[1];
  eptr1 = elem;
  d->dnux_dJx = d->dnux_dJy = d->dnuy_dJy = 0;
  while (eptr1) {
    b3L1 = 0;
    switch (eptr1->type) {
    case T_SEXT:
      b3L1 = ((SEXT*)eptr1->p_elem)->k2 * ((SEXT*)eptr1->p_elem)->length/2;
      break;
    case T_KSEXT:
      b3L1 = ((KSEXT*)eptr1->p_elem)->k2 * ((KSEXT*)eptr1->p_elem)->length/2;
      break;
    case T_KQUSE:
      b3L1 = ((KQUSE*)eptr1->p_elem)->k2 * ((KQUSE*)eptr1->p_elem)->length/2;
      break;
    case T_RBEN:
    case T_SBEN:
      b3L1 = ((BEND*)eptr1->p_elem)->k2 * ((BEND*)eptr1->p_elem)->length/2;
      break;
    case T_CSBEND:
      b3L1 = ((CSBEND*)eptr1->p_elem)->k2 * ((CSBEND*)eptr1->p_elem)->length/2;
      break;
    case T_CSRCSBEND:
      b3L1 = ((CSRCSBEND*)eptr1->p_elem)->k2 * ((CSRCSBEND*)eptr1->p_elem)->length/2;
      break;
    default:
      break;
    }      
    if (b3L1) {
      if (eptr1->pred) {
        betax1 = (eptr1->twiss->betax + eptr1->pred->twiss->betax)/2;
        phix1  = (eptr1->twiss->phix + eptr1->pred->twiss->phix)/2;
        betay1 = (eptr1->twiss->betay + eptr1->pred->twiss->betay)/2;
        phiy1  = (eptr1->twiss->phiy + eptr1->pred->twiss->phiy)/2;
      } else {
        betax1 = (eptr1->twiss->betax + twiss0->betax)/2;
        phix1  = (eptr1->twiss->phix + twiss0->phix)/2;
        betay1 = (eptr1->twiss->betay + twiss0->betay)/2;
        phiy1  = (eptr1->twiss->phiy + twiss0->phiy)/2;
      }
      eptr2 = elem;
      while (eptr2) {
        b3L2 = 0;
        switch (eptr2->type) {
        case T_SEXT:
          b3L2 = ((SEXT*)eptr2->p_elem)->k2 * ((SEXT*)eptr2->p_elem)->length/2;
          break;
        case T_KSEXT:
          b3L2 = ((KSEXT*)eptr2->p_elem)->k2 * ((KSEXT*)eptr2->p_elem)->length/2;
          break;
        case T_KQUSE:
          b3L2 = ((KQUSE*)eptr2->p_elem)->k2 * ((KQUSE*)eptr2->p_elem)->length/2;
          break;
        case T_RBEN:
        case T_SBEN:
          b3L2 = ((BEND*)eptr2->p_elem)->k2 * ((BEND*)eptr2->p_elem)->length/2;
          break;
        case T_CSBEND:
          b3L2 = ((CSBEND*)eptr2->p_elem)->k2 * ((CSBEND*)eptr2->p_elem)->length/2;
          break;
        case T_CSRCSBEND:
          b3L2 = ((CSRCSBEND*)eptr2->p_elem)->k2 * ((CSRCSBEND*)eptr2->p_elem)->length/2;
          break;
        default:
          break;
        }
        if (b3L2) {
          if (eptr2->pred) {
            betax2 = (eptr2->twiss->betax + eptr2->pred->twiss->betax)/2;
            phix2  = (eptr2->twiss->phix + eptr2->pred->twiss->phix)/2;
            betay2 = (eptr2->twiss->betay + eptr2->pred->twiss->betay)/2;
            phiy2  = (eptr2->twiss->phiy + eptr2->pred->twiss->phiy)/2;
          } else {
            betax2 = (eptr2->twiss->betax + twiss0->betax)/2;
            phix2  = (eptr2->twiss->phix + twiss0->phix)/2;
            betay2 = (eptr2->twiss->betay + twiss0->betay)/2;
            phiy2  = (eptr2->twiss->phiy + twiss0->phiy)/2;
          }
          d->dnux_dJx += b3L1*b3L2/(-16*PI)*pow(betax1*betax2, 1.5)*
                          (3*cos(fabs(phix1-phix2)-PI*nux)/sin(PI*nux) + cos(fabs(3*(phix1-phix2))-3*PI*nux)/sin(3*PI*nux));
          d->dnux_dJy += b3L1*b3L2/(8*PI)*sqrt(betax1*betax2)*betay1*
            (2*betax2*cos(fabs(phix1-phix2)-PI*nux)/sin(PI*nux) 
             - betay2*cos(fabs(phix1-phix2)+2*fabs(phiy1-phiy2)-PI*(nux+2*nuy))/sin(PI*(nux+2*nuy))
             + betay2*cos(fabs(phix1-phix2)-2*fabs(phiy1-phiy2)-PI*(nux-2*nuy))/sin(PI*(nux-2*nuy)));
          d->dnuy_dJy += b3L1*b3L2/(-16*PI)*sqrt(betax1*betax2)*betay1*betay2*
            (4*cos(fabs(phix1-phix2)-PI*nux)/sin(PI*nux) 
             + cos(fabs(phix1-phix2)+2*fabs(phiy1-phiy2)-PI*(nux+2*nuy))/sin(PI*(nux+2*nuy)) 
             + cos(fabs(phix1-phix2)-2*fabs(phiy1-phiy2)-PI*(nux-2*nuy))/sin(PI*(nux-2*nuy)));
        }
        eptr2 = eptr2->succ;
      }
    }
    eptr1 = eptr1->succ;
  }
}




