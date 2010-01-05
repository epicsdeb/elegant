/*************************************************************************\
* Copyright (c) 2007 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2007 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: moments.c
 * purpose: computation of beam moments
 *
 * Michael Borland, 2007
 */
#include "mdb.h"
#include "track.h"
#include "moments.h"
#include "matlib.h"
#include <stddef.h>

void determineEquilibriumMoments(double **R, double *D, SIGMA_MATRIX *sigma0);
void propagateBeamMoments(RUN *run, LINE_LIST *beamline, double *traj);
void storeFitpointMomentsParameters(MARK *mark, char *name, long occurence, SIGMA_MATRIX *sigma0, double *centroid);
void prepareMomentsArray(double *data, ELEMENT_LIST *elem, double *sigma);
void computeNaturalEmittances(VMATRIX *M, double *sigma, double *emittance);

static long momentsInitialized = 0;
static long SDDSMomentsInitialized = 0;
static SDDS_TABLE SDDSMoments;
static long momentsCount = 0;

#define IC_ELEMENT 0
#define IC_OCCURENCE 1
#define IC_TYPE 2
#define IC_S 3
#define IC_PCENTRAL 4
#define IC_EMITTANCE (IC_PCENTRAL+1+21+6)
#define IC_SBETA (IC_EMITTANCE+3)
#define IC_EMITBETA (IC_SBETA+10)
#define N_COLUMNS (IC_PCENTRAL+1+21+6+3+10+2)
static SDDS_DEFINITION column_definition[N_COLUMNS] = {
{"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
{"ElementOccurence", "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
{"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    {"s", "&column name=s, type=double, units=m, description=Distance &end"},
    {"pCentral0", "&column name=pCentral0, type=double, units=\"m$be$nc\", symbol=\"p$bcent$n\", description=\"Initial central momentum\" &end"},
    {"s1",    "&column name=s1, symbol=\"$gs$r$b1$n\", units=m, type=double, description=\"sqrt(<x*x>)\" &end"},
    {"s12",    "&column name=s12, symbol=\"$gs$r$b12$n\", units=m, type=double, description=\"<x*xp'>\" &end"},
    {"s13",    "&column name=s13, symbol=\"$gs$r$b13$n\", units=\"m$a2$n\", type=double, description=\"<x*y>\" &end"},
    {"s14",    "&column name=s14, symbol=\"$gs$r$b14$n\", units=m, type=double, description=\"<x*y'>\" &end"},
    {"s15",    "&column name=s15, symbol=\"$gs$r$b15$n\", units=\"m$a2$n\", type=double, description=\"<x*s>\" &end"},
    {"s16",    "&column name=s16, symbol=\"$gs$r$b16$n\", units=m, type=double, description=\"<x*delta>\" &end"},
    {"s2",    "&column name=s2, symbol=\"$gs$r$b2$n\", type=double, description=\"sqrt(<x'*x'>)\" &end"},
    {"s23",    "&column name=s23, symbol=\"$gs$r$b23$n\", units=m, type=double, description=\"<x'*y>\" &end"},
    {"s24",    "&column name=s24, symbol=\"$gs$r$b24$n\", type=double, description=\"<x'*y'>\" &end"},
    {"s25",    "&column name=s25, symbol=\"$gs$r$b25$n\", units=m, type=double, description=\"<x'*s>\" &end"},
    {"s26",    "&column name=s26, symbol=\"$gs$r$b26$n\", type=double, description=\"<x'*delta>\" &end"},
    {"s3",    "&column name=s3, symbol=\"$gs$r$b3$n\", units=m, type=double, description=\"sqrt(<y*y>)\" &end"},
    {"s34",    "&column name=s34, symbol=\"$gs$r$b34$n\", units=m, type=double, description=\"<y*y'>\" &end"},
    {"s35",    "&column name=s35, symbol=\"$gs$r$b35$n\", units=\"m$a2$n\", type=double, description=\"<y*s>\" &end"},
    {"s36",    "&column name=s36, symbol=\"$gs$r$b36$n\", units=m, type=double, description=\"<y*delta>\" &end"},
    {"s4",    "&column name=s4, symbol=\"$gs$r$b4$n\", type=double, description=\"sqrt(<y'*y')>\" &end"},
    {"s45",    "&column name=s45, symbol=\"$gs$r$b45$n\", units=m, type=double, description=\"<y'*s>\" &end"},
    {"s46",    "&column name=s46, symbol=\"$gs$r$b46$n\", type=double, description=\"<s'*delta>\" &end"},
    {"s5",    "&column name=s5, symbol=\"$gs$r$b5$n\", units=m, type=double, description=\"sqrt(<s*s>)\" &end"},
    {"s56",    "&column name=s56, symbol=\"$gs$r$b56$n\", units=m, type=double, description=\"<s*delta>\" &end"},
    {"s6",    "&column name=s6, symbol=\"$gs$r$b6$n\", type=double, description=\"sqrt(<delta*delta>)\" &end"},
    {"c1",    "&column name=c1, symbol=\"c$b1$n\", units=m, type=double, description=\"<x>\" &end"},
    {"c2",    "&column name=c2, symbol=\"c$b2$n\", units=, type=double, description=\"<x'>\" &end"},
    {"c3",    "&column name=c3, symbol=\"c$b3$n\", units=m, type=double, description=\"<y>\" &end"},
    {"c4",    "&column name=c4, symbol=\"c$b4$n\", units=, type=double, description=\"<y'>\" &end"},
    {"c5",    "&column name=c5, symbol=\"c$b5$n\", units=m, type=double, description=\"<s>\" &end"},
    {"c6",    "&column name=c6, symbol=\"c$b6$n\", units=, type=double, description=\"<delta>\" &end"},
    {"ex",    "&column name=ex, symbol=\"$ge$r$bx$n\", units=m, type=double, description=\"Projected horizontal emittance\" &end"},
    {"ey",    "&column name=ey, symbol=\"$ge$r$by$n\", units=m, type=double, description=\"Projected vertical emittance\" &end"},
    {"ez",    "&column name=ez, symbol=\"$ge$r$bz$n\", units=m, type=double, description=\"Projected longitudinal emittance\" &end"},
    {"s1beta",    "&column name=s1beta, symbol=\"$gs$r$b1,$gb$r$n\", units=m, type=double, description=\"sqrt(<x*x>) (betatron)\" &end"},
    {"s12beta",    "&column name=s12beta, symbol=\"$gs$r$b12,$gb$r$n\", units=m, type=double, description=\"<x*xp'> (betatron)\" &end"},
    {"s13beta",    "&column name=s13beta, symbol=\"$gs$r$b13,$gb$r$n\", units=\"m$a2$n\", type=double, description=\"<x*y> (betatron)\" &end"},
    {"s14beta",    "&column name=s14beta, symbol=\"$gs$r$b14,$gb$r$n\", units=m, type=double, description=\"<x*y'> (betatron)\" &end"},
    {"s2beta",    "&column name=s2beta, symbol=\"$gs$r$b2,$gb$r$n\", type=double, description=\"sqrt(<x'*x'>) (betatron)\" &end"},
    {"s23beta",    "&column name=s23beta, symbol=\"$gs$r$b23,$gb$r$n\", units=m, type=double, description=\"<x'*y> (betatron)\" &end"},
    {"s24beta",    "&column name=s24beta, symbol=\"$gs$r$b24,$gb$r$n\", type=double, description=\"<x'*y'> (betatron)\" &end"},
    {"s3beta",    "&column name=s3beta, symbol=\"$gs$r$b3,$gb$r$n\", units=m, type=double, description=\"sqrt(<y*y>) (betatron)\" &end"},
    {"s34beta",    "&column name=s34beta, symbol=\"$gs$r$b34,$gb$r$n\", units=m, type=double, description=\"<y*y'> (betatron)\" &end"},
    {"s4beta",    "&column name=s4beta, symbol=\"$gs$r$b4,$gb$r$n\", type=double, description=\"sqrt(<y'*y')> (betatron)\" &end"},
    {"exbeta",    "&column name=exbeta, symbol=\"$ge$r$bx,$gb$r$n\", units=m, type=double, description=\"Projected horizontal betatron emittance\" &end"},
    {"eybeta",    "&column name=eybeta, symbol=\"$ge$r$by,$gb$r$n\", units=m, type=double, description=\"Projected vertical betatron emittance\" &end"},
};

#define IP_STEP 0
#define IP_STAGE 1
#define IP_PCENTRAL 2
#define IP_E1 3
#define IP_E2 4
#define IP_E3 5
#define N_PARAMETERS IP_E3+1
static SDDS_DEFINITION parameter_definition[N_PARAMETERS] = {
{"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
{"Stage", "&parameter name=Stage, type=string, description=\"Stage of computation\" &end"},
{"pCentral", "&parameter name=pCentral, type=double, units=\"m$be$nc\", description=\"Central momentum\" &end"},
{"e1", "&parameter name=e1, symbol=\"$ge$r$b1$n\", type=double, units=m,  description=\"Emittance of mode 1\" &end"},
{"e2", "&parameter name=e2, symbol=\"$ge$r$b2$n\", type=double, units=m,  description=\"Emittance of mode 2\" &end"},
{"e3", "&parameter name=e3, symbol=\"$ge$r$b3$n\", type=double, units=m,  description=\"Emittance of mode 3\" &end"},
} ;

void dumpBeamMoments(
  LINE_LIST *beamline,
  long n_elem,
  long final_values_only,
  long tune_corrected,
  RUN *run
  )
{
  double data[N_COLUMNS], *emit;
  long j, row_count, elemCheck;
  char *stage;
  ELEMENT_LIST *elem;
  SIGMA_MATRIX *sigma0;
  double eNatural[3] = {0,0,0};
  
  if (tune_corrected==1)
    stage = "tunes corrected";
  else
    stage = "tunes uncorrected";

  if (!SDDS_StartTable(&SDDSMoments, final_values_only?1:n_elem+1)) {
    SDDS_SetError("Problem starting SDDS table (dumpBeamMoments)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  
  if (!SDDS_SetParameters(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                          IP_STEP, momentsCount, IP_STAGE, stage, 
                          IP_PCENTRAL, run->p_central, -1)) {
    SDDS_SetError("Problem setting SDDS parameters (dumpBeamMoments 1)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  elem   = beamline->elem_twiss;
  sigma0 = beamline->sigmaMatrix0;
  emit = data+IC_EMITTANCE;
  if (!final_values_only) {
    row_count = 0;
    data[IC_PCENTRAL] = elem->Pref_input;
    prepareMomentsArray(data, elem, (double*)sigma0->sigma);
    data[IC_S] = 0;     /* position for initial values */
    for (j=IC_S; j<N_COLUMNS; j++)
      if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count++, 
                           IC_ELEMENT, "_BEG_", IC_OCCURENCE, (long)1, IC_TYPE, "MARK", -1)) {
      SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    elemCheck = 0;
    while (elem) {
      data[IC_S] = elem->end_pos;     /* position */
      data[IC_PCENTRAL] = elem->Pref_output;
      if (!elem->sigmaMatrix)
        bomb("Sigma matrix data not computed prior to dumpBeamMoments() call (2)", NULL);
      prepareMomentsArray(data, elem, elem->sigmaMatrix->sigma);
      for (j=IC_S; j<N_COLUMNS; j++)
        if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, j, data[j], -1)) {
          SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row_count, 
                             IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                             IC_TYPE, entity_name[elem->type], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      elemCheck++;
      row_count++;
      elem = elem->succ;
    }
    if (elemCheck!=n_elem)
      bomb("element count error in dumpBeamMoments()", NULL);
  }
  else {
    /* find final element */
    elemCheck = 0;
    while (1) {
      if (!elem->sigmaMatrix)
        bomb("Sigma matrix data not computed prior to dumpBeamMoments() call (3)", NULL);
      elemCheck++;
      if (!elem->succ)
        break;
      elem = elem->succ;
    }
    if (elemCheck!=n_elem)
      bomb("element count error in dumpBeamMoments()", NULL);
    prepareMomentsArray(data, elem, elem->sigmaMatrix->sigma);
    for (j=IC_S; j<N_COLUMNS; j++)
      if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, j, data[j], -1)) {
        SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    if (!SDDS_SetRowValues(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, 0, 
                           IC_ELEMENT, elem->name, IC_OCCURENCE, elem->occurence, 
                           IC_TYPE, entity_name[elem->type], -1)) {
      SDDS_SetError("Problem setting SDDS rows (dumpBeamMoments)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }

  if (equilibrium)
    computeNaturalEmittances(beamline->Mld, beamline->sigmaMatrix0->sigma, eNatural); 
  if (!SDDS_SetParameters(&SDDSMoments, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE,
                          IP_E1, eNatural[0], 
                          IP_E2, eNatural[1], 
                          IP_E3, eNatural[2], -1)) {
    SDDS_SetError("Problem setting SDDS emittance parameters (dumpBeamMoments)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!SDDS_WriteTable(&SDDSMoments)) {
    SDDS_SetError("Unable to write Twiss parameter data (dumpBeamMoments)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(&SDDSMoments);
  if (!SDDS_EraseData(&SDDSMoments)) {
    SDDS_SetError("Unable to erase Twiss parameter data (dumpBeamMoments)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (tune_corrected)
    momentsCount++;
}

void setupMomentsOutput(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, long *doMomentsOutput,
                        long default_order)
{
  /* process namelist input */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&moments_output, nltext);
  if (echoNamelists) print_namelist(stdout, &moments_output);
  
#if USE_MPI
  if (!writePermitted)
    filename = NULL;
#endif 

  if (filename)
    filename = compose_filename(filename, run->rootname);
  *doMomentsOutput = output_at_each_step;
  
  if (reference_file && matched)
    bomb("reference_file and matched=1 are incompatible", NULL);
  if (!matched) {
    if (reference_file) {
      if (reference_element && reference_element_occurrence<0)
        bomb("invalid value of reference_element_occurrence---use 0 for last occurrence, >=1 for specific occurrence.", NULL);
      bomb("reference file feature not implemented yet.", NULL);
    }
    if (beta_x<=0 || beta_y<=0 || beta_z<=0 || emit_x<0 || emit_y<0 || emit_z<0)
      bomb("invalid initial beta-functions given in moments_output namelist", NULL);
  }

  if (filename) {
    SDDS_ElegantOutputSetup(&SDDSMoments, filename, SDDS_BINARY, 1, "Beam moments",
                            run->runfile, run->lattice, 
                            parameter_definition,  N_PARAMETERS, 
                            column_definition, N_COLUMNS, "setupMomentsOutput",
                            SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    SDDSMomentsInitialized = 1;
    momentsCount = 0;
  }
  else
    SDDSMomentsInitialized = 0;
  momentsInitialized = 1;
}

void finishMomentsOutput(void)
{
  if (SDDSMomentsInitialized && !SDDS_Terminate(&SDDSMoments)) {
    SDDS_SetError("Problem terminating SDDS output (finishMomentsOutput)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  SDDSMomentsInitialized = momentsCount = momentsInitialized = 0;
}

long runMomentsOutput(RUN *run, LINE_LIST *beamline, double *startingCoord, long tune_corrected, long writeToFile)
{
  ELEMENT_LIST *eptr, *elast;
  long n_elem, last_n_elem, i;
  
#ifdef DEBUG
  fprintf(stdout, "now in runMomentsOutput\n");
  fflush(stdout);
#endif

  if (!momentsInitialized)
    return 1;
  
  if (tune_corrected==0 && !output_before_tune_correction)
    return 1;

  /* Computations will start at the beginning of the beamline, or at the
   * first recirculation element 
   */
  eptr = beamline->elem_twiss = &(beamline->elem);
  n_elem = last_n_elem = beamline->n_elems;
  while (eptr) {
    if (eptr->type==T_RECIRC) {
      last_n_elem = n_elem;
      beamline->elem_twiss = beamline->elem_recirc = eptr;
    }
    eptr = eptr->succ;
    n_elem --;
  }
  n_elem = last_n_elem;

  if (!beamline->sigmaMatrix0)
    beamline->sigmaMatrix0 = tmalloc(sizeof(*(beamline->sigmaMatrix0)));

  if (verbosity>0) {
    printf("\nPerforming beam moments computation.\n");
    fflush(stdout);
  }

  if (beamline->Mld) {
    free_matrices(beamline->Mld);
    free(beamline->Mld);
    beamline->Mld = NULL;
  }
    
  if (startingCoord) {
    VMATRIX *M1;
    M1 = tmalloc(sizeof(*M1));
    initialize_matrices(M1, 1);
    for (i=0; i<6; i++) {
      M1->C[i] = startingCoord[i];
      M1->R[i][i] = 1;
    }
    beamline->Mld = accumulateRadiationMatrices(beamline->elem_twiss, run, M1, 1, radiation, n_slices);
    free_matrices(M1);
    free(M1);
  }
  else {
    beamline->Mld = accumulateRadiationMatrices(beamline->elem_twiss, run, NULL, 1, radiation, n_slices);
  }

  if (verbosity>0) {
    char text[200];
    sprintf(text, "** One-turn, on-orbit matrix with%sradiation:", radiation?" ":"out ");
    print_matrices(stdout, text, beamline->Mld);
  }
  if (verbosity>1) {
    long j;
    printf("** One-turn diffusion matrix: \n");
    for (i=0; i<6; i++) 
      for (j=0; j<6; j++)
        printf("%13.6e%c", beamline->elast->accumD[sigmaIndex3[i][j]], j==5?'\n':' ');
  }
  
  if (equilibrium) {
    /* Compute equilibrium moments */
    determineEquilibriumMoments(beamline->Mld->R, beamline->elast->accumD, beamline->sigmaMatrix0);
  } else {
    if (matched) {
      /* Use periodic lattice functions */
      if (!(beamline->flags&BEAMLINE_TWISS_DONE))
	bomb("no twiss parameters computed for matched moments propogation", NULL);
      /* Determine starting moments from twiss parameter values */
      setStartingMoments(beamline->sigmaMatrix0,
			 emit_x, beamline->twiss0->betax, beamline->twiss0->alphax, beamline->twiss0->etax, beamline->twiss0->etapx,
			 emit_y, beamline->twiss0->betay, beamline->twiss0->alphay, beamline->twiss0->etay, beamline->twiss0->etapy, 
			 emit_z, beta_z, alpha_z);
    } else 
      /* Determine starting moments from twiss parameter values */
      setStartingMoments(beamline->sigmaMatrix0,
			 emit_x, beta_x, alpha_x, eta_x, etap_x,
			 emit_y, beta_y, alpha_y, eta_y, etap_y, 
			 emit_z, beta_z, alpha_z);
  }
  if (verbosity>1) {
    long j;
    printf("** Starting sigma matrix: \n");
    for (i=0; i<6; i++) 
      for (j=0; j<6; j++)
        printf("%13.6e%c", beamline->sigmaMatrix0->sigma[sigmaIndex3[i][j]], j==5?'\n':' ');
  }
  
  /* Propagate moments to each element */
  propagateBeamMoments(run, beamline, startingCoord);

  elast = beamline->elast;
  if (verbosity>1) {
    long j;
    printf("** Final sigma matrix: \n");
    for (i=0; i<6; i++) 
      for (j=0; j<6; j++)
        printf("%13.6e%c", elast->sigmaMatrix->sigma[sigmaIndex3[i][j]], j==5?'\n':' ');
  }

  if (SDDSMomentsInitialized && writeToFile)
    dumpBeamMoments(beamline, n_elem, final_values_only, tune_corrected, run);

  return 1;
}


void setStartingMoments(SIGMA_MATRIX *sm, 
                        double emit_x, double beta_x, double alpha_x, double eta_x, double etap_x,
                        double emit_y, double beta_y, double alpha_y, double eta_y, double etap_y,
                        double emit_z, double beta_z, double alpha_z)
{
  double sDelta;

  sm->sigma[sigmaIndex3[4][4]] = emit_z*beta_z;
  sm->sigma[sigmaIndex3[5][5]] = emit_z/beta_z*(1+sqr(alpha_z));
  sm->sigma[sigmaIndex3[4][5]] = -emit_z*alpha_z;
  sDelta = sqrt(sm->sigma[sigmaIndex3[5][5]]);
  
  sm->sigma[sigmaIndex3[0][0]] = emit_x*beta_x + sqr(sDelta*eta_x);
  sm->sigma[sigmaIndex3[1][1]] = emit_x/beta_x*(1+sqr(alpha_x))+sqr(sDelta*etap_x);
  sm->sigma[sigmaIndex3[0][1]] = -emit_x*alpha_x + sqr(sDelta)*eta_x*etap_x;
  sm->sigma[sigmaIndex3[0][5]] = sqr(sDelta)*eta_x;
  sm->sigma[sigmaIndex3[1][5]] = sqr(sDelta)*etap_x;

  sm->sigma[sigmaIndex3[2][2]] = emit_y*beta_y + sqr(sDelta*eta_y);
  sm->sigma[sigmaIndex3[3][3]] = emit_y/beta_y*(1+sqr(alpha_y))+sqr(sDelta*etap_y);
  sm->sigma[sigmaIndex3[2][3]] = -emit_y*alpha_y + sqr(sDelta)*eta_y*etap_y;
  sm->sigma[sigmaIndex3[2][5]] = sqr(sDelta)*eta_y;
  sm->sigma[sigmaIndex3[3][5]] = sqr(sDelta)*etap_y;
}

void fillSigmaPropagationMatrix(double **Ms, double **R)
{
  long i, j, k, l, m;
  double Rik;
  
  for (i=0; i<6; i++) {
    for (j=i; j<6; j++) {
      m = sigmaIndex3[i][j];
      for (k=0; k<6; k++) {
        Rik = R[i][k];
        for (l=k; l<6; l++) {
          Ms[m][sigmaIndex3[k][l]] = Rik*R[j][l] + (k!=l ? R[i][l]*R[j][k] : 0);
        }
      }
    }
  }
}

void propagateBeamMoments(RUN *run, LINE_LIST *beamline, double *traj)
{
  long i, j;
  ELEMENT_LIST *elem;
  VMATRIX *M1, *M2, *Me;
  SIGMA_MATRIX *S1, *S2;
  double path[6];
  MATRIX *Ms;
  
  /* Allocate memory to store sigma matrix as we propagate, copy initial matrix */
  S1 = tmalloc(sizeof(*S1));
  S2 = tmalloc(sizeof(*S1));
  memcpy(S1, beamline->sigmaMatrix0, sizeof(*S1));
  
  M1 = tmalloc(sizeof(*M1));
  M2 = tmalloc(sizeof(*M2));
  initialize_matrices(M1, 1);
  initialize_matrices(M2, 1);
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

  elem = beamline->elem_twiss;
  m_alloc(&Ms, 21, 21);
  while (elem) {
    if (!(elem->sigmaMatrix)) 
      elem->sigmaMatrix = tmalloc(sizeof(*(elem->sigmaMatrix)));

    if ((Me = elem->Mld)) {
      fillSigmaPropagationMatrix(Ms->a, Me->R);

      for (i=0; i<21; i++) {
        S2->sigma[i] = 0;
        for (j=0; j<21; j++) {
          S2->sigma[i] += Ms->a[i][j]*S1->sigma[j];
        }
      }
      if (elem->D) 
        for (i=0; i<21; i++)
          S2->sigma[i] += elem->D[i];
      memcpy(elem->sigmaMatrix, S2, sizeof(*S2));
      memcpy(S1, S2, sizeof(*S2));
    } else 
      /* Assume it doesn't modify the sigma matrix */
      memcpy(elem->sigmaMatrix, S1, sizeof(*S1));
    if (elem->type==T_MARK && ((MARK*)elem->p_elem)->fitpoint) 
      storeFitpointMomentsParameters((MARK*)elem->p_elem, elem->name, elem->occurence, elem->sigmaMatrix, elem->Mld->C);
    elem = elem->succ;
  }
  free_matrices(M1); free(M1);
  free_matrices(M2); free(M2);
  free(S1);
  free(S2);
  m_free(&Ms);
}


void determineEquilibriumMoments
  (
   double **R,                   /* revolution matrix (input) */
   double *D,                    /* diffusion matrix (input) */
   SIGMA_MATRIX *sigmaMatrix0    /* sigma matrix (output) */
   )
{
  MATRIX *Ms, *Md, *M1, *M2, *M3;
  long i;
  
  m_alloc(&Ms, 21, 21);   /* sigma matrix propagator */
  m_alloc(&Md, 21,  1);   /* diffiusion matrix */
  m_alloc(&M1, 21, 21);   /* work matrix */
  m_alloc(&M2, 21, 21);   /* work matrix */
  m_alloc(&M3, 21,  1);   /* work matrix */

  fillSigmaPropagationMatrix(Ms->a, R);

  m_zero(Md);
  for (i=0; i<21; i++)
    Md->a[i][0] = D[i];

  /* S = Inv(I-Ms) * D */
  m_identity(M1);
  m_subtract(M2, M1, Ms);
  m_invert(M1, M2);
  m_mult(M3, M1, Md);
  
  for (i=0; i<21; i++)
    sigmaMatrix0->sigma[i] = M3->a[i][0];
  
  m_free(&Ms);
  m_free(&Md);
  m_free(&M1);
  m_free(&M2);
  m_free(&M3);
}

void storeFitpointMomentsParameters(MARK *mark, char *name, long occurence, SIGMA_MATRIX *sigma0, double *centroid)
{
  char s[1000];
  long i;

  if (!(mark->init_flags&32)) {
    mark->moments_mem = tmalloc(sizeof(*(mark->moments_mem))*(21+6));
    mark->init_flags |= 32;
    for (i=0; i<21; i++) {
      sprintf(s, "%s#%ld.s%ld%ldm", name, occurence, 
              sigmaIndex1[i]+1, sigmaIndex2[i]+1);
      mark->moments_mem[i] = rpn_create_mem(s, 0);
    }
    for (i=0; i<6; i++) {
      sprintf(s, "%s#%ld.c%ldm", name, occurence, i+1);
      mark->moments_mem[i+21] = rpn_create_mem(s, 0);
    }
  }
  for (i=0; i<21; i++) 
    rpn_store(sigma0->sigma[i], NULL, mark->moments_mem[i]);
  for (i=0; i<6; i++) 
    rpn_store(centroid[i], NULL, mark->moments_mem[i+21]);
}


void prepareMomentsArray(double *data, ELEMENT_LIST *elem, double *sigma)
{
  long i, j, k, l, plane;
  double *emit;
  double sBeta[4][4];

  for (i=0; i<N_COLUMNS; i++)
    data[i] = -1;
  
  data[IC_S] = elem->end_pos;     /* position */
  data[IC_PCENTRAL] = elem->Pref_output;
  copy_doubles(data+IC_PCENTRAL+1, sigma, 21);
  copy_doubles(data+IC_PCENTRAL+1+21, elem->Mld->C, 6);
  for (i=0; i<6; i++) {
    k = sigmaIndex3[i][i];
    data[IC_PCENTRAL+1+k] = sqrt(data[IC_PCENTRAL+1+k]);
  }

  emit = data+IC_EMITTANCE;
  for (plane=0; plane<3; plane++) {
    emit[plane] = sigma[sigmaIndex3[0+plane*2][0+plane*2]]*sigma[sigmaIndex3[1+plane*2][1+plane*2]] - 
        sqr(sigma[sigmaIndex3[0+plane*2][1+plane*2]]);
    if (emit[plane]>0)
      emit[plane] = sqrt(emit[plane]);
    else
      emit[plane] = -1;
  }

  /* compute betatron quantities */
  for (i=l=0; i<4; i++)
    for (j=i; j<4; j++, l++) {
      data[IC_SBETA+l] = sigma[sigmaIndex3[i][j]] -
        sigma[sigmaIndex3[i][5]]*sigma[sigmaIndex3[j][5]]/sigma[sigmaIndex3[5][5]];
      sBeta[i][j] = sBeta[j][i] = data[IC_SBETA+l];
      if (i==j) {
        if (data[IC_SBETA+l]>0)
          data[IC_SBETA+l] = sqrt(data[IC_SBETA+l]);
        else
          data[IC_SBETA+l] = 0;
      }
    }
  
  emit = data+IC_EMITBETA;
  for (plane=0; plane<2; plane++) {
    emit[plane] = sBeta[2*plane][2*plane]*sBeta[2*plane+1][2*plane+1] - 
      sqr(sBeta[2*plane][2*plane+1]);
    if (emit[plane]>0)
      emit[plane] = sqrt(emit[plane]);
    else
      emit[plane] = -1;
  }
}

/* Following code is from V. Sajaev's calculateEnvelopes.c program, adapted to 
 * elegant by M. Borland and V. Sajaev.
 */

#define MATDIM  6
#define MATDIM2 36

void NormalizeEigenvectors(int dim, double *V, int debug);
double *AddMM1(int sum, int rows, int cols, double *M1, double *M2);
double *MatrixProduct1 (int rows1, int cols1, double *T1, int rows2, int cols2, double *T2);
double *TransposeM(int rows, int cols, double *M);
void MatrixPrintout1(char *string, double *AA, int N, int M);
void dgeev_();

void computeNaturalEmittances(VMATRIX *Mld, double *sigmaMatrix, double *emittance) 
{
  int i, j, k, eigenModesNumber, dim=MATDIM;
  double WR[MATDIM], WI[MATDIM], VL[MATDIM2], VR[MATDIM2], M[MATDIM2], Mcopy[MATDIM2], work[1000];
  double *Rdiag;
  char JOBVL, JOBVR;
  int N, LDA, LDVL, LDVR, lwork, info;
  double ReV[MATDIM2], ImV[MATDIM2];
  double *M1, *M2, *M3, *M4;
  
  eigenModesNumber = 3;
 
  /* Copy the revolution matrix into the working buffer */
  for (i=0; i<MATDIM; i++)
    for (j=0; j<MATDIM; j++) {
      M[i*MATDIM+j] = Mld->R[i][j];
    }
  
  memcpy(Mcopy, M, sizeof(*Mcopy)*MATDIM2);
    
  /*--- Calculating eigenvectors using dgeev ... */
  /* VR is right-hand side eigenvectors such that: VR^transp * M = lamdba * VR^transp
     VR[0 to 5] are real components of vector 1, VR[6 to 11] are imagenary components of vector 1 and so on.
     see description of procedure on the web */
#if defined(SUNPERF) || defined(LAPACK) || defined(CLAPACK)
  JOBVL = 'N';
  JOBVR = 'V';
  N = LDA = LDVR = MATDIM;
  LDVL = 1;
  lwork = 1000;
  dgeev_((char*)&JOBVL, (char*)&JOBVR, (int*)&N, (double*)&Mcopy,
         (int*)&LDA, (double*)&WR, (double*)&WI, (double*)&VL,
         (int*)&LDVL, (double*)&VR, (int*)&LDVR, (double*)&work,
         (int*)&lwork, (int*)&info);
#else
  fprintf(stderr, "Error calling dgeev. You will need to install LAPACK and rebuild elegant\n");
  exit(1);
#endif

  if (info != 0) {
    if (info < 0) { printf("Error calling dgeev, argument %d.\n", abs(info)); }
    if (info > 0) { printf("Error running dgeev, calculation of eigenvalue number %d failed.\n", info); }
    exit(1);
  }

  /*--- Sorting of eigenvalues and eigenvectors according to (x,y,z)... */
  SortEigenvalues(WR, WI, VR, dim, eigenModesNumber, 0);

  /*--- Normalization of eigenvectors... */
  NormalizeEigenvectors(dim, VR, 0);

  /*--- Assembling diagonalizing matrix V ---*/
  for (k=0; k<eigenModesNumber; k++) {
    for (i=0; i<dim; i++) {
      ReV[k*2*dim+i] = VR[k*2*dim+i];
      ReV[(k*2+1)*dim+i] = -VR[(k*2+1)*dim+i];
    }
  }
  for (k=0; k<eigenModesNumber; k++) {
    for (i=0; i<dim; i++) {
      ImV[2*k*dim+i] = VR[(2*k+1)*dim+i];
      ImV[(2*k+1)*dim+i] = -VR[2*k*dim+i];
    }
  }

  /* Copy the sigma matrix into the working buffer */
  for (i=0; i<MATDIM; i++)
    for (j=0; j<MATDIM; j++)
      Mcopy[i*MATDIM+j] = sigmaMatrix[sigmaIndex3[i][j]];

  M3 = MatrixProduct1(dim,dim,ReV,dim,dim,Mcopy);
  M4 = TransposeM(dim,dim,ReV);
  M1 = MatrixProduct1(dim,dim,M3,dim,dim,M4);
  free(M3);
  free(M4);

  M3 = MatrixProduct1(dim,dim,ImV,dim,dim,Mcopy); 
  M4 = TransposeM(dim,dim,ImV);
  M2 = MatrixProduct1(dim,dim,M3,dim,dim,M4);
  free(M3);
  free(M4);

  Rdiag=AddMM1(1,dim,dim,M1,M2);
  free(M1);
  free(M2);

  emittance[0] = Rdiag[0];
  emittance[1] = Rdiag[2*dim+2];
  emittance[2] = Rdiag[4*dim+4];
  free(Rdiag);
}

void NormalizeEigenvectors(int dim, double *V, int debug)
{
  int k, i, eigenModesNumber;
  double Norm[3], Vcopy[MATDIM2];
  eigenModesNumber = 3;

  for (i=0; i<MATDIM2; i++) Vcopy[i] = V[i];
  for (k=0; k<eigenModesNumber; k++) {
    Norm[k]=0;
    for (i=0; i<eigenModesNumber; i++) {
      /* Index = Irow*dim + Icolumn */
      Norm[k]+=(V[2*k*dim+2*i+1]*V[(2*k+1)*dim+2*i]-V[2*k*dim+2*i]*V[(2*k+1)*dim+2*i+1])*2;
    }
    Norm[k]=-1.0/sqrt(fabs(Norm[k]));
    if (debug == 4) printf("Normalization coefficient[%d]= %12.4e \n",k,Norm[k]);
  }
  for (k=0; k<eigenModesNumber; k++) {
    for (i=0; i<dim*2; i++) {
      V[k*2*dim+i]=Vcopy[k*2*dim+i]*Norm[k];
    }
  }
}

double *AddMM1(int sum, int rows, int cols, double *M1, double *M2)
/* Adds two matrices */
{
  int i, j;
  double *M3;
  M3 = malloc(sizeof(*M3)*cols*rows);
  for (i=0; i<rows; i++) {
    for (j=0; j<cols; j++) {
      if (sum == 1) {
	M3[i*cols+j] = M1[i*cols+j] + M2[i*cols+j];
      } else {
	M3[i*cols+j] = M1[i*cols+j] - M2[i*cols+j];
      }
    }
  }
  return M3;
}

double *MatrixProduct1 (int rows1, int cols1, double *T1, int rows2, int cols2, double *T2)
/* Calculates T3=T1*T2 */
{
  double *T3;
  int i, j, k;
  T3 = malloc(sizeof(*T3)*rows1*cols2);
  if (cols1 != rows2) {
    printf("Wrong matrix dimension!\n");
    exit(1);
  }
  for (i=0; i<rows1; i++) {
    for (j=0; j<cols2; j++) {
      T3[i*cols2+j]=0;
      for (k=0; k<cols1; k++) {
        T3[i*cols2+j]+=T1[i*cols1+k]*T2[k*cols2+j];
      }
    }
  }
  return T3;
}

double *TransposeM(int rows, int cols, double *M)
{
  int i, j;
  double *Mt;
  Mt = malloc(sizeof(*Mt)*cols*rows);
  for (i=0; i<rows; i++) 
    for (j=0; j<cols; j++)
      Mt[j*cols+i] = M[i*cols+j];
  return Mt;
}

void MatrixPrintout1(char *string, double *AA, int Nrow, int Ncol)
{
  int i, j;
  printf("%s\n", string);
  for (i=0; i<Nrow; i++) {
    for (j=0; j<Ncol; j++) {
      printf("%16.8e", AA[i*Ncol+j]);
    }
    printf("\n");
  }
  printf("\n");
}

