/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: chrom.c
 * purpose: chromaticity correction by adjustment of sextupoles.
 *
 * Michael Borland, 1992
 */
#if defined(SOLARIS) && !defined(__GNUC__)
#include <sunmath.h>
#endif
#include "mdb.h"
#include "track.h"
#include "chromDefs.h"

double computeChromaticityValue(double dR11, double dR12, double dR22, double R11, double R12, double R22,
				double beta0, double beta1, double alpha0, double alpha1,
				double phi1, short periodic);
double computeChromaticBetaValue(double dR11, double dR12, double dR22, double R11, double R12, double R22,
				 double beta0, double beta1, double alpha0, double alpha1,
				 double phi1, double chrom, short periodic);
double computeChromaticAlphaValue(double dR11, double dR12, double dR22, double R11, double R12, double R22,
				 double beta0, double beta1, double alpha0, double alpha1,
				 double phi1, double chrom, short periodic);

static FILE *fp_sl = NULL;
static long alter_defined_values;
static long verbosityLevel = 2;

void setup_chromaticity_correction(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom)
{
    VMATRIX *M;
    ELEMENT_LIST *eptr, *elast;
#include "chrom.h"
    unsigned long unstable;
    
    log_entry("setup_chromaticity_correction");

    cp_str(&sextupoles, "sf sd");

    /* process namelist input */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&chromaticity, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    str_toupper(sextupoles);
    if (echoNamelists) print_namelist(stdout, &chromaticity);

    if (run->default_order<2)
        bombElegant("default order must be >= 2 for chromaticity correction", NULL);

    if (chrom->name)
        tfree(chrom->name);
    chrom->name = tmalloc(sizeof(*chrom->name)*(chrom->n_families=1));
    while ((chrom->name[chrom->n_families-1]=get_token(sextupoles)))
        chrom->name = trealloc(chrom->name, sizeof(*chrom->name)*(chrom->n_families+=1));
    if ((--chrom->n_families)<1)
        bombElegant("too few sextupoles given for chromaticity correction", NULL);
    chrom->chromx = dnux_dp;
    chrom->chromy = dnuy_dp;
    chrom->n_iterations = n_iterations;
    chrom->correction_fraction = correction_fraction;
    alter_defined_values = change_defined_values;
    chrom->strengthLimit = strength_limit;
    chrom->use_perturbed_matrix = use_perturbed_matrix;
    chrom->sextupole_tweek = sextupole_tweek;
    chrom->tolerance = tolerance;
    verbosityLevel = verbosity;
    chrom->exit_on_failure = exit_on_failure;
    
    if (!use_perturbed_matrix) {
      if (!beamline->twiss0 || !beamline->matrix) {
        double beta_x, alpha_x, eta_x, etap_x;
        double beta_y, alpha_y, eta_y, etap_y;

        fprintf(stdout, "Computing periodic Twiss parameters.\n");
        fflush(stdout);

        if (!beamline->twiss0)
          beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));

        eptr = beamline->elem_twiss = &(beamline->elem);
        elast = eptr;
        while (eptr) {
          if (eptr->type==T_RECIRC)
            beamline->elem_twiss = beamline->elem_recirc = eptr;
          elast = eptr;
          eptr = eptr->succ;
        }
        if (beamline->links) {
          /* rebaseline_element_links(beamline->links, run, beamline); */
          if (assert_element_links(beamline->links, run, beamline, STATIC_LINK+DYNAMIC_LINK)) {
            beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
            beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
            beamline->flags &= ~BEAMLINE_RADINT_CURRENT;
          }
        }
        
        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, NULL, run, 
                                                      &unstable, NULL, NULL);
        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, beamline->waists,
                                   NULL, beamline->elem_twiss, run, NULL,
				   beamline->couplingFactor);
      }

      if (!(M=beamline->matrix) || !M->C || !M->R || !M->T)
        bombElegant("something wrong with transfer map for beamline (setup_chromaticity_correction)", NULL);

      computeChromCorrectionMatrix(run, beamline, chrom);
    }
    
#if USE_MPI
    if (!writePermitted)
       strength_log = NULL;
#endif
    if (strength_log) {
        strength_log = compose_filename(strength_log, run->rootname);
        fp_sl = fopen_e(strength_log, "w", 0);
        fprintf(fp_sl, "SDDS1\n&column name=Step, type=long, description=\"Simulation step\" &end\n");
        fprintf(fp_sl, "&column name=K2, type=double, units=\"1/m$a2$n\" &end\n");
        fprintf(fp_sl, "&column name=SextupoleName, type=string  &end\n");
        fprintf(fp_sl, "&data mode=ascii, no_row_counts=1 &end\n");
        fflush(fp_sl);
        }

    log_exit("setup_chromaticity_correction");
    
}

void computeChromCorrectionMatrix(RUN *run, LINE_LIST *beamline, CHROM_CORRECTION *chrom)
{    
    VMATRIX *M;
    double chromx, chromy;
    double chromx0, chromy0;
    double K2=0.0, *K2ptr;
    ELEMENT_LIST *context;
    long i, count, K2_param=0;
    MATRIX *C, *Ct, *CtC, *inv_CtC;
    
    m_alloc(&C, 2, chrom->n_families);
    m_alloc(&Ct, chrom->n_families, 2);
    m_alloc(&CtC, chrom->n_families, chrom->n_families);
    m_alloc(&inv_CtC, chrom->n_families, chrom->n_families);

    if (chrom->T)
      m_free(&(chrom->T));
    if (chrom->dK2)
      m_free(&(chrom->dK2));
    if (chrom->dchrom)
      m_free(&(chrom->dchrom));
    m_alloc(&(chrom->T), chrom->n_families, 2);
    m_alloc(&(chrom->dK2), chrom->n_families, 1);
    m_alloc(&(chrom->dchrom), 2, 1);

    if (verbosityLevel>2) {
      fprintf(stdout, "Computing chromaticity influence matrix for all named sextupoles.\n");
      fflush(stdout);
    }

    computeChromaticities(&chromx0, &chromy0, 
                          NULL, NULL, NULL, NULL, 
			  beamline->twiss0, beamline->elast->twiss, M=beamline->matrix);
    M = NULL;
    for (i=0; i<chrom->n_families; i++) {
        count = 0;
        context = NULL;
        while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stdout, "error: element %s does not have K2 parameter\n", 
                        context->name);
                fflush(stdout);
                exitElegant(1);
                }
            if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                bombElegant("K2ptr NULL in setup_chromaticity_correction", NULL);
            if (count==0)
                K2 = *K2ptr;
            *K2ptr = K2 + chrom->sextupole_tweek;
            if (context->matrix) {
                free_matrices(context->matrix);
                free(context->matrix);
                }       
            compute_matrix(context, run, NULL);
            count++;
            }
        if (count==0) {
          fprintf(stdout, "error: element %s is not in the beamline.\n",  chrom->name[i]);
          fflush(stdout);
          exitElegant(1);
        }
        if (beamline->links) {
          /* rebaseline_element_links(beamline->links, run, beamline); */
          assert_element_links(beamline->links, run, beamline, STATIC_LINK+DYNAMIC_LINK);
        }
        if (M) {
          free_matrices(M);
          free(M);
        }
        M = full_matrix(beamline->elem_twiss, run, 2);
        computeChromaticities(&chromx, &chromy, 
                              NULL, NULL, NULL, NULL, beamline->twiss0, beamline->elast->twiss, M);

        C->a[0][i] = (chromx-chromx0)/chrom->sextupole_tweek;
        C->a[1][i] = (chromy-chromy0)/chrom->sextupole_tweek;
        if (C->a[0][i]==0 || C->a[1][i]==0) {
            fprintf(stdout, "error: element %s does not change the chromaticity!\n", chrom->name[i]);
            fflush(stdout);
            exitElegant(1);
            }
        count = 0;
        context = NULL;
        while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
            if (!count && !(K2_param=confirm_parameter("K2", context->type))) {
                fprintf(stdout, "error: element %s does not have K2 parameter\n", 
                        context->name);
                fflush(stdout);
                exitElegant(1);
                }
            if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                bombElegant("K2ptr NULL in setup_chromaticity_correction", NULL);
            if (!K2ptr)
                bombElegant("K2ptr NULL in setup_chromaticity_correction", NULL);
            *K2ptr = K2;
            if (context->matrix) {
                free_matrices(context->matrix);
                free(context->matrix);
                }
            compute_matrix(context, run, NULL);
            count++;
            }
        }
    if (M) {
      free_matrices(M);
      free(M);
    }
    if (beamline->matrix)
      free_matrices(beamline->matrix);
    beamline->matrix = full_matrix(beamline->elem_twiss, run, run->default_order);

    if (verbosityLevel>1) {
      fprintf(stdout, "\nfamily           dCHROMx/dK2        dCHROMy/dK2\n");
      fflush(stdout);
      for (i=0; i<chrom->n_families; i++)
        fprintf(stdout, "%10s:    %14.7e     %14.7e\n", chrom->name[i], C->a[0][i], C->a[1][i]);
      fflush(stdout);
    }
    
    m_trans(Ct, C);
    m_mult(CtC, Ct, C);
    m_invert(inv_CtC, CtC);
    m_mult(chrom->T, inv_CtC, Ct);

    if (verbosityLevel>1) {
      fprintf(stdout, "\nfamily           dK2/dCHROMx        dK2/dCHROMy\n");
      fflush(stdout);
      for (i=0; i<chrom->n_families; i++)
        fprintf(stdout, "%10s:    %14.7e     %14.7e\n", chrom->name[i], chrom->T->a[i][0], chrom->T->a[i][1]);
      fflush(stdout);
      fprintf(stdout, "\n");
      fflush(stdout);
    }
    
    m_free(&C);
    m_free(&Ct);
    m_free(&CtC);
    m_free(&inv_CtC);
  }



long do_chromaticity_correction(CHROM_CORRECTION *chrom, RUN *run, LINE_LIST *beamline,
            double *clorb, long step, long last_iteration)
{
    VMATRIX *M;
    double chromx0, chromy0, dchromx=0, dchromy=0;
    double K2=0.0, *K2ptr;
    ELEMENT_LIST *context;
    long i, K2_param=0, type=0, iter, count;
    double beta_x, alpha_x, eta_x, etap_x;
    double beta_y, alpha_y, eta_y, etap_y;
    double K2_min, K2_max;
    unsigned long unstable;
    double lastError, presentError;
    char buffer[256];
    
    log_entry("do_chromaticity_correction");

    if (!beamline->matrix || !beamline->twiss0) {
        if (!beamline->twiss0)
            beamline->twiss0 = tmalloc(sizeof(*beamline->twiss0));
        if (!beamline->elem_twiss) {
            ELEMENT_LIST *eptr;
            eptr = beamline->elem_twiss = &(beamline->elem);
            while (eptr) {
                if (eptr->type==T_RECIRC)
                    beamline->elem_twiss = beamline->elem_recirc = eptr;
                eptr = eptr->succ;
                }
            }
        if (beamline->matrix) {
          free_matrices(beamline->matrix);
          free(beamline->matrix);
        }
        beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
						  &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
						  beamline->elem_twiss, clorb, run, &unstable, NULL, NULL);

        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        propagate_twiss_parameters(beamline->twiss0, beamline->tune, beamline->waists,
                                   NULL, beamline->elem_twiss, run, clorb, 
				   beamline->couplingFactor);
        }
    else if (beamline->matrix->order<2) {
      if (beamline->matrix) {
          free_matrices(beamline->matrix);
          free(beamline->matrix);
        }
      beamline->matrix = full_matrix(beamline->elem_twiss, run, 2);
    }

    if (!(M = beamline->matrix) || !M->C || !M->R || !M->T)
        bombElegant("something wrong with transfer map for beamline (do_chromaticity_correction.1)", NULL);

    computeChromaticities(&chromx0, &chromy0, 
			  NULL, NULL, NULL, NULL, beamline->twiss0,
			  beamline->elast->twiss, M);

    if (verbosityLevel>0) {
      fprintf(stdout, "\nAdjusting chromaticities:\n");
      fflush(stdout);
      fprintf(stdout, "initial chromaticities:  %e  %e\n", chromx0, chromy0);
      fflush(stdout);
    }
    
    presentError = DBL_MAX;
    for (iter=0; iter<chrom->n_iterations; iter++) {
        K2_max = -(K2_min = DBL_MAX);
        dchromx = chrom->chromx - chromx0;
        dchromy = chrom->chromy - chromy0;
        if (chrom->tolerance>0 &&
            chrom->tolerance>fabs(dchromx) &&
            chrom->tolerance>fabs(dchromy) )
          break;

        lastError = presentError;
        presentError = sqr(dchromx)+sqr(dchromy);
        if (iter && presentError>lastError) {
          fprintf(stdout, "Error increasing---reducing gain\n");
          fflush(stdout);
          chrom->correction_fraction /= 10;
        }

        if (chrom->use_perturbed_matrix)
          computeChromCorrectionMatrix(run, beamline, chrom);
        chrom->dchrom->a[0][0] = dchromx;
        chrom->dchrom->a[1][0] = dchromy;

        m_mult(chrom->dK2, chrom->T, chrom->dchrom);
        for (i=0; i<chrom->n_families; i++) {
          if (isnan(chrom->correction_fraction*chrom->dK2->a[i][0]) ||
              isinf(chrom->correction_fraction*chrom->dK2->a[i][0]))
            break;
        }
        if (i!=chrom->n_families) {
          fprintf(stdout, "Unable to correct chromaticity---diverged\n");
          fflush(stdout);
          return 0;
        }
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            count = 0;
            while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
                if (count==0 && (K2_param = confirm_parameter("K2", context->type))<0) {
                    fprintf(stdout, "error: element %s doesn't have K2 parameter\n",
                            context->name);
                    fflush(stdout);
                    exitElegant(1);
                    }
                if (!(K2ptr = (double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)))
                    bombElegant("K2ptr NULL in setup_chromaticity_correction", NULL);
                K2 = (*K2ptr += chrom->correction_fraction*chrom->dK2->a[i][0]);
                if (chrom->strengthLimit>0 && chrom->strengthLimit<fabs(K2)) {
                  K2 = *K2ptr = SIGN(K2)*chrom->strengthLimit;
                }
                sprintf(buffer, "%s#%ld.K2", context->name, context->occurence);
                rpn_store(K2, NULL, rpn_create_mem(buffer, 0));
                if (K2>K2_max)
                  K2_max = K2;
                if (K2<K2_min)
                  K2_min = K2;
                if (context->matrix) {
                  free_matrices(context->matrix);
                  free(context->matrix);
                  }
                compute_matrix(context, run, NULL);
                type = context->type;
                count++;
                /* fprintf(stdout, "new value of %s#%ld[K2] is  %.15lg 1/m^3\n", 
                       chrom->name[i], context->occurence, K2);
                   fflush(stdout);
                       */
              }
            if (verbosityLevel>1) {
              fprintf(stdout, "Change for family %ld (%ld sextupoles): %e\n",
                      i, count, chrom->correction_fraction*chrom->dK2->a[i][0]);
              fflush(stdout);
            }
            if (alter_defined_values)
              change_defined_parameter(chrom->name[i], K2_param, type, K2, NULL, LOAD_FLAG_ABSOLUTE);
            }    

        if (beamline->links) {
          /* rebaseline_element_links(beamline->links, run, beamline); */
          assert_element_links(beamline->links, run, beamline, 
                               STATIC_LINK+DYNAMIC_LINK+(alter_defined_values?LINK_ELEMENT_DEFINITION:0));
        }
        
        if (beamline->matrix) {
          free_matrices(beamline->matrix);
          free(beamline->matrix);
        }
        M = beamline->matrix = compute_periodic_twiss(&beta_x, &alpha_x, &eta_x, &etap_x, beamline->tune,
                                                      &beta_y, &alpha_y, &eta_y, &etap_y, beamline->tune+1, 
                                                      beamline->elem_twiss, clorb, run, &unstable, NULL, NULL);
        beamline->twiss0->betax  = beta_x;
        beamline->twiss0->alphax = alpha_x;
        beamline->twiss0->phix   = 0;
        beamline->twiss0->etax   = eta_x;
        beamline->twiss0->etapx  = etap_x;
        beamline->twiss0->betay  = beta_y;
        beamline->twiss0->alphay = alpha_y;
        beamline->twiss0->phiy   = 0;
        beamline->twiss0->etay   = eta_y;
        beamline->twiss0->etapy  = etap_y;
        
        if (!M || !M->C || !M->R || !M->T)
            bombElegant("something wrong with transfer map for beamline (do_chromaticity_correction.2)", NULL);
        computeChromaticities(&chromx0, &chromy0, 
			      NULL, NULL, NULL, NULL, beamline->twiss0, 
			      beamline->elast->twiss, M);
        beamline->chromaticity[0] = chromx0;
        beamline->chromaticity[1] = chromy0;
        if (verbosityLevel>0) {
          fprintf(stdout, "resulting chromaticities:  %e  %e\n", chromx0, chromy0);
          fprintf(stdout, "min, max sextupole strength:  %e  %e  1/m^2\n", K2_min, K2_max);
          fflush(stdout);
        }
      }
    
    if (fp_sl && last_iteration) {
        for (i=0; i<chrom->n_families; i++) {
            context = NULL;
            while ((context=find_element(chrom->name[i], &context, beamline->elem_twiss))) {
              if (( K2_param = confirm_parameter("K2", context->type))<0)
                bombElegant("confirm_parameter doesn't return offset for K2 parameter.\n", NULL);
              fprintf(fp_sl, "%ld %e %s\n", 
                      step,
                      *((double*)(context->p_elem + entity_description[context->type].parameter[K2_param].offset)),
                      chrom->name[i]);
            }
          }
        fflush(fp_sl);
      }
    

    propagate_twiss_parameters(beamline->twiss0, beamline->tune, 
                               beamline->waists, NULL, beamline->elem_twiss, run, clorb,
			       beamline->couplingFactor);
    log_exit("do_chromaticity_correction");

    if (chrom->tolerance>0 &&
        (chrom->tolerance<fabs(dchromx) || chrom->tolerance<fabs(dchromy))
        && chrom->exit_on_failure) {
      fprintf(stdout, "Chromaticity correction failure---exiting!\n");
      exitElegant(1);
    }

    return 1;
  }

void computeChromaticities(double *chromx, double *chromy, 
                           double *dbetax, double *dbetay,
                           double *dalphax, double *dalphay,
                           TWISS *twiss0, TWISS *twiss1, VMATRIX *M)
{
  double computeChromaticDerivRElem(long i, long j, TWISS *twiss, VMATRIX *M);
  double dR11, dR22, dR12, dR33, dR34, dR44;

  dR11 = computeChromaticDerivRElem(1, 1, twiss1, M);
  dR12 = computeChromaticDerivRElem(1, 2, twiss1, M);
  dR22 = computeChromaticDerivRElem(2, 2, twiss1, M);
  dR33 = computeChromaticDerivRElem(3, 3, twiss1, M);
  dR34 = computeChromaticDerivRElem(3, 4, twiss1, M);
  dR44 = computeChromaticDerivRElem(4, 4, twiss1, M);

  *chromx = computeChromaticityValue(dR11, dR12, dR22, 
				     M->R[0][0], M->R[0][1], M->R[1][1],
				     twiss0->betax, twiss1->betax,
				     twiss0->alphax, twiss1->alphax,
				     twiss1->phix, twiss1->periodic);
  *chromy = computeChromaticityValue(dR33, dR34, dR44, 
				     M->R[2][2], M->R[2][3], M->R[3][3],
				     twiss0->betay, twiss1->betay,
				     twiss0->alphay, twiss1->alphay,
				     twiss1->phiy, twiss1->periodic);

  if (dbetax) {
    *dbetax = computeChromaticBetaValue(dR11, dR12, dR22, 
					M->R[0][0], M->R[0][1], M->R[1][1],
					twiss0->betax, twiss1->betax,
					twiss0->alphax, twiss1->alphax,
					twiss1->phix, *chromx, twiss1->periodic);
#ifdef DEBUG
    fprintf(stdout, "dbetax/dp = %e\n", *dbetax);
#endif
  }
  if (dalphax) {
    *dalphax = computeChromaticBetaValue(dR11, dR12, dR22, 
					M->R[0][0], M->R[0][1], M->R[1][1],
					twiss0->betax, twiss1->betax,
					twiss0->alphax, twiss1->alphax,
					twiss1->phix, *chromx, twiss1->periodic);
#ifdef DEBUG
    fprintf(stdout, "dalphax/dp = %e\n", *dalphax);
#endif
  }
  if (dbetay) {
    *dbetay = computeChromaticBetaValue(dR33, dR34, dR44,
					M->R[2][2], M->R[2][3], M->R[3][3],
					twiss0->betay, twiss1->betay,
					twiss0->alphay, twiss1->alphay,
					twiss1->phiy, *chromy, twiss1->periodic);
#ifdef DEBUG
    fprintf(stdout, "dbetay/dp = %e\n", *dbetay);
#endif
  }
  if (dalphay) {
    *dalphay = computeChromaticAlphaValue(dR33, dR34, dR44,
					M->R[2][2], M->R[2][3], M->R[3][3],
					twiss0->betay, twiss1->betay,
					twiss0->alphay, twiss1->alphay,
					twiss1->phiy, *chromy, twiss1->periodic);
#ifdef DEBUG
    fprintf(stdout, "dalphay/dp = %e\n", *dalphay);
#endif
  }
}

void computeHigherOrderChromaticities(LINE_LIST *beamline, double *clorb, RUN *run,
				      long concatOrder, double deltaStep, long deltaPoints,
                                      long quickMode)
{
#define MAX_NDELTA_VALUES 101   /* must be at least 5 */
  double trace[2][MAX_NDELTA_VALUES], delta[MAX_NDELTA_VALUES];
  double eta[6];
  double coef[MAX_NDELTA_VALUES], sCoef[MAX_NDELTA_VALUES], chi;
  long i, p;
  double c1;
  VMATRIX M1, M0, *Mp;

  if (deltaPoints>MAX_NDELTA_VALUES)
    bombElegant("too many points for higher-order chromaticity", NULL);
  if (deltaPoints<5)
    deltaPoints = 5;
  if (!(beamline->matrix))
    bombElegant("no matrix for beamline (computeHigherOrderChromaticities)", NULL);
  
  beamline->chrom2[0] = beamline->chrom2[1] = 0;
  beamline->chrom3[0] = beamline->chrom3[1] = 0;
  initialize_matrices(&M0, 1);
  initialize_matrices(&M1, concatOrder);
  
  eta[0] = beamline->twiss0->etax;
  eta[1] = beamline->twiss0->etapx;
  eta[2] = beamline->twiss0->etay;
  eta[3] = beamline->twiss0->etapy;
  eta[4] = 0;
  eta[5] = 1;
  for (p=0; p<deltaPoints; p++) {
    delta[p] = (p-(deltaPoints/2+1))*deltaStep;
    for (i=0; i<6; i++) {
      M0.C[i] = (clorb?clorb[i]:0)+delta[p]*eta[i] +
	(i<4? sqr(delta[p])*beamline->eta2[i] + pow3(delta[p])*beamline->eta3[i] : 0);
      M0.R[i][i] = 1;
    }
    if (quickMode)
      concat_matrices(Mp=&M1, beamline->matrix, &M0, 0);
    else
      Mp = append_full_matrix(beamline->elem_twiss, run, &M0, concatOrder);
    for (i=0; i<2; i++)
      /* Tr[0,1][p] is the trace for x,y plane for point p */
      trace[i][p] = Mp->R[2*i][2*i] + Mp->R[2*i+1][2*i+1];
    if (!quickMode) {
      free_matrices(Mp);
      free(Mp);
    }
  }
  for (i=0; i<2; i++) {
    lsfn(delta, trace[i], NULL, deltaPoints, (long)(MIN(deltaPoints-2,5)), 
         coef, sCoef, &chi, NULL);
    c1 = -coef[1]/(4*PI*sin(PIx2*beamline->tune[i]));
    beamline->chrom2[i] = (
                           2*coef[2]
                           + 8*sqr(PI)*sqr(c1)*cos(PIx2*beamline->tune[i])
                           )/(-4*PI*sin(PIx2*beamline->tune[i]));
    if (quickMode)
      beamline->chrom3[i] = 0;
    else
      beamline->chrom3[i] = (
                             6*coef[3] 
                             - 16*pow3(PI*c1)*sin(PIx2*beamline->tune[i]) 
                             + 24*sqr(PI)*c1*beamline->chrom2[i]*cos(PIx2*beamline->tune[i])
                             )/(-4*PI*sin(PIx2*beamline->tune[i]));
  }
  free_matrices(&M0);
  free_matrices(&M1);
}

double computeChromaticDerivRElem(long i, long j, TWISS *twiss, VMATRIX *M)
{
  long k;
  double sum, eta[6] = {0,0,0,0,0,0};

  if (!(M->T))
    return 0.0;
  
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  i--;
  j--;
  for (k=sum=0; k<6; k++) {
    if (k>j) {
      sum += eta[k]*M->T[i][k][j];
    } else if (k==j) {
      sum += eta[k]*M->T[i][k][k]*2;
    } else {
      sum += eta[k]*M->T[i][j][k];
    }
  }
#ifdef DEBUG
  fprintf(stdout, "dR%ld%ld/ddelta = %e\n",
          i+1, j+1, sum);
#endif
  return sum;
}  

double computeChromaticDeriv2RElem(long i, long m, TWISS *twiss, VMATRIX *M,
                                   double *eta2q)
{
  long j, k, l, count;
  double sum, eta[6] = {0,0,0,0,0,0};
  double eta2[6] = {0,0,0,0,0,0};
   
  fprintf(stdout, "Computing d2R%ld%ld:\n", i, m);
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  for (j=0; j<4; j++)
    eta2[j] = eta2q[j];

  i--;
  m--;
  sum = 0;
  count = 0;
  if (M->T)
    for (k=0; k<6; k++) {
      if (k>m) {
        sum += eta2[k]*M->T[i][k][m];
      } else if (k==m) {
        sum += eta2[k]*M->T[i][k][k]*2;
      } else {
        sum += eta2[k]*M->T[i][m][k];
      }
      count ++;
    }
  if (M->Q) {
    for (j=0; j<6; j++)
      for (k=0; k<=j; k++)
	for (l=0; l<=k; l++) {
	  if (j==m) {
	    sum += M->Q[i][m][k][l]*eta[k]*eta[l];
	    count ++;
	  }
	  if (k==m) {
	    sum += M->Q[i][j][m][l]*eta[j]*eta[l];
	    count ++;
	  }
	  if (l==m) {
	    sum += M->Q[i][j][k][l]*eta[j]*eta[k];
	    count ++;
	  }
	}
  }
  fprintf(stdout, "sum is %e from %ld terms\n",
	  sum, count);
  return sum;
}  

double computeChromaticDeriv3RElem(long i, long m, TWISS *twiss, VMATRIX *M,
                                   double *eta2q, double *eta3q)
{
  long j, k, l, count;
  double sum, eta[6] = {0,0,0,0,0,0};
  double eta2[6] = {0,0,0,0,0,0};
  double eta3[6] = {0,0,0,0,0,0};
   
  fprintf(stdout, "Computing d3R%ld%ld:\n", i, m);
  eta[0] = twiss->etax;
  eta[1] = twiss->etapx;
  eta[2] = twiss->etay;
  eta[3] = twiss->etapy;
  eta[5] = 1;
  for (j=0; j<4; j++) {
    eta2[j] = eta2q[j];
    eta3[j] = eta3q[j];
  }
  for (j=0; j<6; j++) {
    fprintf(stdout, "eta[%ld] = %e, %e, ", j, eta[j], eta2[j]);
    fprintf(stdout, "%e\n", eta3[j]);
  }
  i--;
  m--;
  sum = count = 0;
  if (M->T)
    for (k=0; k<6; k++) {
      if (k>m) {
	count ++;
        sum += eta3[k]*M->T[i][k][m];
	fprintf(stdout, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, k, m, eta3[k], M->T[i][k][m]);
      } else if (k==m) {
	count ++;
        sum += eta3[k]*M->T[i][k][k]*2;
	fprintf(stdout, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, k, k, eta3[k], M->T[i][k][k]);
      } else {
	count ++;
        sum += eta3[k]*M->T[i][m][k];
	fprintf(stdout, "term eta3[%ld]*T[%ld][%ld][%ld] = %e*%e\n",
		k, i, m, k, eta3[k], M->T[i][m][k]);
      }
    }
  if (M->Q) {
    for (j=0; j<6; j++)
      for (k=0; k<=j; k++)
	for (l=0; l<=k; l++) {
	  if (j==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[k]*eta[l]+eta[k]*eta2[l]);
	    fprintf(stdout, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[k]*eta[l]+eta[k]*eta2[l]), 
		    M->Q[i][j][k][l]);
	  }
	  if (k==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[j]*eta[l]+eta[j]*eta2[l]);
	    fprintf(stdout, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[j]*eta[l]+eta[j]*eta2[l]),
		    M->Q[i][j][k][l]);
	  }
	  if (l==m) {
	    count ++;
	    sum += M->Q[i][j][k][l]*(eta2[j]*eta[k]+eta[j]*eta2[k]);
	    fprintf(stdout, "term eta*U[%ld][%ld][%ld][%ld] = %e*%e\n",
		    i, j, k, l, (eta2[j]*eta[k]+eta[j]*eta2[k]),
		    M->Q[i][j][k][l]);
	  }
	}
  }
  fprintf(stdout, "sum of %ld terms is %e\n", count, sum);
  return sum;
}  

void computeChromaticTuneLimits(LINE_LIST *beamline)
{
  long i, j, n, p;
  double c1, c2, c3, tuneValue[5], solution[2];
  
  for (i=0; i<2; i++) {
    if (beamline->chromDeltaHalfRange<=0) {
      beamline->tuneChromUpper[i] = beamline->tuneChromLower[i] = beamline->tune[i];
    } else {
      tuneValue[0] = beamline->tune[i];
      /* polynomial coefficients */
      c1 = beamline->chromaticity[i];
      c2 = beamline->chrom2[i]/2.0;
      c3 = beamline->chrom3[i]/6.0;
      tuneValue[1] = beamline->tune[i] + 
        beamline->chromDeltaHalfRange*c1 + 
          sqr(beamline->chromDeltaHalfRange)*c2 +
            ipow(beamline->chromDeltaHalfRange, 3)*c3;
      tuneValue[2] = beamline->tune[i] -
        beamline->chromDeltaHalfRange*c1 + 
          sqr(beamline->chromDeltaHalfRange)*c2 -
            ipow(beamline->chromDeltaHalfRange, 3)*c3;
      p = 3;
      /* find extrema */
      n = solveQuadratic(3*c3, 2*c2, c1, solution);
      for (j=0; j<n; j++) {
        if (fabs(solution[j])>beamline->chromDeltaHalfRange)
          continue;
        tuneValue[p] = beamline->tune[i] +
          solution[j]*c1 + sqr(solution[j])*c2 +
            ipow(solution[j], 3)*c3;
        p += 1;
      }
      find_min_max(beamline->tuneChromLower+i, beamline->tuneChromUpper+i, 
                   tuneValue, p);
    }
  }
}

double computeChromaticityValue(double dR11, double dR12, double dR22,
				double R11, double R12, double R22,
				double beta0, double beta1,
				double alpha0, double alpha1,
				double phi1, short periodic)
{
  if (periodic) {
    if (R12==0)
      return DBL_MAX;
    return -(dR11+dR22)/R12*beta0/(2*PIx2);
  }
  else
    return 
      ((dR12*(cos(phi1)+alpha0*sin(phi1)))/sqrt(beta0*beta1) 
       - dR11*sin(phi1)*sqrt(beta0/beta1))/(PIx2) ;
}

double computeChromaticBetaValue(double dR11, double dR12, double dR22,
				 double R11, double R12, double R22,
				 double beta0, double beta1,
				 double alpha0, double alpha1,
				 double phi1, double chrom, short periodic)
{
  if (periodic) {
    if (R12==0)
      return DBL_MAX;
    return beta0*(dR12 - PI*beta0*chrom*(R11+R22))/R12;
  } else 
    return 
      2*(
	 dR11*cos(phi1)*sqrt(beta0*beta1) +
	 dR12*sqrt(beta1/beta0)*(sin(phi1) - cos(phi1)*alpha0)
       );
}

double computeChromaticAlphaValue(double dR11, double dR12, double dR22,
				 double R11, double R12, double R22,
				 double beta0, double beta1,
				 double alpha0, double alpha1,
				 double phi1, double chrom, short periodic)
{
  if (periodic) {
    if (R12==0)
      return DBL_MAX;
    return beta0/(2*R12)*(-alpha0*(R11+R22)*PIx2*chrom + dR11-dR22);
  } else {
    if (sin(phi1)==0)
      return DBL_MAX;
    return 
      /* this ugly, inefficient expression was created by Mathematica(TM) */
      -((-(sqrt(beta0*beta1)*dR12*cos(phi1)) 
	 + sqr(beta0)*sqrt(beta1/beta0)*dR11*sin(phi1) 
	 - alpha0*sqrt(beta0*beta1)*dR12*sin(phi1))/
	(beta0*beta1));
  }
}
