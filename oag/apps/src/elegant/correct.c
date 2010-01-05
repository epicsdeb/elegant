/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/


/* file: correct.c
 * purpose: trajectory/orbit correction for elegant
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "match_string.h"
#include "correctDefs.h"

#define N_CORRECTION_MODES 2
char *correction_mode[N_CORRECTION_MODES] = {
    "trajectory", "orbit"
    } ;

#define GLOBAL_CORRECTION 0
#define ONE_TO_ONE_CORRECTION 1
#define THREAD_CORRECTION 2
#define ONE_TO_BEST_CORRECTION 3
#define ONE_TO_NEXT_CORRECTION 4
#define N_CORRECTION_METHODS 5
char *correction_method[N_CORRECTION_METHODS] = {
    "global", "one-to-one", "thread", "one-to-best", "one-to-next",
    } ;

/* For trajectory correction:
 *   C      :     NMON x NCOR matrix of dQ_mon/dkick_cor 
 *   Qo     :     NMON x 1    matrix of initial monitor positions, Q stands for either x or y 
 *   Q      :     NMON x 1    matrix of final monitor positions
 *   dK     :     NCOR x 1    matrix of changes to correctors 
 *
 *   Q = Qo + C*dK
 *   Want to minimize S = SUM((Qi*wi)^2), where wi is the weight for the ith monitor.
 *   Solution is dK = -Inverse(Trans(C) W C) Trans(C) W Qo, or dK = T Qo.
 *   W is NMON x NMON, and W(i,j) = delta(i,j)*wi.  A monitor with a zero or negative weight is
 *       ignored.
 *   T is NCOR x NMON, and is computed before any errors are added and before any parameters are varied
 */

long global_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, 
                          RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam);
void one_to_one_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run, 
            LINE_LIST *beamline, double *starting_coord, BEAM *beam, long method);
void thread_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run, 
            LINE_LIST *beamline, double *starting_coord, BEAM *beam);
long orbcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **orbit, 
                  long n_iterations, double accuracy,
                  long clorb_iter, double clorb_iter_frac,
                  RUN *run, LINE_LIST *beamline, double *closed_orbit, double *Cdp);
ELEMENT_LIST *find_useable_moni_corr(int32_t *nmon, int32_t *ncor, long **mon_index,
            ELEMENT_LIST ***umoni, ELEMENT_LIST ***ucorr, double **kick_coef, long **sl_index,
            long plane, STEERING_LIST *SL, RUN *run, LINE_LIST *beamline, long recircs);
ELEMENT_LIST *next_element_of_type(ELEMENT_LIST *elem, long type);
ELEMENT_LIST *next_element_of_type2(ELEMENT_LIST *elem, long type1, long type2);
ELEMENT_LIST *next_element_of_types(ELEMENT_LIST *elem, long *type, long n_types, long *index);
long find_parameter_offset(char *param_name, long elem_type);
long zero_correctors_one_plane(ELEMENT_LIST *elem, RUN *run, STEERING_LIST *SL, long plane);
long zero_correctors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
long zero_hcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
long zero_vcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct);
double rms_value(double *data, long n_data);
long steering_corrector(ELEMENT_LIST *eptr, STEERING_LIST *SL, long plane);
void zero_closed_orbit(TRAJECTORY *clorb, long n);
long find_index(long key, long *list, long n_listed);
long add_steer_elem_to_lists(STEERING_LIST *SL, long plane, char *name, char *item, 
                             char *element_type, double tweek, double limit,
                             LINE_LIST *beamline, RUN *run, long forceQuads);
long add_steer_type_to_lists(STEERING_LIST *SL, long plane, long type, char *item, double tweek, double limit,
    LINE_LIST *beamline, RUN *run, long forceQuads);
double compute_kick_coefficient(ELEMENT_LIST *elem, long plane, long type, double corr_tweek, char *name, char *item, RUN *run);
double noise_value(double xamplitude, double xcutoff, long xerror_type);
void do_response_matrix_output(char *filename, char *type, RUN *run, char *beamline_name, CORMON_DATA *CM, 
                               STEERING_LIST *SL, long plane);


static long rpn_x_mem= -1, rpn_y_mem= -1;
static long usePerturbedMatrix = 0, fixedLengthMatrix = 0;

double getMonitorWeight(ELEMENT_LIST *elem);
double getMonitorCalibration(ELEMENT_LIST *elem, long coord);
double getCorrectorCalibration(ELEMENT_LIST *elem, long coord);

#define UNIFORM_ERRORS 0
#define GAUSSIAN_ERRORS 1
#define PLUS_OR_MINUS_ERRORS 2
#define N_ERROR_TYPES 3
static char *known_error_type[N_ERROR_TYPES] = {
    "uniform", "gaussian", "plus_or_minus"
    };


void correction_setup(
    CORRECTION *_correct,     /* the underscore is to avoid conflicts with the namelist "correct" */
    NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline
    )
{
#include "correct.h"
    char *item;

    log_entry("correction_setup");

    if (_correct->traj) {
        free_czarray_2d((void**)_correct->traj, 3, beamline->n_elems+1);
        _correct->traj = NULL;
        }

    if (_correct->CMx) {
        matrix_free(_correct->CMx->T);
        matrix_free(_correct->CMx->C);
        free_czarray_2d((void**)_correct->CMx->kick, _correct->n_iterations+1, _correct->CMx->ncor);
        free_czarray_2d((void**)_correct->CMx->posi, _correct->n_iterations+1, _correct->CMx->nmon);
        tfree(_correct->CMx->mon_index); _correct->CMx->mon_index = NULL;
        tfree(_correct->CMx); _correct->CMx = NULL;
        }

    if (_correct->CMy) {
        matrix_free(_correct->CMy->T);
        matrix_free(_correct->CMy->C);
        free_czarray_2d((void**)_correct->CMy->kick, _correct->n_iterations+1, _correct->CMy->ncor);
        free_czarray_2d((void**)_correct->CMy->posi, _correct->n_iterations+1, _correct->CMy->nmon);
        tfree(_correct->CMy->mon_index); _correct->CMy->mon_index = NULL;
        tfree(_correct->CMy); _correct->CMy = NULL;
        }

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&correct, nltext);
#if USE_MPI 
    if (isSlave) {
       trajectory_output = NULL;
       corrector_output = NULL;
       statistics = NULL;
    }
#endif
    if (echoNamelists) print_namelist(stdout, &correct);
    usePerturbedMatrix = use_perturbed_matrix;
    fixedLengthMatrix = fixed_length_matrix;
    
    /* check for valid input data */
    if ((_correct->mode=match_string(mode, correction_mode, N_CORRECTION_MODES, 0))<0)
        bomb("invalid correction mode", NULL);
    if ((_correct->method=match_string(method, correction_method, N_CORRECTION_METHODS, 0))<0)
        bomb("invalid correction method", NULL);
    if (corrector_tweek[0]==0 || corrector_tweek[0]==0)
        bomb("invalid corrector tweek(s)", NULL);
    if (corrector_limit[0]<0 || corrector_limit[1]<0 ||
        (corrector_limit[0] && corrector_limit[0]<corrector_tweek[0]) ||
        (corrector_limit[1] && corrector_limit[1]<corrector_tweek[1]) )
        bomb("invalid corrector_limit(s)", NULL);
    if (correction_fraction[0]<=0 || correction_fraction[1]<=0)
        bomb("invalid correction_fraction(s)", NULL);
    if (correction_accuracy[0]<0 || correction_accuracy[1]<0)
        bomb("invalid correction_accuracy(s)", NULL);
    if (trajectory_output)
        setup_orb_traj_output(trajectory_output=compose_filename(trajectory_output, run->rootname),
                    correction_mode[_correct->mode], run);
    if (statistics)
        setup_cormon_stats(statistics=compose_filename(statistics, run->rootname), run);
    if (corrector_output)
        setup_corrector_output(corrector_output=compose_filename(corrector_output, run->rootname), run);
    if (closed_orbit_accuracy<=0)
        bomb("closed_orbit_accuracy must be > 0", NULL);
    if (closed_orbit_iteration_fraction<=0 ||
        closed_orbit_iteration_fraction>1)
      bomb("closed_orbit_iteration_fraction must be on (0, 1]", NULL);
    _correct->clorb_accuracy = closed_orbit_accuracy;
    _correct->clorb_iter_fraction = closed_orbit_iteration_fraction;
    _correct->verbose = verbose;
    _correct->track_before_and_after = track_before_and_after;
    _correct->prezero_correctors = prezero_correctors;
    _correct->start_from_centroid = start_from_centroid;
    _correct->use_actual_beam = use_actual_beam;
    if ((_correct->clorb_iterations=closed_orbit_iterations)<=0)
        bomb("closed_orbit_iterations <= 0", NULL);
    if ((_correct->n_iterations = n_iterations)<0)
        bomb("n_iterations < 0", NULL);
    if ((_correct->n_xy_cycles = n_xy_cycles)<0)
        bomb("n_xy_cycles < 0", NULL);
    _correct->minimum_cycles = minimum_cycles;
    if (threading_divisor[0]<=1 ||  threading_divisor[1]<=1)
      bomb("threading_divisors must be >1", NULL);

    _correct->CMx = tmalloc(sizeof(*_correct->CMx));
    _correct->CMy = tmalloc(sizeof(*_correct->CMy));

    _correct->CMx->default_tweek = corrector_tweek[0];
    _correct->CMy->default_tweek = corrector_tweek[1];
    _correct->CMx->default_threading_divisor = threading_divisor[0];
    _correct->CMy->default_threading_divisor = threading_divisor[1];
    _correct->CMx->corr_limit = corrector_limit[0];
    _correct->CMy->corr_limit = corrector_limit[1];
    _correct->CMx->corr_fraction = correction_fraction[0];
    _correct->CMy->corr_fraction = correction_fraction[1];
    _correct->CMx->corr_accuracy = correction_accuracy[0];
    _correct->CMy->corr_accuracy = correction_accuracy[1];
    _correct->CMx->bpm_noise = bpm_noise[0];
    _correct->CMy->bpm_noise = bpm_noise[1];
    if ((_correct->CMx->bpm_noise_distribution 
         = match_string(bpm_noise_distribution[0], known_error_type, N_ERROR_TYPES, 0))<0)
      bomb("unknown noise distribution type", NULL);
    if ((_correct->CMy->bpm_noise_distribution 
         = match_string(bpm_noise_distribution[1], known_error_type, N_ERROR_TYPES, 0))<0)
      bomb("unknown noise distribution type", NULL);
    _correct->CMx->bpm_noise_cutoff = bpm_noise_cutoff[0];
    _correct->CMy->bpm_noise_cutoff = bpm_noise_cutoff[1];
    _correct->CMx->fixed_length = _correct->CMy->fixed_length = fixed_length;
    _correct->response_only = n_iterations==0;
    _correct->CMx->T = _correct->CMy->T = NULL;
    _correct->CMx->C = _correct->CMy->C = NULL;
    
    _correct->CMx->remove_smallest_SVs = remove_smallest_SVs[0];
    _correct->CMx->auto_limit_SVs = auto_limit_SVs[0];
    _correct->CMx->keep_largest_SVs = keep_largest_SVs[0];
    if ((_correct->CMx->minimum_SV_ratio = minimum_SV_ratio[0])>=1)
      bomb("minimum_SV_ratio should be less than 1 to be meaningful", NULL);

    _correct->CMy->remove_smallest_SVs = remove_smallest_SVs[1];
    _correct->CMy->auto_limit_SVs = auto_limit_SVs[1];
    _correct->CMy->keep_largest_SVs = keep_largest_SVs[1];
    if ((_correct->CMy->minimum_SV_ratio = minimum_SV_ratio[1])>=1) 
      bomb("minimum_SV_ratio should be less than 1 to be meaningful", NULL);

    if (verbose)
      fputs("finding correctors/monitors and/or computing correction matrices\n", stdout);
    
    if (_correct->SLx.n_corr_types==0) {
      long found = 0;
      cp_str(&item, "KICK");
      found += add_steer_type_to_lists(&_correct->SLx, 0, T_HCOR, item, _correct->CMx->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run, 0);
      cp_str(&item, "HKICK");
      found += add_steer_type_to_lists(&_correct->SLx, 0, T_HVCOR, item, _correct->CMx->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run, 0);
      cp_str(&item, "HKICK");
      found += add_steer_type_to_lists(&_correct->SLx, 0, T_QUAD, item, _correct->CMx->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run, 0);
      found += add_steer_type_to_lists(&_correct->SLx, 0, T_KQUAD, item, _correct->CMx->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run, 0);
      if (!found)
        bomb("no horizontal steering elements found", NULL);
    }
    if (_correct->SLy.n_corr_types==0) {
      long found = 0;
      cp_str(&item, "KICK");
      found += add_steer_type_to_lists(&_correct->SLy, 2, T_VCOR, item, _correct->CMy->default_tweek, 
                              _correct->CMx->corr_limit, beamline, run, 0);
      cp_str(&item, "VKICK");
      found += add_steer_type_to_lists(&_correct->SLy, 2, T_HVCOR, item, _correct->CMy->default_tweek, 
                              _correct->CMy->corr_limit, beamline, run, 0);
      cp_str(&item, "VKICK");
      found += add_steer_type_to_lists(&_correct->SLy, 2, T_QUAD, item, _correct->CMy->default_tweek, 
                              _correct->CMy->corr_limit, beamline, run, 0);
      found += add_steer_type_to_lists(&_correct->SLy, 2, T_KQUAD, item, _correct->CMy->default_tweek, 
                              _correct->CMy->corr_limit, beamline, run, 0);
      if (!found)
        bomb("no horizontal steering elements found", NULL);
    }

    if (_correct->mode==TRAJECTORY_CORRECTION) {
      compute_trajcor_matrices(_correct->CMx, &_correct->SLx, 0, run, beamline, 
                               _correct->method==THREAD_CORRECTION,
                               !_correct->response_only && _correct->method==GLOBAL_CORRECTION);
      compute_trajcor_matrices(_correct->CMy, &_correct->SLy, 2, run, beamline, 
                               _correct->method==THREAD_CORRECTION,
                               !_correct->response_only && _correct->method==GLOBAL_CORRECTION);
    }
    else if (_correct->mode==ORBIT_CORRECTION) {
      compute_orbcor_matrices(_correct->CMx, &_correct->SLx, 0, run, beamline, 0, 
                              !_correct->response_only, fixed_length_matrix, verbose);
      compute_orbcor_matrices(_correct->CMy, &_correct->SLy, 2, run, beamline, 0, 
                              !_correct->response_only, fixed_length_matrix, verbose);
    }
    else
      bomb("something impossible happened (correction_setup)", NULL);

    if (n_iterations!=0) {
      /* allocate space to store before/after data for correctors and monitors */
      _correct->CMx->kick = (double**)czarray_2d(sizeof(**(_correct->CMx->kick)), _correct->n_iterations+1, _correct->CMx->ncor);
      _correct->CMx->posi = (double**)czarray_2d(sizeof(**(_correct->CMx->posi)), _correct->n_iterations+1, _correct->CMx->nmon);
      _correct->CMy->kick = (double**)czarray_2d(sizeof(**(_correct->CMy->kick)), _correct->n_iterations+1, _correct->CMy->ncor);
      _correct->CMy->posi = (double**)czarray_2d(sizeof(**(_correct->CMy->posi)), _correct->n_iterations+1, _correct->CMy->nmon);
      
      /* Allocate space to store before/after trajectories/closed orbits.
       * After each correction pass through x and y, the first trajectory is the initial one,
       * the second is after x correction, and the third is after x and y correction.
       */
      _correct->traj = (TRAJECTORY**)czarray_2d(sizeof(**_correct->traj), 3, beamline->n_elems+1);
      
    }
    if (verbose) {
      fprintf(stdout, "there are %" PRId32 " useable horizontal monitors and %" PRId32 " useable horizontal correctors\n",
              _correct->CMx->nmon, _correct->CMx->ncor);
      fflush(stdout);
      fprintf(stdout, "there are %" PRId32 " useable   vertical monitors and %" PRId32 " useable   vertical correctors\n",
              _correct->CMy->nmon, _correct->CMy->ncor);
      fflush(stdout);
    }
    rpn_x_mem = rpn_create_mem("x", 0);
    rpn_y_mem = rpn_create_mem("y", 0);

    log_exit("correction_setup");
  }

void add_steering_element(CORRECTION *correct, LINE_LIST *beamline, RUN *run, NAMELIST_TEXT *nltext)
{
#include "steer_elem.h"

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&steering_element, nltext);
  if (echoNamelists) print_namelist(stdout, &steering_element);

  if (limit && (limit<tweek || limit<0))
    bomb("invalid limit specified for steering element", NULL);

  if (plane[0]=='h' || plane[0]=='H')  {
    if (!add_steer_elem_to_lists(&correct->SLx, 0, name, item, element_type, tweek, limit, beamline, run, 1))
      bomb("no match to given element name or type", NULL);
  }
  else if (plane[0]=='v' || plane[0]=='V') {
    if (!add_steer_elem_to_lists(&correct->SLy, 2, name, item, element_type, tweek, limit, beamline, run, 1))
      bomb("no match to given element name or type", NULL);
  }
  else
    bomb("invalid plane specified for steering element", NULL);
}

long add_steer_type_to_lists(STEERING_LIST *SL, long plane, long type, char *item, double tweek, double limit,
                             LINE_LIST *beamline, RUN *run, long forceQuads)
{
  long found = 0;
  ELEMENT_LIST *context;
  context = &(beamline->elem);
  while (context && (context=next_element_of_type(context, type))) {
    found += add_steer_elem_to_lists(SL, plane, context->name, item, NULL, tweek, limit, beamline, run, forceQuads);
    context = context->succ;
  }
  return found;
}

long add_steer_elem_to_lists(STEERING_LIST *SL, long plane, char *name, char *item, 
                             char *element_type, double tweek, double limit, 
                             LINE_LIST *beamline, RUN *run, long forceQuads)
{
  ELEMENT_LIST *context;
  long param_number, i, found;

  if (SL->n_corr_types==0) {
    if (SL->corr_name)    tfree(SL->corr_name);
    SL->corr_name = NULL;
    if (SL->corr_type)    tfree(SL->corr_type);
    SL->corr_type = NULL;
    if (SL->corr_param)   tfree(SL->corr_param);
    SL->corr_param = NULL;
    if (SL->corr_tweek)   tfree(SL->corr_tweek);
    SL->corr_tweek = NULL;
    if (SL->corr_limit)   tfree(SL->corr_limit);
    SL->corr_limit = NULL;
    if (SL->param_offset) tfree(SL->param_offset);
    SL->param_offset = NULL;
    if (SL->param_index) tfree(SL->param_index);
    SL->param_index = NULL;
  }

  if (!name && !element_type)
    bomb("NULL name and element_type passed to add_steer_elem_to_list", NULL);
  if (name) {
    str_toupper(name);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
  }
  else 
    cp_str(&name, "*");
  if (element_type) {
    str_toupper(element_type);
    if (has_wildcards(element_type) && strchr(element_type, '-'))
      element_type = expand_ranges(element_type);
  }

  if (!item)
    bomb("NULL item passed to add_steer_elem_to_list", NULL);
  str_toupper(item);

  context = NULL;
  found = 0;
  
  while ((context=wfind_element(name, &context, &(beamline->elem)))) {
    if (element_type &&
        !wild_match(entity_name[context->type], element_type))
      continue;

    switch (context->type) {
    case T_QUAD:
      if (plane) {
        if (!((QUAD*)(context->p_elem))->xSteering)
          ((QUAD*)(context->p_elem))->xSteering = forceQuads;
      }
      else {
        if (!((QUAD*)(context->p_elem))->ySteering)
          ((QUAD*)(context->p_elem))->ySteering = forceQuads;
      }
      break;
    case T_KQUAD:
      if (plane) {
        if (!((KQUAD*)(context->p_elem))->xSteering)
          ((KQUAD*)(context->p_elem))->xSteering = forceQuads;
      }
      else {
        if (!((KQUAD*)(context->p_elem))->ySteering)
          ((KQUAD*)(context->p_elem))->ySteering = forceQuads;
      }
      break;
    default:
      break;
    }
    
    for (i=0; i<SL->n_corr_types; i++)
      if (strcmp(context->name, SL->corr_name[i])==0)
        break;
    if (i!=SL->n_corr_types)
      continue;

#ifdef DEBUG
  printf("Adding %s to %c plane steering list\n",
         context->name, plane?'y':'x');
#endif

    SL->corr_name    = trealloc(SL->corr_name, (SL->n_corr_types+1)*sizeof(*SL->corr_name));
    SL->corr_type    = trealloc(SL->corr_type, (SL->n_corr_types+1)*sizeof(*SL->corr_type));
    SL->corr_param   = trealloc(SL->corr_param, (SL->n_corr_types+1)*sizeof(*SL->corr_param));
    SL->corr_tweek   = trealloc(SL->corr_tweek, (SL->n_corr_types+1)*sizeof(*SL->corr_tweek));
    SL->corr_limit   = trealloc(SL->corr_limit, (SL->n_corr_types+1)*sizeof(*SL->corr_limit));
    SL->param_offset = trealloc(SL->param_offset, (SL->n_corr_types+1)*sizeof(*SL->param_offset));
    SL->param_index  = trealloc(SL->param_index , (SL->n_corr_types+1)*sizeof(*SL->param_index ));
    
    SL->corr_type[SL->n_corr_types] = context->type;
    cp_str(SL->corr_name+SL->n_corr_types, context->name);

    cp_str(SL->corr_param+SL->n_corr_types, item);
    SL->corr_tweek[SL->n_corr_types] = tweek;
    SL->corr_limit[SL->n_corr_types] = limit;
    
    if ((SL->param_index[SL->n_corr_types]=param_number=confirm_parameter(item, context->type))<0 ||
        entity_description[context->type].parameter[param_number].type!=IS_DOUBLE ||
        (SL->param_offset[SL->n_corr_types]=find_parameter_offset(item, SL->corr_type[SL->n_corr_types]))<0) {
      fprintf(stderr, "No such floating-point parameter (%s) for %s (add_steer_elem_to_lists)\n", 
              item, context->name);
      exit(1);
    }
    
    SL->n_corr_types += 1;    
    found = 1;
  }
  return found;
}


double compute_kick_coefficient(ELEMENT_LIST *elem, long plane, long type, double corr_tweek, char *name, char *item, RUN *run)
{
  double value, coef;
  long param_offset=0, param_number;

  if ((param_number=confirm_parameter(item, type))<0 || (param_offset=find_parameter_offset(item, type))<0)
    bomb("invalid parameter or element type (compute_kick_coefficient)", NULL);

  switch (type) {
  case T_HCOR:
    if (plane==0 && param_offset==find_parameter_offset("KICK", type))
      return(1.0);
    break;
  case T_VCOR:
    if (plane!=0 && param_offset==find_parameter_offset("KICK", type))
      return(1.0);
    break;
  case T_HVCOR:
    if (plane==0 && param_offset==find_parameter_offset("HKICK", type))
      return(1.0);
    if (plane!=0 && param_offset==find_parameter_offset("VKICK", type))
      return(1.0);
    break;
  default:
    if (entity_description[elem->type].flags&HAS_MATRIX) {
      VMATRIX *M1, *M2, *M0;
      M0 = elem->matrix;
      value = *((double*)(elem->p_elem+param_offset));
      *((double*)(elem->p_elem+param_offset)) += corr_tweek;
      M1 = compute_matrix(elem, run, NULL);
      *((double*)(elem->p_elem+param_offset)) -= 2*corr_tweek;
      M2 = compute_matrix(elem, run, NULL);
      if (plane==0) 
        coef = (M1->C[1]-M2->C[1])/(2*corr_tweek);
      else
        coef = (M1->C[3]-M2->C[3])/(2*corr_tweek);
      /*
        fprintf(stdout, "computed kick coefficient for %s.%s: %g rad/%s\n",
        name, item, coef, entity_description[type].parameter[param_number].unit);
        fflush(stdout);
        */
      free_matrices(M1); tfree(M1); M1 = NULL;
      free_matrices(M2); tfree(M2); M2 = NULL;
      *((double*)(elem->p_elem+param_offset)) = value;
      elem->matrix = M0;
      return(coef);
    }
    else {
      fprintf(stdout, "error: element %s does not have a matrix--can't be used for steering.\n",
              elem->name);
      fflush(stdout);
      exit(1);
    }
    break;
  }
  return(0.0);
}

long do_correction(CORRECTION *correct, RUN *run, LINE_LIST *beamline, double *starting_coord, 
                   BEAM *beam, long sim_step, long initial_correction)
{
  long i, i_cycle, x_failed, y_failed, n_iter_taken, bombed=0, final_traj;
  double *closed_orbit, rms_before, rms_after, *Cdp;

  log_entry("do_correction");

  if (correct->response_only)
    return 1;

  closed_orbit = starting_coord;   /* for return of closed orbit at starting point */

  if (correct->verbose && correct->n_iterations>=1 && correct->n_xy_cycles>0) {
    if (correct->CMx->fixed_length && correct->mode==ORBIT_CORRECTION) {
      fputs(" plane   iter     rms kick     rms pos.    max kick     max pos.   mom.offset\n", stdout);
      fputs("                    mrad          mm         mrad          mm           %\n", stdout);
      fputs("-------  ----    ----------   ---------   ----------   ----------  ----------\n", stdout);
    }
    else {
      fputs(" plane   iter     rms kick     rms pos.    max kick     max pos.\n", stdout);
      fputs("                    mrad          mm         mrad          mm\n", stdout);
      fputs("-------  ----    ----------   ---------   ----------   ----------\n", stdout);
    }
  }

  if (correct->prezero_correctors && initial_correction) {
    long changed;
    if (beamline->elem_recirc)
      changed = zero_correctors(beamline->elem_recirc, run, correct);
    else
      changed = zero_correctors(&(beamline->elem), run, correct);
    if (changed && beamline->matrix) {
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
    }
  }

  correct->CMx->n_cycles_done = correct->CMy->n_cycles_done = 0;
  final_traj = 1;
  switch (correct->mode) {
  case TRAJECTORY_CORRECTION:
    x_failed = y_failed = 0;
    if (usePerturbedMatrix) {
      if (!(correct->CMx->nmon==0 || correct->CMx->ncor==0))
        compute_trajcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 0,
                               !correct->response_only);
      if (!(correct->CMy->nmon==0 || correct->CMy->ncor==0))
        compute_trajcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, 0,
                                 !correct->response_only);
    }
    for (i_cycle=0; i_cycle<correct->n_xy_cycles; i_cycle++) {
      final_traj = 1;
      if (!x_failed && correct->CMx->ncor && correct->CMx->nmon) {
        switch (correct->method) {
        case GLOBAL_CORRECTION:
          if (!global_trajcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, correct->n_iterations, 
                                    run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL)))
            return 0;
          break;
        case ONE_TO_ONE_CORRECTION:
        case ONE_TO_BEST_CORRECTION:
        case ONE_TO_NEXT_CORRECTION:
          one_to_one_trajcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, correct->n_iterations, 
                                   run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL),
                                   correct->method);
          break;
        case THREAD_CORRECTION:
          thread_trajcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, correct->n_iterations, 
                               run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL));
          break;
        default:
          bomb("Invalid x trajectory correction mode---this should never happen!", NULL);
          break;
        }
        correct->CMx->n_cycles_done = i_cycle+1;
        if (correct->n_iterations<1)
          break;
        rms_before = rms_value(correct->CMx->posi[0], correct->CMx->nmon);
        rms_after  = rms_value(correct->CMx->posi[correct->n_iterations], correct->CMx->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          x_failed = 1;
          if (correct->verbose)
            fputs("horizontal trajectory diverged--setting correctors to zero\n", stdout);
          zero_hcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy && i_cycle>correct->minimum_cycles &&
              correct->method!=THREAD_CORRECTION) {
            x_failed = 1;
            if (correct->verbose)
              fputs("trajectory not improved--discontinuing horizontal correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 0, correct->CMx->kick, 
                          correct->CMx->ncor, correct->CMx->posi, correct->CMx->nmon, NULL, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || x_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || x_failed))
          dump_corrector_data(correct->CMx, &correct->SLx, correct->n_iterations, "horizontal", sim_step);
      }
      if (!y_failed && correct->CMy->ncor && correct->CMy->nmon) {                    
        final_traj = 2;

        switch (correct->method) {
        case GLOBAL_CORRECTION:
          if (!global_trajcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, correct->n_iterations, 
                                    run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL)))
            return 0;
          break;
        case ONE_TO_ONE_CORRECTION:
        case ONE_TO_BEST_CORRECTION:
        case ONE_TO_NEXT_CORRECTION:
          one_to_one_trajcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, correct->n_iterations, 
                                   run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL),
                                   correct->method);
          break;
        case THREAD_CORRECTION:
          thread_trajcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, correct->n_iterations, 
                               run, beamline, starting_coord, (correct->use_actual_beam?beam:NULL));
          break;
        default:
          bomb("Invalid y trajectory correction mode---this should never happen!", NULL);
          break;
        }

        correct->CMy->n_cycles_done = i_cycle+1;
        rms_before = rms_value(correct->CMy->posi[0], correct->CMy->nmon);
        rms_after  = rms_value(correct->CMy->posi[correct->n_iterations], correct->CMy->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          y_failed = 1;
          if (correct->verbose)
            fputs("vertical trajectory diverged--setting correctors to zero\n", stdout);
          zero_vcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMy->corr_accuracy && i_cycle>correct->minimum_cycles &&
              correct->method!=THREAD_CORRECTION) {
            y_failed = 1;
            if (correct->verbose)
              fputs("trajectory not improved--discontinuing vertical correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 2, correct->CMy->kick, 
                          correct->CMy->ncor, correct->CMy->posi, correct->CMy->nmon, NULL, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || y_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || y_failed))
          dump_corrector_data(correct->CMy, &correct->SLy, correct->n_iterations, "vertical", sim_step);
      }
      if (initial_correction && i_cycle==0 && ((correct->CMx->ncor && correct->CMx->nmon) || (correct->CMy->ncor && correct->CMy->nmon)))
        dump_orb_traj(correct->traj[0], beamline->n_elems, "uncorrected", sim_step);
      if (x_failed && y_failed) {
        if (correct->verbose)
          fputs("trajectory correction discontinued\n", stdout);
        break;
      }
    }
    if (!initial_correction && correct->n_iterations>=1 && ((correct->CMx->ncor && correct->CMx->nmon) || (correct->CMy->ncor && correct->CMy->nmon)))
      dump_orb_traj(correct->traj[final_traj], beamline->n_elems, "corrected", sim_step);
    if (starting_coord)
      for (i=0; i<6; i++)
        starting_coord[i] = 0;  /* don't want to seem to be returning a closed orbit here */
    bombed = 0;
    break;
  case ORBIT_CORRECTION:
    if (correct->CMx->fixed_length)
      Cdp = tmalloc(sizeof(*Cdp)*(correct->n_iterations+1));
    else
      Cdp = NULL;
    x_failed = y_failed = bombed = 0;
    if (usePerturbedMatrix) {
      if (correct->verbose)
        fprintf(stdout, "Computing orbit correction matrices\n"); 
      if (!(correct->CMx->nmon==0 || correct->CMx->ncor==0))
        compute_orbcor_matrices(correct->CMx, &correct->SLx, 0, run, beamline, 0, 
                                !correct->response_only, fixedLengthMatrix, correct->verbose);
      if (!(correct->CMy->nmon==0 || correct->CMy->ncor==0))
        compute_orbcor_matrices(correct->CMy, &correct->SLy, 2, run, beamline, 0, 
                                !correct->response_only, fixedLengthMatrix, correct->verbose);
    }
    for (i_cycle=0; i_cycle<correct->n_xy_cycles; i_cycle++) {
      final_traj = 1;
      if (!x_failed && correct->CMx->ncor && correct->CMx->nmon) {
        if ((n_iter_taken = orbcor_plane(correct->CMx, &correct->SLx, 0, correct->traj, 
                                         correct->n_iterations, correct->clorb_accuracy, 
                                         correct->clorb_iterations, 
                                         correct->clorb_iter_fraction,
                                         run, beamline, 
                                         closed_orbit, Cdp))<0) {
          fprintf(stdout, "Horizontal correction has failed.\n");
          fflush(stdout);
          bombed = 1; 
          break;
        }
        correct->CMx->n_cycles_done = i_cycle+1;
        if (correct->n_iterations<1)
          break;
        rms_before = rms_value(correct->CMx->posi[0], correct->CMx->nmon);
        rms_after  = rms_value(correct->CMx->posi[n_iter_taken], correct->CMx->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          x_failed = 1;
          if (correct->verbose)
            fputs("horizontal orbit diverged--setting correctors to zero\n", stdout);
          zero_hcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy) {
            x_failed = 1;
            if (correct->verbose)
              fputs("orbit not improved--discontinuing horizontal correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 0, correct->CMx->kick, 
                          correct->CMx->ncor, correct->CMx->posi, correct->CMx->nmon, Cdp, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || x_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || x_failed)) 
          dump_corrector_data(correct->CMx, &correct->SLx, correct->n_iterations, "horizontal", sim_step);
      }
      if (!y_failed && correct->CMy->ncor && correct->CMy->nmon) {
        final_traj = 2;
        if ((n_iter_taken = orbcor_plane(correct->CMy, &correct->SLy, 2, correct->traj+1, 
                                         correct->n_iterations, correct->clorb_accuracy, 
                                         correct->clorb_iterations, 
                                         correct->clorb_iter_fraction,
                                         run, beamline, 
                                         closed_orbit, Cdp))<0) {
          bombed = 1;
          fprintf(stdout, "Vertical correction has failed.\n");
          fflush(stdout);
          break;
        }
        correct->CMy->n_cycles_done = i_cycle+1;
        rms_before = rms_value(correct->CMy->posi[0], correct->CMy->nmon);
        rms_after  = rms_value(correct->CMy->posi[n_iter_taken], correct->CMy->nmon);
#if defined(IEEE_MATH)
        if (isnan(rms_before) || isnan(rms_after) || isinf(rms_before) || isinf(rms_after)) {
          y_failed = 1;
          if (correct->verbose)
            fputs("vertical orbit diverged--setting correctors to zero\n", stdout);
          zero_vcorrectors(&(beamline->elem), run, correct);
        }
        else
#endif
          if (rms_before<=rms_after+correct->CMx->corr_accuracy) {
            y_failed = 1;
            if (correct->verbose)
              fputs("orbit not improved--discontinuing vertical correction\n", stdout);
          }
        dump_cormon_stats(correct->verbose, 2, correct->CMy->kick, 
                          correct->CMy->ncor, correct->CMy->posi, correct->CMy->nmon, Cdp, 
                          correct->n_iterations, i_cycle, i_cycle==correct->n_xy_cycles-1 || y_failed,
                          sim_step, initial_correction);
        if (!initial_correction && (i_cycle==correct->n_xy_cycles-1 || y_failed))
          dump_corrector_data(correct->CMy, &correct->SLy, correct->n_iterations, "vertical", sim_step);
      }
      if (initial_correction && i_cycle==0 && ((correct->CMx->ncor && correct->CMx->nmon) || (correct->CMy->ncor && correct->CMy->nmon))) 
        dump_orb_traj(correct->traj[0], beamline->n_elems, "uncorrected", sim_step);
      if (x_failed && y_failed) {
        if (correct->verbose)
          fputs("orbit correction discontinued\n", stdout);
        break;
      }
    }
    if (!initial_correction && !bombed && correct->n_iterations>=1 && ((correct->CMx->ncor && correct->CMx->nmon) || (correct->CMy->ncor && correct->CMy->nmon)))
      dump_orb_traj(correct->traj[final_traj], beamline->n_elems, "corrected", sim_step);
    break;
  }

  beamline->closed_orbit = correct->traj[final_traj];

  log_exit("do_correction");
  return(!bombed);
}

void compute_trajcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, long find_only, long invert)
{
#ifdef DEBUG
  long i_debug;
  FILE *fpdeb;
  char s[100];
#endif
  ELEMENT_LIST *corr, *start;
  TRAJECTORY *traj0, *traj1;
  long kick_offset, i_corr, i_moni, i, equalW;
  long n_part;
  double **one_part, p, p0, kick0, kick1, corr_tweek, corrCalibration, *moniCalibration, W0=0.0;
  double *weight=NULL, conditionNumber;
  VMATRIX *save;
  long i_type;

  start = find_useable_moni_corr(&CM->nmon, &CM->ncor, &CM->mon_index,
                                 &CM->umoni, &CM->ucorr, &CM->kick_coef, &CM->sl_index, coord, SL, run, beamline, 0);
  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction may be unstable (use SV controls).\n");
    fflush(stdout);
  }
  if (CM->ncor==0) {
    fprintf(stdout, "Warning: no correctors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    return;
  }
  if (CM->nmon==0) {
    fprintf(stdout, "Warning: no monitors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    CM->ncor = 0;
    return;
  }

  if (find_only)
    return;

  fprintf(stdout, "computing response matrix...\n");
  fflush(stdout);
  report_stats(stdout, "start");

  /* allocate matrices for this plane */
  if (CM->C)
    matrix_free(CM->C);
  CM->C  = matrix_get(CM->nmon, CM->ncor);   /* Response matrix */
  if (CM->T)
    matrix_free(CM->T);
  CM->T  = NULL;
  
  /* arrays for trajectory data */
  traj0 = tmalloc(sizeof(*traj0)*beamline->n_elems);
  traj1 = tmalloc(sizeof(*traj1)*beamline->n_elems);
  one_part = (double**)czarray_2d(sizeof(**one_part), 1, 7);

  /* find initial trajectory */
  p = p0 = sqrt(sqr(run->ideal_gamma)-1);
  n_part = 1;
  fill_double_array(*one_part, 7, 0.0);
  if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                   traj0, run, 0, 
                   TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL, NULL, NULL))
    bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

#ifdef  DEBUG
  i_debug = 0;
  sprintf(s, "traj%c-%ld.deb", (coord==0?'x':'y'), i_debug++);
  fpdeb = fopen_e(s, "w", 0);
  fprintf(fpdeb, "z (m)\n%c (m)\ninitial trajectory computed in compute_trajcor_matrices\n\n%ld\n", 
          (coord==0?'x':'y'), beamline->n_elems);
  for (i=0; i<beamline->n_elems; i++)
    fprintf(fpdeb, "%e\t%e\n", traj0[i].elem->end_pos, traj0[i].centroid[coord]);
  fclose(fpdeb);
#endif

  /* set up weight matrix and monitor calibration array */
  weight = tmalloc(sizeof(*weight)*CM->nmon);
  moniCalibration = tmalloc(sizeof(*moniCalibration)*CM->nmon);
  equalW = 1;
  for (i_moni=0; i_moni<CM->nmon; i_moni++) {
    weight[i_moni] = getMonitorWeight(CM->umoni[i_moni]);
    if (!i_moni)
      W0 = weight[i_moni];
    else if (W0!=weight[i_moni])
      equalW = 0;
    moniCalibration[i_moni] = getMonitorCalibration(CM->umoni[i_moni], coord);
  }

  for (i_corr = 0; i_corr<CM->ncor; i_corr++) {
    corr = CM->ucorr[i_corr];

    if ((i_type = find_index(corr->type, SL->corr_type, SL->n_corr_types))<0)
      bomb("failed to find corrector type in type list", NULL);

    kick_offset = SL->param_offset[i_type];
    corr_tweek  = SL->corr_tweek[i_type];

    /* record value of corrector */
    kick0 = *((double*)(corr->p_elem+kick_offset));

    /* change the corrector by corr_tweek and compute the new matrix for the corrector */
    kick1 = *((double*)(corr->p_elem+kick_offset)) = kick0 + corr_tweek;
#ifdef DEBUG
    fprintf(stdout, "corrector %s tweeked to %e (type=%ld, offset=%ld)\n", corr->name, *((double*)(corr->p_elem+kick_offset)),
            i_type, kick_offset);
    fflush(stdout);
#endif

    if (corr->matrix)
      save = corr->matrix;
    else
      save = NULL;
    compute_matrix(corr, run, NULL);
#ifdef DEBUG
    print_matrices(stdout, "*** corrector matrix:", corr->matrix);
#endif

    /* track with positively-tweeked corrector */
    p = p0;
    n_part = 1;
    fill_double_array(*one_part, 7, 0.0);
    if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                     traj1, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL, NULL, NULL))
      bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

#ifdef DEBUG
    sprintf(s, "traj%c-%ld.deb", (coord==0?'x':'y'), i_debug++);
    fpdeb = fopen_e(s, "w", 0);
    fprintf(fpdeb, "z (m)\n%c (m)\ntrajectory computed in compute_trajcor_matrices\n%s = %e\n%ld\n", 
            (coord==0?'x':'y'), corr->name, kick1, beamline->n_elems);
    for (i=0; i<beamline->n_elems; i++)
      fprintf(fpdeb, "%e\t%e\n", traj1[i].elem->end_pos, traj1[i].centroid[coord]);
    fclose(fpdeb);
#endif


#if TWO_POINT_TRAJRESPONSE
    /* change the corrector by -corr_tweek and compute the new matrix for the corrector */
    kick1 = *((double*)(corr->p_elem+kick_offset)) = kick0 - corr_tweek;
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
#ifdef DEBUG
    fprintf(stdout, "corrector %s tweeked to %e\n", corr->name, *((double*)(corr->p_elem+kick_offset)));
    fflush(stdout);
#endif
    free_matrices(corr->matrix); corr->matrix = NULL;
    compute_matrix(corr, run, NULL);
#ifdef DEBUG
    print_matrices(stdout, "*** corrector matrix:", corr->matrix);
#endif

    /* track with tweeked corrector */
    p = p0;
    n_part = 1;
    fill_double_array(*one_part, 7, 0.0);
    if (!do_tracking(NULL, one_part, n_part, NULL, beamline, &p, (double**)NULL, (BEAM_SUMS**)NULL, (long*)NULL,
                     traj0, run, 0, TEST_PARTICLES+TIME_DEPENDENCE_OFF, 1, 0, NULL, NULL, NULL, NULL, NULL))
      bomb("tracking failed for test particle (compute_trajcor_matrices())", NULL);

    /* compute coefficients of array C that are driven by this corrector */
    corrCalibration = getCorrectorCalibration(CM->ucorr[i_corr], coord)/(2*corr_tweek);
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      i = CM->mon_index[i_moni];
      Mij(CM->C, i_moni, i_corr) = moniCalibration[i_moni]*corrCalibration*
        (traj1[i].centroid[coord] - traj0[i].centroid[coord]);
    }
#else

    /* compute coefficients of array C that are driven by this corrector */
    corrCalibration = getCorrectorCalibration(CM->ucorr[i_corr], coord)/corr_tweek;
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      i = CM->mon_index[i_moni];
      Mij(CM->C, i_moni, i_corr) = moniCalibration[i_moni]*corrCalibration*
        (traj1[i].centroid[coord] - traj0[i].centroid[coord]);
    }
#endif

    /* change the corrector back */
    *((double*)(corr->p_elem+kick_offset)) = kick0;
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
    if (corr->matrix) {
      free_matrices(corr->matrix);
      corr->matrix = NULL;
    }
    if (save)
      corr->matrix = save;
    else 
      compute_matrix(corr, run, NULL);
  } 
  free(moniCalibration);
  tfree(traj0); traj0 = NULL;
  tfree(traj1); traj1 = NULL;
  free_czarray_2d((void**)one_part, 1, 7); one_part = NULL;
#ifdef DEBUG
  matrix_show(CM->C    , "%13.6le ", "influence matrix\n", stdout);
#endif
  report_stats(stdout, "done");

  if (invert) {
    report_stats(stdout, "Computing correction matrix ");
    fflush(stdout);

    /* compute correction matrix T */
    if (CM->auto_limit_SVs && (CM->C->m < CM->C->n) && CM->remove_smallest_SVs < (CM->C->n - CM->C->m)) {
      CM->remove_smallest_SVs = CM->C->n - CM->C->m;
      printf("Removing %ld smallest singular values to prevent instability\n", (long)CM->remove_smallest_SVs);
    }
    
    CM->T = matrix_invert(CM->C, equalW?NULL:weight, (int32_t)CM->keep_largest_SVs, (int32_t)CM->remove_smallest_SVs,
                          CM->minimum_SV_ratio, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &conditionNumber);
    matrix_scmul(CM->T, -1);

    if (weight)
      free(weight);

    report_stats(stdout, "\ndone.");
    printf("Condition number is %e\n", conditionNumber);
    fflush(stdout);
  }
  
  
#ifdef DEBUG
    matrix_show(CM->T, "%13.6le ", "correction matrix\n", stdout);
#endif
}

long global_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, 
                          RUN *run, LINE_LIST *beamline, double *starting_coord, BEAM *beam)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *traj;
  long iteration, kick_offset;
  long i_moni, i_corr;
  long n_part, i, tracking_flags, sl_index;
  double **particle;
  double p, x, y, reading, fraction, minFraction, param, change;
  MAT *Qo, *dK;
  
  log_entry("global_trajcor_plane");

  if (!matrix_check(CM->C) || !matrix_check(CM->T))
    bomb("corrupted response matrix detected (global_trajcor_plane)", NULL);
  if (!CM->mon_index)
    bomb("monitor index array is NULL (global_trajcor_plane)", NULL);
  if (!CM->posi)
    bomb("monitor readout array is NULL (global_trajcor_plane)", NULL);
  if (!CM->kick)
    bomb("corrector value array is NULL (global_trajcor_plane)", NULL);
  if (!traject)
    bomb("no trajectory arrays supplied (global_trajcor_plane)", NULL);
  if (!CM->T)
    bomb("no inverse matrix computed (global_trajcor_plane)", NULL);

  Qo = matrix_get(CM->nmon, 1);   /* Vector of BPM errors */
   
  if (!beam) {
    particle = (double**)czarray_2d(sizeof(**particle), 1, 7);
    tracking_flags = TEST_PARTICLES;
  }
  else {
    if (beam->n_to_track==0)
      bomb("no particles to track in global_trajcor_plane()", NULL);
    particle = (double**)czarray_2d(sizeof(**particle), beam->n_to_track, 7);
    tracking_flags = TEST_PARTICLES+TEST_PARTICLE_LOSSES;
  }

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction may be unstable (use SV controls)\n");
    fflush(stdout);
  }
  for (iteration=0; iteration<=n_iterations; iteration++) {
    if (!CM->posi[iteration])
      bomb("monitor readout array for this iteration is NULL (global_trajcor_plane)", NULL);
    if (!CM->kick[iteration])
      bomb("corrector value array for this iteration is NULL (global_trajcor_plane)", NULL);
    if (!(traj = traject[iteration?1:0]))
      bomb("trajectory array for this iteration is NULL (global_trajcor_plane)", NULL);

    /* find trajectory */
    p = sqrt(sqr(run->ideal_gamma)-1);

    if (!beam) {
      if (!starting_coord)
        fill_double_array(*particle, 7, 0.0);
      else {
        for (i=0; i<6; i++)
          particle[0][i] = starting_coord[i];
      }
      n_part = 1;
    }
    else
      copy_particles(particle, beam->particle, n_part=beam->n_to_track);

    n_part = do_tracking(NULL, particle, n_part, NULL, beamline, &p, (double**)NULL, 
                         (BEAM_SUMS**)NULL, (long*)NULL,
                         traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL, NULL, NULL);
    if (beam) {
      fprintf(stdout, "%ld particles survived tracking", n_part);
      fflush(stdout);
      if (n_part==0) {
        for (i=0; i<beamline->n_elems+1; i++)
          if (traj[i].n_part==0)
            break;
        if (i!=0 && i<beamline->n_elems+1)
          fprintf(stdout, "---all beam lost before z=%em (element %s)",
                  traj[i].elem->end_pos, traj[i].elem->name);
        fflush(stdout);
        fputc('\n', stdout);
        return 0;
      }
      fputc('\n', stdout);
    }

    /* find readings at monitors and add in reading errors */
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      if (!(eptr=traj[CM->mon_index[i_moni]].elem))
        bomb("invalid element pointer in trajectory array (global_trajcor_plane)", NULL);
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      reading = computeMonitorReading(eptr, coord, x, y, 0);
      if (isnan(reading) || isinf(reading)) 
        return 0;
      CM->posi[iteration][i_moni] = Mij(Qo, i_moni, 0) =  
        reading + (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0);
    }
    
    if (iteration==n_iterations)
      break;

    /* solve for the corrector changes */
    dK = matrix_mult(CM->T, Qo);

#ifdef DEBUG
    matrix_show(Qo, "%13.6le ", "traj matrix\n", stdout);
    matrix_show(dK, "%13.6le ", "kick matrix\n", stdout);
#endif

    /* step through beamline find any kicks that are over their limits */
    minFraction = 1;
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      param = fabs(*((double*)(corr->p_elem+kick_offset)) +
                   (change=Mij(dK, i_corr, 0)/CM->kick_coef[i_corr]*CM->corr_fraction));
      if (SL->corr_limit[sl_index] && param>SL->corr_limit[sl_index]) {
        fraction = fabs((SL->corr_limit[sl_index]-fabs(*((double*)(corr->p_elem+kick_offset))))/change);
        if (fraction<minFraction)
          minFraction = fraction;
      }
    }
    fraction = minFraction*CM->corr_fraction;

#if defined(DEBUG)
    fprintf(stdout, "Changing correctors:");
    fflush(stdout);
#endif
    /* step through beamline and change correctors */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
#if defined(DEBUG)
      fprintf(stdout, "before = %e, ", *((double*)(corr->p_elem+kick_offset)));
      fflush(stdout);
#endif
      if (iteration==0)
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      *((double*)(corr->p_elem+kick_offset)) += Mij(dK, i_corr, 0)*fraction;
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
#if defined(DEBUG)
      fprintf(stdout, "after = %e\n", *((double*)(corr->p_elem+kick_offset)));
      fflush(stdout);
#endif
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    matrix_free(dK);
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
  }

  /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
  beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
  beamline->flags &= ~BEAMLINE_TWISS_CURRENT;

  if (!beam)
    free_czarray_2d((void**)particle, 1, 7);
  else
    free_czarray_2d((void**)particle, beam->n_to_track, 7);
  particle = NULL;
  log_exit("global_trajcor_plane");

  matrix_free(Qo);

  return(1);
}

void one_to_one_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run,
                              LINE_LIST *beamline, double *starting_coord, BEAM *beam, long method)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *traj;
  long iteration, kick_offset;
  long i_moni, i_corr, sl_index;
  long n_part, i, tracking_flags;
  double **particle, param, fraction;
  double p, x, y, reading;
  double response;
  
  log_entry("one_to_one_trajcor_plane");
  
  if (!matrix_check(CM->C))
    bomb("corrupted correction matrix detected (one_to_one_trajcor_plane)", NULL);
  if (!CM->mon_index)
    bomb("monitor index array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!CM->posi)
    bomb("monitor readout array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!CM->kick)
    bomb("corrector value array is NULL (one_to_one_trajcor_plane)", NULL);
  if (!traject)
    bomb("no trajectory arrays supplied (one_to_one_trajcor_plane)", NULL);
  
  if (!beam) {
    particle = (double**)czarray_2d(sizeof(**particle), 1, 7);
    tracking_flags = TEST_PARTICLES;
  }
  else {
    if (beam->n_to_track==0)
      bomb("no particles to track in one_to_one_trajcor_plane()", NULL);
    particle = (double**)czarray_2d(sizeof(**particle), beam->n_to_track, 7);
    tracking_flags = TEST_PARTICLES+TEST_PARTICLE_LOSSES;
  }
  
  for (iteration=0; iteration<=n_iterations; iteration++) {
    if (!CM->posi[iteration])
      bomb("monitor readout array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);
    if (!CM->kick[iteration])
      bomb("corrector value array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);
    if (!(traj = traject[iteration?1:0]))
      bomb("trajectory array for this iteration is NULL (one_to_one_trajcor_plane)", NULL);

    /* record the starting trajectory */
    p = sqrt(sqr(run->ideal_gamma)-1);
    if (!beam) {
      if (!starting_coord)
        fill_double_array(*particle, 7, 0.0);
      else {
        for (i=0; i<7; i++)
          particle[0][i] = starting_coord[i];
      }
      n_part = 1;
    }
    else
      copy_particles(particle, beam->particle, n_part=beam->n_to_track);
    n_part = do_tracking(NULL, particle, n_part, NULL, beamline, &p, (double**)NULL, 
                         (BEAM_SUMS**)NULL, (long*)NULL,
                         traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL, NULL, NULL);
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      if (!(eptr=traj[CM->mon_index[i_moni]].elem))
        bomb("invalid element pointer in trajectory array (one_to_one_trajcor_plane)", NULL);
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      CM->posi[iteration][i_moni] = computeMonitorReading(eptr, coord, x, y, 0);
    }
    if (iteration==n_iterations)
      break;
    
    i_moni = 0;
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      switch (method) {
      case ONE_TO_NEXT_CORRECTION:
      case ONE_TO_ONE_CORRECTION:
        /* Find next BPM downstream of this corrector */
        for ( ; i_moni<CM->nmon; i_moni++) 
          if (CM->ucorr[i_corr]->end_pos < CM->umoni[i_moni]->end_pos)
            break;
        break;
      case ONE_TO_BEST_CORRECTION:
        /* Find BPM with larger response than its neighbors */
        for ( ; i_moni<CM->nmon; i_moni++) 
          if (CM->ucorr[i_corr]->end_pos < CM->umoni[i_moni]->end_pos)
            break;
        if (i_moni!=CM->nmon) {
          response = fabs(Mij(CM->C, i_moni, i_corr));
          for ( ; i_moni<CM->nmon-1; i_moni++) {
            if (response>fabs(Mij(CM->C, i_moni+1, i_corr)))
              break;
            response = fabs(Mij(CM->C, i_moni+1, i_corr));
          }
        }
        break;
      default:
        printf("Error: Unknown method in one_to_one_trajcor_plane: %ld---This shouldn't happen!\n", method);
        exit(1);
        break;
      }
      
      if (i_moni==CM->nmon)
        break;
      
      /* find trajectory */
      p = sqrt(sqr(run->ideal_gamma)-1);
      
      if (!beam) {
        if (!starting_coord)
          fill_double_array(*particle, 7, 0.0);
        else {
          for (i=0; i<7; i++)
            particle[0][i] = starting_coord[i];
        }
        n_part = 1;
      }
      else
        copy_particles(particle, beam->particle, n_part=beam->n_to_track);
      
      n_part = do_tracking(NULL, particle, n_part, NULL, beamline, &p, (double**)NULL, 
                           (BEAM_SUMS**)NULL, (long*)NULL,
                           traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL, NULL, NULL);
      if (traj[CM->mon_index[i_moni]].n_part==0) {
        printf("beam lost at position %e m\n", traj[CM->mon_index[i_moni]].elem->end_pos);
        break;
      }
      
      if (!(eptr=traj[CM->mon_index[i_moni]].elem))
        bomb("invalid element pointer in trajectory array (one_to_one_trajcor_plane)", NULL);
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      reading = computeMonitorReading(eptr, coord, x, y, 0);
      reading += (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0);

#if defined(DEBUG)
      fprintf(stdout, "i_moni = %ld, i_corr = %ld, reading = %e", i_moni, i_corr, reading);
      fflush(stdout);
#endif
      if (iteration==n_iterations || (i_corr>=CM->ncor || Mij(CM->C, i_moni, i_corr)==0)) {
#if defined(DEBUG)
        fprintf(stdout, "\n");
        fflush(stdout);
#endif
        continue;
      }

      /* change corrector */
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];           /* steering list index of this corrector */
      kick_offset = SL->param_offset[sl_index];  /* offset of steering parameter in element structure */
      if (iteration==0)
        /* Record initial value */
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      /* Compute new value of the parameter
       * NewValue = OldValue - BpmReading/ResponseCoefficient*CorrectionFraction
       */
      param = *((double*)(corr->p_elem+kick_offset)) - reading/Mij(CM->C, i_moni, i_corr)*CM->corr_fraction;
      /* Check that we haven't exceeded allowed strength */
      fraction = 1;
      if (param && SL->corr_limit[sl_index])
        if ((fraction = SL->corr_limit[sl_index]/fabs(param))>1)
          fraction = 1;
      *((double*)(corr->p_elem+kick_offset)) = param*fraction;
      /* Record the new kick strength */
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
#if defined(DEBUG)
      fprintf(stdout, ", param = %e, fraction=%e\n", param, fraction);
      fflush(stdout);
#endif
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
      if (beamline->links)
        assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
      
      /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
    }
  }

  if (!beam)
    free_czarray_2d((void**)particle, 1, 7);
  else
    free_czarray_2d((void**)particle, beam->n_to_track, 7);
  particle = NULL;
  log_exit("one_to_one_trajcor_plane");
}

void thread_trajcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **traject, long n_iterations, RUN *run,
                          LINE_LIST *beamline, double *starting_coord, BEAM *beam)
{
  ELEMENT_LIST *corr;
  TRAJECTORY *traj;
  long iteration, kick_offset;
  long i_corr, i_moni, sl_index, iElem, nElems, iBest;
  long n_part, n_left, i, tracking_flags, done, direction, improved, worsened;
  double **particle, param, fraction, bestValue, origValue, lastValue;
  double p, scanStep;
  long iScan, nScan;

  if (!CM->mon_index)
    bomb("monitor index array is NULL (thread__trajcor_plane)", NULL);
  if (!CM->posi)
    bomb("monitor readout array is NULL (thread__trajcor_plane)", NULL);
  if (!CM->kick)
    bomb("corrector value array is NULL (thread__trajcor_plane)", NULL);
  if (!traject)
    bomb("no trajectory arrays supplied (thread__trajcor_plane)", NULL);
  
  if (!beam) {
    particle = (double**)czarray_2d(sizeof(**particle), 1, 7);
  }
  else {
    if (beam->n_to_track==0)
      bomb("no particles to track in thread__trajcor_plane()", NULL);
    particle = (double**)czarray_2d(sizeof(**particle), beam->n_to_track, 7);
  }
  tracking_flags = TEST_PARTICLES+TEST_PARTICLE_LOSSES;
  nElems = beamline->n_elems+1;
  done = 0;

  nScan = CM->default_threading_divisor+1.5;
  scanStep = 1.0/(nScan-1);
  
  for (iteration=0; iteration<=n_iterations; iteration++) {
    if (!CM->posi[iteration])
      bomb("monitor readout array for this iteration is NULL (thread__trajcor_plane)", NULL);
    if (!CM->kick[iteration])
      bomb("corrector value array for this iteration is NULL (thread__trajcor_plane)", NULL);
    if (!(traj = traject[iteration?1:0]))
      bomb("trajectory array for this iteration is NULL (thread__trajcor_plane)", NULL);

    /* Establish baseline distance traveled before loss */
    p = sqrt(sqr(run->ideal_gamma)-1);
    if (!beam) {
      if (!starting_coord)
        fill_double_array(*particle, 7, 0.0);
      else {
        for (i=0; i<7; i++)
          particle[0][i] = starting_coord[i];
      }
      n_part = 1;
    }
    else
      copy_particles(particle, beam->particle, n_part=beam->n_to_track);
    
    n_left = do_tracking(NULL, particle, n_part, NULL, beamline, &p, (double**)NULL, 
                         (BEAM_SUMS**)NULL, (long*)NULL,
                         traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL, NULL, NULL);
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      double x, y;
      x = traj[CM->mon_index[i_moni]].centroid[0];
      y = traj[CM->mon_index[i_moni]].centroid[2];
      CM->posi[iteration][i_moni] = computeMonitorReading(traj[CM->mon_index[i_moni]].elem, coord, x, y, 0);
    }
    if (iteration==n_iterations)
      break;

    if (n_left!=n_part) {
      for (iElem=iBest=0; iElem<nElems; iElem++)
        if (traj[iElem].n_part==0)
          break;
      iBest = iElem;
      
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        corr = CM->ucorr[i_corr];
        sl_index = CM->sl_index[i_corr];           /* steering list index of this corrector */
        kick_offset = SL->param_offset[sl_index];  /* offset of steering parameter in element structure */
        if (iteration==0)
          /* Record initial value */
          CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];

        lastValue = bestValue = origValue = *((double*)(corr->p_elem+kick_offset));
        if (!done && corr->end_pos < traj[iBest].elem->end_pos) {
          improved = 0;
          for (direction=-1; !improved && !done && direction<2 ; direction+=2) {
            worsened = 0;
            for (iScan=0; !worsened && !done && iScan<nScan; iScan++) {
              /* -- Compute new value of the parameter
               *    NewValue = OldValue + CorrectorTweek 
               */
              param = origValue 
                + SL->corr_tweek[sl_index]/CM->kick_coef[i_corr]*iScan*scanStep*direction;
              /* Check that we haven't exceeded allowed strength */
              fraction = 1;
              if (param && SL->corr_limit[sl_index])
                if ((fraction = SL->corr_limit[sl_index]/fabs(param))>1)
                  fraction = 1;
              lastValue = *((double*)(corr->p_elem+kick_offset)) = param*fraction;
              if (corr->matrix) {
                free_matrices(corr->matrix);
                tfree(corr->matrix);
                corr->matrix = NULL;
              }
              compute_matrix(corr, run, NULL);
              if (beamline->links)
                assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
              
              /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
              beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
              beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
              
              /* -- Find trajectory */
              p = sqrt(sqr(run->ideal_gamma)-1);
              if (!beam) {
                if (!starting_coord)
                  fill_double_array(*particle, 7, 0.0);
                else {
                  for (i=0; i<7; i++)
                    particle[0][i] = starting_coord[i];
                }
                n_part = 1;
              }
              else
                copy_particles(particle, beam->particle, n_part=beam->n_to_track);
              n_left = do_tracking(NULL, particle, n_part, NULL, beamline, &p, (double**)NULL, 
                                   (BEAM_SUMS**)NULL, (long*)NULL,
                                   traj, run, 0, tracking_flags, 1, 0, NULL, NULL, NULL, NULL, NULL);
              /* Determine if this is better than the previous best */
              if (n_left!=n_part) {
                for (iElem=0; iElem<nElems; iElem++) {
                  if (traj[iElem].n_part==0)
                    break;
                }
                if (iElem>iBest) {
                  iBest = iElem;
                  bestValue = *((double*)(corr->p_elem+kick_offset));
                  improved = 1;
                } else if (iElem<iBest) {
                  /* getting worse, so quit this scan */
                  worsened = 1;
                }
              } else {
                /* Made it to the end, so quit */
                iBest = nElems;
                bestValue = *((double*)(corr->p_elem+kick_offset));
                done = 1;
                break;
              }
            }
          }
          
          if (lastValue!=bestValue) {
            *((double*)(corr->p_elem+kick_offset)) = bestValue;
            if (corr->matrix) {
              free_matrices(corr->matrix);
              tfree(corr->matrix);
              corr->matrix = NULL;
            }
            compute_matrix(corr, run, NULL);
            if (beamline->links)
              assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);
          }
          if (bestValue!=origValue) {
            printf("Beam now advanced to element %ld\n", iBest);
            fflush(stdout);
          }
          
          /* indicate that beamline concatenation and Twiss parameter computation (if wanted) are not current */
          beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
          beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
        }
        
        /* Record the kick strength */
        CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      }
    } else {
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        corr = CM->ucorr[i_corr];
        sl_index = CM->sl_index[i_corr];           /* steering list index of this corrector */
        kick_offset = SL->param_offset[sl_index];  /* offset of steering parameter in element structure */
        if (iteration==0) 
          /* Record initial value */
          CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];

        /* Record the kick strength */
        CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      }
    }
  }
  
  if (!beam)
    free_czarray_2d((void**)particle, 1, 7);
  else
    free_czarray_2d((void**)particle, beam->n_to_track, 7);
  particle = NULL;
}


ELEMENT_LIST *find_useable_moni_corr(int32_t *nmon, int32_t *ncor, long **mon_index,
                                     ELEMENT_LIST ***umoni, ELEMENT_LIST ***ucorr,
                                     double **kick_coef, long **sl_index, long plane, STEERING_LIST *SL, 
                                     RUN *run, LINE_LIST *beamline, long recircs)
{
  ELEMENT_LIST *corr, *moni, *start;
  long i_elem, moni_follows, corr_seen, index;
  long moni_type_1, moni_type_2;

  log_entry("find_useable_moni_corr");

  switch (plane) {
  case 0:
    moni_type_1 = T_HMON;
    moni_type_2 = T_MONI;
    break;
  case 2:
    moni_type_1 = T_VMON;
    moni_type_2 = T_MONI;
    break;
  default:
    fprintf(stdout, "error: invalid coordinate for correction: %ld\n", plane);
    fflush(stdout);
    exit(1);
    break;
  }

  start = &(beamline->elem);
  if (recircs) {
    /* advance to recirculation point, if there is one */
    while (start) {
      if (start->type==T_RECIRC)
        break;
      start = start->succ;
    }
    if (!start)
      start = &(beamline->elem);
  }
  if (*ncor && *nmon) {
    log_exit("find_useable_moni_corr");
    return start;
  }

  
  *ncor = *nmon = 0;
  *mon_index = NULL;
  *umoni = NULL;
  *ucorr = NULL;
  *kick_coef = NULL;
  *sl_index = NULL;

  /* first count correctors */
  corr = start;
  do {
    /* advance to position of next corrector */
    if (!(corr = next_element_of_types(corr, SL->corr_type, SL->n_corr_types, &index)))
      break;
    if (steering_corrector(corr, SL, plane) && match_string(corr->name, SL->corr_name, SL->n_corr_types, EXACT_MATCH)>=0) {
      *ucorr = trealloc(*ucorr, sizeof(**ucorr)*(*ncor+1));
      (*ucorr)[*ncor] = corr;
      *sl_index = trealloc(*sl_index, sizeof(**sl_index)*(*ncor+1));
      (*sl_index)[*ncor] = index;
      *kick_coef = trealloc(*kick_coef, sizeof(**kick_coef)*(*ncor+1));
      if (!((*kick_coef)[*ncor] = compute_kick_coefficient(corr, plane, corr->type, 
                                                           SL->corr_tweek[index], corr->name, SL->corr_param[index], run))) {
        fprintf(stdout, "error: changing %s.%s does not kick the beam--can't use for steering.\n",
                corr->name, SL->corr_param[index]);
        fflush(stdout);
        exit(1);
      }
      
      if (!recircs) {
        /* Since the beamline doesn't recirculate, advance through remainder of beamline
         * to make sure there are subsequent monitors */
        moni = corr->succ;
        moni_follows = 0;
        do {
          while (moni && !((moni->type==moni_type_1 || moni->type==moni_type_2) &&
                           ((MONI*)moni->p_elem)->weight>0) )
            moni = moni->succ;
          if (moni)
            moni_follows = 1;
        } while (moni && (moni=moni->succ) && !moni_follows);
        if (!moni_follows)
          break;        /* ignore correctors with no monitors following */
      }
      *ncor += 1;
    }
  } while ((corr=corr->succ));

  if (!recircs) {
    /* count all monitors with one or more correctors upstream and non-negative weight */
    moni = start;
    corr_seen = 0;
    i_elem = 0;
    while (moni) {
      if (find_index(moni->type, SL->corr_type, SL->n_corr_types)!=-1 && 
          (steering_corrector(moni, SL, plane) || match_string(moni->name, SL->corr_name, SL->n_corr_types, EXACT_MATCH)>=0))
        corr_seen = 1;
      if ((moni->type==moni_type_1 || moni->type==moni_type_2) 
          && ((MONI*)moni->p_elem)->weight>0 && corr_seen) {
        *nmon += 1;
        *umoni = trealloc(*umoni, *nmon*sizeof(**umoni));
        (*umoni)[*nmon-1] = moni;
        *mon_index = trealloc(*mon_index, *nmon*sizeof(**mon_index));
        (*mon_index)[*nmon-1] = i_elem;
      }
      moni = moni->succ;
      i_elem++;
    }
  }
  else {
    /* find all monitors */
    moni = start;
    i_elem = 0;
    while (moni) {
      if ((moni->type==moni_type_1 || moni->type==moni_type_2) 
          && ((MONI*)moni->p_elem)->weight>0) {
        *nmon += 1;
        *umoni = trealloc(*umoni, *nmon*sizeof(**umoni));
        (*umoni)[*nmon-1] = moni;
        *mon_index = trealloc(*mon_index, *nmon*sizeof(**mon_index));
        (*mon_index)[*nmon-1] = i_elem;
      }
      moni = moni->succ;
      i_elem++;
    }
  } 
  log_exit("find_useable_moni_corr");
  return(start);
}

void compute_orbcor_matrices(CORMON_DATA *CM, STEERING_LIST *SL, long coord, RUN *run, LINE_LIST *beamline, 
                             long find_only, long invert, long fixed_length, long verbose)
{
  ELEMENT_LIST *start;
  long i_corr, i_moni, equalW;
  double coef, htune, moniFactor, *corrFactor, *corrFactorFL, coefFL, W0=0.0;
  double *weight=NULL, conditionNumber;
  char memName[1024];

  start = find_useable_moni_corr(&CM->nmon, &CM->ncor, &CM->mon_index, &CM->umoni, &CM->ucorr, 
                                 &CM->kick_coef, &CM->sl_index, coord, SL, run, beamline, 1);

#ifdef DEBUG
  fprintf(stdout, "finding twiss parameters beginning at %s.\n", start->name);
  fflush(stdout);
#endif
  if (!(beamline->flags&BEAMLINE_TWISS_CURRENT)) {
    if (verbose) {
      fprintf(stdout, "updating twiss parameters...");
      fflush(stdout);
    }
    update_twiss_parameters(run, beamline, NULL);

#ifdef DEBUG
    fprintf(stdout, "Tunes: %e, %e\n", beamline->tune[0], beamline->tune[1]);
    fprintf(stdout, "Initial eta: %e, %e\n", 
            beamline->elem.twiss->etax, 
            beamline->elem.twiss->etay);
    fprintf(stdout, "Final eta: %e, %e\n", 
            beamline->elast->twiss->etax, 
            beamline->elast->twiss->etay);
    fflush(stdout);
#endif
    if (verbose)
      report_stats(stdout, "\ndone: ");
  }

#ifdef DEBUG
  fprintf(stdout, "monitors: ");
  fflush(stdout);
  for (i_moni=0; i_moni<CM->nmon; i_moni++)
    fprintf(stdout, "%s ", CM->umoni[i_moni]->name);
  fflush(stdout);
  fprintf(stdout, "\ncorrectors: ");
  fflush(stdout);
  for (i_corr=0; i_corr<CM->ncor; i_corr++)
    fprintf(stdout, "%s ", CM->ucorr[i_corr]->name);
  fflush(stdout);
  fputc('\n', stdout);
#endif

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction may be unstable (use SV controls).\n");
    fflush(stdout);
  }
  if (CM->ncor==0) {
    fprintf(stdout, "Warning: no correctors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    return;
  }
  if (CM->nmon==0) {
    fprintf(stdout, "Warning: no monitors for %c plane.  No correction done.\n",  (coord==0?'x':'y'));
    fflush(stdout);
    CM->ncor = 0;
    return;
  }

  /* allocate matrices for this plane */
  if (CM->C)
    matrix_free(CM->C);
  CM->C  = matrix_get(CM->nmon, CM->ncor);   /* Response matrix */
  if (CM->T)
    matrix_free(CM->T);
  CM->T  = NULL;

  if (find_only)
    return;

  /* set up weight matrix */
  equalW = 1;
  weight = tmalloc(sizeof(*weight)*CM->nmon);
  for (i_moni=0; i_moni<CM->nmon; i_moni++) {
    weight[i_moni] = getMonitorWeight(CM->umoni[i_moni]);
    if (!i_moni)
      W0 = weight[i_moni];
    else if (weight[i_moni]!=W0)
      equalW = 0;
  }

  /* find transfer matrix from correctors to monitors */
  if ((coef = 2*sin(htune=PIx2*beamline->tune[coord?1:0]/2))==0)
    bomb("can't compute response matrix--beamline unstable", NULL);
  coefFL = beamline->matrix->R[4][5];
  
  corrFactor   = tmalloc(sizeof(*corrFactor)*CM->ncor);
  corrFactorFL = tmalloc(sizeof(*corrFactorFL)*CM->ncor);
  if (verbose) {
    fprintf(stdout, "computing orbit response matrix...");
    fflush(stdout);
  }
  
  switch (coord) {
  case 0:
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      double betax, etax;
      betax = CM->ucorr[i_corr]->twiss->betax;
      etax = CM->ucorr[i_corr]->twiss->etax;
      if (CM->ucorr[i_corr]->pred) {
        betax = (betax + CM->ucorr[i_corr]->pred->twiss->betax)/2;
        etax = (etax + CM->ucorr[i_corr]->pred->twiss->etax)/2;
      }
      corrFactor[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*sqrt(betax)/coef;
      corrFactorFL[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*etax/coefFL;
    }
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      moniFactor = 
        getMonitorCalibration(CM->umoni[i_moni], coord)*sqrt(CM->umoni[i_moni]->twiss->betax);
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        double phi;
        phi = CM->ucorr[i_corr]->twiss->phix;
        if (CM->ucorr[i_corr]->pred)
          phi = (phi + CM->ucorr[i_corr]->pred->twiss->phix)/2;
        Mij(CM->C, i_moni, i_corr)
          = moniFactor*corrFactor[i_corr]*
            cos(htune-fabs(CM->umoni[i_moni]->twiss->phix - phi));
        if (fixed_length)
          Mij(CM->C, i_moni, i_corr) -= CM->umoni[i_moni]->twiss->etax*corrFactorFL[i_corr];
        sprintf(memName, "HR_%s#%ld_%s#%ld.%s",
                CM->umoni[i_moni]->name, CM->umoni[i_moni]->occurence,
                CM->ucorr[i_corr]->name, CM->ucorr[i_corr]->occurence, 
                SL->corr_param[CM->sl_index[i_corr]]);
        rpn_store(Mij(CM->C, i_moni, i_corr), NULL, rpn_create_mem(memName, 0));
      }
    }
    break;
  default:
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      double betay;
      betay = CM->ucorr[i_corr]->twiss->betay;
      if (CM->ucorr[i_corr]->pred)
        betay = (betay + CM->ucorr[i_corr]->pred->twiss->betay)/2;
      corrFactor[i_corr] = 
        getCorrectorCalibration(CM->ucorr[i_corr], coord)*sqrt(betay)/coef;
    }
    for (i_moni=0; i_moni<CM->nmon; i_moni++) {
      moniFactor = 
        getMonitorCalibration(CM->umoni[i_moni], coord)*sqrt(CM->umoni[i_moni]->twiss->betay);
      for (i_corr=0; i_corr<CM->ncor; i_corr++) {
        double phi;
        phi = CM->ucorr[i_corr]->twiss->phiy;
        if (CM->ucorr[i_corr]->pred)
          phi = (phi + CM->ucorr[i_corr]->pred->twiss->phiy)/2;
        Mij(CM->C, i_moni, i_corr) 
          = moniFactor*corrFactor[i_corr]*
            cos(htune-fabs(CM->umoni[i_moni]->twiss->phiy - phi));
        sprintf(memName, "VR_%s#%ld_%s#%ld.%s",
                CM->umoni[i_moni]->name, CM->umoni[i_moni]->occurence,
                CM->ucorr[i_corr]->name, CM->ucorr[i_corr]->occurence,
                SL->corr_param[CM->sl_index[i_corr]]);
        rpn_store(Mij(CM->C, i_moni, i_corr), NULL, rpn_create_mem(memName, 0));
      }
    }
  }
  free(corrFactor);
  if (verbose) {
    report_stats(stdout, "\ndone");
    fflush(stdout);
  }
#ifdef DEBUG
  matrix_show(CM->C    , "%13.6le ", "influence matrix\n", stdout);
#endif

  if (invert) {
    /* compute correction matrix T */
    if (verbose) {
      fprintf(stdout, "computing correction matrix...");
      fflush(stdout);
    }
    if (CM->auto_limit_SVs && (CM->C->m < CM->C->n) && CM->remove_smallest_SVs < (CM->C->n - CM->C->m)) {
      CM->remove_smallest_SVs = CM->C->n - CM->C->m;
      printf("Removing %ld smallest singular values to prevent instability\n", (long)CM->remove_smallest_SVs);
    }
    CM->T = matrix_invert(CM->C, equalW?NULL:weight, (int32_t)CM->keep_largest_SVs, (int32_t)CM->remove_smallest_SVs,
                          CM->minimum_SV_ratio, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &conditionNumber);
    matrix_scmul(CM->T, -1);

    if (weight)
      free(weight);
    if (verbose) {
      report_stats(stdout, "\ndone.");
      printf("Condition number is %e\n", conditionNumber);
      fflush(stdout);
    }
#ifdef DEBUG
    matrix_show(CM->T, "%13.6le ", "correction matrix\n", stdout);
#endif
  }

}

long orbcor_plane(CORMON_DATA *CM, STEERING_LIST *SL, long coord, TRAJECTORY **orbit, long n_iterations, 
                  double clorb_acc, long clorb_iter, double clorb_iter_frac, RUN *run, LINE_LIST *beamline, double *closed_orbit, double *Cdp)
{
  ELEMENT_LIST *corr, *eptr;
  TRAJECTORY *clorb;
  VMATRIX *M;
  long iteration, kick_offset;
  long i_moni, i_corr, i, sl_index;
  double dp, x, y, reading;
  double last_rms_pos, best_rms_pos, rms_pos, corr_fraction;
  double fraction, minFraction, param, change;
  MAT *Qo, *dK;
  
  log_entry("orbcor_plane");

  if (!CM)
    bomb("NULL CORMON_DATA pointer passed to orbcor_plane", NULL);
  if (!orbit)
    bomb("NULL TRAJECTORY pointer passed to orbcor_plane", NULL);
  if (!run)
    bomb("NULL RUN pointer passed to orbcor_plane", NULL);
  if (!beamline)
    bomb("NULL LINE_LIST pointer passed to orbcor_plane", NULL);
  if (CM->ncor && CM->nmon && !CM->T)
    bomb("No inverse matrix computed prior to orbcor_plane", NULL);

  if (closed_orbit)
    dp = closed_orbit[5];
  else
    dp = 0;

  clorb = orbit[0];
  if (closed_orbit) {
    for (i=0; i<4; i++)
      clorb[0].centroid[i] = closed_orbit[i];
  }
  else {
    for (i=0; i<4; i++)
      clorb[0].centroid[i] = 0;
  }

  if (CM->nmon<CM->ncor) {
    fprintf(stdout, "*** Warning: more correctors than monitors for %c plane.\n",  (coord==0?'x':'y'));
    fprintf(stdout, "*** Correction may be unstable (use SV controls)\n");
    fflush(stdout);
  }

  Qo = matrix_get(CM->nmon, 1);   /* Vector of BPM errors */
   
  best_rms_pos = rms_pos = DBL_MAX/4;
  corr_fraction = CM->corr_fraction;
  for (iteration=0; iteration<=n_iterations; iteration++) {
    clorb = orbit[iteration?1:0];
    delete_phase_references();
    if (!(M = beamline->matrix)) {
      if (beamline->elem_twiss)
        M = beamline->matrix = full_matrix(beamline->elem_twiss, run, 1);
      else
        M = beamline->matrix = full_matrix(&(beamline->elem), run, 1);
    }
    if (!M || !M->R)
      bomb("problem calculating full matrix (orbcor_plane)", NULL);
    if (iteration==1)
      for (i=0; i<6; i++)
        orbit[1][0].centroid[i] = orbit[0][0].centroid[i];
    if (!find_closed_orbit(clorb, clorb_acc, clorb_iter, beamline, M, run, dp, 1, CM->fixed_length, NULL, 
                           clorb_iter_frac, NULL)) {
      fprintf(stdout, "Failed to find closed orbit.\n");
      fflush(stdout);
      return(-1);
    }
    if (Cdp)
      Cdp[iteration] = clorb[0].centroid[5];
    if (closed_orbit) {
      for (i=0; i<6; i++) 
        closed_orbit[i] = clorb[0].centroid[i];
    }

    /* find readings at monitors and add in reading errors */
    i = 1;
    last_rms_pos = rms_pos;
    if (best_rms_pos>rms_pos)
      best_rms_pos = rms_pos;
    rms_pos = 0;
    for (i_moni=0; i_moni<CM->nmon; i_moni++, i++) {
      while (CM->umoni[i_moni]!=clorb[i].elem) {
        if (clorb[i].elem->succ)
          i++;
        else
          bomb("missing monitor in closed orbit", NULL);
      }

      x = clorb[i].centroid[0];
      y = clorb[i].centroid[2];
      eptr = clorb[i].elem;
      reading = computeMonitorReading(eptr, coord, x, y, 0);

      Mij(Qo, i_moni, 0) = CM->posi[iteration][i_moni] = reading + 
        (CM->bpm_noise?noise_value(CM->bpm_noise, CM->bpm_noise_cutoff, CM->bpm_noise_distribution):0.0);
      rms_pos += sqr(Mij(Qo, i_moni, 0));
      if (!clorb[i].elem->succ)
        break;
    }
    rms_pos = sqrt(rms_pos/CM->nmon);
    if (iteration==0 && rms_pos>1e9) {
      /* if the closed orbit has RMS > 1e9m, I assume correction won't work and routine bombs */
      fprintf(stdout, "Orbit beyond 10^9 m.  Aborting correction.\n");
      fflush(stdout);
      return(-1);
    }
    if (rms_pos>best_rms_pos*1.01) {
      if (corr_fraction==1 || corr_fraction==0)
        break;            
      fprintf(stdout, "orbit diverging on iteration %ld: RMS=%e m (was %e m)--redoing with correction fraction of %e\n", 
              iteration, rms_pos, best_rms_pos, corr_fraction = sqr(corr_fraction));
      best_rms_pos = rms_pos;
      fflush(stdout);
#if 0
      if (iteration>=1) {
        /* step through beamline and change correctors back to last values */
        for (i_corr=0; i_corr<CM->ncor; i_corr++) {
          corr = CM->ucorr[i_corr];
          sl_index = CM->sl_index[i_corr];
          kick_offset = SL->param_offset[sl_index];
          *((double*)(corr->p_elem+kick_offset)) = CM->kick[iteration-1][i_corr]/CM->kick_coef[i_corr];
          if (corr->matrix) {
            free_matrices(corr->matrix);
            tfree(corr->matrix);
            corr->matrix = NULL;
          }
          compute_matrix(corr, run, NULL);
        }
      }
      if (beamline->links)
        assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

      /* free concatenated matrix to force it to be recomputed with new corrector matrices */
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
      iteration -= 1;
      rms_pos = last_rms_pos;
      continue;
#endif
    }

    if (CM->fixed_length)
      dp = clorb[0].centroid[5];

    if (iteration==n_iterations)
      break;

    if (CM->nmon && CM->ncor)
      /* solve for the corrector kicks */
      dK = matrix_mult(CM->T, Qo);

#ifdef DEBUG
    matrix_show(Qo, "%13.6le ", "traj matrix\n", stdout);
    matrix_show(dK, "%13.6le ", "kick matrix\n", stdout);
#endif

    /* see if any correctors are over their limit */
    minFraction = 1;
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      param = fabs(*((double*)(corr->p_elem+kick_offset)) +
                   (change=Mij(dK, i_corr, 0)/CM->kick_coef[i_corr]*corr_fraction));
      if (SL->corr_limit[sl_index] && param>SL->corr_limit[sl_index]) {
        fraction = fabs((SL->corr_limit[sl_index]-fabs(*((double*)(corr->p_elem+kick_offset))))/change);
        if (fraction<minFraction)
          minFraction = fraction;
      }
    }
    fraction = minFraction*corr_fraction;
    
    /* step through beamline and change correctors */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      if (iteration==0) 
        CM->kick[iteration][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      *((double*)(corr->p_elem+kick_offset)) += Mij(dK, i_corr, 0)/CM->kick_coef[i_corr]*fraction;
      CM->kick[iteration+1][i_corr] = *((double*)(corr->p_elem+kick_offset))*CM->kick_coef[i_corr];
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

    /* free concatenated matrix to force it to be recomputed with new corrector matrices */
    free_matrices(beamline->matrix);
    tfree(beamline->matrix);
    beamline->matrix = NULL;
    matrix_free(dK);
  }
  

  if (rms_pos>1e9) {
    /* if the final closed orbit has RMS > 1e9m, I assume correction didn't work and routine bombs */
    fprintf(stdout, "Orbit beyond 1e9 m.  Aborting correction.\n");
    fflush(stdout);
    return(-1);
  }

  if (rms_pos>best_rms_pos+CM->corr_accuracy && (corr_fraction==1 || corr_fraction==0)) {
    fprintf(stdout, "orbit not improving--iteration terminated at iteration %ld with last result of %e m RMS\n",
            iteration, rms_pos);
    fflush(stdout);
    /* step through beamline and change correctors back to last values */
    for (i_corr=0; i_corr<CM->ncor; i_corr++) {
      corr = CM->ucorr[i_corr];
      sl_index = CM->sl_index[i_corr];
      kick_offset = SL->param_offset[sl_index];
      *((double*)(corr->p_elem+kick_offset)) = CM->kick[iteration][i_corr]/CM->kick_coef[i_corr];
      if (corr->matrix) {
        free_matrices(corr->matrix);
        tfree(corr->matrix);
        corr->matrix = NULL;
      }
      compute_matrix(corr, run, NULL);
    }
    if (beamline->links)
      assert_element_links(beamline->links, run, beamline, DYNAMIC_LINK);

    /* free concatenated matrix to force it to be recomputed with new corrector matrices */
    free_matrices(beamline->matrix);
    tfree(beamline->matrix);
    beamline->matrix = NULL;
  }

  /* indicate that beamline concatenation and Twiss computation (if wanted) are not current */
  beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
  beamline->flags &= ~BEAMLINE_TWISS_CURRENT;

  matrix_free(Qo);
  
  log_exit("orbcor_plane");
  return(iteration);
}


ELEMENT_LIST *next_element_of_type(ELEMENT_LIST *elem, long type)
{
  while (elem) {
    if (elem->type==type)
      return(elem);
    elem = elem->succ;
  }
  return(elem);
}

ELEMENT_LIST *next_element_of_type2(ELEMENT_LIST *elem, long type1, long type2)
{
  while (elem) {
    if (elem->type==type1 || elem->type==type2)
      return(elem);
    elem = elem->succ;
  }
  return(elem);
}

ELEMENT_LIST *next_element_of_types(ELEMENT_LIST *elem, long *type, long n_types, long *index)
{
  register long i;
  while (elem) {
    for (i=0; i<n_types; i++) {
      if (elem->type==type[i]) {
        *index = i;
        return(elem);
      }
    }
    elem = elem->succ;
  }
  return(elem);
}

long find_parameter_offset(char *param_name, long elem_type)
{
  long param;
  if ((param=confirm_parameter(param_name, elem_type))<0)
    return(-1);
  return(entity_description[elem_type].parameter[param].offset);
}


long zero_correctors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  return 
    zero_hcorrectors(elem, run, correct) +
      zero_vcorrectors(elem, run, correct);
}

long zero_correctors_one_plane(ELEMENT_LIST *elem, RUN *run, STEERING_LIST *SL, long plane)
{  
  long n_zeroed = 0, i;
  long paramOffset;

  while (elem) {
    paramOffset = -1;
    if (steering_corrector(elem, SL, plane)) {
      for (i=0; i<SL->n_corr_types; i++)
	if (strcmp(elem->name, SL->corr_name[i])==0)
	  break;
      if (i!=SL->n_corr_types) {
	paramOffset = SL->param_offset[i];
	*((double*)(elem->p_elem+paramOffset)) = 0;
	if (elem->matrix) {
	  free_matrices(elem->matrix);
	  tfree(elem->matrix);
	  elem->matrix = NULL;
	}
	compute_matrix(elem, run, NULL);
	n_zeroed++;
      }
    }
    elem = elem->succ;
  }
  return(n_zeroed);
}

long zero_hcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  long nz;
  nz = zero_correctors_one_plane(elem, run, &(correct->SLx), 0);
  fprintf(stderr, "%ld H correctors set to zero\n", nz);
  return nz;
}

long zero_vcorrectors(ELEMENT_LIST *elem, RUN *run, CORRECTION *correct)
{
  long nz;
  nz = zero_correctors_one_plane(elem, run, &(correct->SLy), 1);
  fprintf(stderr, "%ld V correctors set to zero\n", nz);
  return nz;
}

void rotate_xy(double *x, double *y, double angle)
{
    static double x0, y0, ca, sa;

    x0 = *x;
    y0 = *y;
    *x  =  x0*(ca=cos(angle)) + y0*(sa=sin(angle));
    *y  = -x0*sa              + y0*ca;
    }

double rms_value(double *data, long n_data)
{
    double sum2;
    long i;

    for (i=sum2=0; i<n_data; i++) 
        sum2 += sqr(data[i]);
    if (n_data)
        return(sqrt(sum2/n_data));
    else
        return(0.0);
    }

long steering_corrector(ELEMENT_LIST *eptr, STEERING_LIST *SL, long plane)
{
  long i;
  for (i=0; i<SL->n_corr_types; i++)
    if (strcmp(eptr->name, SL->corr_name[i])==0) {
      switch (eptr->type) {
      case T_HVCOR:
        return ((HVCOR*)(eptr->p_elem))->steering;
      case T_HCOR:
        return ((HCOR*)(eptr->p_elem))->steering;
      case T_VCOR:
        return ((VCOR*)(eptr->p_elem))->steering;
      case T_QUAD:
        if (plane==0)
          return ((QUAD*)(eptr->p_elem))->xSteering;
        return ((QUAD*)(eptr->p_elem))->ySteering;
      case T_KQUAD:
        if (plane==0)
          return ((KQUAD*)(eptr->p_elem))->xSteering;
        return ((KQUAD*)(eptr->p_elem))->ySteering;
      default:
        return 1;
      }
    }
  return 0;
}

long find_index(long key, long *list, long n_listed)
{
    register long i;
    for (i=0; i<n_listed; i++)
        if (key==list[i])
            return(i);
    return(-1);
    }

double noise_value(double xamplitude, double xcutoff, long xerror_type)
{
    switch (xerror_type) {
        case UNIFORM_ERRORS:
            return(2*xamplitude*(random_3(0)-0.5));
        case GAUSSIAN_ERRORS:
            return(gauss_rn_lim(0.0, xamplitude, xcutoff, random_3));
        case PLUS_OR_MINUS_ERRORS:
            /* return a number on [-x-1, x-1]), which is added to 1 in the calling routine
             * (since these are implemented as fractional errors)
             */
            return(xamplitude*(random_3(0)>0.5?1.0:-1.0) - 1);
        default:
            bomb("unknown error type in perturbation()", NULL);
            exit(1);
            break;
        }
    return(0.0);
    }

double computeMonitorReading(ELEMENT_LIST *elem, long coord, double x, double y,
                             unsigned long flags)
/* coord = 0 is x, otherwise y */
{
  double calibration, tilt, reading;
  char *equation;

  switch (elem->type) {
  case T_MONI:  
    x -= ((MONI*)(elem->p_elem))->dx;
    y -= ((MONI*)(elem->p_elem))->dy;
    if (coord==0)
      calibration = ((MONI*)(elem->p_elem))->xcalibration;
    else 
      calibration = ((MONI*)(elem->p_elem))->ycalibration;
    tilt = ((MONI*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    if (coord==0)
      equation = ((MONI*)(elem->p_elem))->x_readout; 
    else
      equation = ((MONI*)(elem->p_elem))->y_readout;
    break;
  case T_HMON: 
    x -= ((HMON*)(elem->p_elem))->dx;
    y -= ((HMON*)(elem->p_elem))->dy;
    calibration = ((HMON*)(elem->p_elem))->calibration;
    tilt = ((HMON*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    equation = ((HMON*)(elem->p_elem))->readout; 
    if (coord!=0)
      bomb("element in horizontal monitor list is not a vertical monitor--internal logic error", NULL);
    break;
  case T_VMON:
    x -= ((VMON*)(elem->p_elem))->dx;
    y -= ((VMON*)(elem->p_elem))->dy;
    calibration = ((VMON*)(elem->p_elem))->calibration;
    tilt = ((VMON*)(elem->p_elem))->tilt;
    if (flags&COMPUTEMONITORREADING_TILT_0)
      tilt = 0;
    if (tilt)
      rotate_xy(&x, &y, tilt);   
    equation = ((VMON*)(elem->p_elem))->readout; 
    if (!coord)
      bomb("element in vertical monitor list is not a vertical monitor--internal logic error", NULL);
    break;
  default:
    fprintf(stdout, "error: element %s found in monitor list--internal logic error\n", 
            elem->name);
    fflush(stdout);
    abort();
    break;
  }
  
  if (flags&COMPUTEMONITORREADING_CAL_1)
    calibration = 1;

  if (equation) {
    rpn_clear();
    rpn_store(x, NULL, rpn_x_mem);
    rpn_store(y, NULL, rpn_y_mem);
    reading = rpn(equation)*calibration;
    if (rpn_check_error()) exit(1);
  }
  else {
    switch (coord) {
    case 0: 
      reading = x*calibration;
      break;
    default:
      reading = y*calibration;
      break;
    }
  }
  return reading;
}

double getMonitorWeight(ELEMENT_LIST *elem)
{
    switch (elem->type) {
      case T_MONI:
        return ((MONI*)(elem->p_elem))->weight;
      case T_HMON:
        return ((HMON*)(elem->p_elem))->weight;
      case T_VMON:
        return ((VMON*)(elem->p_elem))->weight;
        }
    bomb("invalid element type in getMonitorWeight()", NULL);
	return(1);
    }

double getMonitorCalibration(ELEMENT_LIST *elem, long coord)
{
    switch (elem->type) {
      case T_HMON:
        return ((HMON*)(elem->p_elem))->calibration;
      case T_VMON:
        return ((VMON*)(elem->p_elem))->calibration;
      case T_MONI:
        if (coord)
            return ((MONI*)(elem->p_elem))->ycalibration;
        return ((MONI*)(elem->p_elem))->xcalibration;
        }
    bomb("invalid element type in getMonitorCalibration()", NULL);
	return(1);
    }

void setMonitorCalibration(ELEMENT_LIST *elem, double calib, long coord)
{
  switch (elem->type) {
  case T_HMON:
    ((HMON*)(elem->p_elem))->calibration = calib;
    break;
  case T_VMON:
    ((VMON*)(elem->p_elem))->calibration = calib;
    break;
  case T_MONI:
    if (coord)
      ((MONI*)(elem->p_elem))->ycalibration = calib;
    else 
      ((MONI*)(elem->p_elem))->xcalibration = calib;
    break;
  default:
    bomb("invalid element type in setMonitorCalibration()", NULL);
    break;
  }
}


double getCorrectorCalibration(ELEMENT_LIST *elem, long coord)
{
    switch (elem->type) {
      case T_HCOR:
        return ((HCOR*)(elem->p_elem))->calibration;
      case T_VCOR:
        return ((VCOR*)(elem->p_elem))->calibration;
      case T_HVCOR:
        if (coord)
            return ((HVCOR*)(elem->p_elem))->ycalibration;
        return ((HVCOR*)(elem->p_elem))->xcalibration;
      case T_QUAD:
        if (coord) 
          return ((QUAD*)(elem->p_elem))->yKickCalibration;
        return ((QUAD*)(elem->p_elem))->xKickCalibration;
      case T_KQUAD:
        if (coord) 
          return ((KQUAD*)(elem->p_elem))->yKickCalibration;
        return ((KQUAD*)(elem->p_elem))->xKickCalibration;
      default:
        return 1;
        }
    }
