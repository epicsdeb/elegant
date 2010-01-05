/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: vary.c
 * contents: vary_setup(), vary_beamline()
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include "vary.h"

long load_enumerated_values(double **value, char *file, char *column);
void reset_parameter_values(char **elem_name, long *param_number, long *type, long n_elems,
                            LINE_LIST *beamline);

void vary_setup(VARY *_control, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{

    log_entry("vary_setup");

    /* assert defaults (necessary for multiple step runs) */
/*
    n_indices = bunch_frequency = 0;
    n_passes = n_steps = 1;
*/

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&run_control, nltext);
    if (echoNamelists) print_namelist(stdout, &run_control);

    /* check validity of input values */
    if (n_steps<=0 && n_indices<=0)
        bomb("n_steps <= 0  and  n_indices <= 0", NULL);
    if (n_indices<=0) {
        if (bunch_frequency) {
            if (bunch_frequency<0)
                bomb("bunch_frequency<0", NULL);
            }
        }
        
    /* reset flags on elements (necessary for multiple step runs) */
    if (_control->n_elements_to_vary)
        set_element_flags(beamline, _control->element, NULL, NULL, NULL, _control->n_elements_to_vary, 
                PARAMETERS_ARE_STATIC, 0, 1, 0);
    
    /* copy data into run control structure */
    if ((_control->n_indices = n_indices)) {
        _control->index_limit = tmalloc(sizeof(*_control->index_limit)*_control->n_indices);
        fill_long_array(_control->index_limit, _control->n_indices, 0);
        _control->index       = tmalloc(sizeof(*_control->index)*_control->n_indices);
        fill_long_array(_control->index, _control->n_indices, 0);
        }
    _control->i_step  = 0;
    _control->n_steps = n_steps;
    _control->bunch_frequency = bunch_frequency;
    _control->n_passes = n_passes;
    _control->reset_rf_each_step = reset_rf_for_each_step;
    if (first_is_fiducial && n_passes!=1)
      bomb("can't have fiducial beam and multiple passes", NULL);
    _control->fiducial_flag = 0;
    if (first_is_fiducial)
      _control->fiducial_flag = FIRST_BEAM_IS_FIDUCIAL |
        (restrict_fiducialization?RESTRICT_FIDUCIALIZATION:0); 
    
    /* reset flags for elements that may have been varied previously */
    if (_control->n_elements_to_vary) {
        set_element_flags(beamline, _control->element, NULL, NULL, NULL, _control->n_elements_to_vary,
                PARAMETERS_ARE_STATIC, 0, 1, 0);
        if (_control->cell)
            set_element_flags(_control->cell,_control->element, NULL, NULL, NULL, _control->n_elements_to_vary,
                    PARAMETERS_ARE_STATIC, 0, 1, 0);
        }

    _control->cell = NULL;
    _control->n_elements_to_vary = 0;
    _control->new_data_read = 0;
    _control->at_start = 1;
    _control->i_vary = 0;
    log_exit("vary_setup");
    }

void add_varied_element(VARY *_control, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long n_elements_to_vary;
    ELEMENT_LIST *context;
    double value;

    log_entry("add_varied_element");

    if (_control->n_indices<=0)
        bomb("can't vary an element if n_indices==0 in run_control namelist", NULL);
    
    if ((n_elements_to_vary = _control->n_elements_to_vary)==0) {
        if (_control->new_data_read)
            bomb("improper sequencing of variation and tracking", NULL);
        _control->new_data_read = 1;
        }

    _control->element_index = trealloc(_control->element_index, sizeof(*_control->element_index)*(n_elements_to_vary+1));
    _control->element       = trealloc(_control->element, sizeof(*_control->element)*(n_elements_to_vary+1));
    _control->item          = trealloc(_control->item, sizeof(*_control->item)*(n_elements_to_vary+1));
    _control->initial       = trealloc(_control->initial, sizeof(*_control->initial)*(n_elements_to_vary+1));
    _control->final         = trealloc(_control->final, sizeof(*_control->final)*(n_elements_to_vary+1));
    _control->step          = trealloc(_control->step, sizeof(*_control->step)*(n_elements_to_vary+1));
    _control->enumerated_value = trealloc(_control->enumerated_value, sizeof(*_control->enumerated_value)*(n_elements_to_vary+1));
    _control->varied_quan_name = 
            trealloc(_control->varied_quan_name, sizeof(*_control->varied_quan_name)*(n_elements_to_vary+1));
    _control->varied_quan_unit = 
            trealloc(_control->varied_quan_unit, sizeof(*_control->varied_quan_unit)*(n_elements_to_vary+1));
    _control->varied_type   = trealloc(_control->varied_type, sizeof(*_control->varied_type)*(n_elements_to_vary+1));
    _control->varied_quan_value = 
            trealloc(_control->varied_quan_value, sizeof(*_control->varied_quan_value)*(n_elements_to_vary+1));
    _control->varied_param  = trealloc(_control->varied_param, sizeof(*_control->varied_param)*(n_elements_to_vary+1));
    _control->flags  = trealloc(_control->flags, sizeof(*_control->flags)*(n_elements_to_vary+1));

    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    process_namelist(&vary_element, nltext);
    if (echoNamelists) print_namelist(stdout, &vary_element);

    /* check for valid input */
    if (index_number<0 || index_number>=_control->n_indices)
        bomb("invalid index_number", NULL);
    _control->element_index[n_elements_to_vary] = index_number;
    _control->enumerated_value[n_elements_to_vary] = NULL;
    if (enumeration_file) {
        if (!enumeration_column)
            bomb("must supply enumeration_column with enumeration_file", NULL);
        if (!(index_limit = load_enumerated_values(_control->enumerated_value+n_elements_to_vary, 
                                   enumeration_file, enumeration_column)))
            bomb("enumerated_values_file contains no valid values", NULL);
        fprintf(stdout, "%ld values of %s loaded from file %s\n", index_limit, enumeration_column,
               enumeration_file);
        fflush(stdout);
        if (_control->index_limit[index_number]) {
            fprintf(stdout, "Warning: the limit for index %ld is specified more than once.\nThe lowest-valued specification is used.\n", index_number);
            fflush(stdout);
            if (index_limit>_control->index_limit[index_number])
                index_limit = _control->index_limit[index_number];
            }
        initial = _control->enumerated_value[n_elements_to_vary][0];
        _control->index_limit[index_number] = index_limit;
        }
    else {
        if (!_control->index_limit[index_number]) {
            if (index_limit>0)
                _control->index_limit[index_number] = index_limit;
            }
        else if (index_limit>0) {
            fprintf(stdout, "Warning: the limit for index %ld is specified more than once.\nThe first specification is used.\n",
                   index_number);
            fflush(stdout);
            }
        }

    if (name==NULL)
        bomb("element name missing in vary_element namelist", NULL);
    str_toupper(name);
    context = NULL;
    if (!find_element(name, &context, &(beamline->elem))) {
        fprintf(stdout, "error: cannot vary element %s--not in beamline\n", name);
        fflush(stdout);
        exit(1);
        }
    cp_str(&_control->element[n_elements_to_vary], name);
    _control->varied_type[n_elements_to_vary] = context->type;
    str_toupper(item);
    if ((_control->varied_param[n_elements_to_vary] = confirm_parameter(item, context->type))<0) {
        fprintf(stdout, "error: cannot vary %s--no such parameter for %s\n",item, name);
        fflush(stdout);
        exit(1);
        }
    cp_str(&_control->item[n_elements_to_vary], item);
    cp_str(&_control->varied_quan_unit[n_elements_to_vary], 
        entity_description[_control->varied_type[n_elements_to_vary]].parameter[_control->varied_param[n_elements_to_vary]].unit);
    if (geometric && (multiplicative || differential)) {
      fprintf(stdout, "error: geometric is incompatible with multiplicative and differential modes\n");
      exit(1);
    }

    if (multiplicative) {
        if (!get_parameter_value(&value, name, _control->varied_param[n_elements_to_vary], 
                                 _control->varied_type[n_elements_to_vary], beamline))
            bomb("unable to get preset parameter value for item", NULL);
        if (enumeration_file) {
            long i;
            for (i=0; i<index_limit; i++)
                _control->enumerated_value[n_elements_to_vary][i] *= value;
            initial *= value;
            }
        else {
            initial *= value;
            final   *= value;
            }
        }
    if (differential) {
        if (!get_parameter_value(&value, name, _control->varied_param[n_elements_to_vary], 
                                 _control->varied_type[n_elements_to_vary], beamline))
            bomb("unable to get preset parameter value for item", NULL);
        if (enumeration_file) {
            long i;
            for (i=0; i<index_limit; i++)
                _control->enumerated_value[n_elements_to_vary][i] += value;
            initial += value;
            }
        else {
            initial += value;
            final   += value;
            }
        }
    
    _control->initial[n_elements_to_vary] = initial;
    _control->varied_quan_value[n_elements_to_vary] = initial;
    _control->final[n_elements_to_vary]   = final;
    _control->step[n_elements_to_vary]    = 0;        /* will be set after all vary data is read in */
    _control->varied_quan_name[n_elements_to_vary] = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+3));
    sprintf(_control->varied_quan_name[n_elements_to_vary], "%s.%s", name, item);

    _control->flags[n_elements_to_vary] = 0;
    if (geometric && !enumeration_file) {
        _control->flags[n_elements_to_vary] |= VARY_GEOMETRIC;
        if (initial==0 || SIGN(initial)!=SIGN(final))
            bomb("must have initial!=0 and SIGN(initial)=SIGN(final) for geometric variation", NULL);
        }

    _control->n_elements_to_vary = n_elements_to_vary+1;

    log_exit("add_varied_element");
    }


long vary_beamline(VARY *_control, ERRORVAL *errcon, RUN *run, LINE_LIST *beamline)
{
  long i, ret_val, do_perturbations, step_incremented, parameters_loaded;
  ELEMENT_LINKS *links;

  log_entry("vary_beamline");
  links = beamline->links;

#if DEBUG
  fprintf(stdout, "vary_beamline called: i_step = %ld, i_vary = %ld\n", 
          _control->i_step, _control->i_vary);
  fflush(stdout);
#endif

  if ((_control->bunch_frequency==0 && _control->reset_rf_each_step) || 
      _control->i_step==0)
    delete_phase_references();
  reset_special_elements(beamline, _control->reset_rf_each_step);

  do_perturbations = step_incremented = 0;

  if (links && links->n_links) 
    reset_element_links(links, run, beamline);

  for (i=0, _control->indexLimitProduct=1; i<_control->n_indices; i++) {
    if ((_control->indexLimitProduct *= _control->index_limit[i])<=0) {
      fprintf(stdout, "index %ld has a limit of <= 0", i);
      fflush(stdout);
      exit(1);
    }
  }
  
  if (_control->n_indices==0 || _control->at_start) {
    do_perturbations = 1;
    if (errcon->n_items) {
#if DEBUG
      fputs("(re)asserting unperturbed values", stdout);
#endif
      log_entry("vary_beamline.1");
      /* assert unperturbed values */
      reset_parameter_values(errcon->name, errcon->param_number, errcon->elem_type,
                             errcon->n_items, beamline);
      if (errcon->new_data_read) {
        /* set element flags to indicate perturbation of parameters that change the matrix */
        set_element_flags(beamline, errcon->name, errcon->flags, errcon->elem_type, errcon->param_number,
                          errcon->n_items, PARAMETERS_ARE_PERTURBED, VMATRIX_IS_PERTURBED, 0, PRE_CORRECTION);
        errcon->new_data_read = 0;
      }
      log_exit("vary_beamline.1");
    }
  }
  if (_control->n_elements_to_vary) {
#if DEBUG
    fputs("inside vary-elements section", stdout);
#endif
    log_entry("vary_beamline.2");
    check_VARY_structure(_control, "vary_beamline");

    if (_control->new_data_read || _control->at_start) {
      /* assert initial values of varied parameters */
      assert_parameter_values(_control->element, _control->varied_param, _control->varied_type,
                              _control->initial, _control->n_elements_to_vary, beamline);
      if (_control->cell) 
        assert_parameter_values(_control->element, _control->varied_param, _control->varied_type,
                                _control->initial, _control->n_elements_to_vary, _control->cell);
      /* check for zero index_limits */
      for (i=0; i<_control->n_indices; i++) {
        if (_control->index_limit[i]==0)
          bomb("index limit must be given for each index", NULL);
        _control->index[i] = 0;
        _control->varied_quan_value[i] = _control->initial[i];
      }
      /* calculate step sizes after checking that data is okay for all indices */
      for (i=0; i<_control->n_elements_to_vary; i++) {
        if (_control->index_limit[_control->element_index[i]]>1) {
          if (_control->flags[i]&VARY_GEOMETRIC)
            _control->step[i] = pow(_control->final[i]/_control->initial[i], 
                                    1./(_control->index_limit[_control->element_index[i]]-1));
          else
            _control->step[i] = (_control->final[i]-_control->initial[i])/
              (_control->index_limit[_control->element_index[i]]-1);
        }
        else if (_control->index_limit[_control->element_index[i]]==1)
          _control->step[i] = 0;
        else
          bomb("index limits must be set for all indices", NULL);
      }
      /* set element flags to indicate variation of parameters */
      set_element_flags(beamline, _control->element, NULL, _control->varied_type, _control->varied_param, 
                        _control->n_elements_to_vary, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
      if (_control->cell)
        set_element_flags(_control->cell, _control->element, NULL, _control->varied_type, _control->varied_param, 
                          _control->n_elements_to_vary, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
      _control->new_data_read = _control->at_start = 0;
      _control->i_vary = 1;
      _control->i_step++;
      step_incremented = 1;
      fprintf(stdout, "vary counter reset\n");
      fflush(stdout);
    }
    else {
#if DEBUG
      fputs("calling advance_values", stdout);
#endif
      if (advance_values1(_control->varied_quan_value, _control->n_elements_to_vary, _control->element_index,
                          _control->initial, _control->step, _control->enumerated_value, 
                          _control->index, _control->index_limit, _control->flags, _control->n_indices)<0)  {
        if (_control->n_steps && _control->i_step>=_control->n_steps) {
          log_exit("vary_beamline.2");
          log_exit("vary_beamline");
          return(0);
        }
        _control->at_start = 1;
        ret_val = vary_beamline(_control, errcon, run, beamline);
        log_exit("vary_beamline.2");
        log_exit("vary_beamline");
        return(ret_val);
      }
      fprintf(stdout, "counter advanced: ");
      fflush(stdout);
      for (i=0; i<_control->n_indices; i++)
        fprintf(stdout, "%4ld ", _control->index[i]);
      fflush(stdout);
      fprintf(stdout, "\nvalues advanced: ");
      fflush(stdout);
      for (i=0; i<_control->n_elements_to_vary; i++)
        fprintf(stdout, "%e ", _control->varied_quan_value[i]);
      fflush(stdout);
      fprintf(stdout, "\n");
      fflush(stdout);
      assert_parameter_values(_control->element, _control->varied_param, _control->varied_type,
                              _control->varied_quan_value, _control->n_elements_to_vary, beamline);
      if (_control->cell)
        assert_parameter_values(_control->element, _control->varied_param, _control->varied_type,
                                _control->varied_quan_value, _control->n_elements_to_vary, _control->cell);
      /* set element flags to indicate variation of parameters */
      set_element_flags(beamline, _control->element, NULL, _control->varied_type, _control->varied_param, 
                        _control->n_elements_to_vary, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
      if (_control->cell)
        set_element_flags(_control->cell, _control->element, NULL, _control->varied_type, _control->varied_param, 
                          _control->n_elements_to_vary, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
      _control->i_vary++;
    }
    log_exit("vary_beamline.2");
  }

  parameters_loaded = 0;
  if (!(_control->i_step>=_control->n_steps && !_control->n_indices))
    parameters_loaded = do_load_parameters(beamline, 0);

  if (errcon->n_items && do_perturbations) {
#if DEBUG
    fputs("doing perturbation", stdout);
#endif
    log_entry("vary_beamline.3");
    /* calculate random errors and add them to the existing value (which may not be the unperturbed value
       if the parameter is both varied and perturbed) */
    if (_control->n_steps && (_control->i_step-(step_incremented?1:0))>=_control->n_steps) {
      log_exit("vary_beamline.3");
      log_exit("vary_beamline");
      return(0);
    }
    assert_perturbations(errcon->name, errcon->param_number, errcon->elem_type,
                         errcon->n_items, errcon->error_level, errcon->error_cutoff, errcon->error_type, 
                         errcon->error_value, errcon->flags, errcon->bind_number,
                         errcon->boundTo, errcon->sMin, errcon->sMax,
                         errcon->fp_log, _control->i_step, beamline, 
                         PRE_CORRECTION+
                         ((_control->i_step==0 && errcon->no_errors_first_step)?FORCE_ZERO_ERRORS:0));
    /* set element flags to indicate perturbation of parameters that change the matrix */
    set_element_flags(beamline, errcon->name, errcon->flags, errcon->elem_type, errcon->param_number,
                      errcon->n_items, PARAMETERS_ARE_PERTURBED, VMATRIX_IS_PERTURBED, 0, PRE_CORRECTION);
    if (!step_incremented) {
      _control->i_step++;
      step_incremented = 1;
    }
    log_exit("vary_beamline.3");
  }        

  if (_control->n_elements_to_vary || errcon->n_items || parameters_loaded==PARAMETERS_LOADED) {
#if DEBUG
    fputs("computing matrices", stdout);
#endif
    log_entry("vary_beamline.4");
    /* compute matrices for perturbed elements */
    rebaseline_element_links(links, run, beamline);
    _control->i_step -= step_incremented;
    fprintf(stdout, "%ld matrices (re)computed\n", 
            (i=compute_changed_matrices(beamline, run)
             + assert_element_links(links, run, beamline, STATIC_LINK+DYNAMIC_LINK))
            + (_control->cell?compute_changed_matrices(_control->cell, run):0) );
    _control->i_step += step_incremented;
    fflush(stdout);
    if (i) {
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
      beamline->flags &= ~BEAMLINE_RADINT_CURRENT;
      
    }
    if (i && beamline->matrix) {
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
    }
    if (!step_incremented) {
      _control->i_step++;
      step_incremented = 1;
    }
    fprintf(stdout, "tracking step %ld.%ld\n", _control->i_step, _control->i_vary);
    fflush(stdout);
    log_exit("vary_beamline.4");
    log_exit("vary_beamline");
    return(1);
  }

  if ((_control->i_step>=_control->n_steps || parameters_loaded==PARAMETERS_ENDED) && !_control->n_indices) {
    log_exit("vary_beamline");
    return(0);
  }

  if (links && links->n_links) {
    rebaseline_element_links(links, run, beamline);
    fprintf(stdout, "%ld matrices (re)computed\n", i=assert_element_links(links, run, beamline, STATIC_LINK+DYNAMIC_LINK));
    fflush(stdout);
    if (i) {
      beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
      beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
      beamline->flags &= ~BEAMLINE_RADINT_CURRENT;
    }
    if (i && beamline->matrix) {
      free_matrices(beamline->matrix);
      tfree(beamline->matrix);
      beamline->matrix = NULL;
    }
  }

  fprintf(stdout, "tracking step %ld\n", ++_control->i_step);
  fflush(stdout);

  log_exit("vary_beamline");
  return(1);
}

long perturb_beamline(VARY *_control, ERRORVAL *errcon, RUN *run, LINE_LIST *beamline)
{
    long i;
    ELEMENT_LINKS *links;

    log_entry("perturb_beamline");
    links = beamline->links;

    if (errcon->n_items) {
        /* calculate random errors and add them to the existing value (which may not be the unperturbed value
           if the parameter is both varied and perturbed) */
        assert_perturbations(errcon->name, errcon->param_number, errcon->elem_type,
            errcon->n_items, errcon->error_level, errcon->error_cutoff, errcon->error_type, 
            errcon->error_value, errcon->flags, errcon->bind_number, errcon->boundTo,
            errcon->sMin, errcon->sMax,                             
            errcon->fp_log, _control->i_step-1, beamline, POST_CORRECTION+
                             ((_control->i_step==0 && errcon->no_errors_first_step)?FORCE_ZERO_ERRORS:0));
        /* set element flags to indicate perturbation of parameters that change the matrix */
        set_element_flags(beamline, errcon->name, errcon->flags, errcon->elem_type, errcon->param_number,
                    errcon->n_items, PARAMETERS_ARE_PERTURBED, VMATRIX_IS_PERTURBED, 0, POST_CORRECTION);
        }        
    if (_control->n_elements_to_vary || errcon->n_items) {
        /* compute matrices for perturbed elements */
        fprintf(stdout, "%ld matrices (re)computed after correction\n", 
               (i=compute_changed_matrices(beamline, run)
                + assert_element_links(links, run, beamline, DYNAMIC_LINK+POST_CORRECTION_LINK))
               + (_control->cell?compute_changed_matrices(_control->cell, run):0) );
        fflush(stdout);
        if (i) {
            beamline->flags &= ~BEAMLINE_CONCAT_CURRENT;
            beamline->flags &= ~BEAMLINE_TWISS_CURRENT;
            beamline->flags &= ~BEAMLINE_RADINT_CURRENT;
            }
        if (i && beamline->matrix) {
            free_matrices(beamline->matrix);
            tfree(beamline->matrix);
            beamline->matrix = NULL;
            }
        log_exit("perturb_beamline");
        return(1);
        }
    log_exit("perturb_beamline");
    return(0);
    }

void set_element_flags(LINE_LIST *beamline, char **elem_name, long *elem_perturb_flags,
        long *type, long *param, long n_elems, 
        long pflag, long mflag, long overwrite, long permit_flags)
{
    long i_elem;
    ELEMENT_LIST *eptr;

    log_entry("set_element_flags");
    if (!elem_name)
        bomb("elem_name array is NULL (set_element_flags)", NULL);
    if (type && !param)
        bomb("param array is NULL when type array is not (set_element_flags)", NULL);
    if (!beamline)
        bomb("beamline pointer is NULL (set_element_flags)", NULL);

    for (i_elem=0; i_elem<n_elems; i_elem++) {
        eptr = NULL;
        if (!elem_name[i_elem]) {
            fprintf(stdout, "error: name missing for element %ld (set_element_flags)\n", i_elem);
            fflush(stdout);
            if (i_elem!=0)
                fprintf(stdout, "preceeding element is %s\n", elem_name[i_elem-1]);
                fflush(stdout);
            if (i_elem!=n_elems-1)
                fprintf(stdout, "suceeding element is %s\n", elem_name[i_elem+1]);
                fflush(stdout);
            }
        while (find_element(elem_name[i_elem], &eptr, &(beamline->elem))) {
            if (elem_perturb_flags && !(elem_perturb_flags[i_elem]&permit_flags))
                continue;
            if (overwrite) 
                eptr->flags = pflag;
            else 
                eptr->flags |= pflag;
            if (type && entity_description[type[i_elem]].flags&HAS_MATRIX &&
                (entity_description[type[i_elem]].parameter[param[i_elem]].flags&PARAM_CHANGES_MATRIX))
                eptr->flags |= mflag;
            }
        }
    log_exit("set_element_flags");
    }
        
void reset_parameter_values(char **elem_name, long *param_number, long *type, long n_elems,
            LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    char *p_elem, *p_elem0;
    long i_elem, elem_type, data_type, param;

    log_entry("reset_parameter_values");
    if (!elem_name)
        bomb("elem_name array is NULL (reset_parameter_values)", NULL);
    if (!param_number)
        bomb("param_number array is NULL (reset_parameter_values)", NULL);
    if (!type)
        bomb("type array is NULL (reset_parameter_values)", NULL);
    if (!beamline)
        bomb("beamline pointer is NULL (reset_parameter_values)", NULL);

    for (i_elem=0; i_elem<n_elems; i_elem++) {
        eptr = NULL;
        elem_type = type[i_elem];
        param     = param_number[i_elem];
        data_type = entity_description[elem_type].parameter[param].type;
        if (!elem_name[i_elem]) {
            fprintf(stdout, "error: name missing for element %ld (reset_parameter_values)\n", i_elem);
            fflush(stdout);
            if (i_elem!=0)
                fprintf(stdout, "preceeding element is %s\n", elem_name[i_elem-1]);
                fflush(stdout);
            if (i_elem!=n_elems-1)
                fprintf(stdout, "suceeding element is %s\n", elem_name[i_elem+1]);
                fflush(stdout);
            }
        while (find_element(elem_name[i_elem], &eptr, &(beamline->elem))) {
            p_elem = eptr->p_elem;
            p_elem0 = eptr->p_elem0;
            switch (data_type) {
                case IS_DOUBLE:
                    *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
		      *((double*)(p_elem0+entity_description[elem_type].parameter[param].offset)) ;
                    break;
                case IS_LONG:
                    *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
		      *((long*)(p_elem0+entity_description[elem_type].parameter[param].offset)) ;
                    break;
                case IS_STRING:
                default:
                    bomb("unknown/invalid variable quantity", NULL);
                    exit(1);
                }
            }
        }
    log_exit("reset_parameter_values");
    }

void assert_parameter_values(char **elem_name, long *param_number, long *type, double *value, long n_elems,
            LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    char *p_elem;
    long i_elem, elem_type, data_type, param;

    log_entry("assert_parameter_values");
    if (!elem_name)
        bomb("elem_name array is NULL (assert_parameter_values)", NULL);
    if (!param_number)
        bomb("param_number array is NULL (assert_parameter_values)", NULL);
    if (!type)
        bomb("type array is NULL (assert_parameter_values)", NULL);
    if (!value)
        bomb("value array is NULL (assert_parameter_values)", NULL);
    if (!beamline)
        bomb("beamline pointer is NULL (assert_parameter_values)", NULL);

    for (i_elem=0; i_elem<n_elems; i_elem++) {
        eptr = NULL;
        elem_type = type[i_elem];
        param     = param_number[i_elem];
        data_type = entity_description[elem_type].parameter[param].type;
        if (!elem_name[i_elem]) {
            fprintf(stdout, "error: name missing for element %ld (assert_parameter_values)\n", i_elem);
            fflush(stdout);
            if (i_elem!=0)
                fprintf(stdout, "preceeding element is %s\n", elem_name[i_elem-1]);
                fflush(stdout);
            if (i_elem!=n_elems-1)
                fprintf(stdout, "suceeding element is %s\n", elem_name[i_elem+1]);
                fflush(stdout);
            }
        while (find_element(elem_name[i_elem], &eptr, &(beamline->elem))) {
            p_elem = eptr->p_elem;
            switch (data_type) {
                case IS_DOUBLE:
                    *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = value[i_elem];
                    break;
                case IS_LONG:
                    *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                      nearestInteger(value[i_elem]);
                    break;
                case IS_STRING:
                default:
                    bomb("unknown/invalid variable quantity", NULL);
                    exit(1);
                }
            }
        }
    log_exit("assert_parameter_values");
    }

long get_parameter_value(double *value, char *elem_name, long param_number, long type, LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    char *p_elem;
    long data_type;

    log_entry("get_parameter_value");
    if (!value)
        bomb("value array is NULL (get_parameter_value)", NULL);
    if (!elem_name)
        bomb("elem_name array is NULL (get_parameter_value)", NULL);
    if (!beamline)
        bomb("beamline pointer is NULL (get_parameter_value)", NULL);

    eptr = NULL;
    data_type = entity_description[type].parameter[param_number].type;
    while (find_element(elem_name, &eptr, &(beamline->elem))) {
        p_elem = eptr->p_elem;
        switch (data_type) {
            case IS_DOUBLE:
                *value = *((double*)(p_elem+entity_description[type].parameter[param_number].offset));
                log_exit("get_parameter_value");
                return(1);
            case IS_LONG:
                *value = *((long*)(p_elem+entity_description[type].parameter[param_number].offset));
                log_exit("get_parameter_value");
                return(1);
            case IS_STRING:
                log_exit("get_parameter_value");
                return(0);
            default:
                bomb("unknown/invalid variable quantity", NULL);
                exit(1);
                break;
            }
        }
    log_exit("get_parameter_value");
    return(0);
    }

void assert_perturbations(char **elem_name, long *param_number, long *type, long n_elems,
        double *amplitude, double *cutoff, long *error_type, double *perturb, long *elem_perturb_flags, 
        long *bind_number, long *bound_to, double *sMin, double *sMax,
        FILE *fp_log, long step, LINE_LIST *beamline, long permit_flags)
{
    ELEMENT_LIST *eptr;
    char *p_elem;
    long i_elem, elem_type, data_type, param, i_group, first_output;
    double delta=0.0;

    log_entry("assert_perturbations");
    if (!elem_name)
        bomb("element name array is NULL (assert_perturbations)", NULL);
    if (!param_number)
        bomb("param_number array is NULL (assert_perturbations)", NULL);
    if (!type)
        bomb("type array is NULL (assert_perturbations)", NULL);
    if (!amplitude)
        bomb("amplitude array is NULL (assert_perturbations)", NULL);
    if (!cutoff)
        bomb("cutoff array is NULL (assert_perturbations)", NULL);
    if (!error_type)
        bomb("error_type array is NULL (assert_perturbations)", NULL);
    if (!perturb)
        bomb("perturb array is NULL (assert_perturbations)", NULL);
    if (!elem_perturb_flags)
        bomb("elem_perturb_flags array is NULL (assert_perturbations)", NULL);
    if (!bind_number)
        bomb("bind_number array is NULL (assert_perturbations)", NULL);
    if (!beamline)
        bomb("beamline structure pointer is NULL (assert_perturbations)", NULL);

    first_output = 1;
    
    for (i_elem=0; i_elem<n_elems; i_elem++) {
        eptr = NULL;
        elem_type = type[i_elem];
        param     = param_number[i_elem];
        data_type = entity_description[elem_type].parameter[param].type;
        i_group = 0;
        delta = DBL_MAX;
        if (!elem_name[i_elem]) {
            fprintf(stdout, "error: name missing for element %ld (assert_perturbations)\n", i_elem);
            fflush(stdout);
            if (i_elem!=0)
                fprintf(stdout, "preceeding element is %s\n", elem_name[i_elem-1]);
                fflush(stdout);
            if (i_elem!=n_elems-1)
                fprintf(stdout, "suceeding element is %s\n", elem_name[i_elem+1]);
                fflush(stdout);
            }
        while (find_element(elem_name[i_elem], &eptr, &(beamline->elem))) {
            p_elem = eptr->p_elem;
            if (!(elem_perturb_flags[i_elem]&permit_flags))
                continue;
            if ((sMin[i_elem]>=0 && eptr->end_pos<sMin[i_elem]) ||
                (sMax[i_elem]>=0 && eptr->end_pos>sMax[i_elem]))
                continue;
            switch (data_type) {
                case IS_DOUBLE:
                    if (permit_flags&FORCE_ZERO_ERRORS) {
                      delta = 0;
                      if (elem_perturb_flags[i_elem]&NONADDITIVE_ERRORS)
                        delta = *((double*)(p_elem+entity_description[elem_type].parameter[param].offset));
                    }
                    else  {
                      if (!(elem_perturb_flags[i_elem]&BIND_ERRORS_MASK) || 
                          (bind_number[i_elem]>=1 && i_group%bind_number[i_elem]==0) ||
                          i_group==0) 
                        delta = perturbation(amplitude[i_elem], cutoff[i_elem], error_type[i_elem]);
                      if (bound_to[i_elem]>=0)
                        delta = perturb[bound_to[i_elem]];
                      if (elem_perturb_flags[i_elem]&FRACTIONAL_ERRORS)
                        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= (1+delta);
                      else {
                        if (elem_perturb_flags[i_elem]&NONADDITIVE_ERRORS)
                          *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = delta;
                        else
                          *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) += delta;
                      }
                    }
                    if (fp_log) {
                        if (first_output) {
                            first_output = 0;
                            fprintf(fp_log, "%ld              ! simulation step\n", step);
                            if (permit_flags&POST_CORRECTION)
                                fprintf(fp_log, "post-correction\n");
                            else
                                fprintf(fp_log, "pre-correction\n");
                            }
                        fprintf(fp_log, "%21.15e %21.15e %10s %10s %ld %10s\n",
                                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)), delta,
                                entity_description[elem_type].parameter[param].name, eptr->name,
                                eptr->occurence, entity_name[elem_type]);
                        }
                    break;
                case IS_LONG:
                    if (permit_flags&FORCE_ZERO_ERRORS) {
                      delta = 0;
                      if (elem_perturb_flags[i_elem]&NONADDITIVE_ERRORS)
                        delta = *((long*)(p_elem+entity_description[elem_type].parameter[param].offset));
                    }
                    else {
                      if (!(elem_perturb_flags[i_elem]&BIND_ERRORS_MASK) || 
                          (bind_number[i_elem]>=1 && i_group%bind_number[i_elem]==0) ||
                          i_group==0) 
                        delta = perturbation(amplitude[i_elem], cutoff[i_elem], error_type[i_elem]);
                      else if (bound_to[i_elem]>=0)
                        delta = perturb[bound_to[i_elem]];
                      if (elem_perturb_flags[i_elem]&FRACTIONAL_ERRORS)
                        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= (1+delta);
                      else {
                        if (elem_perturb_flags[i_elem]&NONADDITIVE_ERRORS)
                          *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                            nearestInteger(delta);
                        else
                          *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) += 
                            nearestInteger(delta);
                      }
                    }
                    break;
                case IS_STRING:
                default:
                    bomb("unknown/invalid variable quantity", NULL);
                    exit(1);
                }
            if (i_group++==0) 
              /* save value if this is the first from a group */
              perturb[i_elem] = delta;
            if (elem_perturb_flags[i_elem]&ANTIBIND_ERRORS)
                delta = -delta;
            }
        }
    if (fp_log) {
        if (!first_output)
            fputc('\n', fp_log);    /* end-of-table indicator */
        fflush(fp_log);
        }
    log_exit("assert_perturbations");
    }

long compute_changed_matrices(LINE_LIST *beamline, RUN *run)
{
  ELEMENT_LIST *eptr;
  long n_changed;
  double Pref_input;
  
  log_entry("compute_changed_matrices");

  n_changed = 0;
  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->pred)
      Pref_input = eptr->pred->Pref_output;
    else
      Pref_input = run->p_central;
    if (eptr->flags&VMATRIX_IS_VARIED || eptr->flags&VMATRIX_IS_PERTURBED || 
        (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix==NULL) ||
        (Pref_input!=eptr->Pref_input && entity_description[eptr->type].flags&MAT_CHW_ENERGY)) {
      if (eptr->matrix) {
        free_matrices(eptr->matrix);
        tfree(eptr->matrix);
        eptr->matrix = NULL;
      }
      compute_matrix(eptr, run, NULL);
      n_changed++;
      eptr->flags &= ~VMATRIX_IS_PERTURBED;
      eptr->flags &= ~VMATRIX_IS_VARIED;
    } else if (Pref_input!=eptr->Pref_input)
        eptr->Pref_input = eptr->Pref_output = Pref_input;
    eptr = eptr->succ;
  }
  log_exit("compute_changed_matrices");
  return(n_changed);
}

void check_VARY_structure(VARY *vary, char *caller)
{
    if (!vary->index) {
        fprintf(stdout, "VARY structure index array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->index_limit) {
        fprintf(stdout, "VARY structure index_limit array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->element_index) {
        fprintf(stdout, "VARY structure element_index array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->element) {
        fprintf(stdout, "VARY structure element array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->item) {
        fprintf(stdout, "VARY structure item array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->initial) {
        fprintf(stdout, "VARY structure initial array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->final) {
        fprintf(stdout, "VARY structure final array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->step) {
        fprintf(stdout, "VARY structure step array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->varied_quan_name) {
        fprintf(stdout, "VARY structure varied_quan_name array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->varied_quan_unit) {
        fprintf(stdout, "VARY structure varied_quan_unit array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->varied_type) {
        fprintf(stdout, "VARY structure varied_type array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->varied_quan_value) {
        fprintf(stdout, "VARY structure varied_quan_value array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->varied_param) {
        fprintf(stdout, "VARY structure varied_param array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    if (!vary->flags) {
        fprintf(stdout, "VARY structure flags array is NULL (%s)", caller);
        fflush(stdout);
        abort();
        }
    }

long load_enumerated_values(double **value, char *file, char *column)
{
    SDDS_TABLE SDDS_table;
    long count=0;

    if (!value)
        bomb("NULL value pointer passed (load_enumerated_values)", NULL);
    if (!file)
        bomb("NULL filename passed (load_enumerated_values)", NULL);
    if (!column)
        bomb("NULL column name passed (load_enumerated_values)", NULL);
    
    if (!SDDS_InitializeInputFromSearchPath(&SDDS_table, file) || 
        SDDS_ReadTable(&SDDS_table)!=1 ||
        !(count=SDDS_CountRowsOfInterest(&SDDS_table)) ||
        !(*value = SDDS_GetColumnInDoubles(&SDDS_table, column)) ||
        !SDDS_Terminate(&SDDS_table))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    return(count);
    }
