/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: link_elements.c
 * contents: element_link_control(), add_element_links(), assert_element_links()
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "link_elements.h"

#define DEBUG 0

void element_link_control(ELEMENT_LINKS *links, NAMELIST_TEXT *nltext, RUN *run_cond, LINE_LIST *beamline)
{
    long i, j, flag;

    log_entry("element_link_control");

    /* reset namelist variables to defaults */
/*
    clear_links = summarize_links = verbosity = 0;
 */

    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&link_control, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &link_control);
    
    if (summarize_links) {
        fprintf(stdout, "\nsummary of element links:\n");
        fflush(stdout);
        if (!links->n_links)
            fprintf(stdout, "    no links defined\n");
            fflush(stdout);
        if (!links->target_name || !links->item || !links->source_name || !links->equation ||
                !links->n_targets || !links->target_elem || !links->source_elem)
            bombElegant("link structure has null pointers", NULL);
        for (i=0; i<links->n_links; i++) {
            j = 0;
            flag = links->flags[i]&LINK_FLAG_MASK;
            while ((flag = flag/2))
                j ++;
            if (j>=N_LINK_MODES)
                bombElegant("unknown link mode detected during link summary", NULL);
            fprintf(stdout, "%s.%s linked (%s) to %s within [%e, %e] using equation \"%s\"  --  %ld occurences:\n",
                    links->target_name[i], links->item[i], 
                    link_mode[j],
                    links->source_name[i],
                    links->minimum[i], links->maximum[i], links->equation[i],
                    links->n_targets[i]);
            fflush(stdout);
            for (j=0; j<links->n_targets[i]; j++)
                fprintf(stdout, "   %s#%ld at z=%.15gm linked to %s#%ld at z=%.15gm\n", 
                    links->target_elem[i][j]->name, links->target_elem[i][j]->occurence, links->target_elem[i][j]->end_pos,
                    links->source_elem[i][j]->name, links->source_elem[i][j]->occurence, links->source_elem[i][j]->end_pos);
                fflush(stdout);
            }
        fputc('\n', stdout);
        log_exit("element_link_control");
        return;
        }

    if (clear_links) {
        links->n_links = 0;
        }

    log_exit("element_link_control");
    }

void add_element_links(ELEMENT_LINKS *links, NAMELIST_TEXT *nltext, LINE_LIST *beamline)
{
    long n_links, src_position_code=0, n_targets, n_sources, mode_code=0;
    long targets, iTarget, j;
    char **targetList;
    ELEMENT_LIST *t_context, *s_context, **eptr, *eptr1;
    double dz_min, dz;
#if DEBUG
    long i;
#endif

    log_entry("add_element_links");

    /* set namelist variables to defaults */
    target = item = source = equation = exclude = NULL;
    /* must initialize these hear rather than in the .nl file
     * to avoid problems with str_tolower() and other operations
     */
    cp_str(&source_position, "before");
    cp_str(&mode, "dynamic");

    /* process namelist text */
    if (processNamelist(&link_elements, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (target)          str_toupper(target);
    if (exclude)         str_toupper(exclude);
    if (item)            str_toupper(item);
    if (source)          str_toupper(source);
    if (source_position) 
      str_tolower(source_position);
    else
      cp_str(&source_position, "nearest");
    if (mode) 
      str_tolower(mode);
    else
      cp_str(&mode, "dynamic");
    if (echoNamelists) print_namelist(stdout, &link_elements);

    /* check for valid input */
    if (!target)
        bombElegant("link target not named", NULL);
    if (!item)
        bombElegant("link item not named", NULL);
    if (!source)
        bombElegant("link source not named", NULL);
    if (!equation)
        bombElegant("link equation not given", NULL);
    if (!source_position || (src_position_code=match_string(source_position, src_position_name, N_SRC_POSITIONS, 0))<0)
        bombElegant("source_position not given/unknown", NULL);
    if (!mode || (mode_code=match_string(mode, link_mode, N_LINK_MODES, 0))<0)
        bombElegant("link mode not known", NULL);
    if (minimum>maximum)
      bombElegant("minimum>maximum", NULL);
    
    t_context = s_context = NULL;

    if (has_wildcards(target) && strchr(target, '-'))
      target = expand_ranges(target);
    if (exclude && strlen(exclude) && has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
    
    if (!(t_context=wfind_element(target, &t_context, &(beamline->elem)))) {
      fprintf(stdout, "error: cannot make link with target element %s--not in beamline\n", target);
      fflush(stdout);
      exitElegant(1);
    }
    if (!(s_context=find_element(source, &s_context, &(beamline->elem)))) {
      fprintf(stdout, "error: cannot make link with source element %s--not in beamline\n", source);
      fflush(stdout);
      exitElegant(1);
    }

    targets = 0;
    targetList = NULL;
    /* make a list of all the unique element names that match this (possibly wildcard) target */
    do {
      int32_t duplic;
      if (!exclude || !strlen(exclude) || !wild_match(t_context->name, exclude)) {
        targetList = SDDS_Realloc(targetList, sizeof(*targetList)*(targets+1));
        binaryInsert((void**)targetList, targets, t_context->name, strcmp, &duplic);
        if (!duplic)
          targets++;
      }
    } while ((t_context=wfind_element(target, &t_context, &(beamline->elem))));
    if (!targets)
      bombElegant("cannot make link--no targets found\n", NULL);
      
    /* note that targets==1 if all the targets have the same name ! */
    for (iTarget=0; iTarget<targets; iTarget++) {
      n_links = links->n_links;
      target = targetList[iTarget];
      t_context = NULL;
      t_context = find_element(target, &t_context, &(beamline->elem));

      /* expand the arrays */
      links->target_name     = trealloc(links->target_name, sizeof(*links->target_name)*(n_links+1));
      links->target_elem     = trealloc(links->target_elem, sizeof(*links->target_elem)*(n_links+1));
      links->item            = trealloc(links->item, sizeof(*links->item)*(n_links+1));
      links->target_param    = trealloc(links->target_param, sizeof(*links->target_param)*(n_links+1));
      links->source_name     = trealloc(links->source_name, sizeof(*links->source_name)*(n_links+1));
      links->source_position = trealloc(links->source_position, sizeof(*links->source_position)*(n_links+1));
      links->flags           = trealloc(links->flags, sizeof(*links->flags)*(n_links+1));
      links->source_elem     = trealloc(links->source_elem, sizeof(*links->source_elem)*(n_links+1));
      links->equation        = trealloc(links->equation, sizeof(*links->equation)*(n_links+1));
      links->n_targets       = trealloc(links->n_targets, sizeof(*links->n_targets)*(n_links+1));
      links->initial_value   = trealloc(links->initial_value, sizeof(*links->initial_value)*(n_links+1));
      links->baseline_value  = trealloc(links->baseline_value, sizeof(*links->baseline_value)*(n_links+1));
      links->minimum   = trealloc(links->minimum, sizeof(*links->minimum)*(n_links+1));
      links->maximum   = trealloc(links->maximum, sizeof(*links->maximum)*(n_links+1));

      /* copy the basic data */
      cp_str(links->target_name+n_links, target);
      cp_str(links->item+n_links, item);
      cp_str(links->source_name+n_links, source);
      cp_str(links->equation+n_links, equation);
      links->source_position[n_links] = src_position_code;
      links->flags[n_links] = link_mode_flag[mode_code];
      links->minimum[n_links] = minimum;
      links->maximum[n_links] = maximum;
      
      /* make the list of pointers to targets */
      eptr = tmalloc(sizeof(*eptr));
      eptr[0] = t_context;
      if ((links->target_param[n_links] = confirm_parameter(item, t_context->type))<0) {
        fprintf(stdout, "error: element %s does not have a parameter %s\n", target, item);
        fflush(stdout);
        exitElegant(1);
      }
      n_targets = 1;
      while ((t_context=find_element(target, &t_context, &(beamline->elem)))) {
        eptr = trealloc(eptr, sizeof(*eptr)*(n_targets+1));
        eptr[n_targets] = t_context;
        n_targets++;
      }
      links->baseline_value[n_links] = tmalloc(sizeof(*links->baseline_value[n_links])*n_targets);
      links->n_targets[n_links] = n_targets;
      links->target_elem[n_links] = eptr;
      t_context = links->target_elem[n_links][0];
      switch (entity_description[eptr[0]->type].parameter[links->target_param[n_links]].type) {
      case IS_DOUBLE:
        links->initial_value[n_links] = 
          *((double*)(eptr[0]->p_elem+entity_description[eptr[0]->type].parameter[links->target_param[n_links]].offset));
        break;
      case IS_LONG:
        links->initial_value[n_links] = 
          *((long*)(eptr[0]->p_elem+entity_description[eptr[0]->type].parameter[links->target_param[n_links]].offset));
        break;
      default:
        bombElegant("invalid type of item for target of link", NULL);
        break;
      }
      for (j=0; j<n_targets; j++)
        links->baseline_value[n_links][j] = links->initial_value[n_links];

      /* make the list of pointers to sources */
      if (iTarget) {
        s_context = NULL;
        if (!(s_context=find_element(source, &s_context, &(beamline->elem)))) {
          fprintf(stdout, "error: cannot make link with source element %s--not in beamline\n", source);
          fflush(stdout);
          exitElegant(1);
        }
      }
      eptr = tmalloc(sizeof(*eptr)*(n_targets));
      if (src_position_code==SRC_POSITION_SAME_OCCURENCE) {
        n_sources = 0;
        while (n_sources<n_targets) {
          eptr1 = NULL;
          s_context = NULL;
          while (find_element(source, &s_context, &(beamline->elem))) {
            if (s_context->occurence==links->target_elem[n_links][n_sources]->occurence) {
              eptr1 = s_context;
              break;
            }
          }
          if (!eptr1) {
            fprintf(stdout, "error: no %s element is found with the same occurence number as the %ld-th %s element--can't link as requested\n",
                    source, n_sources, target);
            fflush(stdout);
            exitElegant(1);
          }
          eptr[n_sources++] = eptr1;
        }
      }
      else if (src_position_code==SRC_POSITION_NEAREST) {
        n_sources = 0;
        while (n_sources<n_targets) {
          dz_min = DBL_MAX;
          eptr1 = NULL;
          s_context = NULL;
          while (find_element(source, &s_context, &(beamline->elem))) {
            if ((dz = fabs(s_context->end_pos-links->target_elem[n_links][n_sources]->end_pos))<dz_min) {
              eptr1 = s_context;
              dz_min = dz;
            }
          }
          if (!eptr1) {
            fprintf(stdout, "error: no %s element is found near the %ld-th %s element--can't link as requested\n",
                    source, n_sources, target);
            fflush(stdout);
            exitElegant(1);
          }
          eptr[n_sources++] = eptr1;
        }
      }
      else if (src_position_code==SRC_POSITION_ADJACENT) {
        n_sources = 0;
        while (n_sources<n_targets) {
          eptr1 = NULL;
          if ((eptr1=links->target_elem[n_links][n_sources]->pred)) {
            if (strcmp(eptr1->name, source)!=0)
              eptr1 = NULL;
          }
          if (!eptr1 && (eptr1=links->target_elem[n_links][n_sources]->succ)) {
            if (strcmp(eptr1->name, source)!=0)
              eptr1 = NULL;
          }
          if (!eptr1) {
            fprintf(stdout, "error: no %s element is found adjacent to the %ld-th %s element--can't link as requested\n",
                    source, n_sources, target);
            fflush(stdout);
            exitElegant(1);
          }
          eptr[n_sources++] = eptr1;
        }
      }
      else if (src_position_code==SRC_POSITION_BEFORE) {
        if (links->target_elem[n_links][0]->end_pos<s_context->end_pos) {
          fprintf(stdout, "error: there is no %s element before the first %s element--can't link as requested\n",
                  source, target);
          fflush(stdout);
          exitElegant(1);
        }
        eptr[0] = s_context;
        n_sources = 0;
        while (n_sources<n_targets) {
          eptr1 = NULL;
          do {
            if (s_context->end_pos<links->target_elem[n_links][n_sources]->end_pos)
              eptr1 = s_context;
            else if (s_context->end_pos==links->target_elem[n_links][n_sources]->end_pos) {
              eptr1 = s_context;
              break;
            }
            else
              break;
          } while (find_element(source, &s_context, &(beamline->elem)));
          if (!eptr1) {
            fprintf(stdout, "error: no %s element is found before the %ld-th %s element--can't link as requested\n",
                    source, n_sources, target);
            fflush(stdout);
            exitElegant(1);
          }
          eptr[n_sources++] = eptr1;
          s_context = eptr[n_sources-1];
        }
      }
      else if (src_position_code==SRC_POSITION_AFTER) {
        if (links->target_elem[n_links][0]->end_pos>=s_context->end_pos) {
          /* search for first source element after first target element */
          while (find_element(source, &s_context, &(beamline->elem))) {
            if (links->target_elem[n_links][0]->end_pos<s_context->end_pos)
              break;
          }
          if (!s_context) {
            fprintf(stdout, "error: no %s element after the first %s element--can't link as requested\n",
                    source, target);
            fflush(stdout);
            exitElegant(1);
          }
        }
        eptr[0] = s_context;
        n_sources = 1;
        while (n_sources<n_targets) {
          s_context = links->target_elem[n_links][n_sources-1];
          while (find_element(source, &s_context, &(beamline->elem))) {
            if (s_context->end_pos>links->target_elem[n_links][n_sources]->end_pos)
              break;
          }
          if (!s_context) {
            fprintf(stdout, "error: no %s element is found after the %ld-th %s element--can't link as requested\n",
                    source, n_sources, target);
            fflush(stdout);
            exitElegant(1);
          }
          eptr[n_sources++] = s_context;
        }
      }

      links->source_elem[n_links] = eptr;

#if DEBUG
      fprintf(stdout, "list of targets and sources:\n");
      fflush(stdout);
      for (i=0; i<n_targets; i++)
        fprintf(stdout, "%s at z=%em linked to %s at z=%em\n", 
                links->target_elem[n_links][i]->name, links->target_elem[n_links][i]->end_pos,
                links->source_elem[n_links][i]->name, links->source_elem[n_links][i]->end_pos);
      fflush(stdout);
#endif

      links->n_links += 1;
    }
    
    log_exit("add_element_links");
    }

long assert_element_links(ELEMENT_LINKS *links, RUN *run, LINE_LIST *beamline, long flags)
{
    long i_link, i_elem, i_item, matrices_changed;
    long elem_type, data_type, param;
    double value;
    ELEMENT_LIST **targ, **sour;
    char *p_elem;
    short lengthChanged = 0;
    
    log_entry("assert_element_links");
    if (!links || links->n_links==0) {
        log_exit("assert_element_links");
        return(0);
        }

    if (!links->target_name || !links->item || !links->source_name || !links->equation ||
            !links->n_targets || !links->target_elem || !links->source_elem) {
        fputs("error: link structure has null pointers (assert_element_links)", stdout);
        abort();
        }

    matrices_changed = 0;
    for (i_link=0; i_link<links->n_links; i_link++) {
        if (!(flags&links->flags[i_link]))
            continue;
        targ = links->target_elem[i_link];
        sour = links->source_elem[i_link];
        elem_type = targ[0]->type;
        param     = links->target_param[i_link];
        data_type = entity_description[elem_type].parameter[param].type;
#if DEBUG
        fprintf(stdout, "asserting %ld links of %s.%s to %s\n", links->n_targets[i_link],
            links->target_name[i_link], links->item[i_link], links->source_name[i_link]);
        fflush(stdout);
        fprintf(stdout, "source type is %ld, with %ld parameters\n", sour[0]->type, 
                    entity_description[sour[0]->type].n_params);
        fflush(stdout);
#endif
        for (i_elem=0; i_elem<links->n_targets[i_link]; i_elem++) {
#if DEBUG
            fprintf(stdout, "  working on element %ld\n", i_elem);
            fflush(stdout);
#endif
            p_elem = sour[i_elem]->p_elem;
            for (i_item=0; i_item<entity_description[sour[0]->type].n_params; i_item++) {
  	        char s[1024];
		double value0;
                switch (entity_description[sour[0]->type].parameter[i_item].type) {
                    case IS_DOUBLE:
                        value = *((double*)(p_elem+entity_description[sour[0]->type].parameter[i_item].offset));
                        rpn_store(value, NULL, rpn_create_mem(entity_description[sour[0]->type].parameter[i_item].name, 0));
                        value0 = *((double*)(sour[i_elem]->p_elem0+entity_description[sour[0]->type].parameter[i_item].offset));
			sprintf(s, "%s0", entity_description[sour[0]->type].parameter[i_item].name);
                        rpn_store(value0, NULL, rpn_create_mem(s, 0));
                        break;
                    case IS_LONG:
                        value = *((long*)(p_elem+entity_description[sour[0]->type].parameter[i_item].offset));
                        rpn_store(value, NULL, rpn_create_mem(entity_description[sour[0]->type].parameter[i_item].name, 0));
                        value0 = *((long*)(sour[i_elem]->p_elem0+entity_description[sour[0]->type].parameter[i_item].offset));
			sprintf(s, "%s0", entity_description[sour[0]->type].parameter[i_item].name);
                        rpn_store(value0, NULL, rpn_create_mem(s, 0));
                        break;
                    default:
                        break;
                    }
#if DEBUG 
                fprintf(stdout, "    asserting value %e for %s\n", value, entity_description[sour[0]->type].parameter[i_item].name);
                fprintf(stdout, "    asserting value %e for %s0\n", value0, entity_description[sour[0]->type].parameter[i_item].name);
                fflush(stdout);
#endif
                }
            p_elem = targ[i_elem]->p_elem;

            rpn_clear();
            /* push original value onto stack */
            if (verbosity>1)
                fprintf(stdout, "prior value is %.15g for %s#%ld.%s at z=%.15gm\n",
                        links->baseline_value[i_link][i_elem], 
                        links->target_name[i_link], targ[i_elem]->occurence, links->item[i_link], 
                        targ[i_elem]->end_pos);
            push_num(links->baseline_value[i_link][i_elem]);
            value = rpn(links->equation[i_link]);
            if (rpn_check_error()) exitElegant(1);
            rpn_clear();
            if (value>links->maximum[i_link])
              value = links->maximum[i_link];
            if (value<links->minimum[i_link])
              value = links->minimum[i_link];
            if (verbosity>0)
                fprintf(stdout, "asserting value %.15g for %s#%ld.%s at z=%.15gm\n",
                    value, links->target_name[i_link], targ[i_elem]->occurence, links->item[i_link], targ[i_elem]->end_pos);
                fflush(stdout);
            if (entity_description[elem_type].flags&HAS_LENGTH && 
                entity_description[elem_type].parameter[param].offset==0)
              lengthChanged = 1;
            switch (data_type) {
                case IS_DOUBLE:
                    *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = value;
                    break;
                case IS_LONG:
                    *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                      nearestInteger(value);
                    break;
                case IS_STRING:
                default:
                    bombElegant("unknown/invalid variable quantity (assert_element_links)", NULL);
                    exitElegant(1);
                }
            if (flags&LINK_ELEMENT_DEFINITION) 
                change_defined_parameter_values(&targ[i_elem]->name, &param, &targ[i_elem]->type, &value, 1);
            if ((entity_description[targ[0]->type].parameter[param].flags&PARAM_CHANGES_MATRIX) && targ[i_elem]->matrix) {
                free_matrices(targ[i_elem]->matrix);
                tfree(targ[i_elem]->matrix);
                compute_matrix(targ[i_elem], run, NULL);
                matrices_changed++;
                }
            }
        }
#if DEBUG
    print_line(stdout, beamline);
#endif
    if (lengthChanged)
      compute_end_positions(beamline);
    log_exit("assert_element_links");
    return(matrices_changed);
    }

void reset_element_links(ELEMENT_LINKS *links, RUN *run, LINE_LIST *beamline)
{
    long i_link, i_elem;
    long elem_type, data_type, param;
    ELEMENT_LIST **targ;
    char *p_elem;

    log_entry("reset_element_links");
    if (!links || links->n_links==0) {
        log_exit("reset_element_links");
        return;
        }

    if (!links->target_name || !links->item || !links->equation ||
            !links->n_targets || !links->target_elem) {
        fputs("error: link structure has null pointers (reset_element_links)", stdout);
        abort();
        }
    
    for (i_link=0; i_link<links->n_links; i_link++) {
        targ = links->target_elem[i_link];
        elem_type = targ[0]->type;
        param     = links->target_param[i_link];
        data_type = entity_description[elem_type].parameter[param].type;
        for (i_elem=0; i_elem<links->n_targets[i_link]; i_elem++) {
            p_elem = targ[i_elem]->p_elem;
            switch (data_type) {
                case IS_DOUBLE:
                    *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                      links->baseline_value[i_link][i_elem] = links->initial_value[i_link];
                    break;
                case IS_LONG:
                    *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
                      links->baseline_value[i_link][i_elem] = nearestInteger(links->initial_value[i_link]);
                    break;
                case IS_STRING:
                default:
                    bombElegant("unknown/invalid variable quantity (reset_element_links)", NULL);
                    exitElegant(1);
                }
            if ((entity_description[targ[0]->type].parameter[param].flags&PARAM_CHANGES_MATRIX) && targ[i_elem]->matrix) {
                free_matrices(targ[i_elem]->matrix);
                tfree(targ[i_elem]->matrix);
                targ[i_elem]->matrix = NULL;
                }
            }
        }
    log_exit("reset_element_links");
    }


void rebaseline_element_links(ELEMENT_LINKS *links, RUN *run, LINE_LIST *beamline)
{
    long i_link, i_elem;
    long elem_type, data_type, param;
    ELEMENT_LIST **targ;
    char *p_elem;

    if (!links || links->n_links==0) {
      return;
    }

    if (!links->target_name || !links->item || !links->equation ||
            !links->n_targets || !links->target_elem) {
        fputs("error: link structure has null pointers (rebaseline_element_links)", stdout);
        abort();
        }
    
    for (i_link=0; i_link<links->n_links; i_link++) {
        targ = links->target_elem[i_link];
        elem_type = targ[0]->type;
        param     = links->target_param[i_link];
        data_type = entity_description[elem_type].parameter[param].type;
        for (i_elem=0; i_elem<links->n_targets[i_link]; i_elem++) {
            p_elem = targ[i_elem]->p_elem;
            switch (data_type) {
                case IS_DOUBLE:
                    links->baseline_value[i_link][i_elem] = *((double*)(p_elem+entity_description[elem_type].parameter[param].offset));
                    break;
                case IS_LONG:
                    links->baseline_value[i_link][i_elem] = *((long*)(p_elem+entity_description[elem_type].parameter[param].offset));
                    break;
                case IS_STRING:
                default:
                    bombElegant("unknown/invalid variable quantity (reset_element_links)", NULL);
                    exitElegant(1);
                }
            }
        }
    }


