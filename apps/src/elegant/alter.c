/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "alter.h"
#include "match_string.h"

#define DEBUG 0

void do_alter_element(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long thisType, lastType, iParam=0, nMatches;
    ELEMENT_LIST *context, *eptr;
    char *p_elem;
    char *p_elem0;
    char **changedDefinedParameter = NULL;
    long nChangedDefinedParameter = 0;

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&alter_elements, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &alter_elements);

    if (!name || !strlen(name))
      bombElegant("no name given", NULL);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
    if (!item || !strlen(item))
      bombElegant("no item given", NULL);
    if (multiplicative) {
      if (!differential)
	/* convert to fractional change */
	value = value-1;
      differential = 0;
    }
    if ((((s_start>=0 && s_end>=0) ? 1 : 0) +
         ((start_occurence!=0 && end_occurence!=0) ? 1 : 0 ) +
         ((after || before) ? 1 : 0 ))>1)
      bombElegant("can't combine start_occurence/end_occurence, s_start/s_end, and after/before---use one method only", NULL);
    if (start_occurence>end_occurence) 
      bombElegant("start_occurence > end_occurence", NULL);
    if (after || before) {
      ELEMENT_LIST *context;
      context = NULL;
      s_start = -DBL_MAX;
      s_end = DBL_MAX;
      if (after && strlen(after)) {
        if (!(context=find_element(after, &context, &(beamline->elem)))) {
          fprintf(stdout, "Element %s not found in beamline.\n", after);
          exitElegant(1);
        }
        s_start = context->end_pos;
        if (find_element(after, &context, &(beamline->elem))) {
          fprintf(stdout, "Element %s found in beamline more than once.\n", after);
          exitElegant(1);
        }
        fprintf(stdout, "%s found at s = %le m\n", after, s_start);
        fflush(stdout);
      }
      context = NULL;
      if (before && strlen(before)) {
        if (!(context=find_element(before, &context, &(beamline->elem)))) {
          fprintf(stdout, "Element %s not found in beamline.\n", before);
          exitElegant(1);
        }
        s_end = context->end_pos;
        if (find_element(before, &context, &(beamline->elem))) {
          fprintf(stdout, "Element %s found in beamline more than once.\n", before);
          exitElegant(1);
        }
        fprintf(stdout, "%s found at s = %le m\n", before, s_end);
        fflush(stdout);
      }
      if (s_start>s_end) 
        bombElegant("'after' element follows 'before' element!", NULL);
    }
    if (s_start>s_end)
      bombElegant("s_start > s_end", NULL);
    if (type) {
      long i;
      str_toupper(type);
      if (has_wildcards(type) && strchr(type, '-'))
        type = expand_ranges(type);
      for (i=0; i<N_TYPES; i++)
	if (wild_match(entity_name[i], type))
	  break;
      if (i==N_TYPES)
	bombElegant("type pattern does not match any known type", NULL);
    }    
    if (exclude && has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
      
    context = NULL;
    lastType = -1;
    nMatches = 0;
    while ((eptr=wfind_element(name, &context, &(beamline->elem)))) {
      if (exclude && strlen(exclude) && wild_match(eptr->name, exclude))
        continue;
      if (start_occurence!=0 && end_occurence!=0 && 
          (eptr->occurence<start_occurence || eptr->occurence>end_occurence))
        continue;
      if (s_start>=0 && s_end>=0 &&
          (eptr->end_pos<s_start || eptr->end_pos>s_end))
        continue;
      if (type && !wild_match(entity_name[context->type], type))
        continue;
      if ((thisType = eptr->type)!=lastType) {
        lastType = thisType;
        iParam = confirm_parameter(item, thisType);
      }
      if (iParam<0) {
	fprintf(stderr, "%s: element %s does not have parameter %s\n", 
		allow_missing_parameters?"Warning":"Error",
		eptr->name, item);
	if (!allow_missing_parameters)
	  exitElegant(1);
	continue;
      }
      nMatches++;
      p_elem = eptr->p_elem;
      p_elem0 = eptr->p_elem0;
      switch (entity_description[eptr->type].parameter[iParam].type) {
      case IS_DOUBLE:
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %21.15e to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bombElegant("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (differential)
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += value;
        else if (multiplicative)
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+value;
        else
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) = value;
        *((double*)(p_elem0+entity_description[eptr->type].parameter[iParam].offset)) = 
          *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) ;
        if (verbose) {
          fprintf(stdout, "%21.15e\n",
                  *((double*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
            fflush(stdout);
        }
        break;
      case IS_LONG:
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %ld to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bombElegant("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (differential)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) += 
            nearestInteger(value);
        else if (multiplicative)
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) *= 1+value;
        else
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) =
            nearestInteger(value);
        *((long*)(p_elem0+entity_description[eptr->type].parameter[iParam].offset)) = 
          *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)) ;
        if (verbose) {
          fprintf(stdout, "%ld\n",
                  *((long*)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      case IS_STRING:
        if (string_value==NULL) {
          fprintf(stderr, "Error: string_value is NULL for alter_elements, but parameter %s of %s is a string parameter\n",
                  entity_description[eptr->type].parameter[iParam].name, eptr->name);
          exitElegant(1);
        }
        /* unfortunately, can't free the existing pointer as I can't be sure that it isn't
         * pointing to static memory
         */
        if (verbose)
          fprintf(stdout, "Changing %s.%s from %s to ",
                  eptr->name, 
                  entity_description[eptr->type].parameter[iParam].name, 
                  *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
        cp_str((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset),
               string_value);
        cp_str((char**)(p_elem0+entity_description[eptr->type].parameter[iParam].offset),
               string_value);
        /* this step could be very inefficient */
        if ((nMatches==1 || has_wildcards(name)) &&
	    (!nChangedDefinedParameter ||
	    match_string(eptr->name, changedDefinedParameter, 
			 nChangedDefinedParameter, EXACT_MATCH)==-1)) {
          change_defined_parameter(eptr->name, iParam, thisType, value, NULL, 
                                   differential?LOAD_FLAG_DIFFERENTIAL:
                                   (multiplicative?LOAD_FLAG_FRACTIONAL:LOAD_FLAG_ABSOLUTE));
	  if (!(changedDefinedParameter=SDDS_Realloc(changedDefinedParameter,
						    sizeof(*changedDefinedParameter)*
						     (nChangedDefinedParameter+1)))) {
	    bombElegant("memory allocation failure (alter_elements)", NULL);
	  }
	  changedDefinedParameter[nChangedDefinedParameter++] = eptr->name;
	}
        if (verbose) {
          fprintf(stdout, "%s\n",
                  *((char**)(p_elem+entity_description[eptr->type].parameter[iParam].offset)));
          fflush(stdout);
        }
        break;
      }
      eptr->flags |= 
        PARAMETERS_ARE_PERTURBED |
          ((entity_description[eptr->type].parameter[iParam].flags&PARAM_CHANGES_MATRIX)?VMATRIX_IS_PERTURBED:0);
    }
    if (nMatches==0)  {
      if (allow_missing_elements)
        fprintf(stdout, "Warning: no matches for %s\n", name);
      else {
        fprintf(stdout, "Error: no matches for %s\n", name);
        exitElegant(1);
      }
    } else 
      compute_end_positions(beamline);
    if (changedDefinedParameter)
      free(changedDefinedParameter);
  }


