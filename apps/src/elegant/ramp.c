/*************************************************************************\
* Copyright (c) 2010 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2010 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "ramp.h"

void addRampElements(RAMP_DATA *rampData, NAMELIST_TEXT *nltext, LINE_LIST *beamline)
{
  long n_items, n_added, firstIndexInGroup;
  ELEMENT_LIST *context;
  double sMin = -DBL_MAX, sMax = DBL_MAX;

  /* process namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&ramp_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (name==NULL) {
    if (!type)
      bombElegant("element name missing in ramp_elements namelist", NULL);
    SDDS_CopyString(&name, "*");
  }
  if (item==NULL)
    bombElegant("item name missing in ramp_elements namelist", NULL);
  if (echoNamelists) print_namelist(stdout, &ramp_elements);
  if (start_pass>=end_pass)
    bombElegant("start_pass >= end_pass", NULL);
  if (exponent<=0)
    bombElegant("exponent <= 0", NULL);
  
  n_added = 0;
  n_items = rampData->nItems;
  context = NULL;
  if (after && strlen(after)) {
    if (!(context=find_element(after, &context, &(beamline->elem)))) {
      fprintf(stdout, "Element %s not found in beamline.\n", after);
      exitElegant(1);
    }
    sMin = context->end_pos;
    if (find_element(after, &context, &(beamline->elem))) {
      fprintf(stdout, "Element %s found in beamline more than once.\n", after);
      exitElegant(1);
    }
    fprintf(stdout, "%s found at s = %le m\n", after, sMin);
    fflush(stdout);
  }
  context = NULL;
  if (before && strlen(before)) {
    if (!(context=find_element(before, &context, &(beamline->elem)))) {
      fprintf(stdout, "Element %s not found in beamline.\n", before);
      exitElegant(1);
    }
    sMax = context->end_pos;
    if (find_element(before, &context, &(beamline->elem))) {
      fprintf(stdout, "Element %s found in beamline more than once.\n", after);
      exitElegant(1);
    }
    fprintf(stdout, "%s found at s = %le m\n", before, sMax);
    fflush(stdout);
  }
  if (after && before && sMin>sMax) {
    fprintf(stdout, "Element %s is not upstream of %s!\n",
            before, after);
    exitElegant(1);
  }
  if (type && has_wildcards(type) && strchr(type, '-'))
    type = expand_ranges(type);
  if (has_wildcards(name)) {
    if (strchr(name, '-'))
      name = expand_ranges(name);
    str_toupper(name);
    firstIndexInGroup = -1;
    while ((context=wfind_element(name, &context, &(beamline->elem)))) {
      if (type && !wild_match(entity_name[context->type], type))
        continue;
      if ((sMin>=0 && context->end_pos<sMin) ||
          (sMax>=0 && context->end_pos>sMax) ||
          (s_start>=0 && context->end_pos<s_start) ||
          (s_end>=0 && context->end_pos>s_end) ||
          (start_occurence && context->occurence<start_occurence) ||
          (end_occurence && context->occurence>end_occurence) )
        continue;

      rampData->element          = SDDS_Realloc(rampData->element, sizeof(*rampData->element)*(n_items+1));
      rampData->parameterNumber  = SDDS_Realloc(rampData->parameterNumber, sizeof(*rampData->parameterNumber)*(n_items+1));
      rampData->flags            = SDDS_Realloc(rampData->flags, sizeof(*rampData->flags)*(n_items+1));
      rampData->unperturbedValue = SDDS_Realloc(rampData->unperturbedValue, sizeof(*rampData->unperturbedValue)*(n_items+1));
      rampData->startPass        = SDDS_Realloc(rampData->startPass, sizeof(*rampData->startPass)*(n_items+1));
      rampData->endPass          = SDDS_Realloc(rampData->endPass, sizeof(*rampData->endPass)*(n_items+1));
      rampData->startValue       = SDDS_Realloc(rampData->startValue, sizeof(*rampData->startValue)*(n_items+1));
      rampData->endValue         = SDDS_Realloc(rampData->endValue, sizeof(*rampData->endValue)*(n_items+1));
      rampData->exponent         = SDDS_Realloc(rampData->exponent, sizeof(*rampData->exponent)*(n_items+1));

      rampData->element[n_items] = context;
      rampData->flags[n_items] = (multiplicative?MULTIPLICATIVE_RAMP:0) + (differential?DIFFERENTIAL_RAMP:0) 
        + (verbose?VERBOSE_RAMP:0);
      rampData->startPass[n_items] = start_pass;
      rampData->endPass[n_items] = end_pass;
      rampData->startValue[n_items] = start_value;
      rampData->endValue[n_items] = end_value;
      rampData->exponent[n_items] = exponent;
      
      if ((rampData->parameterNumber[n_items] = confirm_parameter(item, context->type))<0) {
        fprintf(stdout, "error: cannot ramp %s---no such parameter for %s (wildcard name: %s)\n",item, context->name, name);
        fflush(stdout);
        exitElegant(1);
      }

      rampData->unperturbedValue[n_items] 
        = parameter_value(context->name, context->type, rampData->parameterNumber[n_items], beamline);

      if (rampData->unperturbedValue[n_items]==0 && rampData->flags[n_items]&MULTIPLICATIVE_RAMP) {
        fprintf(stdout, "***\7\7\7 warning: you've specified multiplicative modulation for %s.%s, but the unperturbed value is zero.\nThis may be an error.\n", 
                context->name, item);
        fflush(stdout);
      }
      rampData->nItems = ++n_items;
      n_added++;
      if (firstIndexInGroup==-1)
        firstIndexInGroup = n_items-1;
    }
  }
  else {
    str_toupper(name);
    if (!(context=find_element(name, &context, &(beamline->elem)))) {
      fprintf(stdout, "error: cannot ramp element %s--not in beamline\n", name);
      fflush(stdout);
      exitElegant(1);
    }
    rampData->element          = SDDS_Realloc(rampData->element, sizeof(*rampData->element)*(n_items+1));
    rampData->parameterNumber  = SDDS_Realloc(rampData->parameterNumber, sizeof(*rampData->parameterNumber)*(n_items+1));
    rampData->flags            = SDDS_Realloc(rampData->flags, sizeof(*rampData->flags)*(n_items+1));
    rampData->unperturbedValue = SDDS_Realloc(rampData->unperturbedValue, sizeof(*rampData->unperturbedValue)*(n_items+1));
    rampData->startPass        = SDDS_Realloc(rampData->startPass, sizeof(*rampData->startPass)*(n_items+1));
    rampData->endPass          = SDDS_Realloc(rampData->endPass, sizeof(*rampData->endPass)*(n_items+1));
    rampData->startValue       = SDDS_Realloc(rampData->startValue, sizeof(*rampData->startValue)*(n_items+1));
    rampData->endValue         = SDDS_Realloc(rampData->endValue, sizeof(*rampData->endValue)*(n_items+1));
    rampData->exponent         = SDDS_Realloc(rampData->exponent, sizeof(*rampData->exponent)*(n_items+1));
    
    rampData->flags[n_items] = (multiplicative?MULTIPLICATIVE_RAMP:0) + (differential?DIFFERENTIAL_RAMP:0) 
      + (verbose?VERBOSE_RAMP:0);
    rampData->startPass[n_items] = start_pass;
    rampData->endPass[n_items] = end_pass;
    rampData->startValue[n_items] = start_value;
    rampData->endValue[n_items] = end_value;
    rampData->exponent[n_items] = exponent;

    rampData->element[n_items] = context;
    
    if ((rampData->parameterNumber[n_items] = confirm_parameter(item, context->type))<0) {
      fprintf(stdout, "error: cannot ramp %s--no such parameter for %s (wildcard name: %s)\n",item, context->name, name);
      fflush(stdout);
      exitElegant(1);
    }
    rampData->unperturbedValue[n_items] 
      = parameter_value(context->name, context->type, rampData->parameterNumber[n_items], beamline);
    if (rampData->unperturbedValue[n_items]==0 && rampData->flags[n_items]&MULTIPLICATIVE_RAMP) {
      fprintf(stdout, "***\7\7\7 warning: you've specified multiplicative modulation for %s.%s, but the unperturbed value is zero.\nThis may be an error.\n", 
              context->name,  item);
      fflush(stdout);
    }
    rampData->nItems += 1;
    n_added++;
  }
  

  if (!n_added) {
    fprintf(stdout, "error: no match given modulation\n");
    fflush(stdout);
    exitElegant(1);
  }
}

long applyElementRamps(RAMP_DATA *rampData, double pCentral, RUN *run, long iPass)
{
  long iMod, matricesUpdated;
  double modulation, value;
  long type, param;
  char *p_elem;
  
  matricesUpdated = 0;
  
  for (iMod=0; iMod<rampData->nItems; iMod++) {

    type = rampData->element[iMod]->type;
    param = rampData->parameterNumber[iMod];
    p_elem = (char*)(rampData->element[iMod]->p_elem);
    value = rampData->unperturbedValue[iMod];

    if (iPass<=rampData->startPass[iMod])
      modulation = rampData->startValue[iMod];
    else if (iPass>=rampData->endPass[iMod])
      modulation = rampData->endValue[iMod];
    else
      modulation = rampData->startValue[iMod] 
        + (rampData->endValue[iMod]-rampData->startValue[iMod])*
          pow((1.0*iPass-rampData->startPass[iMod])/(1.0*rampData->endPass[iMod]-rampData->startPass[iMod]), 
              rampData->exponent[iMod]);
    
    if (rampData->flags[iMod]&DIFFERENTIAL_RAMP) {
      if (rampData->flags[iMod]&MULTIPLICATIVE_RAMP)
        value = (1+modulation)*value;
      else
        value = value + modulation;
    } else {
      if (rampData->flags[iMod]&MULTIPLICATIVE_RAMP)
        value = value*modulation;
      else
        value = modulation;
    }

    switch (entity_description[type].parameter[param].type)  {
    case IS_DOUBLE:
      *((double*)(p_elem+entity_description[type].parameter[param].offset)) = value;
      if (rampData->flags[iMod]&VERBOSE_RAMP) 
        fprintf(stdout, "Ramp value for element %s, parameter %s is %le at pass %ld (originally %le)\n",
                rampData->element[iMod]->name,
                entity_description[type].parameter[param].name, value, iPass, rampData->unperturbedValue[iMod]);
      break;
    case IS_LONG:
      *((long*)(p_elem+entity_description[type].parameter[param].offset)) = value + 0.5;
      if (rampData->flags[iMod]&VERBOSE_RAMP) 
        fprintf(stdout, "Ramp value for element %s, parameter %s is %ld at pass %ld (originally %ld)\n",
                rampData->element[iMod]->name,
                entity_description[type].parameter[param].name, (long)(value+0.5), iPass, (long)(rampData->unperturbedValue[iMod]));
      break;
    default:
      break;
    }
    
    if (entity_description[type].flags&HAS_MATRIX && 
        entity_description[type].parameter[param].flags&PARAM_CHANGES_MATRIX) {
      /* update the matrix */
      if (rampData->element[iMod]->matrix) {
        free_matrices(rampData->element[iMod]->matrix);
        tfree(rampData->element[iMod]->matrix);
        compute_matrix(rampData->element[iMod], run, NULL);
        matricesUpdated ++;
      }
    }
  }

  return matricesUpdated;
}


      
