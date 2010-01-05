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
#include "match_string.h"
#include "insert_elements.h"

typedef struct {
  char *name, *type, *exclude, *elemDef;
  long nskip, add_end, total, occur[100];
} ADD_ELEM;

static ADD_ELEM addElem;
static long add_elem_flag = 0;

long getAddElemFlag() 
{
  return (add_elem_flag);
}

char *getElemDefinition()
{
  return (addElem.elemDef);
} 

long getAddEndFlag()
{
  return (addElem.add_end);
} 

void do_insert_elements(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&insert_elements, nltext);
  if (echoNamelists) print_namelist(stdout, &insert_elements);

  if (disable)
    return;
  if (!skip && !add_at_end && !total_occurrences)
    bomb("skip, add_at_end and total_occurrences can not be zero at the same time", NULL);

  if (total_occurrences) {
    if (skip)
      bomb("skip and total_occurrences can not be used together. One must have to be set to zero", NULL);
    if (has_wildcards(name))
      bomb("element name has to be specified if you use the occurrence feature", NULL);
  }
  add_elem_flag = 0;
  
  if ((!name || !strlen(name)) && !type)
    bomb("name or type needs to be given", NULL);
  if (!element_def || !strlen(element_def))
    bomb("element's definition is not given", NULL);
 
  if (name) {
    str_toupper(name);
    if (has_wildcards(name) && strchr(name, '-'))
      name = expand_ranges(name);
  }
  if (type) {
    str_toupper(type);
    if (has_wildcards(type) && strchr(type, '-'))
      type = expand_ranges(type);
    for (i=0; i<N_TYPES; i++)
      if (wild_match(entity_name[i], type))
	break;
    if (i==N_TYPES) {
      fprintf(stderr, "type pattern %s does not match any known type", type);
      exit(1);
    }
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }
  str_to_upper_quotes(element_def);

  add_elem_flag = 1;
  addElem.nskip = skip;
  addElem.add_end =0;
  if (add_at_end)
    addElem.add_end = add_at_end;
  addElem.name = name;
  addElem.type = type;
  addElem.exclude = exclude;
  addElem.elemDef = element_def;
  delete_spaces(addElem.elemDef);

  addElem.total = total_occurrences;
  for (i=0; i< addElem.total; i++) {
    addElem.occur[i] =  occurrence[i];
  }

  beamline = get_beamline(NULL, beamline->name, run->p_central, 0);
  add_elem_flag = 0;

  return;
}

long insertElem(char *name, long type, long *skip, long occurPosition) 
{
  long i;

  if (addElem.exclude && wild_match(name, addElem.exclude))
    return(0);
  if (addElem.name && !wild_match(name, addElem.name))
    return(0);
  if (addElem.type && !wild_match(entity_name[type], addElem.type))
    return(0);

  if (addElem.total) {
    for (i=0; i<addElem.total; i++) {
      if (occurPosition == addElem.occur[i])
        return(1);
    }
    return(0);
  }

  if (addElem.nskip) {
    (*skip)++;
    if (*skip < addElem.nskip)
      return(0);
    *skip = 0;
    return(1);
  }
  return(0);
}

