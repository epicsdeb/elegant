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
#include "replace_elements.h"

typedef struct {
  char *name, *type, *exclude, *elemDef;
  long nskip, total, occur[100];
} DEL_ELEM;

static DEL_ELEM delElem;
/* flag=0: do nothing; 
   flag=-1: delete element; 
   flag=1: replace element with new element. */
static long  delete_elem_flag = 0;

long getDelElemFlag() 
{
  return (delete_elem_flag);
}

char *getElemDefinition1()
{
  return (delElem.elemDef);
} 

void do_replace_elements(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&replace_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &replace_elements);

  if (disable)
    return;
  if (!skip && !total_occurrences)
    bombElegant("skip and total_occurrences can not be zero at the same time", NULL);

  if (total_occurrences) {
    if (skip)
      bombElegant("skip and total_occurrences can not be used together. One must have to be set to zero", NULL);
    if (has_wildcards(name))
      bombElegant("element name has to be specified if you use the occurrence feature", NULL);
  }
  delete_elem_flag = 0;
  
  if ((!name || !strlen(name)) && !type)
    bombElegant("name or type needs to be given", NULL);
 
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
      exitElegant(1);
    }
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }

  if (!element_def || !strlen(element_def)) {
    delete_elem_flag = -1;
  } else {
    delete_elem_flag = 1;
    str_toupper(element_def);
    delElem.elemDef = element_def;
    delete_spaces(delElem.elemDef);
  }

  delElem.nskip = skip;
  delElem.name = name;
  delElem.type = type;
  delElem.exclude = exclude;

  delElem.total = total_occurrences;
  for (i=0; i< delElem.total; i++) {
    delElem.occur[i] =  occurrence[i];
  }

  beamline = get_beamline(NULL, beamline->name, run->p_central, 0);
  delete_elem_flag = 0;

  return;
}

long replaceElem(char *name, long type, long *skip, long occurPosition) 
{
  long i;

  if (delElem.exclude && wild_match(name, delElem.exclude))
    return(0);
  if (delElem.name && !wild_match(name, delElem.name))
    return(0);
  if (delElem.type && !wild_match(entity_name[type], delElem.type))
    return(0);

  if (delElem.total) {
    for (i=0; i<delElem.total; i++) {
      if (occurPosition == delElem.occur[i])
        return(delete_elem_flag);
    }
    return(0);
  }

  if (delElem.nskip) {
    (*skip)++;
    if (*skip < delElem.nskip)
      return(0);
    *skip = 0;
    return(delete_elem_flag);
  }
  return(0);
}

