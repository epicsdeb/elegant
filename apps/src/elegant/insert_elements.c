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
  char **name, *type, *exclude, *elemDef;
  long nNames, nskip, add_end, total, occur[100];
  double sStart, sEnd;
} ADD_ELEM;

static ADD_ELEM addElem;
static long add_elem_flag = 0;

#define COPYLEN 1024
static char elementDefCopy[COPYLEN+1];
static long insertCount = 0;

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
  long i, nNames;
  char **nameList;
  
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&insert_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &insert_elements);

  if (disable)
    return;
  if (!skip && !add_at_end && !total_occurrences)
    bombElegant("skip, add_at_end and total_occurrences can not be zero at the same time", NULL);

  if (total_occurrences) {
    if (skip)
      bombElegant("skip and total_occurrences can not be used together. One must have to be set to zero", NULL);
    if (has_wildcards(name))
      bombElegant("element name has to be specified if you use the occurrence feature", NULL);
  }
  add_elem_flag = 0;
  
  if ((!name || !strlen(name)) && !type)
    bombElegant("name or type needs to be given", NULL);
  if (!element_def || !strlen(element_def))
    bombElegant("element's definition is not given", NULL);

  nameList = NULL;
  nNames = 0;
  if (name) {
    char *ptr;
    str_toupper(name);
    while (ptr=get_token_t(name, ", ")) {
      nameList = SDDS_Realloc(nameList, sizeof(*nameList)*(nNames+1));
      nameList[nNames] = ptr;
      if (has_wildcards(ptr) && strchr(ptr, '-')) {
        nameList[nNames] = expand_ranges(ptr);
        free(ptr);
      } else
        nameList[nNames] = ptr;
      nNames++;
    }
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
  str_to_upper_quotes(element_def);

  add_elem_flag = 1;
  addElem.nskip = skip;
  addElem.add_end =0;
  if (add_at_end)
    addElem.add_end = add_at_end;
  addElem.name = nameList;
  addElem.nNames = nNames;
  addElem.type = type;
  addElem.exclude = exclude;
  addElem.elemDef = element_def;
  addElem.sStart = s_start;
  addElem.sEnd = s_end;
  delete_spaces(addElem.elemDef);
  strncpy(elementDefCopy, element_def, COPYLEN);

  addElem.total = total_occurrences;
  for (i=0; i< addElem.total; i++) {
    addElem.occur[i] =  occurrence[i];
  }

  insertCount = 0;
  beamline = get_beamline(NULL, beamline->name, run->p_central, 0);
  add_elem_flag = 0;
  if (verbose)
    printf("%ld elements inserted in total\n", insertCount);

  return;
}

long insertElem(char *name, long type, long *skip, long occurPosition, double endPosition) 
{
  long i;

  if (addElem.exclude && wild_match(name, addElem.exclude))
    return(0);
  if (addElem.type && !wild_match(entity_name[type], addElem.type))
    return(0);
  if (addElem.name) {
    for (i=0; i<addElem.nNames; i++) {
      if (addElem.name[i] && wild_match(name, addElem.name[i]))
        break;
    }
    if (i==addElem.nNames)
      return(0);
  }

  if ((addElem.sStart>=0 && endPosition<addElem.sStart) ||
      (addElem.sEnd>=0 && endPosition>addElem.sEnd))
    return 0;

  if (addElem.total) {
    for (i=0; i<addElem.total; i++) {
      if (occurPosition == addElem.occur[i]) {
	if (verbose)
	  printf("Adding \"%s\" after occurrence %ld of %s\n",
		 elementDefCopy, occurPosition, name);
	insertCount ++;
	return(1);
      }
    }
    return(0);
  }

  if (addElem.nskip) {
    (*skip)++;
    if (*skip < addElem.nskip)
      return(0);
    *skip = 0;
    if (verbose)
      printf("Adding \"%s\" after occurrence %ld of %s\n",
	     elementDefCopy, occurPosition, name);
    insertCount ++;
    return(1);
  }
  return(0);
}

