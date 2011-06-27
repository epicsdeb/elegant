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

#define DEBUG 0

typedef struct {
  char *name, *type, *exclude;
  long newType;
} TRANSMUTE_SPEC;

static TRANSMUTE_SPEC *transmutationSpec = NULL;
static long transmutationSpecs = 0;

void addTransmutationSpec(char *name, char *type, char *exclude,
                          long newType)
{
  if (!(transmutationSpec 
	= SDDS_Realloc(transmutationSpec,
		       sizeof(*transmutationSpec)*(transmutationSpecs+1))))
    bombElegant("memory allocation failure", NULL);
  transmutationSpec[transmutationSpecs].name = NULL;
  transmutationSpec[transmutationSpecs].type = NULL;
  transmutationSpec[transmutationSpecs].exclude = NULL;
  transmutationSpec[transmutationSpecs].newType = newType;
  if ((name &&
       !SDDS_CopyString(&transmutationSpec[transmutationSpecs].name, name)) ||
      (type &&
       !SDDS_CopyString(&transmutationSpec[transmutationSpecs].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&transmutationSpec[transmutationSpecs].exclude, exclude)))
    bombElegant("memory allocation failure", NULL);
  
  transmutationSpecs++;
}

void clearTransmutationSpecs() 
{
  while (transmutationSpecs--) {
    if (transmutationSpec[transmutationSpecs].name)
      free(transmutationSpec[transmutationSpecs].name);
    if (transmutationSpec[transmutationSpecs].type)
      free(transmutationSpec[transmutationSpecs].type);
    if (transmutationSpec[transmutationSpecs].exclude)
      free(transmutationSpec[transmutationSpecs].exclude);
  }
  free(transmutationSpec);
  transmutationSpec = NULL;
}

long elementTransmutation(char *name, long type) 
{
  long i;
  for (i=0; i<transmutationSpecs; i++) {
    if (transmutationSpec[i].exclude && wild_match(name, transmutationSpec[i].exclude))
      continue;
    if (transmutationSpec[i].name && !wild_match(name, transmutationSpec[i].name))
      continue;
    if (transmutationSpec[i].type && !wild_match(entity_name[type], transmutationSpec[i].type))
      continue;
    return transmutationSpec[i].newType;
  }
  return -1;
}

#include "transmuteElements.h"

void setupTransmuteElements(NAMELIST_TEXT *nltext, RUN *run, 
			 LINE_LIST *beamline)
{
  long i, j, newType;
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&transmute_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &transmute_elements);

  if (clear_all) {
    clearTransmutationSpecs();
    if (!name && !type)
      return;
  }
  if (disable)
    return;

  if (!new_type)
    bombElegant("new_type must be given", NULL);
  str_toupper(new_type);
  j = -1;
  for (i=0; i<N_TYPES; i++) {
    if (strncmp(entity_name[i], new_type, strlen(new_type))==0) {
      if (j>=0) 
        bombElegant("new_type matches more than one element type", NULL);
      j = i;
    }
  }
  if (j==-1)
    bombElegant("new_type does not match a known element type", NULL);
  newType = j;
  
  if (!name || !strlen(name))
    bombElegant("no name given", NULL);
  str_toupper(name);
  if (has_wildcards(name) && strchr(name, '-'))
    name = expand_ranges(name);
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
  
  addTransmutationSpec(name, type, exclude, newType);
}

