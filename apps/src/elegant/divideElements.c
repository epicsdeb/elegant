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
  long divisions;
  double maximumLength;
} DIVIDE_SPEC;

static DIVIDE_SPEC *divisionSpec = NULL;
static long divisionSpecs = 0;

void addDivisionSpec(char *name, char *type, char *exclude,
		     long divisions,
		     double maximum_length)
{
  if (!(divisionSpec 
	= SDDS_Realloc(divisionSpec,
		       sizeof(*divisionSpec)*(divisionSpecs+1))))
    bombElegant("memory allocation failure", NULL);
  divisionSpec[divisionSpecs].name = NULL;
  divisionSpec[divisionSpecs].type = NULL;
  divisionSpec[divisionSpecs].exclude = NULL;
  divisionSpec[divisionSpecs].divisions = divisions;
  divisionSpec[divisionSpecs].maximumLength = maximum_length;
  if ((name &&
       !SDDS_CopyString(&divisionSpec[divisionSpecs].name, name)) ||
      (type &&
       !SDDS_CopyString(&divisionSpec[divisionSpecs].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&divisionSpec[divisionSpecs].exclude, exclude)))
    bombElegant("memory allocation failure", NULL);
  
  divisionSpecs++;
}

void clearDivisionSpecs() 
{
  while (divisionSpecs--) {
    if (divisionSpec[divisionSpecs].name)
      free(divisionSpec[divisionSpecs].name);
    if (divisionSpec[divisionSpecs].type)
      free(divisionSpec[divisionSpecs].type);
    if (divisionSpec[divisionSpecs].exclude)
      free(divisionSpec[divisionSpecs].exclude);
  }
  free(divisionSpec);
  divisionSpec = NULL;
}

long elementDivisions(char *name, char *type, double length)
{
  long i, div;
  for (i=0; i<divisionSpecs; i++) {
    if (divisionSpec[i].exclude && wild_match(name, divisionSpec[i].exclude))
      continue;
    if (divisionSpec[i].name && !wild_match(name, divisionSpec[i].name))
      continue;
    if (divisionSpec[i].type && !wild_match(type, divisionSpec[i].type))
      continue;
    if (divisionSpec[i].divisions)
      return divisionSpec[i].divisions;
    if (divisionSpec[i].maximumLength) {
      if ((div = (length/divisionSpec[i].maximumLength+0.5))<1)
	return 1;
      return div;
    }
    bombElegant("Invalid division specification seen.  This shouldn't happen.", 
	 NULL);
  }
  return 1;
}

#include "divideElements.h"

void setupDivideElements(NAMELIST_TEXT *nltext, RUN *run, 
			 LINE_LIST *beamline)
{
  long i;
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&divide_elements, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &divide_elements);

  if (clear) {
    clearDivisionSpecs();
    if (!divisions && !maximum_length && !name && !type)
      return;
  }

  if (divisions<=0 && maximum_length<=0)
    bombElegant("either divisions or maximum_length must be positive", NULL);
  if (divisions<0)
    bombElegant("divisions<0 doesn't make sense", NULL);
  if (maximum_length<0)
    bombElegant("maximum_length<0 doesn't make sense", NULL);
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
    if (i==N_TYPES)
      bombElegant("type pattern does not match any known type", NULL);
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }
  addDivisionSpec(name, type, exclude, divisions, maximum_length);
}

