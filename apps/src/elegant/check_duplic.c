/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: check_duplic.c
 * contents: check_duplic_elem(), check_duplic_line()
 * purpose: check for duplicate definitions of elements or beamlines.
 *          Subroutines should be called after each element/line is added
 *          since they only check to see that the most recently added
 *          item is not a duplication of another.
 *
 * Michael Borland, 1988 
 */
#include "mdb.h"
#include "track.h"

long check_duplic_elem(
                       ELEMENT_LIST **elem0,          /* root of element list */
                       ELEMENT_LIST **new_elem,       /* newly added element  */
                       char *nameToCheck,             /* don't add new element, just check for this name */
                       long n_elems
                       )
{
  char *new_name;
  long i, comparison, hi, lo, mid, checkOnly = 0;
  ELEMENT_LIST *elem, *insertionPoint, *elast;
  ELEMENT_LIST **elemArray = NULL;

  if (!elem0 || !*elem0)
    bombElegant("root pointer of element list is null (check_duplic_elem)", NULL);
  if ((!new_elem || !*new_elem) && !nameToCheck)
    bombElegant("new element pointer is null (check_duplic_elem)", NULL);
  
  if (n_elems>=1) {
    if (nameToCheck) {
      checkOnly = 1;
      new_name = nameToCheck;
    }
    else
      new_name = (*new_elem)->name;
#ifdef DEBUG    
    fprintf(stdout, "checking name %s against %ld elements in list\n",
            new_name, n_elems);
    fflush(stdout);
#endif
    
    elemArray = trealloc(elemArray, sizeof(*elemArray)*(n_elems+10));
    elem = *elem0;
    for (i=0; i<n_elems; i++) {
      elemArray[i] = elem;
      elem = elem->succ;
    }
    
    if (nameToCheck) {
      checkOnly = 1;
      new_name = nameToCheck;
    }
    else
      new_name = (*new_elem)->name;
    elem = *elem0;
    insertionPoint = NULL;

    lo = 0; 
    hi = n_elems-1;
    if ((comparison=strcmp(new_name, elemArray[lo]->name))==0) {
      if (checkOnly) {
        free(elemArray);
#ifdef DEBUG
        fprintf(stdout, "**** check reveals problem\n");
#endif
        return 1;
      }
      fprintf(stdout, "error: multiple definitions of element %s\n",
              new_name);
      fflush(stdout);
      exitElegant(1);
    }
    else if (comparison<0)
      insertionPoint = elemArray[lo];
    if (!insertionPoint) {
      if ((comparison=strcmp(new_name, elemArray[hi]->name))==0) {
        if (checkOnly) {
          free(elemArray);
#ifdef DEBUG
          fprintf(stdout, "**** check reveals problem\n");
#endif
          return 1;
        }
        fprintf(stdout, "error: multiple definitions of element %s\n",
                new_name);
        fflush(stdout);
        exitElegant(1);
      }
      else if (comparison<0) {
        do {
          mid = (lo+hi)/2;
          if ((comparison=strcmp(new_name, elemArray[mid]->name))==0) {
            if (checkOnly) {
              free(elemArray);
#ifdef DEBUG
              fprintf(stdout, "**** check reveals problem\n");
#endif
              return 1;
            }
            fprintf(stdout, "error: multiple definitions of element %s\n",
                    new_name);
            fflush(stdout);
            exitElegant(1);
          }
          else if (comparison<0)
            hi = mid;
          else
            lo = mid;
        } while ((hi-lo)>1);
        if (strcmp(new_name, elemArray[hi]->name)<0)
          insertionPoint = elemArray[hi];
      }
    }
    if (!checkOnly && insertionPoint) {
#ifdef DEBUG
      fprintf(stdout, "inserting %s before %s\n", new_name, insertionPoint->name);
      fflush(stdout);
#endif
      elast = (*new_elem)->pred;
      if ((*new_elem)->pred)
        (*new_elem)->pred->succ = NULL;
      (*new_elem)->pred = insertionPoint->pred;
      if (insertionPoint->pred)
        insertionPoint->pred->succ = (*new_elem);
      (*new_elem)->succ = insertionPoint;
      insertionPoint->pred = (*new_elem);
      if (insertionPoint == *elem0)
        *elem0 = (*new_elem);
      *new_elem = elast;
    }
    free(elemArray);
  }

#ifdef DEBUG
  elem  = *elem0;
  i = 0;
  fprintf(stdout, "Elements: \n");
  while (elem && elem->name) {
    fprintf(stdout, "  %s,%c", elem->name, ++i%10==0?'\n':' ');
    fflush(stdout);
    elem = elem->succ;
  }
  fprintf(stdout, "\n");
  fflush(stdout);
#endif

  return 0;
}

long check_duplic_line(
                       LINE_LIST *line,          /* root of line list */
                       char *new_name,
                       long n_lines,
                       long checkOnly
                       )
{
  long i;

  n_lines--;
#ifdef DEBUG
  fprintf(stdout, "Checking name %s against %ld lines\n", new_name, n_lines);
#endif
  for (i=0; i<n_lines; i++) {
#ifdef DEBUG
    fprintf(stdout, "%s,%c", line->name, (i+1)%10==0?'\n':' ');
#endif
    if (strcmp(new_name, line->name)==0) {
      if (checkOnly) {
#ifdef DEBUG
        fprintf(stdout, "+++ check reveals problem\n");
#endif
        return 1;
      }
      *line->name = 0;
      fprintf(stdout, "error: multiple definitions of line %s\n", 
             new_name);
      fflush(stdout);
      exitElegant(1);
    }
    line = line->succ;
  }
#ifdef DEBUG
  fprintf(stdout, "\n");
  fflush(stdout);
#endif
  return 0;
}


