/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file    : mad_parse.c
 * contents: fill_line(), fill_elem(), expand(), copy_element(),
 *           copy_line(), tell_type(), get_param_name(), find_param()
 *
 *           Routines to fill LINE and ELEMENT structures by parsing
 *           a line of text.
 *
 * Michael Borland, 1987, 1988, 1989, 1991.
 */
#include "mdb.h"
#include "track.h"
#include <ctype.h>
#include <string.h>
#include <memory.h>

long is_simple(char *s);
void copy_p_elem(char *target, char *source, long type);

long max_name_length = 100;

long set_max_name_length(long length)
{
    long tmp;
    if (length<=0)
        bombElegant("maximum name length set to nonpositive number", NULL);
    tmp = max_name_length;
    max_name_length = length;
    return(tmp);
    }

/* routine: fill_line()
 * purpose: fill a LINE_LIST structure by parsing a line of text.
 */

void fill_line(
    LINE_LIST *line,         /* linked list of line definitions */
    long nl,                 /* number of line definitions so far */
    ELEMENT_LIST *elem,      /* linked list of element definitions */
    long ne,                 /* number of element definitions so far */
    char *s                  /* text to parse into a new line definition */
    )
{
    register char *ptr;
    register long i;
    ELEMENT_LIST *leptr;
    LINE_LIST *lptr;

    log_entry("fill_line");

#ifdef DEBUG    
    fprintf(stdout, "expanding beamline: %s\n", s);
    fflush(stdout);
#endif

    /* get pointer to empty spot in list of LINE definitions */
    lptr = line;
    for (i=0; i<nl; i++) 
        lptr = lptr->succ;
    if (lptr==NULL) 
        bombElegant("improperly extended linked-list (fill_line)", NULL);

    /* allocate string in which to save original beamline definition for later use */
    lptr->definition = tmalloc(sizeof(char)*(strlen(s)+12));

    /* get beamline name */
    if ((ptr = strchr(s, ','))==NULL) {
        fprintf(stdout, "error parsing LINE: %s\n", s);
        fflush(stdout);
        exitElegant(1);
        }
    *ptr = 0;
    lptr->name = get_token_t(s, " ,\011");
#ifdef DEBUG
    fprintf(stdout, "line name is %s\n", lptr->name);
    fflush(stdout);
#endif
    if (strpbrk(lptr->name, ":.,/-_+abcdefghijklmnopqrstuvwyxz "))
        sprintf(lptr->definition, "\"%s\": LINE=%s", lptr->name, ptr+1);
    else
        sprintf(lptr->definition, "%s: LINE=%s", lptr->name, ptr+1);
    ptr += 2;

    leptr = &(lptr->elem);
    lptr->n_elems = 0;
    lptr->flags = 0;

    expand_line(leptr, lptr, ptr, line, nl, elem, ne, lptr->name);
    while (leptr->succ!=NULL) {
      leptr = leptr->succ;
    }
    if (leptr->pred==NULL)
      bombElegant("programming error: problem in linked list of beamlines", NULL);
    (leptr->pred)->succ = NULL;
    log_exit("fill_line");
    }

extern ELEMENT_LIST *expand_line(
    ELEMENT_LIST *leptr,     /* pointer to list of elements for new line */
    LINE_LIST *lptr,         /* pointer to the line structure new line */
    char *s,                 /* string defining the new line */
    LINE_LIST *line,         /* pointer to linked-list of line structures for existing lines */
    long nl,                 /* number of existing lines */
    ELEMENT_LIST *elem,      /* pointer to linked-list of element structures for defined elements */
    long ne,                 /* number of defined elements */
    char *part_of           /* name of the group that elements will be part of if they are not members of sub-lines */
    )
{
    register char *ptr, *ptr1;
    char *ptrs=NULL, *ptr2=NULL;
    register long i, l;
    long multiplier=1, reverse=0, simple, count1;
    char *line_save=NULL;

    log_entry("expand_line");

    cp_str(&line_save, s);

#ifdef DEBUG
    fprintf(stdout, "line is %s\n", s);
    fflush(stdout);
#endif

    count1 = 0;
    while ((ptr2 = get_token_tq(s, ",", ",", "(\"", ")\""))) {
      ptr1 = ptr2;
      ptr = ptr1+strlen(ptr1)-1;
      while (isspace(*ptr))
        ptr--;
      if (*ptr==')')
        *ptr = 0;
      delete_bounding(ptr1, "\"");
#ifdef DEBUG 
      fprintf(stdout, "ptr1 = %s\n", ptr1);
      fflush(stdout);
#endif
      reverse = 0;
      if (*ptr1=='-') {
        reverse = 1;
        ptr1++;
      }
      
      multiplier = 1;
      if (isdigit(*ptr1)) {
        if ((ptr=strchr(ptr1, '*'))) {
          *ptr = 0;
          if (1!=sscanf(ptr1, "%ld", &multiplier)) {
            fprintf(stdout, "problem with line: %s\n", line_save);
            fflush(stdout);
            fprintf(stdout, "position is %s\n", ptr1);
            fflush(stdout);
            bombElegant("improper element multiplication", NULL);
          }
          ptr1 = ptr+1;
          if (*ptr1=='(' || *ptr1=='"')
            ptr1++;
        }
      }
#ifdef DEBUG
      fprintf(stdout, "reverse = %ld, multiplier = %ld\n", reverse, multiplier);
      fflush(stdout);
      fprintf(stdout, "ptr1 = %s\n",ptr1);
      fflush(stdout);
#endif
      
      simple = is_simple(ptr1);
#ifdef DEBUG
      fprintf(stdout, "simple = %ld\n", simple);
      fflush(stdout);
#endif
      count1 += multiplier;
      if (simple) {
        /* add simple elements to the line's element list */
        lptr->n_elems += 
          (l=expand_phys(leptr, ptr1, elem, ne, line, nl, reverse, multiplier, part_of));
#ifdef DEBUG
        fprintf(stdout, "number of elements now %ld\n", lptr->n_elems);
#endif
        for (i=0; i<l; i++) 
          leptr = leptr->succ;
      }
      else {
        /* add elements from a parenthesized list to the line's element list */
        cp_str(&ptrs, ptr1);
        for (i=0; i<multiplier; i++) {
          strcpy_ss(ptr1, ptrs);
          leptr = expand_line(leptr, lptr, ptr1, line, nl, elem, ne, part_of);
        }
#ifdef DEBUG
        fprintf(stdout, "number of elements now %ld\n", lptr->n_elems);
#endif
        tfree(ptrs);
        ptrs = NULL;
      }
      free(ptr2);
    }

    if (count1>lptr->n_elems) {
      fprintf(stdout, "Problem with line parsing: expected at least %ld elements, but got %ld\n",
              count1, lptr->n_elems);
      print_line(stdout, lptr);
    }
    if (line_save) free(line_save);
    log_exit("fill_line");
    return(leptr);
}

long is_simple(char *s)
{
    register char *ptr;
    log_entry("is_simple");
    ptr = s;
    while (*ptr) {
        if (*ptr=='(' || *ptr==',' || *ptr=='*') {
            log_exit("is_simple");
            return(0);
            }
        ptr++;
        }
    log_exit("is_simple");
    return(1);
    }


/* routine: fill_elem()
 * purpose: fill a single physical element slot in an element list
 */

void fill_elem(ELEMENT_LIST *eptr, char *s, long type, FILE *fp_input)
{
  BEND *bptr;
  MATR *matr;

  log_entry("fill_elem");
  eptr->end_pos = eptr->flags = 0;

    if ((eptr->type = type)>=N_TYPES || type<=0) {
        fprintf(stdout, "unknown element type %ld in fill_elem()\n", type);
        fflush(stdout);
        fprintf(stdout, "remainder of line is: \n%s\n", s);
        fflush(stdout);
        exitElegant(1);
        }

    eptr->name = get_token_t(s, " ,\011");
    if (((long)strlen(eptr->name))>max_name_length) {
        fprintf(stdout, "warning: element name %s truncated to %ld  characters\n",
                eptr->name, max_name_length);
        fflush(stdout);
        eptr->name[max_name_length] = 0;
        }
    cp_str(&eptr->definition_text, s);
    eptr->matrix = NULL;
    eptr->group = NULL;
    eptr->p_elem = tmalloc(entity_description[type].structure_size);
    eptr->p_elem0 = tmalloc(entity_description[type].structure_size);
    zero_memory(eptr->p_elem, entity_description[type].structure_size);
    zero_memory(eptr->p_elem0, entity_description[type].structure_size);

    parse_element((char*)(eptr->p_elem), 
            entity_description[type].parameter, 
            entity_description[type].n_params,
            s, eptr, entity_name[type]);

    switch (type) {
    case T_RBEN:
      bptr = (BEND*)(eptr->p_elem);
      bptr->e1 += bptr->angle/2;
      bptr->e2 += bptr->angle/2;
      type = eptr->type = T_SBEN;
      if (fabs(bptr->angle)>1e-14)
        bptr->length *= (bptr->angle/2)/sin(bptr->angle/2);
      break;
    case T_PEPPOT:
      parse_pepper_pot((PEPPOT*)(eptr->p_elem), fp_input, eptr->name);
      break;
    case T_MATR:
      matr = (MATR*)eptr->p_elem;
      if (!matr->matrix_read) {
        char *filename;
        FILE *fpm;
        if (!(filename=findFileInSearchPath(matr->filename))) {
          fprintf(stderr,"Unable to find MATR file %s\n", matr->filename);
          exitElegant(1);
        }
        fprintf(stdout, "File %s found: %s\n", matr->filename, filename);
        fpm = fopen_e(filename, "r", 0);
        free(filename);
        matr->M.order = matr->order;
        initialize_matrices(&(matr->M), matr->order);
        if (!read_matrices(&(matr->M), fpm)) {
          fprintf(stdout, "error reading matrix from file %s\n", matr->filename);
          fflush(stdout);
          abort();
        }
        fclose(fpm);
        matr->matrix_read = 1;
        matr->length = matr->M.C[4];
      }
      break;
    case T_SCRAPER:
      interpretScraperDirection((SCRAPER*)(eptr->p_elem));
      break;
    default:
      break;
    }
        
    eptr->flags = PARAMETERS_ARE_STATIC;     /* default, to be changed by variation and error settings */

    copy_p_elem(eptr->p_elem0, eptr->p_elem, type);

    log_exit("fill_elem");
    }

/* routine: copy_named_element()
 * purpose: fill a single physical element slot in an element list
 *          by copying a named, previously-defined element
 */

void copy_named_element(ELEMENT_LIST *eptr, char *s, ELEMENT_LIST *elem)
{
    char *name, *match;

    log_entry("copy_named_element");
    eptr->end_pos = eptr->flags = 0;

    if (!(name=get_token_tq(s, "", ":", " \"", " \""))) {
        log_exit("copy_named_element");
        return;
        }
    delete_bounding(name, "\"");

    if (!(match = get_token_tq(s, "", ",=", " \"", " \""))) {
        log_exit("copy_named_element");
        bombElegant("programming error--unable to recover source element name (copy_named_element)", NULL);
        }
    delete_bounding(match, "\"");

    if (((long)strlen(name))>max_name_length) {
        fprintf(stdout, "warning: element name %s truncated to %ld characters\n",
                eptr->name, max_name_length);
        fflush(stdout);
        name[max_name_length] = 0;
        }

#ifdef DEBUG
    fprintf(stdout, "seeking match for name %s to copy for defining name %s\n",
           match, name);
    fflush(stdout);
#endif

    while (elem && elem->name) {
        if (strcmp(elem->name, match)==0)
            break;
        elem = elem->succ;
        }
    if (!elem->name) {
        log_exit("copy_named_element");
        fprintf(stdout, "unable to define %s as copy of %s--source element does not exist\n",
                name, match);
        fflush(stdout);
        exitElegant(1);
        }
#ifdef DEBUG
    fprintf(stdout, "match found for %s: %s\n", name, match);
    fflush(stdout);
#endif

    eptr->p_elem = tmalloc(entity_description[elem->type].structure_size);
    eptr->p_elem0 = tmalloc(entity_description[elem->type].structure_size);
    copy_p_elem(eptr->p_elem, elem->p_elem, elem->type);
    copy_p_elem(eptr->p_elem0, eptr->p_elem, elem->type);
    cp_str(&eptr->name, name);
    eptr->group = NULL;
    if (elem->group)
      cp_str(&eptr->group, elem->group);
    eptr->matrix = NULL;
    eptr->flags = PARAMETERS_ARE_STATIC;
    eptr->type = elem->type;
    log_exit("copy_named_element");
    }

static long *entity_name_length = NULL;

/* routine: expand_phys()
 * purpose: expand a named entity into physical elements 
 */

long expand_phys(
    ELEMENT_LIST *leptr,         /* list being filled with expansion */
    char *entity,                /* entity to expand */
    ELEMENT_LIST *elem_list,     /* list of all physical elements */
    long ne,                     /* number of same */
    LINE_LIST *line_list,        /* list of already-existing LINEs */
    long nl,                     /* number of same */
    long reverse, 
    long multiplier,
    char *part_of                /* name of line this is to be part of */
    )
{
  long ie, il, i, j, comparison, div;
  char trunc_char;
  ELEMENT_LIST *elem0;
  double length;

  log_entry("expand_phys");
  
  /* search for truncated name, for the case that it turns out
   * to be an element */
  trunc_char = 0;
  if (((long)strlen(entity))>max_name_length) {
    trunc_char = entity[max_name_length];
    entity[max_name_length] = 0;
  }
  elem0 = elem_list;
  div = 1;
  for (ie=0; ie<ne; ie++) {
    if ((comparison=strcmp(entity, elem_list->name))==0) {
#ifdef DEBUG
      fprintf(stdout, "adding element %ld*%s\n", multiplier, entity);
      fflush(stdout);
#endif
      if (entity_description[elem_list->type].flags&DIVIDE_OK &&
	  entity_description[elem_list->type].flags&HAS_LENGTH &&
	  (length = *(double*)(elem_list->p_elem))>0) {
	div = elementDivisions(entity, entity_name[elem_list->type], length);
        if (div>1)
          elem_list->divisions = div;
      }
#ifdef DEBUG
      fprintf(stderr,  "Dividing %s %ld times\n",
	      entity, div);
#endif
      for (i=0; i<multiplier; i++) {
	long j;
	for (j=0; j<div; j++) {
	  copy_element(leptr, elem_list, reverse, j, div);
	  leptr->part_of = elem_list->part_of?elem_list->part_of:part_of;
#ifdef DEBUG
          fprintf(stdout, "expand_phys copied: name=%s divisions=%ld first=%hd j=%ld\n",
                  leptr->name, leptr->divisions, leptr->firstOfDivGroup, j);
#endif
	  extend_elem_list(&leptr);
	}
      }
      return(multiplier*div);
    }
    if (comparison<0)
      break;
      elem_list = elem_list->succ;
  }
  
  /* it isn't an element, so search list of beam-lines for occurence
   * of the full name 
   */
  if (trunc_char)
    entity[max_name_length] = trunc_char;
  for (il=0; il<nl; il++) {
    if (strcmp(line_list->name, entity)==0) {
      ie = 0;
      for (i=0; i<multiplier; i++) {
	ie += copy_line(leptr, &(line_list->elem), line_list->n_elems, reverse, entity);
	for (j=0; j<line_list->n_elems; j++) 
	  leptr = leptr->succ;
      }
      log_exit("expand_phys");
      return(ie);
    }
    line_list = line_list->succ;
  }
  
  fprintf(stdout, "no expansion for entity %s\n", entity);
  fflush(stdout);
  fprintf(stdout, "known elements are:\n");
  fflush(stdout);
  print_elem_names(stdout, elem0, 100);
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); /* Make sure the information can be printed before aborting */
#endif
  exitElegant(1);
  return(0);
}

/* routine: copy_element()
 * purpose: copy an element 
 */

void copy_element(ELEMENT_LIST *e1, ELEMENT_LIST *e2, long reverse, long division, long divisions)
{

    log_entry("copy_element");

    cp_str(&e1->name, e2->name);
    e1->group = NULL;
    if (e2->group)
      cp_str(&e1->group, e2->group);
    e1->end_pos = 0;
    e1->flags   = 0;
    e1->matrix  = NULL;
    e1->type    = e2->type;
    e1->p_elem  = tmalloc(entity_description[e1->type].structure_size);
    e1->p_elem0 = tmalloc(entity_description[e1->type].structure_size);
    copy_p_elem(e1->p_elem, e2->p_elem, e1->type);
    e1->divisions = 1;
    if (reverse) {
      BEND *bptr;
      KSBEND *ksbptr;
      CSBEND *csbptr;
      CSRCSBEND *csrbptr;
      NIBEND *nibptr;
      switch (e1->type) {
      case T_SBEN:
      case T_RBEN:
        bptr = (BEND*)e1->p_elem;
        SWAP_DOUBLE(bptr->e1, bptr->e2);
        SWAP_LONG(bptr->edge1_effects, bptr->edge2_effects);
        break;
      case T_KSBEND:
        ksbptr = (KSBEND*)e1->p_elem;
        SWAP_DOUBLE(ksbptr->e1, ksbptr->e2);
        break;
      case T_NIBEND:
        nibptr = (NIBEND*)e1->p_elem;
        SWAP_DOUBLE(nibptr->e1, nibptr->e2);
        break;
      case T_CSBEND:
        csbptr = (CSBEND*)e1->p_elem;
        SWAP_DOUBLE(csbptr->e1, csbptr->e2);
        SWAP_LONG(csbptr->edge1_effects, csbptr->edge2_effects);
        break;
      case T_CSRCSBEND:
        csrbptr = (CSRCSBEND*)e1->p_elem;
        SWAP_DOUBLE(csrbptr->e1, csrbptr->e2);
        SWAP_LONG(csrbptr->edge1_effects, csrbptr->edge2_effects);
        break;
      default:
        break;
      }
    }
    e1->firstOfDivGroup = 0;
    if (divisions>1) {
      if (division==0)
        e1->firstOfDivGroup = 1;
      if (entity_description[e1->type].flags&HAS_LENGTH)
        *(double*)(e1->p_elem) /= divisions;
      if (e1->type==T_RFCA) {
        RFCA *rfca;
        rfca = (RFCA*)e1->p_elem;
        rfca->volt /= divisions;
      }
      else if (IS_BEND(e1->type)) {
        BEND *bptr; CSBEND *csbptr;
        switch (e1->type) {
        case T_SBEN:
          bptr = (BEND*)e1->p_elem;
          bptr->angle /= divisions;
          break;
        case T_CSBEND:
          csbptr = (CSBEND*)e1->p_elem;
          csbptr->angle /= divisions;
          csbptr->n_kicks = csbptr->n_kicks / divisions + 1;
          break;
        default:
          printf("Internal error: Attempt to divide angle for element that is not supported (type=%ld)\n",
                 e1->type);
          exitElegant(1);
          break;
        }
      }
      e1->divisions = divisions;
    }
    copy_p_elem(e1->p_elem0, e2->p_elem, e1->type);
    log_exit("copy_element");
  }

/* routine: copy_line
 * purpose: copy a series of elements from a linked-list 
 */

long copy_line(ELEMENT_LIST *e1, ELEMENT_LIST *e2, long ne, long reverse, char *part_of)
{
    register long i;

    log_entry("copy_line");
    
    if (!reverse) {
        for (i=0; i<ne; i++) {
            copy_element(e1, e2, reverse, 0, 0);
            e1->part_of = e2->part_of?e2->part_of:part_of;
            e1->divisions = e2->divisions;
            e1->firstOfDivGroup = e2->firstOfDivGroup;
            extend_elem_list(&e1);
            e2 = e2->succ;
            }                
        }
    else {
        for (i=0; i<(ne-1); i++)
            e2 = e2->succ;
        for (i=0; i<ne; i++) {
            copy_element(e1, e2, reverse, 0, 0);
            e1->part_of = e2->part_of?e2->part_of:part_of;
            e1->divisions = e2->divisions;
            e1->firstOfDivGroup = e2->firstOfDivGroup;
            extend_elem_list(&e1);
            e2 = e2->pred;
            } 
        }
    log_exit("copy_line");
    return(ne);
    }

/* routine: tell_type()
 * purpose: return the integer type code for a physical element 
 *         
 */

long tell_type(char *s, ELEMENT_LIST *elem)
{
    long i, l, match_found, return_value, comparison;
    char *ptr, *ptr1, *name, *name1, *buffer;
    
    match_found = 0;
    return_value = T_NODEF;

    log_entry("tell_type");

    if (!(name1=get_token_tq(s, "", ":", " \"", " \""))) {
        log_exit("tell_type");
        return(T_NODEF);
        }
    name = name1;
    delete_bounding(name, "\"");
#ifdef DEBUG
    fprintf(stdout, "first token on line: >%s<\n", name);
    fflush(stdout);
    fprintf(stdout, "remainder of line: >%s<\n", s);
    fflush(stdout);
#endif
    ptr1 = get_token_tq(s, "", ",=", " \"", " \"");
    ptr = ptr1;
    delete_bounding(ptr, "\"");
#ifdef DEBUG
    fprintf(stdout, "second token: >%s<\n", ptr);
    fflush(stdout);
#endif
    if (!ptr || is_blank(ptr)) {
#ifdef DEBUG
        fprintf(stdout, "second token is blank\n");
        fflush(stdout);
#endif
        if ((ptr=strchr(name, ','))) {
            strcpy_ss(s, ptr+1);
            *ptr = 0;
            }
        ptr = name;
        name = NULL;
        }
    else {
        insert(s, ",");
        insert(s, name);
        /* s now has the form "name,<remainder-of-s>".  The original second token of s is
         * pointed to by ptr, but is not in s. */
        }
#ifdef DEBUG
    fprintf(stdout, "seeking to match %s to entity name\n", ptr);
    fflush(stdout);
#endif
    if (entity_name_length==NULL) {
        entity_name_length = tmalloc(sizeof(*entity_name_length)*N_TYPES);
        for (i=0; i<N_TYPES; i++)
            entity_name_length[i] = strlen(entity_name[i]);
        }
    l = strlen(ptr);
    for (i=0; i<N_TYPES; i++) {
        if (strncmp(ptr, entity_name[i], MIN(l, entity_name_length[i]))==0)  {
            if (match_found) {
                fprintf(stdout, "error: item %s is ambiguous--specify more of the item name\n", ptr);
                fflush(stdout);
                exitElegant(1);
                }
            if (l==entity_name_length[i]) {
                if (name1) free(name1);
                if (ptr1) free(ptr1);
                log_exit("tell_type");
                return(i);
                }
            return_value = i;
            match_found = 1;
#ifdef DEBUG
            fprintf(stdout, "%s matches entity name %s\n", ptr, entity_name[i]);
            fflush(stdout);
#endif
            }
        }
#ifdef DEBUG
    fprintf(stdout, "seeking to match %s to command name\n", ptr);
    fflush(stdout);
#endif
    for (i=1; i<N_MADCOMS; i++) {
        if (strncmp(ptr, madcom_name[i], MIN(l, strlen(madcom_name[i])))==0) {
            if (match_found) {
                fprintf(stdout, "error: item %s is ambiguous--specify more of the item name\n", ptr);
                fflush(stdout);
                exitElegant(1);
                }
            if (l==((long)strlen(madcom_name[i]))) {
                if (name1) free(name1);
                if (ptr1) free(ptr1);
                log_exit("tell_type");
                return(-i-1);
                }
            return_value = -i-1;
            match_found = 1;
#ifdef DEBUG
            fprintf(stdout, "%s matches command name %s\n", ptr, madcom_name[i]);
            fflush(stdout);
#endif
            }
        }
#ifdef DEBUG
    fprintf(stdout, "seeking to match %s to element name\n", ptr);
    fflush(stdout);
#endif
    /* this code depends on the elements being previously sorted into
     * strcmp() order.  This is done by check_duplic_elem
     */
    while (elem && elem->name) {
        if ((comparison=strcmp(elem->name, ptr))==0) {
            if (match_found)
                fprintf(stdout, "warning: reference to item %s is ambiguous--assuming element copy desired\n", ptr);
                fflush(stdout);
            if (!elem->definition_text)
                bombElegant("element copy with no definition_text--internal error", NULL);
            buffer = tmalloc(sizeof(*buffer)*(strlen(elem->definition_text)+strlen(s)+1));
            if ((ptr=strchr(s, ','))) {
#ifdef DEBUG
                fprintf(stdout, "inserting into: %s\n", s);
                fflush(stdout);
                fprintf(stdout, "inserting: %s\n", elem->definition_text);
                fflush(stdout);
#endif
                insert(ptr, elem->definition_text);
                }
            else {
                strcat(s, ",");
                strcat(s, buffer);
                }
            free(buffer);
#ifdef DEBUG
            fprintf(stdout, "modified input text for element %s:\n%s\n", name, s);
            fflush(stdout);
#endif
            if (name1) free(name1);
            if (ptr1) free(ptr1);
            log_exit("tell_type");
            return(elem->type);
            }
        if (comparison>1)
            break;
        elem = elem->succ;
        }
    if (name1) free(name1);
    if (ptr1) free(ptr1);
    log_exit("tell_type");
    return(return_value);
    }


/* routine: get_param_name()
 * purpose: extract the parameter name from a string of the type
 *          "L=3.0331".  A pointer to the name is returned and the
 *          string is altered to bring the number to the beginning.
 */

char *get_param_name(char *s)
{
    char *ptr, *ptr1=NULL;

    log_entry("get_param_name");

    if ((ptr=strchr(s, '='))==NULL) {
        fprintf(stdout, "get_param_name(): no parameter name found in string %s\n", s);
        fflush(stdout);
        exitElegant(1);
        }
    *ptr = 0;
    cp_str(&ptr1, s);
    strcpy_ss(s, ptr+1);
    log_exit("get_param_name");
    return(ptr1);
    }

/* routine: find_param
 * purpose: find the element parameter with a given name, and return a
 *          pointer to the name (or NULL if not found).  The string is
 *          altered so that the value to be scanned appears at the
 *          beginning.
 */

char *find_param(char *s, char *param)
{
    char *ptr;
    long l; 

    log_entry("find_param");

    if (!*s) {
        log_exit("find_param");
        return(NULL);
        }

    l = strlen(param);

    while ((ptr=get_param_name(s))) {
        if (strncmp(param, ptr, l)==0)  {
            log_exit("find_param");
            return(ptr);
            }
        }
    log_exit("find_param");
    return(NULL);
    }

void unknown_parameter(char *parameter, 
    char *element, char *type_name, char *caller)
{
    fprintf(stdout, "error: unknown parameter %s used for %s %s (%s)\n",
            parameter, type_name, element, caller);
    fflush(stdout);
    exitElegant(1);
    }

void parse_element(
    char *p_elem,
    PARAMETER *parameter,
    long n_params,
    char *string,
    ELEMENT_LIST *eptr,
    char *type_name
    )
{
  long i, difference, isGroup, pType=0;
  char *ptr, *ptr1=NULL, *rpn_token;

  log_entry("parse_element");

  if (!(entity_description[eptr->type].flags&OFFSETS_CHECKED)) {
    if (parameter[0].offset<0) {
      fprintf(stdout, "error: bad initial parameter offset for element type %s\n", type_name);
      fflush(stdout);
      exitElegant(1);
    }
    for (i=1; i<n_params; i++) {
      if ((difference=parameter[i].offset-parameter[i-1].offset)<0) {
        fprintf(stdout, "error: bad parameter offset for element type %s, parameter %ld (%s)\n",
                type_name, i, parameter[i].name?parameter[i].name:"NULL");
        fflush(stdout);
        exitElegant(1);
      }
    }
    entity_description[eptr->type].flags |= OFFSETS_CHECKED;
  }

  for (i=0; i<n_params; i++)  {
    switch (parameter[i].type) {
    case IS_DOUBLE:
      *(double*)(p_elem+parameter[i].offset) = parameter[i].number;
      break;
    case IS_LONG:
      *(long*)(p_elem+parameter[i].offset) = parameter[i].integer;
      break;
    case IS_STRING:
      if (parameter[i].string==NULL)
        *(char**)(p_elem+parameter[i].offset) = NULL;
      else
        cp_str((char**)(p_elem+parameter[i].offset), 
               parameter[i].string);
      break;
    }
  }
  
  while ((ptr=get_token(string))) {
#ifdef DEBUG
    fprintf(stdout, "Parsing %s\n", ptr);
#endif
    ptr1 = get_param_name(ptr);
    if (!ptr1) {
      fprintf(stdout, "error getting parameter name for element %s\n", eptr->name);
      exitElegant(1);
    }
    isGroup = 0;
    if (strcmp(ptr1, "GROUP")==0) {
      pType = IS_STRING;
      isGroup = 1;
    } else {
      for (i=0; i<n_params; i++) 
        if (strcmp(parameter[i].name, ptr1)==0)
          break;
      if (i==n_params) {
        /* try again, ignoring underscores */
        for (i=0; i<n_params; i++) {
          if (strcmp_skip(parameter[i].name, ptr1, "_")==0)
            break;
        }
      }
      if (i<n_params)
        pType = parameter[i].type;
    }
    if (i==n_params && !isGroup) {
      fprintf(stdout, "Error: unknown parameter %s used for %s %s (%s)\n",
              ptr1, eptr->name, type_name, "parse_element");
      fflush(stdout);
      fputs("valid parameters are:", stdout);
      for (i=0; i<n_params; i++) {
        switch (parameter[i].type) {
        case IS_DOUBLE:
          fprintf(stdout, "%s (%.16f %s)\n",
                  parameter[i].name, parameter[i].number,
                  parameter[i].unit);
          fflush(stdout);
          break;
        case IS_LONG:
          fprintf(stdout, "%s (%ld %s)\n",
                  parameter[i].name, parameter[i].integer,
                  parameter[i].unit);
          fflush(stdout);
          break;
        case IS_STRING:
          fprintf(stdout, "%s (\"%s\")\n",
                  parameter[i].name, 
                  parameter[i].string==NULL?"{null}":
                  parameter[i].string);
          fflush(stdout);
          break;
        }
      }
      exitElegant(1);
    }
    free(ptr1);
    if (!*ptr) {
      fprintf(stdout, "Error: missing value for parameter %s of %s\n",
              parameter[i].name, eptr->name);
      exitElegant(1);
    }
    switch (pType) {
    case IS_DOUBLE:
      if (!isdigit(*ptr) && *ptr!='.' && *ptr!='-' && *ptr!='+') {
        rpn_token = get_token(ptr);
        SDDS_UnescapeQuotes(rpn_token, '"');
        *((double*)(p_elem+parameter[i].offset)) = rpn(rpn_token);
        if (rpn_check_error()) exitElegant(1);
        fprintf(stdout, "computed value for %s.%s is %.15e\n", eptr->name, parameter[i].name, 
                *((double*)(p_elem+parameter[i].offset)));
        fflush(stdout);
      }
      else {
        if (sscanf(ptr, "%lf", (double*)(p_elem+parameter[i].offset))!=1) {
          printf("Error scanning token %s for double value for parameter %s of %s.  Please check syntax.\n", 
                 ptr, parameter[i].name, eptr->name);
          exitElegant(1);
        }
      }
      break;
    case IS_LONG:
      if (!isdigit(*ptr) && *ptr!='-' && *ptr!='+') {
        rpn_token = get_token(ptr);
        SDDS_UnescapeQuotes(rpn_token, '"');
        *((long*)(p_elem+parameter[i].offset)) = rpn(rpn_token);
        if (rpn_check_error()) exitElegant(1);
        fprintf(stdout, "computed value for %s.%s is %ld\n", 
                eptr->name, parameter[i].name, 
                *((long*)(p_elem+parameter[i].offset)));
        fflush(stdout);
      } 
      else {
        if (sscanf(ptr, "%ld", (long*)(p_elem+parameter[i].offset))!=1) {
          printf("Error scanning token %s for integer value for parameter %s of %s.  Please check syntax.\n", 
                 ptr, parameter[i].name, eptr->name);
          exitElegant(1);
        }
      }
      break;
    case IS_STRING:
      if (!isGroup)
        *(char**)(p_elem+parameter[i].offset) = get_token(ptr);
      else
        eptr->group = get_token(ptr);
      break;
    default:
      bombElegant("unknown data type in parse_element!", NULL);
      break;
    }
    free(ptr);
  }
  log_exit("parse_element");
}

void parse_pepper_pot(
    PEPPOT *peppot,
    FILE *fp,
    char *name
    )
{
    long n_holes, i_hole;
    double *xc, *yc, x, y, cos_tilt, sin_tilt;
    char s[100];

    log_entry("parse_pepper_pot");

    n_holes = peppot->n_holes;
    peppot->x = xc = array_1d(sizeof(double), 0, n_holes-1);
    peppot->y = yc = array_1d(sizeof(double), 0, n_holes-1);
    i_hole = 0;
    sin_tilt = sin(peppot->tilt);
    cos_tilt = cos(peppot->tilt);

    for (i_hole=0; i_hole<n_holes; i_hole++) {
        if (!fgets(s, 100, fp) ||
                !get_double(&x, s) || !get_double(&y, s)) {
            fprintf(stdout, "error: data missing for pepper-pot plate %s\n",
                name);
            fflush(stdout);
            exitElegant(1);
            }
        xc[i_hole] =  x*cos_tilt + y*sin_tilt;
        yc[i_hole] = -x*sin_tilt + y*cos_tilt;
        }

/*
    fprintf(stdout, "pepper pot plate %s:\n", name);
    fflush(stdout);
    fprintf(stdout, "L=%lf  RADII=%lf  N_HOLES=%ld\n", 
            peppot->length, peppot->radii, peppot->n_holes);
    fflush(stdout);
    for (i_hole=0; i_hole<n_holes; i_hole++) 
        fprintf(stdout, "x = %lf    y = %lf\n",
                peppot->x[i_hole], peppot->y[i_hole]);
        fflush(stdout);
 */
    log_exit("parse_pepper_pot");
    }

void interpretScraperDirection(SCRAPER *scraper)
{
  if (scraper->insert_from) {
    switch (toupper(scraper->insert_from[1])) {
    case 'Y': case 'V':
      scraper->direction = 1;
      break;
    case 'X': case 'H':
      scraper->direction = 0;
      break;
    default:
      fprintf(stdout, "Error: invalid scraper insert_from parameter: %s\n",
              scraper->insert_from);
      fflush(stdout);
      bombElegant("scraper insert_from axis letter is not one of x, h, y, or v", NULL);
      break;
    }
    if (scraper->insert_from[0]=='-')
      scraper->direction += 2;
    scraper->insert_from = NULL;
  }
}

void copy_p_elem(char *target, char *source, long type)
{
  long i;

  for (i=0; i<entity_description[type].n_params; i++) {
    switch (entity_description[type].parameter[i].type) {
    case IS_DOUBLE:
      *((double*)(target+entity_description[type].parameter[i].offset)) = 
	*((double*)(source+entity_description[type].parameter[i].offset)) ;
      break;
    case IS_LONG:
      *((long*)(target+entity_description[type].parameter[i].offset)) = 
	*((long*)(source+entity_description[type].parameter[i].offset)) ;
      break;
    case IS_STRING:
      cp_str(((char**)(target+entity_description[type].parameter[i].offset)),
	     *((char**)(source+entity_description[type].parameter[i].offset)) );
      break;
    default:
      bombElegant("invalid parameter type (copy_p_elem).  Seek professional help!", NULL);
      break;
    }
  }

  if (type==T_PEPPOT)  {
    /* Need to modernize pepper-pot element to read data from SDDS file */
    PEPPOT *pps, *ppt;
    pps = (PEPPOT*)source;
    ppt = (PEPPOT*)target;
    ppt->x = array_1d(sizeof(double), 0, ppt->n_holes-1);
    ppt->y = array_1d(sizeof(double), 0, ppt->n_holes-1);
    memcpy(ppt->x, pps->x, sizeof(double)*ppt->n_holes);
    memcpy(ppt->y, pps->y, sizeof(double)*ppt->n_holes);
  }

}
