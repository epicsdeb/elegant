/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: get_beamline()
 * purpose: read a mad-format lattice and return a pointer to a linked
 *          list for the beamline specified. 
 * 
 *	    It is assumed that the MAD-style input file is in
 *          the usual units of meters, radians, etc.  The 
 *          output has the same units.
 *
 * Michael Borland, 1989
 */
#include "mdb.h"
#include "track.h"
#include <ctype.h>
#include "match_string.h"

void show_elem(ELEMENT_LIST *eptr, long type);
void process_rename_request(char *s, char **name, long n_names);

/* elem: root of linked-list of ELEM structures 
 * This list contains the definitions of all elements as supplied in the
 * input file.  An important use of this list is to keep track of the lattice
 * that will be saved with save_lattice.
 */
static ELEMENT_LIST *elem;   

/* line: root of linked-list of LINE structures 
 * This list contains the definitions of all beamlines as supplied in the
 * input file.  Each beamline contains instances of the elements that it
 * contains.  I.e., it does not refer explicitly to the structures in
 * elem
 */
static LINE_LIST *line;      

typedef struct input_object {
  void *ptr;   /* points to an ELEMENT_LIST or LINE_LIST */
  long isLine;
  struct input_object *next;
} INPUT_OBJECT;
static INPUT_OBJECT inputObject, *lastInputObject=NULL;

void addToInputObjectList(void *ptr, long isLine)
{
  if (lastInputObject==NULL)
    lastInputObject = &inputObject;
  else {
    lastInputObject->next = tmalloc(sizeof(INPUT_OBJECT));
    lastInputObject = lastInputObject->next;
  }
  lastInputObject->ptr = ptr;
  lastInputObject->isLine = isLine;
  lastInputObject->next = NULL;
}

void freeInputObjects()
{
  INPUT_OBJECT *ptr;
  lastInputObject = inputObject.next;
  while (lastInputObject) {
    ptr = lastInputObject;
    lastInputObject = lastInputObject->next;
    free(ptr);
  }
  lastInputObject = NULL;
}


#define MAX_LINE_LENGTH 128*16384 
#define MAX_FILE_NESTING 10
LINE_LIST *get_beamline(char *madfile, char *use_beamline, double p_central, long echo)
{
  long type=0, i;
  long occurence, iMad;
  static ELEMENT_LIST *eptr, *eptr1, *eptr_sc, *eptr_add, *eptr_del;
  static LINE_LIST *lptr;
  static long n_elems, n_lines;
  FILE *fp_mad[MAX_FILE_NESTING];
  char *s, *t=NULL, *ptr=NULL;
  char occurence_s[8], eptr_name[1024];
  ntuple *nBx, *nBy, *nBz;
  double ftable_length;
  
  log_entry("get_beamline");

  if (!(s=malloc(sizeof(*s)*MAX_LINE_LENGTH)) ||
      !(t=malloc(sizeof(*s)*MAX_LINE_LENGTH)))
    bombElegant("memory allocation failure (get_beamline)", NULL);
      
  if (madfile) {
    char *filename;

#ifdef DEBUG
    fprintf(stdout, "reading from file %s\n", madfile);
    fflush(stdout);
#endif
    if (!(filename=findFileInSearchPath(madfile))) {
      fprintf(stderr, "Unable to find file %s\n", madfile);
      exitElegant(1);
    }
    fp_mad[0] = fopen_e(filename, "r", 0);
    free(filename);
    
    iMad = 0;
    
    elem = tmalloc(sizeof(*elem));
    line = tmalloc(sizeof(*line));
    
    elem->pred = elem->succ = NULL;
    elem->name = NULL;
    eptr      = elem;
    n_elems   = 0;        /* number of physical elements in linked-list */
    line->pred = line->succ = NULL;
    line->name = NULL;
    lptr      = line;        
    n_lines   = 0;        /* number of line definitions in linked-list  */
    
    /* assemble linked-list of simple elements and a separate linked-list
       of fully expanded line definitions */
    while (iMad>=0) {
      while (cfgets(s, MAX_LINE_LENGTH, fp_mad[iMad])) {
        if (echo) {
          fprintf(stdout, "%s\n", s);
          fflush(stdout);
        }
	if (s[0]=='%') {
	  /* rpn command */
	  chop_nl(s);
	  rpn(s+1);
	}
        if (s[0]=='#' && strncmp(s, "#INCLUDE:", strlen("#INCLUDE:"))==0) {
          char *filename;
          if (++iMad==MAX_FILE_NESTING) 
            bombElegant("files nested too deeply", NULL);
          ptr = get_token(s+strlen("#INCLUDE:"));
          if (echo) {
            fprintf(stdout, "reading file %s\n", ptr);
          }
          if (!(filename=findFileInSearchPath(ptr))) {
            fprintf(stderr, "Unable to find file %s\n", ptr);
            exitElegant(1);
          }
          fp_mad[iMad] = fopen_e(filename, "r", 0);
          free(filename);
          continue;
        }
        strcpy_ss(t, s);
        if ((type = tell_type(s, elem))==T_NODEF) {
          if (!is_blank(s))
            fprintf(stdout, "warning: no recognized statement on line: %s\n", t);
          fflush(stdout);
          continue;
        }
#ifdef DEBUG
        fprintf(stdout, "type code = %ld\n", type);
        fflush(stdout);
#endif
        if (type==T_RENAME)
          process_rename_request(s, entity_name, N_TYPES);
        else if (type==T_TITLE) {
          fgets(s, MAX_LINE_LENGTH, fp_mad[iMad]);
          compressString(s, " ");
          if (s[i=strlen(s)-1]==' ')
            s[i] = 0;
        }    
        else if (type==T_USE || type==T_RETURN) 
          break;
        else if (type==T_LINE) {
#ifdef DEBUG
          fprintf(stdout, "current element list is:\n");
          fflush(stdout);
          print_elem_list(stdout, elem);
#endif
          fill_line(line, n_lines, elem, n_elems, s);
          addToInputObjectList((void*)lptr, 1);
          if (strchr(lptr->name, '#')) {
            fprintf(stdout, "The name %s is invalid for a beamline: # is a reserved character.\n", lptr->name);
            exitElegant(1);
          }
          if (check_duplic_elem(&elem, NULL, lptr->name, n_elems)) {
            fprintf(stdout, "line definition %s conflicts with element of same name\n", lptr->name);
            exitElegant(1);
          }
          check_duplic_line(line, lptr->name, n_lines+1, 0);
#ifdef DEBUG 
          fprintf(stdout, "\n****** expanded line %s:\n", lptr->name);
          fflush(stdout);
          print_line(stdout, lptr); 
#endif
          extend_line_list(&lptr);
          n_lines++;
        }
        else {
          if (type==T_ECOPY) {
#ifdef DEBUG
            fprintf(stdout, "copying existing element\n");
            fflush(stdout);
#endif
            strcpy_ss(s, t);
            copy_named_element(eptr, s, elem);
            if (strchr(eptr->name, '#')) {
              fprintf(stdout, "The name %s is invalid for an element: # is a reserved character.\n", eptr->name);
              exitElegant(1);
            }
          }
          else {
            long newType;
            double length;
#ifdef DEBUG
            fprintf(stdout, "creating new element\n");
            fflush(stdout);
#endif
            fill_elem(eptr, s, type, fp_mad[iMad]);
            addToInputObjectList((void*)eptr, 0);
            if (strchr(eptr->name, '#')) {
              fprintf(stdout, "The name %s is invalid for an element: # is a reserved character.\n", eptr->name);
              exitElegant(1);
            }
            length = 0;
            if ((newType=elementTransmutation(eptr->name, eptr->type))!=eptr->type &&
                newType>=0) {
              if (entity_description[eptr->type].flags&HAS_LENGTH) {
                length = ((DRIFT*)eptr->p_elem)->length;
                if (length && !(entity_description[newType].flags&HAS_LENGTH)) {
                  fprintf(stderr, "error: can't transmute %s %s into %s---would change length of beamline\n",
                          entity_name[eptr->type], eptr->name, 
                          entity_name[newType]);
                  exitElegant(1);
                }
              }
              free(eptr->p_elem);
              eptr->p_elem = tmalloc(entity_description[type].structure_size);
              zero_memory(eptr->p_elem, entity_description[type].structure_size);
              eptr->type = newType;
              if (entity_description[newType].flags&HAS_LENGTH)
                ((DRIFT*)eptr->p_elem)->length = length;
            }
            if (check_duplic_line(line, eptr->name, n_lines+1, 1)) {
              fprintf(stdout, "element %s conflicts with line with same name\n", eptr->name);
              exitElegant(1);
            }
            check_duplic_elem(&elem, &eptr, NULL, n_elems);
          }
#ifdef DEBUG
          print_elem(stdout, elem);
#endif
          extend_elem_list(&eptr);
          n_elems++;
        }
      }
      fclose(fp_mad[iMad]);
      iMad--;
    }
    if (n_elems==0 || n_lines==0) {
      fprintf(stdout, "Insufficient (recognizable) data in file.\n");
      fflush(stdout);
      exitElegant(1);
    }

    if (getSCMULTSpecCount()) {
      fill_elem(eptr, getSCMULTName(), T_SCMULT, NULL);
      eptr_sc = eptr;
      check_duplic_elem(&elem, &eptr, NULL, n_elems);
      extend_elem_list(&eptr);
      n_elems++;  	
    }
    
/* since the lists were being extended before it was known that
       the was another object to put in them, must eliminate references
       to the most recently added nodes. 
       */
    (eptr->pred)->succ = NULL;
    (lptr->pred)->succ = NULL;
    eptr = eptr->pred;
    lptr = lptr->pred;
  }
  else {
    s[0] = 0;
    type = T_NODEF;

    if (getAddElemFlag()) {
      /* go to the last elements in linked-list */
      eptr = elem;
      while (eptr->succ) 
        eptr = eptr->succ;
      /* extend the list for accommodating new element */
      extend_elem_list(&eptr);
      /* add new element to linked-list */
      s = getElemDefinition();
      if ((type = tell_type(s, elem))==T_NODEF) {
        if (!is_blank(s))
          fprintf(stdout, "warning: no recognized statement on line: %s\n", s);
        fflush(stdout);
      }
      fill_elem(eptr, s, type, NULL);
      if (strchr(eptr->name, '#')) {
        fprintf(stdout, "The name %s is invalid for an element: # is a reserved character.\n", eptr->name);
        exitElegant(1);
      }
      eptr_add = eptr;
      if (check_duplic_line(line, eptr->name, n_lines+1, 1)) {
        fprintf(stdout, "element %s conflicts with line with same name\n", eptr->name);
        exitElegant(1);
      }
      check_duplic_elem(&elem, &eptr, NULL, n_elems);
      n_elems++;  	
    }

    if (getDelElemFlag()==1) {
      /* go to the last elements in linked-list */
      eptr = elem;
      while (eptr->succ) 
        eptr = eptr->succ;
      /* extend the list for accommodating new element */
      extend_elem_list(&eptr);
      /* add new element to linked-list */
      s = getElemDefinition1();
      if ((type = tell_type(s, elem))==T_NODEF) {
        if (!is_blank(s))
          fprintf(stdout, "warning: no recognized statement on line: %s\n", s);
        fflush(stdout);
      }
      fill_elem(eptr, s, type, NULL);
      if (strchr(eptr->name, '#')) {
        fprintf(stdout, "The name %s is invalid for an element: # is a reserved character.\n", eptr->name);
        exitElegant(1);
      }
      /* This is mis spelled. Should be eptr_replace. */
      eptr_del = eptr;
      if (check_duplic_line(line, eptr->name, n_lines+1, 1)) {
        fprintf(stdout, "element %s conflicts with line with same name\n", eptr->name);
        exitElegant(1);
      }
      check_duplic_elem(&elem, &eptr, NULL, n_elems);
      n_elems++;  	
    }
  }
  

  if (type!=T_USE && use_beamline==NULL) {
    if (n_lines==0)
      bombElegant("no beam-line defined\n", NULL);
    fprintf(stdout, "no USE statement--will use line %s\n", lptr->name);
    fflush(stdout);
  }
  else {
    if (!use_beamline) {
      if ((ptr=get_token(s))==NULL) 
        bombElegant("no line named in USE statement", NULL);
    }
    else
      ptr = str_toupper(use_beamline);
    lptr = line;
    while (lptr) {
      if (strcmp(lptr->name, ptr)==0) 
        break;
      lptr = lptr->succ;
    }
    if (lptr==NULL) {
      fprintf(stdout, "no definition of beam-line %s\n", ptr);
      fflush(stdout);
      exitElegant(1);
    }
    if (!use_beamline) {
      free(ptr);
    }
  }

  /* these really aren't necessary, since I clear memory upon allocation */
  lptr->elem_recirc = lptr->elem_twiss = lptr->elast = NULL;
  lptr->twiss0 = NULL;
  lptr->matrix = NULL;

  if (getSCMULTSpecCount()) {
  	long skip = 0;
  	long nelem = 0;
  	eptr = &(lptr->elem);
  	while (eptr) {
  		if (eptr->type == T_SCMULT) {					/* the code allow user put scmult explicitly */
  			eptr = eptr->succ;
  			continue;
  		}
  		if (insertSCMULT(eptr->name, eptr->type, &skip)) {
  			add_element(eptr, eptr_sc); 
  			eptr = eptr->succ;							/* move pointer to new added element */
  			nelem++;
   		}
  		if (eptr->succ==NULL && skip!=0) {				/* add space charge element to the end of line */
  			add_element(eptr, eptr_sc);
  			eptr = eptr->succ;							/* this is very impotant to get off the loop */
  			nelem++;
 		}   			
  		eptr = eptr->succ; 
  	}
  	lptr->n_elems += nelem;
  } 

  if (getAddElemFlag()) {
    long skip = 0;
    long nelem = 0;
    eptr = &(lptr->elem);
    while (eptr) {
      /* The end position will have been set in a previous call to get_beamline(), prior to 
         definition of insertions */
      if (insertElem(eptr->name, eptr->type, &skip, eptr->occurence, eptr->end_pos)) {
        add_element(eptr, eptr_add); 
        eptr = eptr->succ;		/* move pointer to new added element */
        nelem++;
      }
      if (eptr->succ==NULL && getAddEndFlag()) {	/* add element to the end of line if request */
        add_element(eptr, eptr_add);
        eptr = eptr->succ;				/* this is very impotant to get off the loop */
        nelem++;
      }   			
      eptr = eptr->succ; 
    }
    lptr->n_elems += nelem;
  } 

  if (getDelElemFlag()) {
    long skip = 0;
    long flag;
    eptr = &(lptr->elem);
    while (eptr) {
      flag = replaceElem(eptr->name, eptr->type, &skip, eptr->occurence);
      if (flag == 1) {
        if (!eptr->pred)
          bombElegant("Can not replace the first element in beamline", NULL);
        eptr = replace_element(eptr, eptr_del);
      }
      if (flag == -1) {
        if (!eptr->pred)
          bombElegant("Can not remove the first element in beamline", NULL);
        eptr = rm_element(eptr);
        lptr->n_elems--;
      } 
      eptr = eptr->succ; 
    }
  } 

  /* go through and give occurence numbers to each element */
  eptr = &(lptr->elem);
  while (eptr) {
    eptr->occurence = 0;
    lptr->elast = eptr;
    eptr = eptr->succ;
  }

  eptr = &(lptr->elem);
  lptr->flags = 0;
  while (eptr) {
    eptr->Pref_input = eptr->Pref_output = p_central;
    if (eptr->occurence==0) {
      /* this is the first occurence of this element--go through and find any others */
      if (eptr->type==T_FTABLE)  {
        initializeFTable((FTABLE*)eptr->p_elem);
        nBx = ((FTABLE*)eptr->p_elem)->Bx;
        nBy = ((FTABLE*)eptr->p_elem)->By;
        nBz = ((FTABLE*)eptr->p_elem)->Bz;
        ftable_length = ((FTABLE*)eptr->p_elem)->length;
      }        

      eptr->occurence = occurence = 1;
      eptr1 = eptr->succ;
      while (eptr1) {
        if (strcmp(eptr->name, eptr1->name)==0) {
          if (eptr1->type==T_FTABLE)  {
            ((FTABLE*)eptr1->p_elem)->initialized=1;
            ((FTABLE*)eptr1->p_elem)->length=ftable_length;
            ((FTABLE*)eptr1->p_elem)->Bx = nBx;
            ((FTABLE*)eptr1->p_elem)->By = nBy;
            ((FTABLE*)eptr1->p_elem)->Bz = nBz;
          }
          eptr1->occurence = ++occurence;
        }
        eptr1 = eptr1->succ;
      }
    }
    if (eptr->type==T_SREFFECTS)
      lptr->flags |= BEAMLINE_TWISS_WANTED;
    eptr = eptr->succ;
  }
  /* create a hash table with the size of 2^12, it can grow automatically if necessary */
  if (!load_hash)
     load_hash = hcreate(12);  
  eptr = &(lptr->elem);
  while (eptr) {
    /* use "eptr->name+eptr->occurence" as the key, and eptr's address as the value for hash table*/
      sprintf (occurence_s, "#%ld", eptr->occurence);
      strcpy (eptr_name, eptr->name);
      strcat (eptr_name, occurence_s);
      hadd (load_hash, eptr_name, strlen(eptr_name), (void*)eptr);
      eptr = eptr->succ;
  }
  compute_end_positions(lptr);
  free(s);
  free(t);
  
  return(lptr);
}

static long negativeLengthWarningsLeft = 100;

double compute_end_positions(LINE_LIST *lptr) 
{
    double z, l, theta, z_recirc;
    static ELEMENT_LIST *eptr;
    long i_elem, recircPresent;
    
    /* use length data to establish z coordinates at end of each element */
    /* also check for duplicate recirculation elements and set occurence numbers to 0 */
    eptr = &(lptr->elem);
    z = z_recirc = 0;
    theta = 0;
    i_elem = 0;
    recircPresent = 0;
    do {
        if (!(entity_description[eptr->type].flags&HAS_LENGTH))
            l = 0;
        else
            l = (*((double*)eptr->p_elem));
        if (eptr->type==T_SBEN || eptr->type==T_RBEN)
            theta += ((BEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_KSBEND)
            theta += ((KSBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_NIBEND)
            theta += ((NIBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_NISEPT)
            theta += ((NISEPT*)eptr->p_elem)->angle;
        else if (eptr->type==T_CSBEND)
            theta += ((CSBEND*)eptr->p_elem)->angle;
        else if (eptr->type==T_EMATRIX)
            theta += ((EMATRIX*)eptr->p_elem)->angle;
        else if (eptr->type==T_FTABLE && ((FTABLE*)eptr->p_elem)->angle) {
            theta += ((FTABLE*)eptr->p_elem)->angle;
            l = ((FTABLE*)eptr->p_elem)->l0/2./sin(((FTABLE*)eptr->p_elem)->angle/2.)*((FTABLE*)eptr->p_elem)->angle;
        } 
        else if (eptr->type==T_RECIRC) {
            if (recircPresent)
                bombElegant("multiple recirculation (RECIRC) elements in beamline--this doesn't make sense", NULL);
            lptr->elem_recirc = eptr;
            lptr->i_recirc = i_elem;
            z_recirc = z;
            recircPresent = 1;
            }
        if (l<0 && negativeLengthWarningsLeft>0) {
            fprintf(stdout, "warning(1): element %s has negative length = %e\n", eptr->name, l);
            if (--negativeLengthWarningsLeft==0)
              fprintf(stdout, "Further negative length warnings will be suppressed.\n");
            fflush(stdout);
          }
        eptr->end_pos = z + l;
        eptr->end_theta = theta ;
        z = eptr->end_pos;
        i_elem++;
        } while ((eptr=eptr->succ));

    lptr->revolution_length = z - z_recirc;
    return lptr->revolution_length;
    }

void show_elem(ELEMENT_LIST *eptr, long type)
{
    long j;
    char *ptr;
    PARAMETER *parameter;

    log_entry("show_elem");

    parameter = entity_description[type].parameter;
    fprintf(stdout,  "%s %s at z=%em:\n", 
        entity_name[type], eptr->name, eptr->end_pos);
    fflush(stdout);
    for (j=0; j<entity_description[type].n_params; j++) {
        switch (parameter[j].type) {
            case IS_DOUBLE:
                fprintf(stdout,  "    %s = %.16e with offset %ld\n", 
                    parameter[j].name, 
                    *(double*)(eptr->p_elem+parameter[j].offset),
                    parameter[j].offset);
                fflush(stdout);
                break;
            case IS_LONG:
                fprintf(stdout,  "    %s = %ld with offset %ld\n", 
                    parameter[j].name, 
                    *(long *)(eptr->p_elem+parameter[j].offset),
                    parameter[j].offset);
                fflush(stdout);
                break;
            case IS_STRING:
                if ((ptr = *(char**)(eptr->p_elem+parameter[j].offset))) {
                    fprintf(stdout,  "    %s = %s\n", parameter[j].name, ptr);
                    fflush(stdout);
                  }
                else {
                    fprintf(stdout,  "    %s = %s\n", parameter[j].name, ptr);
                    fflush(stdout);
                  }
                break;
            }
        }
    log_exit("show_elem");
    }

void free_elements(ELEMENT_LIST *elemlist)
{
    ELEMENT_LIST *eptr;

    log_entry("free_elements");

    if (elemlist) {
        eptr= elemlist;
        }
    else  {
        eptr = elem;
        elem = NULL;
#ifdef DEBUG
        fprintf(stdout, "freeing elements in main list\n");
        fflush(stdout);
#endif
        }
    while (eptr) {
#ifdef DEBUG
        fprintf(stdout, "freeing memory for element %s of type %s\n", 
                eptr->name?eptr->name:"NULL",
                (eptr->type>=0 && eptr->type<N_TYPES)?entity_name[eptr->type]:"NULL");
        fflush(stdout);
#endif
        if (eptr->type==T_MATR) {
          eptr = eptr->succ;
          continue;
        }
        if (eptr->type==T_WATCH) {
	  WATCH *wptr;
	  if ((wptr = (WATCH*)eptr->p_elem)) {
	    if (wptr->initialized && !SDDS_Terminate(&wptr->SDDS_table)) {
	      SDDS_SetError("Problem terminate watch-point SDDS file (free_elements)");
	      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	    }
	  }
	}
#ifdef DEBUG
        fprintf(stdout, "pointers: p_elem = %x   name = %x   matrix = %x\n",
            eptr->p_elem, eptr->name, eptr->matrix);
        fflush(stdout);
        fprintf(stdout, "          pred = %x     succ = %x  \n",
            eptr->pred, eptr->succ);
        fflush(stdout);
#endif
        tfree(eptr->p_elem);
        eptr->p_elem = NULL;
        tfree(eptr->p_elem0);
        eptr->p_elem0 = NULL;
        tfree(eptr->name);
        eptr->name = NULL;
        tfree(eptr->definition_text);
        eptr->definition_text = NULL;
        if (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix) {
            free_matrices(eptr->matrix);
            free(eptr->matrix);
            eptr->matrix = NULL;
            }
        if (eptr->accumMatrix) {
          free_matrices(eptr->accumMatrix);
          free(eptr->accumMatrix);
          eptr->accumMatrix = NULL;
        }
        if (eptr->succ) {
            eptr = eptr->succ;
            free(eptr->pred);
            eptr->pred = NULL;
            }
        else {
            free(eptr);
            break;
            }
        }
    log_exit("free_elements");
    }

void free_beamlines(LINE_LIST *beamline)
{
    LINE_LIST *lptr;

    log_entry("free_beamlines");

    if (beamline) {
        lptr = beamline;
        }
    else {
        lptr = line;
        line = NULL;
#ifdef DEBUG
        fprintf(stdout, "freeing main beamline list\n", NULL);
        fflush(stdout);
#endif
        }
    while (lptr) {
#ifdef DEBUG
        fprintf(stdout, "*************************\nfreeing memory for beamline %s with %ld elements\n",
            lptr->name?lptr->name:"NULL", lptr->n_elems);
        fflush(stdout);
        fprintf(stdout, "pointers:   name = %x    succ = %x   pred = %x\n",
            lptr->name, lptr->succ, lptr->pred);
        fflush(stdout);
#endif
        if (lptr->definition) {
            tfree(lptr->definition);
            lptr->definition = NULL;
            }
        if (lptr->name) {
            tfree(lptr->name);
            lptr->name = NULL;
            }
        if (lptr->n_elems) {
            free_elements((lptr->elem).succ);
            /* should free name etc. for lptr->elem also */
            lptr->n_elems = 0;
            lptr->flags = 0;
            }
        if (lptr->succ) {
            lptr = lptr->succ;
            tfree(lptr->pred);
            lptr->pred = NULL;
            }
        else {
            tfree(lptr);
            break;
            }
        } 
    if (load_hash) { 
       hdestroy(load_hash);                         /* destroy hash table */  
       load_hash = NULL;
    }      
    log_exit("free_beamlines");
}

void delete_matrix_data(LINE_LIST *beamline)
{
    LINE_LIST *lptr;
    ELEMENT_LIST *eptr;

    log_entry("delete_matrix_data");

    if (beamline) {
        lptr = beamline;
        }
    else {
        lptr = line;
        }
    while (lptr) {
        if (lptr->n_elems) {
            eptr = &(lptr->elem);
            while (eptr) {
                if (entity_description[eptr->type].flags&HAS_MATRIX && eptr->matrix) {
                    free_matrices(eptr->matrix);
                    tfree(eptr->matrix);
                    eptr->matrix = NULL;
                    }
                if (eptr->accumMatrix) {
                  free_matrices(eptr->accumMatrix);
                  tfree(eptr->accumMatrix);
                  eptr->accumMatrix = NULL;
                }
                eptr = eptr->succ;
                }
            lptr->n_elems = 0;
            lptr->flags = 0;
            }
        lptr = lptr->succ;
        }
    log_exit("delete_matrix_data");
    }


/* routine: do_save_lattice()
 * purpose: save the element and beamline definitions to a file
 * 
 * Michael Borland, 1991
 */
#include "save_lattice.h"

void do_save_lattice(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
  FILE *fp;
  ELEMENT_LIST *eptr;
  LINE_LIST *lptr;
  long j;
  double dvalue;
  long lvalue;
  char *ptr;
  PARAMETER *parameter;
  char s[16384], t[1024], name[1024];
  INPUT_OBJECT *object;
  
  log_entry("do_save_lattice");

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&save_lattice, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &save_lattice);

  /* check for valid data */
  if (filename==NULL)
    bombElegant("no filename given to save lattice to", NULL);
  if (str_in(filename, "%s"))
    filename = compose_filename(filename, run->rootname);
  fp = fopen_e(filename, "w", FOPEN_INFORM_OF_OPEN);

  if (!output_seq) {
    object = &inputObject;
    do {
      if (!object->isLine) {
        eptr = (ELEMENT_LIST*)(object->ptr);
        parameter = entity_description[eptr->type].parameter;
        if (strpbrk(eptr->name, ":.,/-_+abcdefghijklmnopqrstuvwyxz "))
          sprintf(name, "\"%s\"", eptr->name);
        else
          strcpy_ss(name, eptr->name);
        sprintf(s, "%s: %s,", name, entity_name[eptr->type]);
        if (eptr->group && strlen(eptr->group)) {
          sprintf(t, "GROUP=\"%s\",", eptr->group);
          strcat(s, t);
        }
        for (j=0; j<entity_description[eptr->type].n_params; j++) {
          switch (parameter[j].type) {
          case IS_DOUBLE:
            dvalue = *(double*)(eptr->p_elem+parameter[j].offset);
            /*
              if ((parameter[j].flags&PARAM_DIVISION_RELATED) && 
              eptr->divisions>1) {
              fprintf(stderr, "Multiplying %s by %ld\n",
              parameter[j].name, eptr->divisions);
              dvalue *= eptr->divisions;
              }
            */
            if (!suppress_defaults || dvalue!=parameter[j].number) {
              /* value is not the default, so add to output */
              sprintf(t, "%s=%.16g", parameter[j].name, dvalue);
              strcat(s, t);
              if (j!=entity_description[eptr->type].n_params-1)
                strcat(s, ",");
            }
            break;
          case IS_LONG:
            lvalue = *(long *)(eptr->p_elem+parameter[j].offset);
            if (!suppress_defaults || lvalue!=parameter[j].integer) {
              /* value is not the default, so add to output */
              sprintf(t, "%s=%ld", parameter[j].name, lvalue);
              strcat(s, t);
              if (j!=entity_description[eptr->type].n_params-1)
                strcat(s, ",");
            }
            break;
          case IS_STRING:
            ptr = *(char**)(eptr->p_elem+parameter[j].offset);
            if (ptr &&
                (!suppress_defaults || !parameter[j].string || strcmp(ptr, parameter[j].string)!=0)) {
              sprintf(t, "%s=\"%s\"", parameter[j].name, ptr);
              strcat(s, t);
              if (j!=entity_description[eptr->type].n_params-1)
                strcat(s, ",");
            }
            break;
          }
        }
        if (s[j=strlen(s)-1]==',')
          s[j] = 0;
        print_with_continuation(fp, s, 79);
        eptr = eptr->succ;
      } else {
        lptr = (LINE_LIST*)(object->ptr);
        print_with_continuation(fp, lptr->definition, 79);
      }
    } while ((object=object->next));
    
  } else {
    /* first write element definition */
    long type;
    long nline=1, nelem=0;
    for (type=1; type<N_TYPES; type++) {
      eptr = &(beamline->elem);
      while (eptr) {
        if ((eptr->occurence == 1) && (eptr->type == type)) {
          parameter = entity_description[eptr->type].parameter;
          if (strpbrk(eptr->name, ":.,/-_+abcdefghijklmnopqrstuvwyxz "))
            sprintf(name, "\"%s\"", eptr->name);
          else
            strcpy_ss(name, eptr->name);
          sprintf(s, "%s: %s,", name, entity_name[eptr->type]);
          if (eptr->group && strlen(eptr->group)) {
            sprintf(t, "GROUP=\"%s\",", eptr->group);
            strcat(s, t);
          }
          for (j=0; j<entity_description[eptr->type].n_params; j++) {
            switch (parameter[j].type) {
            case IS_DOUBLE:
              dvalue = *(double*)(eptr->p_elem+parameter[j].offset);
              /*
                if ((parameter[j].flags&PARAM_DIVISION_RELATED) && 
                eptr->divisions>1) {
                fprintf(stderr, "Multiplying %s by %ld\n",
                parameter[j].name, eptr->divisions);
                dvalue *= eptr->divisions;
                }
              */
              if (!suppress_defaults || dvalue!=parameter[j].number) {
                /* value is not the default, so add to output */
                sprintf(t, "%s=%.16g", parameter[j].name, dvalue);
                strcat(s, t);
                if (j!=entity_description[eptr->type].n_params-1)
                  strcat(s, ",");
              }
              break;
            case IS_LONG:
              lvalue = *(long *)(eptr->p_elem+parameter[j].offset);
              if (!suppress_defaults || lvalue!=parameter[j].integer) {
                /* value is not the default, so add to output */
                sprintf(t, "%s=%ld", parameter[j].name, lvalue);
                strcat(s, t);
                if (j!=entity_description[eptr->type].n_params-1)
                  strcat(s, ",");
              }
              break;
            case IS_STRING:
              ptr = *(char**)(eptr->p_elem+parameter[j].offset);
              if (ptr &&
                  (!suppress_defaults || !parameter[j].string || strcmp(ptr, parameter[j].string)!=0)) {
                sprintf(t, "%s=\"%s\"", parameter[j].name, ptr);
                strcat(s, t);
                if (j!=entity_description[eptr->type].n_params-1)
                  strcat(s, ",");
              }
              break;
            }
          }
          if (s[j=strlen(s)-1]==',')
            s[j] = 0;
          print_with_continuation(fp, s, 79);
        }
        eptr = eptr->succ;
      }
    }
    /* Write beamline sequence now, each line has 40 elements limitation */
    eptr = &(beamline->elem);
    sprintf(s, "L%04ld: LINE = (", nline);
    while (eptr) {
      nelem++;
      if (strpbrk(eptr->name, ":.,/-_+abcdefghijklmnopqrstuvwyxz "))
        sprintf(name, "\"%s\"", eptr->name);
      else
        strcpy_ss(name, eptr->name);
      strcat(s, name);
      strcat(s, ",");

      eptr = eptr->succ;
      if (nelem == 40) {
        nline++;
        nelem=0;
        if (s[j=strlen(s)-1]==',')
          s[j] = 0;
        strcat(s, ")");
        print_with_continuation(fp, s, 79);
        sprintf(s, "L%04ld: LINE = (", nline);
      }
    }
    if (nelem) {
      nline++;
      if (s[j=strlen(s)-1]==',')
        s[j] = 0;
      strcat(s, ")");
      print_with_continuation(fp, s, 79);
    }

    sprintf(s, "%s: LINE = (",  beamline->name);
    for (j=1; j<nline; j++) {
      sprintf(t, " L%04ld,",j);
      strcat(s, t);
    }
    if (s[j=strlen(s)-1]==',')
      s[j] = 0;
    strcat(s, ")");
    print_with_continuation(fp, s, 79);
  }
  
  if (beamline && beamline->name)
    fprintf(fp, "USE,%s\n", beamline->name);

  fprintf(fp, "RETURN\n");
  fclose(fp);

  log_exit("do_save_lattice");
}

void print_with_continuation(FILE *fp, char *s, long endcol)
{
  char c, *ptr;
  long l, isContin;

  isContin = 0;
  while ((l=strlen(s))) {
    if (isContin)
      fputc(' ', fp);
    if (l>endcol) {
      ptr = s+endcol-2;
      while (ptr!=s && *ptr!=',')
        ptr--;
      if (ptr==s)
        c = *(ptr = s + endcol - 1);
      else {
        ptr++;
        c = *ptr;
      }
      *ptr = 0;
      fputs(s, fp);
      fputs("&\n", fp);
      isContin = 1;
      s = ptr;
      *ptr = c;
    }
    else {
      fputs(s, fp);
      fputc('\n', fp);
      log_exit("print_with_continuation");
      return;
    }
  }
}

/* Change defined parameter values in the reference list elem.
 * This routine allows changing a number of parameters for a number of differently-named elements
 */
void change_defined_parameter_values(char **elem_name, long *param_number, long *type, 
                                     double *value, long n_elems)
{
  ELEMENT_LIST *eptr;
  char *p_elem;
  long i_elem, elem_type, data_type, param;
  double dValue;
  
  log_entry("change_defined_parameter_values");

  for (i_elem=0; i_elem<n_elems; i_elem++) {
    eptr = NULL;
    elem_type = type[i_elem];
    param     = param_number[i_elem];
    data_type = entity_description[elem_type].parameter[param].type;
    while (find_element(elem_name[i_elem], &eptr, elem)) {
      p_elem = eptr->p_elem;
      switch (data_type) {
      case IS_DOUBLE:
        dValue = value[i_elem];
        if ((entity_description[elem_type].parameter[param].flags&PARAM_DIVISION_RELATED) &&
            eptr->divisions)
          dValue *= eptr->divisions;
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) = dValue;
#if DEBUG
        fprintf(stdout, "   changing parameter %s of %s #%ld to %e\n",
                entity_description[elem_type].parameter[param].name,
                eptr->name, eptr->occurence,
                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
#endif
        break;
      case IS_LONG:
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
          nearestInteger(value[i_elem]);
#if DEBUG
        fprintf(stdout, "   changing parameter %s of %s #%ld to %ld\n",
                entity_description[elem_type].parameter[param].name,
                eptr->name, eptr->occurence,
                *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
        fflush(stdout);
#endif
        break;
      case IS_STRING:
      default:
        bombElegant("unknown/invalid variable quantity", NULL);
        exitElegant(1);
      }
    }
  }
  log_exit("change_defined_parameter_values");
}


/* Change defined parameter values in the reference list elem.
 * This routine allows changing a single parameter for a single element name.
 */
void change_defined_parameter_divopt(char *elem_name, long param, long elem_type, 
                                      double value, char *valueString, unsigned long mode, 
                                      long checkDiv)
{
  ELEMENT_LIST *eptr;
  char *p_elem;
  long data_type;

  log_entry("change_defined_parameter");

  eptr = NULL;
  data_type = entity_description[elem_type].parameter[param].type;
  if (mode&LOAD_FLAG_IGNORE)
    return;
  while (find_element(elem_name, &eptr, elem)) {
    p_elem = eptr->p_elem;
    switch (data_type) {
    case IS_DOUBLE:
      if (valueString) {
        if (!sscanf(valueString, "%lf", &value)) {
          fprintf(stdout, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
          fflush(stdout);
          exitElegant(1);
        }
      }
      if (checkDiv && eptr->divisions>1 &&
          (entity_description[elem_type].parameter[param].flags&PARAM_DIVISION_RELATED))
        value /= eptr->divisions;
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition (mode %s) %s.%s from %e to ", 
                (mode&LOAD_FLAG_ABSOLUTE)?"absolute":
                ((mode&LOAD_FLAG_DIFFERENTIAL)?"differential":
                 (mode&LOAD_FLAG_FRACTIONAL)?"fractional":"unknown"),
                elem_name, entity_description[elem_type].parameter[param].name,
                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
      fflush(stdout);
      if (mode&LOAD_FLAG_ABSOLUTE) {
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) =  value;
      }
      else if (mode&LOAD_FLAG_DIFFERENTIAL) {
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) += value;
      }
      else if (mode&LOAD_FLAG_FRACTIONAL) {
        *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= 1+value;
      }
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "%e\n", 
                *((double*)(p_elem+entity_description[elem_type].parameter[param].offset)));
      fflush(stdout);
      break;
    case IS_LONG:
      if (valueString) {
        if (!sscanf(valueString, "%lf", &value)) {
          fprintf(stdout, "Error (change_defined_parameter): unable to scan double from \"%s\"\n", valueString);
          fflush(stdout);
          exitElegant(1);
        }
      }
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition (mode %s) %s.%s from %ld to ",
                (mode&LOAD_FLAG_ABSOLUTE)?"absolute":
                ((mode&LOAD_FLAG_DIFFERENTIAL)?"differential":
                 (mode&LOAD_FLAG_FRACTIONAL)?"fractional":"unknown"),
                elem_name, entity_description[elem_type].parameter[param].name,
                *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
      fflush(stdout);
      if (mode&LOAD_FLAG_ABSOLUTE) {
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) = 
          nearestInteger(value);
      }
      else if (mode&LOAD_FLAG_DIFFERENTIAL) {
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) += 
          nearestInteger(value);
      }
      else if (mode&LOAD_FLAG_FRACTIONAL) {
        *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)) *= 1+value;
      }
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "%ld\n",
                *((long*)(p_elem+entity_description[elem_type].parameter[param].offset)));
      fflush(stdout);
      break;
    case IS_STRING:
      if (mode&LOAD_FLAG_VERBOSE)
        fprintf(stdout, "Changing definition %s.%s from %s to %s\n",
                elem_name, entity_description[elem_type].parameter[param].name,
                *((char**)(p_elem+entity_description[elem_type].parameter[param].offset)),
                valueString);
      fflush(stdout);
      if (!SDDS_CopyString(((char**)(p_elem+entity_description[elem_type].parameter[param].offset)), 
                           valueString)) {
        fprintf(stdout, "Error (change_defined_parameter): unable to copy string parameter value\n");
        fflush(stdout);
        exitElegant(1);
      }
      break;
    default:
      bombElegant("unknown/invalid variable quantity", NULL);
      exitElegant(1);
    }
  }
  log_exit("change_defined_parameter");
}

void change_defined_parameter(char *elem_name, long param, long elem_type, 
                              double value, char *valueString, unsigned long mode)
{
  change_defined_parameter_divopt(elem_name, param, elem_type, value, valueString, mode, 0);
}


void process_rename_request(char *s, char **name, long n_names)
{
  long i;
  char *old, *new, *ptr;

  log_entry("process_rename_request");
  if (!(ptr=strchr(s, '='))) 
    bombElegant("invalid syntax for RENAME", NULL);
  *ptr++ = 0;
  old = s;
  str_toupper(trim_spaces(old = s));
  while (*old==',' || *old==' ')
    old++;
  str_toupper(trim_spaces(new = ptr));
  if (match_string(new, name, n_names, EXACT_MATCH)>=0) {
    fprintf(stdout, "error: can't rename to name %s--already exists\n", new);
    fflush(stdout);
    exitElegant(1);
  }
  if ((i=match_string(old, name, n_names, EXACT_MATCH))<0) {
    fprintf(stdout, "error: can't rename %s to %s--%s not recognized\n", old, new, old);
    fflush(stdout);
    exitElegant(1);
  }
  fprintf(stdout, "warning: element %s now known as %s\n", old, new);
  fflush(stdout);
  cp_str(name+i, new);
  log_exit("process_rename_request");
}


long nearestInteger(double value)
{
  if (value<0)
    return -1*((long)(-value+0.5));
  return (long)(value+0.5);
}

/* add element "elem1" after "elem0" */
void add_element(ELEMENT_LIST *elem0, ELEMENT_LIST *elem1) 
{
  ELEMENT_LIST *eptr;
  eptr = tmalloc(sizeof(*eptr));
  copy_element(eptr, elem1, 0, 0, 0);

  eptr->pred = elem0;
  eptr->succ = elem0->succ;
  if (elem0->succ)
    (elem0->succ)->pred = eptr;
  elem0->succ = eptr;
}

/* remove element "elem" from the list */
ELEMENT_LIST *rm_element(ELEMENT_LIST *elem) 		
{
  (elem->pred)->succ = elem->succ;
  if (elem->succ)
    (elem->succ)->pred = elem->pred;
  return (elem->pred);
}

/* replace element "elem0" with "elem1" */
ELEMENT_LIST *replace_element(ELEMENT_LIST *elem0, ELEMENT_LIST *elem1) 
{
  ELEMENT_LIST *eptr;
  eptr = tmalloc(sizeof(*eptr));
  copy_element(eptr, elem1, 0, 0, 0);

  (elem0->pred)->succ = eptr;
  if (elem0->succ)
    (elem0->succ)->pred = eptr;
  eptr->pred = elem0->pred;
  eptr->succ = elem0->succ;
  return (eptr);
}
/* This is called at beginning to avoid multiple calls for same element at different locations */
void initializeFTable(FTABLE *ftable)
{
  long i;
  
  ftable->Bx = readbookn(ftable->inputFile, 1);
  ftable->By = readbookn(ftable->inputFile, 2);
  ftable->Bz = readbookn(ftable->inputFile, 3);

  if ((ftable->Bx->nD !=3)||(ftable->By->nD !=3)||(ftable->Bz->nD !=3))
    bombElegant("ND must be 3 for field table %s.", ftable->inputFile);
  ftable->length = ftable->Bz->xmax[2] - ftable->Bz->xmin[2];
  if (fabs((ftable->l0 + ftable->l1 + ftable->l2)-ftable->length)>1e-12)
    bombElegant("L+L1+L2 != field length in file %s.", ftable->inputFile);

  for (i=0; i<ftable->Bz->nD; i++) {
    ftable->Bx->xmin[i] -= ftable->Bx->dx[i]/2;
    ftable->By->xmin[i] -= ftable->By->dx[i]/2;
    ftable->Bz->xmin[i] -= ftable->Bz->dx[i]/2;
    ftable->Bx->xmax[i] += ftable->Bx->dx[i]/2;
    ftable->By->xmax[i] += ftable->By->dx[i]/2;
    ftable->Bz->xmax[i] += ftable->Bz->dx[i]/2;
  }
  ftable->initialized = 1;
  return;
}

