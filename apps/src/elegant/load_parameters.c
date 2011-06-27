/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: load_parameters.c
 * contents: setup_load_parameters(), do_load_parameters()
 *
 * Michael Borland, 1993
 */
#include "mdb.h"
#include "track.h"
#include "load_parameters.h"
#include "SDDS.h"
#include "match_string.h"

#define DEBUG 0

/* structure to store and manage a load_parameters request */
typedef struct {
    SDDS_TABLE table;
    char *filename;
    unsigned long flags;
#define COMMAND_FLAG_CHANGE_DEFINITIONS 0x0001UL
#define COMMAND_FLAG_IGNORE             0x0002UL
#define COMMAND_FLAG_IGNORE_OCCURENCE   0x0004UL
#define COMMAND_FLAG_USE_FIRST          0x0008UL
    char *includeNamePattern, *includeItemPattern, *includeTypePattern;
    char *excludeNamePattern, *excludeItemPattern, *excludeTypePattern;
    long last_code;          /* return code from SDDS_ReadTable */
    short string_data;       /* if non-zero, indicates data stored as strings */   
    double *starting_value;  /* only for numerical data */
    void **reset_address;
    short *value_type;
    ELEMENT_LIST **element;
    long *element_flags;
    long values;
    } LOAD_PARAMETERS;

/* variables to store and manage load_parameters requests */
static LOAD_PARAMETERS *load_request = NULL;
static long load_requests = 0;
static long load_parameters_setup = 0;

/* names of the SDDS columns that will be used */
static char *Element_ColumnName = "ElementName";
static char *Occurence_ColumnName = "ElementOccurence";
static char *Parameter_ColumnName = "ElementParameter";
static char *Value_ColumnName = "ParameterValue";
static char *ValueString_ColumnName = "ParameterValueString";
static char *Mode_ColumnName = "ParameterMode";
static char *ElementType_ColumnName = "ElementType";

/* the order here must match the order of the #define's in track.h */
#define LOAD_MODES             4
static char *load_mode[LOAD_MODES] = {
    "absolute", "differential", "ignore", "fractional"
    } ;

long setup_load_parameters_for_file(char *filename, RUN *run, LINE_LIST *beamline);

static long missingElementWarningsLeft = 100, missingParameterWarningsLeft = 100;

long setup_load_parameters(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
  long i=0;
  
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&load_parameters, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &load_parameters);

  if (force_occurence_data && use_first)
    bombElegant("Error: force_occurence_data=1 and use_first=1 not meaningful\n", NULL);
  
  if (filename_list && strlen(filename_list)) {
    fprintf(stdout, "Using list of filenames for parameter loading\n");
  }
  else {
    if (!filename && !clear_settings)
      bombElegant("filename or filename_list must be given for load_parameters unless you are clearing settings", NULL);
    fprintf(stdout, "Using single filename for parameter loading\n");
  }

  if (clear_settings && load_requests) {
    for (i=0; i<load_requests; i++) {
      if (!SDDS_Terminate(&load_request[i].table))
        bombElegant("problem terminating load_parameters table", NULL);
    }
    load_requests = 0;
  }
  if (clear_settings)
    return 0;
  load_parameters_setup = 1;

  if (filename_list) {
    char *filename0;
    while ((filename0 = get_token(filename_list)) != NULL)
      i = setup_load_parameters_for_file(filename0, run, beamline);
    return i;
  }
  else
    return setup_load_parameters_for_file(filename, run, beamline);
}

long setup_load_parameters_for_file(char *filename, RUN *run, LINE_LIST *beamline)
{
  long index;
  
  load_request = trealloc(load_request, sizeof(*load_request)*(load_requests+1));
  load_request[load_requests].flags = change_defined_values?COMMAND_FLAG_CHANGE_DEFINITIONS:0;
  load_request[load_requests].filename = compose_filename(filename, run->rootname);
  if (change_defined_values && !force_occurence_data)
    load_request[load_requests].flags |= COMMAND_FLAG_IGNORE_OCCURENCE;
  if (use_first)
    load_request[load_requests].flags |= COMMAND_FLAG_USE_FIRST;
  load_request[load_requests].includeNamePattern = include_name_pattern;
  load_request[load_requests].includeItemPattern = include_item_pattern;
  load_request[load_requests].includeTypePattern = include_type_pattern;
  load_request[load_requests].excludeNamePattern = exclude_name_pattern;
  load_request[load_requests].excludeItemPattern = exclude_item_pattern;
  load_request[load_requests].excludeTypePattern = exclude_type_pattern;

#ifdef USE_MPE /* use the MPE library */    
  int event1a, event1b;
  event1a = MPE_Log_get_event_number();
  event1b = MPE_Log_get_event_number();
  if(isMaster) 
    MPE_Describe_state(event1a, event1b, "load_parameters", "blue");
  MPE_Log_event(event1a, 0, "start load_parameters"); /* record time spent on reading input */ 
#endif
  
  SDDS_ClearErrors();
#if SDDS_MPI_IO 
  /* All the processes will read the wake file, but not in parallel.
     Zero the Memory when call  SDDS_InitializeInput */
  load_request[load_requests].table.parallel_io = 0; 
#endif

  if (!SDDS_InitializeInputFromSearchPath(&load_request[load_requests].table, 
                                          load_request[load_requests].filename)) {
    fprintf(stdout, "%s: couldn't initialize SDDS input for load_parameters file %s\n", 
            load_request[load_requests].filename,
            allow_missing_files?"warning":"error");
    fflush(stdout);
    if (!allow_missing_files) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    return 1;
  }
  if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Element_ColumnName))<0 ||
      SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
    fprintf(stdout, "Column \"%s\" is not in file %s or is not of string type.\n", Element_ColumnName,
            load_request[load_requests].filename);
    fflush(stdout);
    exitElegant(1);
  }
  if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Parameter_ColumnName))<0 ||
      SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
    fprintf(stdout, "Column \"%s\" is not in file %s or is not of string type.\n", Parameter_ColumnName, 
            load_request[load_requests].filename);
    fflush(stdout);
    exitElegant(1);
  }
  load_request[load_requests].string_data = 0;
  if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Value_ColumnName))>=0) {
    if (SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_DOUBLE) {
      fprintf(stdout, "Column \"%s\" is not in file %s or is not of double-precision type.\n", 
              Value_ColumnName, load_request[load_requests].filename);
      fflush(stdout);
      exitElegant(1);
    } 
  } else {
    if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, ValueString_ColumnName))<0 ||
        SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
      fprintf(stdout, "Column \"%s\" is not in file %s or is not of string type.\n",
              ValueString_ColumnName, load_request[load_requests].filename);
      fflush(stdout);
      exitElegant(1);
    }
    load_request[load_requests].string_data = 1;
  }
  if ((include_type_pattern || exclude_type_pattern) &&
      ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, ElementType_ColumnName))<0 ||
       SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING)) {
    fprintf(stdout, "include_type_pattern and/or exclude_type_pattern given, but\n");
    fprintf(stdout, "column \"%s\" is not in file %s or is not of string type.\n",
            ElementType_ColumnName, load_request[load_requests].filename);
    fflush(stdout);
    exitElegant(1);
  }

  /* The Mode column is optional: */
  if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Mode_ColumnName))>=0 &&
      SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_STRING) {
    fprintf(stdout, "Column \"%s\" is in file %s but is not of string type.\n", 
            Mode_ColumnName, load_request[load_requests].filename);
    fflush(stdout);
    exitElegant(1);
  }
  /* The Occurence column is optional: */ 
  if ((index=SDDS_GetColumnIndex(&load_request[load_requests].table, Occurence_ColumnName))>=0 && 
      SDDS_GetColumnType(&load_request[load_requests].table, index)!=SDDS_LONG) {
    fprintf(stdout, "Column \"%s\" is in file %s but is not of long-integer type.\n", 
            Occurence_ColumnName, load_request[load_requests].filename);
    fflush(stdout);
    exitElegant(1);
  }
  
  if (SDDS_NumberOfErrors()) {
    fprintf(stdout, "error: an uncaught error occured in SDDS routines (setup_load_parameters):\n");
    fflush(stdout);
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }

  load_request[load_requests].last_code = 0;
  load_request[load_requests].starting_value = NULL;
  load_request[load_requests].element = NULL;
  load_request[load_requests].reset_address = NULL;
  load_request[load_requests].value_type = NULL;
  load_request[load_requests].element_flags = NULL;
  load_request[load_requests].values = 0;
  load_requests++;
  if (load_request[load_requests-1].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) {
    /* do this right away so that it gets propagated into error and vary operations */
    do_load_parameters(beamline, 1);
    fprintf(stdout, "New length per pass: %21.15e m\n",
            compute_end_positions(beamline));
#ifdef  USE_MPE
  MPE_Log_event(event1b, 0, "end load_parameters");
#endif
    /* No reason to keep this, so just decrement the counter */
    load_requests --;
    return 1;
  }
#ifdef  USE_MPE
  MPE_Log_event(event1b, 0, "end load_parameters");
#endif
  return 0;
}


long do_load_parameters(LINE_LIST *beamline, long change_definitions)
{
  long i, j, mode_flags, code, rows,  param, allFilesRead, allFilesIgnored;
  char **element, **parameter, **type, **mode, *p_elem, *p_elem0, *ptr, lastMissingElement[100];
  double *value, newValue;
  char **valueString;
  ELEMENT_LIST *eptr;
  long element_missing;
  int32_t numberChanged, totalNumberChanged = 0;
  long lastMissingOccurence = 0;
  int32_t *occurence;
  char elem_param[1024];
  htab *hash_table;

  if (!load_requests || !load_parameters_setup)
    return NO_LOAD_PARAMETERS;
  allFilesRead = 1;
  allFilesIgnored = 1;
  
  for (i=0; i<load_requests; i++) {
    if (load_request[i].flags&COMMAND_FLAG_IGNORE)
      continue;

    hash_table = hcreate(12); /* create a hash table with the size of 2^12, it can grow automatically if necessary */

    allFilesIgnored = 0;
    if (load_request[i].last_code) {
      for (j=0; j<load_request[i].values; j++) {
        load_request[i].element[j]->flags = load_request[i].element_flags[j];
        switch (load_request[i].value_type[j]) {
        case IS_DOUBLE:
          *((double*)(load_request[i].reset_address[j])) = load_request[i].starting_value[j];
          break;
        case IS_LONG:
          *((long*)(load_request[i].reset_address[j])) = nearestInteger(load_request[i].starting_value[j]);
          break;
        default:
          fprintf(stdout, "internal error: invalid value type for load request value restoration\n");
          fflush(stdout);
          exitElegant(1);
          break;
        }
      }
    }

    if ((code=load_request[i].last_code=SDDS_ReadTable(&load_request[i].table))<1) {
      free(load_request[i].reset_address);
      load_request[i].reset_address = NULL;
      free(load_request[i].value_type);
      load_request[i].value_type = NULL;
      free(load_request[i].element_flags);
      load_request[i].element_flags = NULL;
      free(load_request[i].starting_value);
      load_request[i].starting_value = NULL;
      free(load_request[i].element);
      load_request[i].element = NULL;
      if (code<0) {
        fprintf(stdout, "warning: file %s ends unexpectedly (code=%ld)\n", load_request[i].filename,
                code);
        fflush(stdout);
        load_request[i].flags |= COMMAND_FLAG_IGNORE;
        continue;
      }
      fprintf(stdout, "Error: problem reading data from load_parameters file %s\n", load_request[i].filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    allFilesRead = 0;
    SDDS_SetRowFlags(&load_request[i].table, 1);
    if ((rows=SDDS_CountRowsOfInterest(&load_request[i].table))<=0) {
      load_request[i].last_code = 0;
      load_request[i].flags |= COMMAND_FLAG_IGNORE;
      continue;
    }
    mode = NULL;
    if (!(element  =(char **)SDDS_GetColumn(&load_request[i].table, Element_ColumnName)) ||
        !(parameter=(char **)SDDS_GetColumn(&load_request[i].table, Parameter_ColumnName)) || 
        (SDDS_GetColumnIndex(&load_request[i].table, Mode_ColumnName)>=0 &&
         !(mode     =(char **)SDDS_GetColumn(&load_request[i].table, Mode_ColumnName)))) {
      fprintf(stdout, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    type = NULL;
    if ((load_request[i].includeTypePattern || load_request[i].excludeTypePattern) &&
	!(type = (char**)SDDS_GetColumn(&load_request[i].table, ElementType_ColumnName))) {
      fprintf(stdout, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    valueString = NULL;
    value = NULL;
    if ((!load_request[i].string_data  &&
         !(value    =(double*)SDDS_GetColumn(&load_request[i].table, Value_ColumnName))) ||
        (load_request[i].string_data &&
         !(valueString = (char **)SDDS_GetColumn(&load_request[i].table, ValueString_ColumnName)))) {
      fprintf(stdout, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }

    occurence = NULL;
    if (!(load_request[i].flags&COMMAND_FLAG_IGNORE_OCCURENCE)) {
      if (verbose)
	fprintf(stdout, "Using occurence data.\n");
      if (SDDS_GetColumnIndex(&load_request[i].table, Occurence_ColumnName)>=0) {
        if (!(occurence = (int32_t *)SDDS_GetColumn(&load_request[i].table, Occurence_ColumnName))) {
          fprintf(stdout, "Error: problem accessing data from load_parameters file %s\n", load_request[i].filename);
          fflush(stdout);
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
      }
    }

    load_request[i].values = 0;
    element_missing = 0;
    lastMissingElement[0] = 0;
    lastMissingOccurence = 0;
    for (j=0; j<rows; j++) {
      /* If the user gives the use_first flag, then we load only the first instance for
       * any parameter.   If occurence data is present, then we load the first instance
       * for each occurence.  Otherwise, we load the first instance ignoring occurrence
       * (this would happen in change_defined_values mode).
       * If use_first is not given, then we load all values, which is the slowest option
       * but also the original default behavior.
       */
      if (load_request[i].flags&COMMAND_FLAG_USE_FIRST){
        if (!occurence) {
          strcpy(elem_param, element[j]);
          strcat(elem_param, parameter[j]);
        } else {
          sprintf(elem_param, "%s%s%"PRId32, element[j], parameter[j], occurence[j]);
        }
        if (!hadd(hash_table, elem_param, strlen(elem_param), NULL))
          continue;
      }
      eptr = NULL;
      if (occurence) {
        if (occurence[j]<=lastMissingOccurence &&
            strcmp(lastMissingElement, element[j])==0)
          continue;
      } else {
        if (strcmp(lastMissingElement, element[j])==0)
          continue;
      }
      if (load_request[i].includeNamePattern &&
          !wild_match(element[j], load_request[i].includeNamePattern))
        continue;
      if (load_request[i].includeItemPattern &&
          !wild_match(parameter[j], load_request[i].includeItemPattern))
        continue;
      if (load_request[i].includeTypePattern &&
          !wild_match(type[j], load_request[i].includeTypePattern))
        continue;
      if (load_request[i].excludeNamePattern &&
          wild_match(element[j], load_request[i].excludeNamePattern))
        continue;
      if (load_request[i].excludeItemPattern &&
          wild_match(parameter[j], load_request[i].excludeItemPattern))
        continue;
      if (load_request[i].excludeTypePattern &&
          wild_match(type[j], load_request[i].excludeTypePattern))
        continue;
      element_missing = 0;

      /* if occurence is available, we can take advantage of hash table */
      if ((occurence && (!find_element_hash(element[j], occurence[j], &eptr, &(beamline->elem)))) ||
	  (!occurence && (!find_element_hash(element[j], 1,  &eptr, &(beamline->elem)))))  {
	if (occurence) {
	  if (missingElementWarningsLeft || !allow_missing_elements) {
	    fprintf(stdout, "%s: unable to find occurence %" PRId32 " of element %s (do_load_parameters)\n", 
		    allow_missing_elements?"Warning":"Error",
		    occurence[j], element[j]);
	    if (allow_missing_elements && --missingElementWarningsLeft==0)
	      fprintf(stdout, "Further missing elements warnings suppressed\n");
	    fflush(stdout);
	  }
	} else {
	  if (missingElementWarningsLeft || !allow_missing_elements) {
	    fprintf(stdout, "%s: unable to find element %s (do_load_parameters)\n", 
		    allow_missing_elements?"Warning":"Error",
		    element[j]);
	    if (allow_missing_elements && --missingElementWarningsLeft==0)
	      fprintf(stdout, "Further missing elements warnings suppressed\n");
	    fflush(stdout);
	  }
	}
	if (!allow_missing_elements)
	  exitElegant(1);
	element_missing = 1;
      }
  
      if (element_missing) {
        if (occurence)
          lastMissingOccurence = occurence[j];
        strncpy(lastMissingElement, element[j], 99);
        continue;
      }
      lastMissingElement[0] = 0;
      lastMissingOccurence = 0;
      if ((param = confirm_parameter(parameter[j], eptr->type))<0) {
        if (missingParameterWarningsLeft || !allow_missing_parameters) {
#if !USE_MPI
          fprintf(stdout, "%s: element %s does not have a parameter %s (do_load_parameters)\n",
                  allow_missing_parameters?"Warning":"Error",
                  eptr->name, parameter[j]);
#else
	  if (allow_missing_parameters)
	    fprintf(stdout, "%s: element %s does not have a parameter %s (do_load_parameters)\n",
		    "Warning", eptr->name, parameter[j]);
	  else
	    fprintf(stderr, "%s: element %s does not have a parameter %s (do_load_parameters)\n",
		    "Error", eptr->name, parameter[j]);
#endif
          if (allow_missing_parameters && --missingParameterWarningsLeft==0)
            fprintf(stdout, "Further missing parameters warnings suppressed\n");
          fflush(stdout);
        }
        if (!allow_missing_parameters)
          exitElegant(1);
        continue;
      }
      mode_flags = 0;
      if (mode) 
        while ((ptr=get_token_t(mode[j], " \t,+"))) {
          long k;
          if ((k=match_string(ptr, load_mode, LOAD_MODES, UNIQUE_MATCH))<0) {
            fprintf(stdout, "Error: unknown/ambiguous mode specifier %s (do_load_parameters)\nKnown specifiers are:\n",
                    ptr);
            fflush(stdout);
            for (k=0; k<LOAD_MODES; k++)
              fprintf(stdout, "    %s\n", load_mode[k]);
              fflush(stdout);
            exitElegant(1);
          }
          mode_flags |= 1<<k;
          free(ptr);
        }
      if (mode_flags==0)
        mode_flags = LOAD_FLAG_ABSOLUTE;
      if (mode_flags&LOAD_FLAG_IGNORE)
        continue;
      if (verbose)
        fprintf(stdout, "Working on row %ld of file\n", j);
      
      if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) {
        change_defined_parameter(element[j], param, eptr->type, value?value[j]:0, 
                                 valueString?valueString[j]:NULL, 
                                 mode_flags+(verbose?LOAD_FLAG_VERBOSE:0));
      }
      numberChanged = 0;
      do {
        numberChanged++;
        p_elem = eptr->p_elem;
        p_elem0 = eptr->p_elem0;
        load_request[i].reset_address 
          = trealloc(load_request[i].reset_address,
                     sizeof(*load_request[i].reset_address)*(load_request[i].values+1));
        load_request[i].value_type 
          = trealloc(load_request[i].value_type,
                     sizeof(*load_request[i].value_type)*(load_request[i].values+1));
        load_request[i].element_flags 
          = trealloc(load_request[i].element_flags,
                     sizeof(*load_request[i].element_flags)*(load_request[i].values+1));
        load_request[i].starting_value 
          = trealloc(load_request[i].starting_value,
                     sizeof(*load_request[i].starting_value)*(load_request[i].values+1));
        load_request[i].element 
          = trealloc(load_request[i].element,
                     sizeof(*load_request[i].element)*(load_request[i].values+1));
        load_request[i].reset_address[load_request[i].values]
          = ((double*)(p_elem+entity_description[eptr->type].parameter[param].offset));
        load_request[i].element[load_request[i].values] = eptr;
        switch (entity_description[eptr->type].parameter[param].type) {
        case IS_DOUBLE:
          if (valueString) {
            if (!sscanf(valueString[j], "%lf", &newValue)) {
              fprintf(stdout, "Error: unable to scan double from \"%s\"\n", valueString[j]);
              fflush(stdout);
              exitElegant(1);
            }
          } else {
            newValue = value[j];
          }
          if (eptr->divisions>1 && (entity_description[eptr->type].parameter[param].flags&PARAM_DIVISION_RELATED))
            newValue /= eptr->divisions;
          if (verbose)
            fprintf(stdout, "Changing %s.%s #%" PRId32 "  from %21.15e to ",
                    eptr->name, 
                    entity_description[eptr->type].parameter[param].name, 
                    occurence?occurence[j]:numberChanged,
                    *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          load_request[i].starting_value[load_request[i].values]
            = *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset));
          load_request[i].value_type[load_request[i].values] = IS_DOUBLE;
          if (mode_flags&LOAD_FLAG_ABSOLUTE) {
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) = newValue;
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((double*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) = newValue;
	      
	  } else if (mode_flags&LOAD_FLAG_DIFFERENTIAL) {
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) += newValue;
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((double*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) += newValue;
	  } else if (mode_flags&LOAD_FLAG_FRACTIONAL) {
            *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((double*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
	  }
          if (verbose)
            fprintf(stdout, "%21.15e (%21.15e)\n",
                    *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)),
                    *((double*)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          break;
        case IS_LONG:
          if (valueString) {
            if (!sscanf(valueString[j], "%lf", &newValue)) {
              fprintf(stdout, "Error: unable to scan double from \"%s\"\n", valueString[j]);
              fflush(stdout);
              exitElegant(1);
            }
          } else {
            newValue = value[j];
          }
          if (verbose)
            fprintf(stdout, "Changing %s.%s #%" PRId32 "  from %ld  to ",
                    eptr->name,
                    entity_description[eptr->type].parameter[param].name, numberChanged,
                    *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          load_request[i].starting_value[load_request[i].values]
            = *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset));
          load_request[i].value_type[load_request[i].values] = IS_LONG;
          if (mode_flags&LOAD_FLAG_ABSOLUTE) {
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) = 
              nearestInteger(newValue);
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((long*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) = 
		nearestInteger(newValue);
	  } else if (mode_flags&LOAD_FLAG_DIFFERENTIAL) {
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) += 
              nearestInteger(newValue);
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((long*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) = 
		nearestInteger(newValue);
	  } else if (mode_flags&LOAD_FLAG_FRACTIONAL) {
            *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
	    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS)
	      *((long*)(p_elem0+entity_description[eptr->type].parameter[param].offset)) *= 1+newValue;
	  }
          if (verbose)
            fprintf(stdout, "%ld \n",
                    *((long*)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          break;
        case IS_STRING:
          load_request[i].value_type[load_request[i].values] = IS_STRING;
          if (verbose)
            fprintf(stdout, "Changing %s.%s #%" PRId32 "  from %s to ",
                    eptr->name,
                    entity_description[eptr->type].parameter[param].name, numberChanged,
                    *((char**)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          if (!SDDS_CopyString((char**)(p_elem+entity_description[eptr->type].parameter[param].offset),
                               valueString[j])) {
            fprintf(stdout, "Error (do_load_parameters): unable to copy value string\n");
            fflush(stdout);
            exitElegant(1);
          }
	  if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) {
	    if (!SDDS_CopyString((char**)(p_elem0+entity_description[eptr->type].parameter[param].offset),
				 valueString[j])) {
	      fprintf(stdout, "Error (do_load_parameters): unable to copy value string\n");
	      fflush(stdout);
	      exitElegant(1);
	    }
	  }
          if (verbose)
            fprintf(stdout, "%s\n", 
                    *((char**)(p_elem+entity_description[eptr->type].parameter[param].offset)));
            fflush(stdout);
          break;
        default:
          fprintf(stdout,
                  "Error: can't load parameter value for parameter %s of %s--not numeric parameter (do_load_parameters)\n",
                  parameter[j], element[j]);
          fflush(stdout);
          break;
        }
        eptr->flags |= 
          PARAMETERS_ARE_PERTURBED |
            ((entity_description[eptr->type].parameter[param].flags&PARAM_CHANGES_MATRIX)?VMATRIX_IS_PERTURBED:0);
        load_request[i].element_flags[load_request[i].values] = eptr->flags;
        load_request[i].values++;
      } while (!occurence && find_element_hash(element[j], numberChanged+1, &eptr, &(beamline->elem)));
      free(element[j]);
      free(parameter[j]);
      if (mode)
        free(mode[j]);
      totalNumberChanged += numberChanged;
    }
    
    free(element);
    free(parameter);
    if (mode)
      free(mode);
    free(value);
    if (occurence)
      free(occurence);
    if (load_request[i].flags&COMMAND_FLAG_CHANGE_DEFINITIONS) {
      free(load_request[i].reset_address);
      load_request[i].reset_address = NULL;
      free(load_request[i].value_type);
      load_request[i].value_type = NULL;
      free(load_request[i].element_flags);
      load_request[i].element_flags = NULL;
      free(load_request[i].starting_value);
      load_request[i].starting_value = NULL;
      free(load_request[i].element);
      load_request[i].element = NULL;
      load_request[i].flags |= COMMAND_FLAG_IGNORE;   /* ignore hereafter */
    }

    if (hash_table) { 
      hdestroy(hash_table);                         /* destroy hash table */  
      hash_table = NULL;
    } 
  }
  
  if (!allFilesRead || allFilesIgnored) {
    compute_end_positions(beamline);
    fprintf(stdout, "%" PRId32 "  parameter values loaded\n", totalNumberChanged);
    return PARAMETERS_LOADED;
  }
  return PARAMETERS_ENDED;
}

void finish_load_parameters()
{
  long i;
  for (i=0; i<load_requests; i++) {
    if (!SDDS_Terminate(&load_request[i].table)) {
      fprintf(stdout, "Error: unable to terminate load_parameters SDDS file %s\n", load_request[i].filename);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    free(load_request[i].filename);
  }
  if (load_requests)
    free(load_request);
  load_request = NULL;
  load_requests = 0;
}

static long dumpingLatticeParameters = 0;
static long iElementName, iElementParameter, iParameterValue, iElementType, iOccurence, iElementGroup;
static SDDS_DATASET SDDS_dumpLattice;
void dumpLatticeParameters(char *filename, RUN *run, LINE_LIST *beamline)
{
  SDDS_DATASET *SDDSout;
  long iElem, iParam;
  ELEMENT_LIST *eptr;
  PARAMETER *parameter;
  long row, maxRows, doSave;
  double value=0.0;
  
  SDDSout = &SDDS_dumpLattice;
  if (!dumpingLatticeParameters) {
    if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 0, NULL, NULL, filename) ||
        (iElementName=SDDS_DefineColumn(SDDSout, "ElementName", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        (iElementParameter=SDDS_DefineColumn(SDDSout, "ElementParameter", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        (iParameterValue=SDDS_DefineColumn(SDDSout, "ParameterValue", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
        (iElementType=SDDS_DefineColumn(SDDSout, "ElementType", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        (iOccurence=SDDS_DefineColumn(SDDSout, "ElementOccurence", NULL, NULL, NULL, NULL, SDDS_LONG, 0))<0 ||
        (iElementGroup=SDDS_DefineColumn(SDDSout, "ElementGroup", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
        !SDDS_WriteLayout(SDDSout)) {
      fprintf(stdout, "Problem setting up parameter output file\n");
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    dumpingLatticeParameters = 1;
  }

  if (!SDDS_StartPage(SDDSout, maxRows=beamline->n_elems*10))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  row = 0;
  
  eptr = &(beamline->elem);
  for (iElem=0; iElem<beamline->n_elems; iElem++) {
    /*
      fprintf(stdout, "name=%s, divisions=%" PRId32 "   first=%hd\n",
      eptr->name, eptr->divisions, eptr->firstOfDivGroup);
      */
    if (eptr->divisions>1 && !eptr->firstOfDivGroup) {
      /* don't emit data for every member of a divided group */
      eptr = eptr->succ;
      continue;
    }
    if (!(parameter = entity_description[eptr->type].parameter))
      SDDS_Bomb("parameter entry is NULL for entity description (dumpLatticeParameters)");
    if (!(eptr->name))
      SDDS_Bomb("element name is NULL (dumpLatticeParameters)");
    for (iParam=0; iParam<entity_description[eptr->type].n_params; iParam++) {
      doSave = 1;
      switch (parameter[iParam].type) {
      case IS_DOUBLE: 
        value = *(double*)(eptr->p_elem+parameter[iParam].offset);
        if (parameter[iParam].flags&PARAM_DIVISION_RELATED &&
            eptr->divisions>1)
          value *= eptr->divisions;
        break;
      case IS_LONG:
        value = *(long*)(eptr->p_elem+parameter[iParam].offset);
        break;
      default:
        doSave = 0;
        break;
      }
      /* some kludges to avoid saving things that shouldn't be saved as the user didn't
         set them in the first place
         */
      if (strcmp(parameter[iParam].name, "PHASE_REFERENCE")==0 && value>LONG_MAX/2) {
        doSave = 0;
      }

      if (doSave) {
        if (!(parameter[iParam].name)) 
          SDDS_Bomb("parameter name is NULL (dumpLatticeParameters)");
        if (row>=maxRows) {
          if (!SDDS_LengthenTable(SDDSout, 1000))
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          maxRows += 1000;
        }
        if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                               iOccurence, eptr->occurence,
                               iElementName, eptr->name,
                               iElementParameter, parameter[iParam].name,
                               iParameterValue, value, 
                               iElementType, entity_name[eptr->type],
                               iElementGroup, eptr->group?eptr->group:"",
                               -1)) {
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        row++;
      }
    }
    eptr = eptr->succ;
  }
  if (!SDDS_WritePage(SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!inhibitFileSync)
    SDDS_DoFSync(SDDSout);
}

void finishLatticeParametersFile() 
{
  if (dumpingLatticeParameters && !SDDS_Terminate(&SDDS_dumpLattice))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  dumpingLatticeParameters = 0;
}

