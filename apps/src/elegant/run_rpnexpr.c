/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: run_rpn.c
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#include "run_rpnexpr.h"

void run_rpn_expression(NAMELIST_TEXT *nltext)
{
/*    expression = NULL; 
 */

    log_entry("run_rpn_expression");
   
    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&rpn_expression, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &rpn_expression);

    if (expression) {
        rpn(expression);
        if (rpn_check_error()) exitElegant(1);
        }
    log_exit("run_rpn_expression");
    }

void run_rpn_load(NAMELIST_TEXT *nltext, RUN *run)
{
  SDDS_DATASET SDDSin;
  long code, foundPage, iColumn, matchRow, rows, iParameter;
  int32_t columns, parameters;
  char *parameterValue = NULL;
  double *data, data1;
  char **columnName, **matchColumnData, *memName = NULL;
  char **parameterName;
  
  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&rpn_load, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &rpn_load);
  
  if (match_column && strlen(match_column)) {
    if (use_row!=-1) {
      fprintf(stdout, "Error: you asked to match a column and also gave use_row.\n");
      exitElegant(1);
    } 
    if (!match_column_value || !strlen(match_column_value)) {
      fprintf(stdout, "Error: you must give match_column_value with match_column\n");
      exitElegant(1);
    }
  }
  if (match_parameter && strlen(match_parameter)) {
    if (use_page!=-1) {
      fprintf(stdout, "Error: you asked to match a parameter and also gave use_page.\n");
      exitElegant(1);
    }
    if (!match_parameter_value || !strlen(match_parameter_value)) {
      fprintf(stdout, "Error: you must give match_parameter_value with match_parameter\n");
      exitElegant(1);
    }
  }
    
  if (!filename || !strlen(filename)) {
    fprintf(stdout, "Error: no filename given for rpn_load.\n");
    exitElegant(1);
  }

  filename = compose_filename(filename, run->rootname);
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, filename)) {
    fprintf(stdout, "Error: couldn't initialize SDDS input for %s\n",
            filename);
    exitElegant(1);
  }

  foundPage = 0;
  while ((code=SDDS_ReadPage(&SDDSin))>0) {
    if (use_page>0) {
      if (code==use_page) {
        foundPage = 1;
        break;
      }
      continue;
    }
    if (match_parameter && strlen(match_parameter)) {
      if (!(parameterValue=SDDS_GetParameterAsString(&SDDSin, match_parameter, NULL)))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!wild_match(parameterValue, match_parameter_value))
        continue;
      foundPage = 1;
      break;
    }
    if (use_page==-1 && SDDS_CheckEndOfFile(&SDDSin)==1) {
      foundPage = 1;
      break;
    }
  }

  if (!foundPage) {
    fprintf(stdout, "Error: no appropriate page found\n");
    exitElegant(1);
  }

  if (!load_parameters) {
    if ((columnName = SDDS_GetColumnNames(&SDDSin, &columns))==NULL) {
      fprintf(stdout, "Warning: No columns in file!\n");
      return;
    }

    rows = SDDS_RowCount(&SDDSin);
    matchRow = rows-1;
    if (use_row!=-1) {
      if (use_row>=rows) {
        fprintf(stdout, "Error: number of rows in file (%ld) less than needed for use_row=%ld\n",
                rows, use_row);
        exitElegant(1);
      }
      matchRow = use_row;
    } 

    if (match_column) {
      if (SDDS_GetNamedColumnType(&SDDSin, match_column)!=SDDS_STRING) {
        fprintf(stdout, "Error: column %s nonexistent or not string type.\n",
                match_column);
        exitElegant(1);
      }
      if (!(matchColumnData=SDDS_GetColumn(&SDDSin, match_column))) {
        fprintf(stdout, "Error: unable to get data for column %s\n", match_column);
        exitElegant(1);
      }
      if (matching_row_number<0) {
        /* use last match */
        for (matchRow=rows-1; matchRow>=0; matchRow--)
          if (wild_match(matchColumnData[matchRow], match_column_value))
            break;
      } else {
        /* use nth match */
        for (matchRow=0; matchRow<rows; matchRow++)
          if (wild_match(matchColumnData[matchRow], match_column_value) &&
              matching_row_number-- == 0)
            break;
      }
      
      if (matchRow<0 || matchRow>=rows) {
        fprintf(stdout, "Error: unable to find match for %s in column %s\n",
                match_column_value, match_column);
        exitElegant(1);
      }
      SDDS_FreeStringArray(matchColumnData, rows);
    }
    
    for (iColumn=0; iColumn<columns; iColumn++) {
      switch (SDDS_GetNamedColumnType(&SDDSin, columnName[iColumn])) {
      case SDDS_CHARACTER:
      case SDDS_STRING:
        break;
      default:
        if (!(data=SDDS_GetColumnInDoubles(&SDDSin, columnName[iColumn]))) {
          fprintf(stdout, "Error: unable to get data for column %s as numerical data.\n",
                  columnName[iColumn]);
          exitElegant(1);
        }
        if (!(memName=SDDS_Realloc(memName, sizeof(*memName)*((tag?strlen(tag):0)+strlen(columnName[iColumn])+2)))) {
          fprintf(stdout, "Memory allocation failure trying to create memory name for loaded data\n");
          exitElegant(1);
        }
        if (tag && strlen(tag))
          sprintf(memName, "%s.%s", tag, columnName[iColumn]);
        else
          sprintf(memName, "%s", columnName[iColumn]);
        rpn_store(data[matchRow], NULL, rpn_create_mem(memName, 0));
        fprintf(stdout, "%le --> %s\n", data[matchRow], memName);
        free(columnName[iColumn]);
        free(data);
      }
    }
    if (memName)
      free(memName);
    if (columnName)
      free(columnName);
  } else {
    /* load data from parameters */
    if ((parameterName = SDDS_GetParameterNames(&SDDSin, &parameters))==NULL) {
      fprintf(stdout, "Warning: No parameters in file!\n");
      return;
    }

    for (iParameter=0; iParameter<parameters; iParameter++) {
      switch (SDDS_GetNamedParameterType(&SDDSin, parameterName[iParameter])) {
      case SDDS_CHARACTER:
      case SDDS_STRING:
        break;
      default:
        if (!SDDS_GetParameterAsDouble(&SDDSin, parameterName[iParameter], &data1)) {
          fprintf(stdout, "Error: unable to get data for parameter %s as numerical data.\n",
                  parameterName[iParameter]);
          exitElegant(1);
        }
        if (!(memName=SDDS_Realloc(memName, sizeof(*memName)*((tag?strlen(tag):0)+strlen(parameterName[iParameter])+2)))) {
          fprintf(stdout, "Memory allocation failure trying to create memory name for loaded data\n");
          exitElegant(1);
        }
        if (tag && strlen(tag))
          sprintf(memName, "%s.%s", tag, parameterName[iParameter]);
        else
          sprintf(memName, "%s", parameterName[iParameter]);
        rpn_store(data1, NULL, rpn_create_mem(memName, 0));
        fprintf(stdout, "%le --> %s\n", data1,  memName);
        free(parameterName[iParameter]);
      }
    }
    if (memName)
      free(memName);
    if (parameterName)
      free(parameterName);
  }
}

