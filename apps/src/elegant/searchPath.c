/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* Copyright Argonne National Lab and Michael Borland, 2001 */

#include "track.h"
#include "mdb.h"

static char *search_path = NULL;
void setSearchPath(char *input)
{
  if (search_path)
    free(search_path);
  if (input)
    cp_str(&search_path, input);
  else
    search_path = NULL;
}

char *findFileInSearchPath(char *filename)
{
  char *path, *pathList, *tmpName;
  char *sddsTags = NULL;
  
  if (!filename || !strlen(filename))
    return NULL;
  if ((sddsTags=strchr(filename, '='))) {
    /* <filename>=<x>+<y> form ? */
    if (!strchr(sddsTags+1, '+'))
      sddsTags = NULL;
    else 
      /* yes */
      *sddsTags++ = 0;
  }
  if (search_path && strlen(search_path)) {
    cp_str(&pathList, search_path);
    while ((path=get_token(pathList))) {
      tmpName = malloc(strlen(filename)+strlen(path)+2+(sddsTags?strlen(sddsTags)+2:0));
      sprintf(tmpName, "%s/%s", path, filename);
      if (fexists(tmpName)) {
        if (sddsTags) {
          /* put the sddsTags back on the end */
          strcat(tmpName, "=");
          strcat(tmpName, sddsTags);
        }
        return tmpName;
      }
    }
  }
  if (fexists(filename)) {
    if (sddsTags)
      *(sddsTags-1) = '=';
    cp_str(&tmpName, filename);
    return tmpName;
  }
  return NULL;
}
  

long getTableFromSearchPath(TABLE *tab, char *file, long sampleInterval, long flags)
{
  char *filename;
  long value;
  if (!(filename=findFileInSearchPath(file)))
    return 0;
  value = get_table(tab, filename, sampleInterval, flags);
  free(filename);
  return value;
}

long SDDS_InitializeInputFromSearchPath(SDDS_DATASET *SDDSin, char *file)
{
  char *filename;
  long value;
  if (!(filename=findFileInSearchPath(file)))
    return 0;
  fprintf(stdout, "File %s found: %s\n", file, filename);
  value = SDDS_InitializeInput(SDDSin, filename);
  free(filename);
  return value;
}

