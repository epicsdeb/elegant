/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: compose_fn.c
 *
 * Michael Borland, 1994.
 */
#include "mdb.h"
#include "track.h"

char *compose_filename(char *template, char *root_name)
{
    char *ptr;
    
    if (str_in(template, "%s")) {
        ptr = tmalloc(sizeof(char)*(strlen(template)+strlen(root_name)+1));
        sprintf(ptr, template, root_name);
        return(ptr);
        }
    else
        return(template);
    }

char *compose_filename_occurence(char *template, char *root_name, long occurence)
{
  char *ptr_s, *ptr_ld;
  char *ptr, *filename;
  char format[10];
  long i;
  
  ptr_s = str_in(template, "%s");
  ptr_ld = str_in(template, "%ld");
  for (i=1; i<20; i++) {
    sprintf(format, "%%%ldld", i);
    if ((ptr=str_in(template, format)) && (!ptr_ld || ptr_ld>ptr))
      ptr_ld = ptr;
    sprintf(format, "%%0%ldld", i);
    if ((ptr=str_in(template, format)) && (!ptr_ld || ptr_ld>ptr))
      ptr_ld = ptr;
  }
  if (ptr_ld != NULL) {
    sprintf(format, "-%s", ptr_ld);
    if ((ptr=str_in(template, format)))
      ptr_ld = ptr;
  }
  if (!occurence && ptr_ld) {
    filename = tmalloc(sizeof(char)*(strlen(template)));
    ptr = template;
    for (i=0; i<strlen(template); i++) {
      if (ptr<ptr_ld || ptr>ptr_ld+strlen(ptr_ld))
        strncat(filename, ptr, 1);
      ptr++;
    }
    ptr = tmalloc(sizeof(char)*(strlen(filename)+strlen(root_name)+100));
    if (ptr_s)
      sprintf(ptr, filename, root_name);
    else 
      ptr = filename;
    return ptr;
  }
  
  ptr = tmalloc(sizeof(char)*(strlen(template)+strlen(root_name)+100));
  if (ptr_s && ptr_ld) {
    if (ptr_s>ptr_ld)
      sprintf(ptr, template, occurence, root_name);
    else
      sprintf(ptr, template, root_name, occurence);
  } else if (ptr_s)
    sprintf(ptr, template, root_name);
  else if (ptr_ld) 
    sprintf(ptr, template, occurence);
  else 
    ptr = template;
  
  return ptr;
}
  
