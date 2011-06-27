/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* prototypes for this file are in prot.out */
/* routine: cfgets()
 * purpose: read from a file, skipping MAD-style comments.
 *          Returns the next non-comment line, including continuations 
 *          (marked with an & in the last column).  The routine also 
 *          deletes all spaces from the line, as a kludge for the parsing 
 *          routines.
 * 
 * Michael Borland, 1987
 */
#include "mdb.h"
#include "track.h"
#include <ctype.h>
#if defined(_WIN32)
#include <stdlib.h>
#else
#if defined(UNIX) && defined(GNU_C)
long toupper(char c);
long tolower(char c);
#endif
#endif

void delete_spaces(char *s);
void str_to_upper_quotes(char *s);
char *cfgets1(char *s, long n, FILE *fpin, long strip, long start);

char *cfgets(char *s, long n, FILE *fpin)
{
  s[0] = 0;
  if (!cfgets1(s, n, fpin, 1, 1))
    return NULL;
  if (s[0]!='%') 
    str_to_upper_quotes(s);
  return s;
}

char *cfgets1(char *s, long n, FILE *fpin, long strip, long start)
{
  register long l;

  while (fgets(s, n, fpin)) {
    if (start && s[0]=='%')
      strip = 0;
    if (s[0]=='!')
      continue;
    start = 0;
    chop_nl(s);
    if (strip)
      delete_spaces(s);
    l = strlen(s);
    while (l!=0 && s[l]<27)
      l--;
    if (s[l]=='&') {
      s[l] = 0;
      cfgets1(s+l, n-l, fpin, strip, 0);
      return(s);
    }
    else {
      return(s);
    }
  }
  return NULL;
}

void delete_spaces(char *s)
{
    char *ptr, *ptr0;
    ptr0 = s;
    while (*s) {
        if (*s=='"' && (ptr0==s || *(s-1)!='\\')) {
            s++;
            while (*s && (*s!='"' || *(s-1)=='\\'))
                s++;
            if (*s=='"' && *(s-1)!='\\')
                s++;
            }
        else if (*s==' ' || *s=='\011') {
            ptr = s++;
            while (*s==' ' || *s=='\011')
                s++;
            strcpy_ss(ptr, s);
            s = ptr;
            }
        else
            s++;
        }
    }

void str_to_upper_quotes(char *s)
{
  char *ptr0;
  ptr0 = s;
  if (!s)
    return;
  while (*s) {
    if (*s=='"' && (ptr0==s || *(s-1)!='\\')) {
      while (*++s && (*s!='"' || *(s-1)=='\\'))
        ;
    }
    else 
      *s = toupper(*s);
    s++;
  }
}
