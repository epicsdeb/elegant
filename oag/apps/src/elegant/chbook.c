/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include <stdio.h>
#include "mdb.h"
#include "SDDS.h"
#include "constants.h"
#include "chbook.h"

book1 *chbook1(char *vName, double xmin, double xmax, int xbins)
{
  book1 *book;
  book = (book1 *)malloc(sizeof(*book));

  book->vname = vName;

  book->xmin = xmin;
  book->xmax = xmax;
  book->dx = (xmax-xmin)/(double)xbins;
  book->xbins = xbins;  
  book->length = xbins;
  book->count = 0;

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfill1(book1 *bName, double x, double weight)
{
  int index;

  index = (int)((x - bName->xmin)/bName->dx);
  if (index < 0) index =0;
  if (index > bName->xbins-1) index = bName->xbins-1;
  
  bName->value[index] += weight;
  bName->count++;

  return;
}

void chprint1(book1 *bName, char *filename, char *description, int verbosity)
{
  SDDS_DATASET outPage;
  double *v;
  int i;

  v = calloc(sizeof(*v), bName->length);
  for(i=0; i<bName->length; i++) {
    v[i] = ((double)i+0.5)*bName->dx+bName->xmin;
  }

  /* Open file for writting */
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

  if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                             description, description, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&outPage, "VariableName", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xInterval", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xMinimum", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xDimension", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineColumn(&outPage, "x", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "weight", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&outPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "VariableName", bName->vname,
                          "xInterval", bName->dx,
                          "xMinimum", bName->xmin,
                          "xDimension", bName->length,
                          "Total_count", bName->count, NULL) ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v,            bName->length, "x") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, bName->value, bName->length, "weight") ||
      !SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free(v);
  return;
}

void free_hbook1 (book1 *x)
{
  free(x->value);
  free(x);
  return;
}

book2 *chbook2(char *xName, char *yName, double xmin, double xmax, double ymin, double ymax, int xbins, int ybins)
{
  book2 *book;

  book = (book2 *)malloc(sizeof(*book));

  book->xname = xName;
  book->yname = yName;

  book->xmin = xmin;
  book->xmax = xmax;
  book->dx = (xmax-xmin)/(double)xbins;
  book->xbins = xbins;  
  book->ymin = ymin;
  book->ymax = ymax;
  book->dy = (ymax-ymin)/(double)ybins;
  book->ybins = ybins;  

  book->length = xbins*ybins;
  book->count = 0;

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfill2(book2 *bName, double x, double y, double weight)
{
  int index[2], i;

  index[0] = (int)((x - bName->xmin)/bName->dx);
  if (index[0] < 0) index[0] =0;
  if (index[0] > (bName->xbins-1)) index[0] = bName->xbins-1;
  
  index[1] = (int)((y - bName->ymin)/bName->dy);
  if (index[1] < 0) index[1] =0;
  if (index[1] > (bName->ybins-1)) index[1] = bName->ybins-1;
  
  i = index[0]*bName->xbins + index[1];
  bName->value[i] += weight;
  bName->count++;

  return;
}

void chprint2(book2 *bName, char *filename, char *description, int verbosity)
{
  SDDS_DATASET outPage;

  /* Open file for writting */
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

  if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                             description, description, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&outPage, "Variable1Name", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Variable2Name", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xInterval", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xMinimum", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "xDimension", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "yInterval", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "yMinimum", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "yDimension", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineColumn(&outPage, "weight", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&outPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, bName->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Variable1Name", "x",
                          "Variable2Name", "y",
                          "xInterval", bName->dx,
                          "xMinimum", bName->xmin,
                          "xDimension", bName->xbins,
                          "yInterval", bName->dy,
                          "yMinimum", bName->ymin,
                          "yDimension", bName->ybins,
                          "Total_count", bName->count, NULL) ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, bName->value, bName->length, "weight") ||
      !SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void free_hbook2 (book2 *x)
{
  free(x->value);
  free(x);
  return;
}

ntuple *chbookn(char **vName, int ND, double *xmin, double *xmax, int *xbins)
{
  ntuple *book;
  int i;

  book = (ntuple *)malloc(sizeof(*book));

  book->nD = ND;

  book->vname = calloc(sizeof(*book->vname), ND);
  book->xmin = calloc(sizeof(*book->xmin), ND);
  book->xmax = calloc(sizeof(*book->xmax), ND);
  book->dx = calloc(sizeof(*book->dx), ND);
  book->xbins = calloc(sizeof(*book->xbins), ND);
  
  book->count = 0;
  book->length = 1;
  for (i=0; i<ND; i++) {
    book->vname[i] = vName[i];
    book->length *= xbins[i];
    book->xbins[i] = xbins[i];
    book->xmin[i] = xmin[i];
    book->xmax[i] = xmax[i];
    book->dx[i] = (xmax[i]-xmin[i])/(double)xbins[i];
  }

  book->value = calloc(sizeof(*book->value), book->length);

  return book;
}

void chfilln(ntuple *bName, double *x, double weight)
{
  int index;
  int i, temp;

  index = 0;
  for (i=0; i<bName->nD; i++) {
    temp = (int) ((x[i] - bName->xmin[i])/bName->dx[i]);
    if (temp < 0) temp = 0;
    if (temp > bName->xbins[i]-1) temp = bName->xbins[i]-1;
    if(i>0) index = index*bName->xbins[i-1]+temp;
    if(i==0) index = temp;
  }
  
  bName->value[index] += weight;
  bName->count++;
  /*
char a[100]="";    
  sprintf( a, "%d", num );
  for(i=0;i<6;i++) {
    strcpy(a + strlen(a), vName[i]);
    strcpy(a + strlen(a), ",");
  }
  printf("%s\n",a);                 
*/

  return;
}

void free_hbookn (ntuple *x)
{
  free(x->vname);
  free(x->xmin);
  free(x->xmax);
  free(x->dx);
  free(x->xbins);
  free(x->value);
  free(x);
  return;
}

void checkbook2 (book2 *bName) {
  printf("xname=%s, yname=%s\n", bName->xname, bName->yname);
  printf("xmin=%g, xmax=%g, dx=%g, xbins=%d\n", bName->xmin, bName->xmax, bName->dx, bName->xbins);
  printf("ymin=%g, ymax=%g, dy=%g, ybins=%d\n", bName->ymin, bName->ymax, bName->dy, bName->ybins);
  printf("length=%ld, count=%ld\n", bName->length, bName->count);

}
