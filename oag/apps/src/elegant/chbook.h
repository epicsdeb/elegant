/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* prototypes for chbook.c */

typedef struct {
  char *vname;
  double xmin, xmax, dx;
  int xbins;
  double *value;
  long length, count;
} book1;

typedef struct {
  char *xname, *yname;
  double xmin, xmax, dx;
  int xbins;
  double ymin, ymax, dy;
  int ybins;
  double *value;
  long length, count;
} book2;

typedef struct {
  char **vname;
  int nD;
  double *xmin, *xmax, *dx;
  int *xbins;
  double *value;
  long length, count;
} ntuple;

book1 *chbook1(char *vName, double xmin, double xmax, int xbins);
void chfill1(book1 *bName, double x, double weight);
void chprint1(book1 *bName, char *filename, char *description, int verbosity);
void free_hbook1(book1 *x);

book2 *chbook2(char *xName, char *yName, double xmin, double xmax, double ymin, double ymax, int xbins, int ybins);
void chfill2(book2 *bName, double x, double y, double weight);
void chprint2(book2 *bName, char *filename, char *description, int verbosity);
void free_hbook2(book2 *x);

ntuple *chbookn(char **vName, int NDimension, double *xmin, double *xmax, int *xbins);
void chfilln(ntuple *bName, double *x, double weight);
void free_hbookn(ntuple *x);

void checkbook2 (book2 *bName);
