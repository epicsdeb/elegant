/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: taylorSeries() (derived from polynomial_kicks())
 * purpose: apply a taylor series map to coordinates.
 * 
 * Louis Emery 2004
 */

/*
 $Log: not supported by cvs2svn $
 Revision 1.3  2010/02/05 17:48:07  soliday
 Fixed minor issues related to compiler warnings.

 Revision 1.2  2010/02/04 15:17:27  borland
 No longer use the bomb() routine.  Instead, use bombElegant(), which allows
 better control of what happens when exiting.  Added "failed" semaphore option.
 Switched from process_namelist() to processNamelist() for better response
 to errors.
 Includes Y. Wang's changes to parallelize shell and line beam types.

 Revision 1.1  2004/03/26 15:59:52  borland
 First version in repository, by L. Emery.

*/

#include "mdb.h"
#include "track.h"

#define N_MAPCOORDINATES 6
/* qx and qy is px/Po and py/Po */
char *mapCoordinateList[N_MAPCOORDINATES] = {
  "x", "qx", "y", "qy", "cdt", "delta"
   };

void initialize_taylorSeries(TAYLORSERIES *taylorSeries)
{
  SDDS_DATASET SDDSin;
  char buffer[1024];
  TAYLORSERIES_DATA *data;
  char *Coordinate;
  long readCode, problem, index;
  
  if (taylorSeries->elementInitialized)
    return;
  if (!taylorSeries->filename)
    bombElegant("TAYLORSERIES element doesn't have filename", NULL);
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, taylorSeries->filename)) {
    sprintf(buffer, "Problem opening file %s (TAYLORSERIES)\n", taylorSeries->filename);
    SDDS_SetError(buffer);
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  if (SDDS_CheckParameter(&SDDSin, "Coordinate", NULL, SDDS_STRING, stdout)!=SDDS_CHECK_OK )
    bombElegant("problems with Coordinate parameter data in TAYLORSERIES input file", NULL);
  if (SDDS_CheckColumn(&SDDSin, "Ix", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "Iqx", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "Iy", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "Iqy", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "Icdt", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "Idelta", NULL, SDDS_ANY_INTEGER_TYPE, stdout)!=SDDS_CHECK_OK )
    bombElegant("problems with integer power data in TAYLORSERIES input file", NULL);
  if (SDDS_CheckColumn(&SDDSin, "Coefficient", NULL, SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK )
    bombElegant("problems with Coefficient data in TAYLORSERIES input file", NULL);
  while ((readCode=SDDS_ReadPage(&SDDSin))>0)  {
    /* determine the coordinate that is mapped */
    if (!(SDDS_GetParameter(&SDDSin, "Coordinate", &Coordinate))) {
      sprintf(buffer, "problem reading Coordinate parameter data for TAYLORSERIES file %s\n", taylorSeries->filename);
      SDDS_SetError(buffer);
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    /* index is the index of the map array in element definition */
    if (0>(index = match_string(Coordinate, mapCoordinateList, N_MAPCOORDINATES, 0))) {
      sprintf(buffer, "invalid Coordinate parameter data for TAYLORSERIES file %s\n", taylorSeries->filename);
      SDDS_SetError(buffer);
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    data = &taylorSeries->coord[index];
    if ((data->terms = SDDS_RowCount(&SDDSin))<=0) {
      /* presently make having zero rows an error 
	 later we can make the default behavior of 0 rows be
	 equivalent to unity map. */
      fprintf(stdout, "Warning: no data in TAYLORSERIES file %s\n", taylorSeries->filename);
      fflush(stdout);
      SDDS_Terminate(&SDDSin);
      return;
    }

    if (!(data->Ix=SDDS_GetColumnInLong(&SDDSin, "Ix")) ||
	!(data->Iqx=SDDS_GetColumnInLong(&SDDSin, "Iqx")) ||
	!(data->Iy=SDDS_GetColumnInLong(&SDDSin, "Iy")) ||
	!(data->Iqy=SDDS_GetColumnInLong(&SDDSin, "Iqy")) ||
	!(data->Icdt=SDDS_GetColumnInLong(&SDDSin, "Icdt")) ||
	!(data->Idelta=SDDS_GetColumnInLong(&SDDSin, "Idelta")) ||
	!(data->Coefficient=SDDS_GetColumnInDoubles(&SDDSin, "Coefficient"))) {
      sprintf(buffer, "Unable to read data for TAYLORSERIES file %s\n", taylorSeries->filename);
      SDDS_SetError(buffer);
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
      exitElegant(1);
    }
    data->mapInitialized = 1;
  }
  SDDS_Terminate(&SDDSin);
  problem = 0;
  for ( index=0 ; index < 6 ; index++ ) {
    data = &taylorSeries->coord[index];
    if (!data->mapInitialized) {
      problem = 1;
      sprintf(buffer, "problem with finding map for %s for TAYLORSERIES file %s\n", mapCoordinateList[index],taylorSeries->filename);
      SDDS_SetError(buffer);
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
    }
  }
  if (problem) {
    exitElegant(1);
  }
  taylorSeries->elementInitialized = 1;
}

long taylorSeries_tracking(
    double **particle,  /* initial/final phase-space coordinates */
    long n_part,        /* number of particles */
    TAYLORSERIES *taylorSeries,       /* taylorSeries element structure */
    double p_error,     /* p_nominal/p_central */
    double Po,          /* p_central */
    double **accepted,
    double z_start
    )
{
  TAYLORSERIES_DATA data;
  double dx, dy, dz;  /* transverse offsets of the element center */
  long i_part, i_top, is_lost=0;
  double *coord;
  double cos_tilt, sin_tilt;
  double x, xp, y, yp;
  double qx, qy, denom, cdt, dp;
  long i, j;
  double p, beta0, outputCoord[6];
  
    if (!particle)
        bombElegant("particle array is null (taylorSeries)", NULL);

    if (!taylorSeries)
        bombElegant("null TAYLORSERIES pointer (taylorSeries)", NULL);

    if (!taylorSeries->elementInitialized) 
      initialize_taylorSeries(taylorSeries);

    cos_tilt = cos(taylorSeries->tilt);
    sin_tilt = sin(taylorSeries->tilt);
    dx = taylorSeries->dx;
    dy = taylorSeries->dy;
    dz = taylorSeries->dz;

    i_top = n_part-1;
    for (i_part=0; i_part<=i_top; i_part++) {
        if (!(coord = particle[i_part])) {
            fprintf(stdout, "null coordinate pointer for particle %ld (taylorSeries)", i_part);
            fflush(stdout);
            abort();
            }
        if (accepted && !accepted[i_part]) {
            fprintf(stdout, "null accepted coordinates pointer for particle %ld (taylorSeries)", i_part);
            fflush(stdout);
            abort();
            }

        /* calculate coordinates in rotated and offset frame */
        coord[4] += dz*sqrt(1 + sqr(coord[1]) + sqr(coord[3]));
        coord[0] += -dx + dz*coord[1];
        coord[2] += -dy + dz*coord[3];

        x  =   cos_tilt*coord[0] + sin_tilt*coord[2];
        y  = - sin_tilt*coord[0] + cos_tilt*coord[2];
        xp =   cos_tilt*coord[1] + sin_tilt*coord[3];
        yp = - sin_tilt*coord[1] + cos_tilt*coord[3];

#if defined(IEEE_MATH)
        if (isnan(x) || isnan(xp) || isnan(y) || isnan(yp)) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }
#endif
        if (FABS(x)>COORD_LIMIT || FABS(y)>COORD_LIMIT ||
            FABS(xp)>SLOPE_LIMIT || FABS(yp)>SLOPE_LIMIT) {
            SWAP_PTR(particle[i_part], particle[i_top]);
            if (accepted)
                SWAP_PTR(accepted[i_part], accepted[i_top]);
            particle[i_top][4] = z_start;
            particle[i_top][5] = Po*(1+particle[i_top][5]);
            i_top--;
            i_part--;
            continue;
            }

	/* convert xp yp (angles) to qx=px/Po and qy=py/Po */
	dp = coord[5];
	p = Po*(1+dp);
        qx = (1+dp)*xp/(denom=sqrt(1+sqr(xp)+sqr(yp)));
        qy = (1+dp)*yp/denom;
	cdt = 0; /* I would imaging that the cdt coordinate doesn't appear in maps */
        beta0 = p/sqrt(sqr(p)+1);

        /* apply map */
	for ( i=0; i < 6; i++ ) {
	  data = taylorSeries->coord[i];
	  outputCoord[i] = 0.0;
	  for ( j=0; j < data.terms; j++) {
	    outputCoord[i] += data.Coefficient[j] * ipow(x,data.Ix[j]) * ipow(qx,data.Iqx[j])  
	      * ipow(y,data.Iy[j]) * ipow(qy,data.Iqy[j]) 
	      * ipow(cdt,data.Icdt[j]) * ipow(dp,data.Idelta[j]);
	  }
	}
	x = outputCoord[0];
	qx = outputCoord[1];
	y = outputCoord[2];
	qy = outputCoord[3];
	cdt = outputCoord[4];
	dp = outputCoord[5];

	/* convert qx qy to xp yp */
	if ((denom=sqr(1+dp)-sqr(qx)-sqr(qy))<=0) {
	  is_lost = 1;
	  break;
	}
	xp = qx/(denom=sqrt(denom));
	yp = qy/denom;
	
        /* undo the rotation and store in place of initial coordinates */
        /* don't need to change coord[0] or coord[2] since x and y are unchanged */
        coord[0] = cos_tilt*x - sin_tilt*y;
        coord[1] = cos_tilt*xp - sin_tilt*yp;
        coord[2] = sin_tilt*x + cos_tilt*y;
        coord[3] = sin_tilt*xp + cos_tilt*yp;
	/* the map gives time coordinates as c dt while elegant
	   gives time coordinates as beta c dt */
	coord[4] += outputCoord[4] * beta0;
        coord[5] = dp;

        /* remove the coordinate offsets */
        coord[0] += dx - coord[1]*dz;
        coord[2] += dy - coord[3]*dz;
        coord[4] -= dz*sqrt(1+ sqr(coord[1]) + sqr(coord[3]));
        }
    return(i_top+1);
    }
