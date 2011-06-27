/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#if USE_MPI 
#include "mpi.h" /* Defines the interface to MPI allowing the use of all MPI functions. */
#if USE_MPE
#include "mpe.h" /* Defines the MPE library */ 
#endif
#endif
#include <complex>
#include "mdb.h"
#include "track.h"
#include "match_string.h"

#define DEBUG 0

typedef struct {
  char *name, *type, *exclude;
} SC_SPEC;

static SC_SPEC *scSpec = NULL;
static long No_scSpec = 0;
static SPACE_CHARGE *sc = NULL;

void linearSCKick(double *coord, ELEMENT_LIST *eptr, double *center);
int  nonlinearSCKick(double *coord, ELEMENT_LIST *eptr, double *center, 
                     double sigmax, double sigmay, double *kick);

long getSCMULTSpecCount() 
{
  return (No_scSpec);
}
char *getSCMULTName()
{
  return (sc->name);
}

void addSCSpec(char *name, char *type, char *exclude)
{
  if (!(scSpec 
	= (SC_SPEC*)SDDS_Realloc(scSpec,
		       sizeof(*scSpec)*(No_scSpec+1))))
    bombElegant((char*)"memory allocation failure", NULL);
  scSpec[No_scSpec].name = NULL;
  scSpec[No_scSpec].type = NULL;
  scSpec[No_scSpec].exclude = NULL;
  if ((name &&
       !SDDS_CopyString(&scSpec[No_scSpec].name, name)) ||
      (type &&
       !SDDS_CopyString(&scSpec[No_scSpec].type, type)) ||
      (exclude &&
       !SDDS_CopyString(&scSpec[No_scSpec].exclude, exclude)))
    bombElegant((char*)"memory allocation failure", NULL);
  
  No_scSpec++;
}

void clearSCSpecs() 
{
  while (No_scSpec--) {
    if (scSpec[No_scSpec].name)
      free(scSpec[No_scSpec].name);
    if (scSpec[No_scSpec].type)
      free(scSpec[No_scSpec].type);
    if (scSpec[No_scSpec].exclude)
      free(scSpec[No_scSpec].exclude);
  }
  free(scSpec);
  scSpec = NULL;
}

long insertSCMULT(char *name, long type, long *occurrence) 
{
  long i;
  for (i=0; i<No_scSpec; i++) {
    if (scSpec[i].exclude && wild_match(name, scSpec[i].exclude))
      continue;
    if (scSpec[i].name && !wild_match(name, scSpec[i].name))
      continue;
    if (scSpec[i].type && !wild_match(entity_name[type], scSpec[i].type))
      continue;
    (*occurrence)++;
    break;
  }

  if (*occurrence < sc->nskip || sc->nskip==0)
    return(0);

  *occurrence = 0;
  return(1);
}

#include "insertSCeffects.h"
void setupSCEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  long i;
  
  if (!No_scSpec && !(sc = (SPACE_CHARGE*)SDDS_Realloc(sc, sizeof(*sc))))
    bombElegant((char*)"memory allocation failure", NULL);                

  /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&insert_sceffects, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &insert_sceffects);

  if (clear) {
    clearSCSpecs();
    if (!name && !type)
      return;
  }
  if (disable)
    return;
  
  if (!name || !strlen(name))
    bombElegant((char*)"no name given", NULL);
  str_toupper(name);
  if (has_wildcards(name) && strchr(name, '-'))
    name = expand_ranges(name);
  if (type) {
    str_toupper(type);
    if (has_wildcards(type) && strchr(type, '-'))
      type = expand_ranges(type);
    for (i=0; i<N_TYPES; i++)
      if (wild_match(entity_name[i], type))
	break;
    if (i==N_TYPES) {
      fprintf(stderr, (char*)"type pattern %s does not match any known type", type);
      exitElegant(1);
    }
  }
  if (exclude) {
    str_toupper(exclude);
    if (has_wildcards(exclude) && strchr(exclude, '-'))
      exclude = expand_ranges(exclude);
  }
  
  addSCSpec(name, type, exclude);
  cp_str(&(sc->name),element_prefix);

  sc->nskip = 0;
  sc->nonlinear = 0; 
  sc->horizontal = sc->vertical = sc->longitudinal =0;
  if (skip)
    sc->nskip = skip;

  if (vertical)
    sc->vertical = vertical;

  if (horizontal)
    sc->horizontal = horizontal;

  if (longitudinal)
    sc->longitudinal =longitudinal;

  if (nonlinear)
    sc->nonlinear = nonlinear;
}

/* track through space charge element */
void trackThroughSCMULT(double **part, long np, ELEMENT_LIST *eptr)
{
  long i;
#if USE_MPI
  long np_total;
#endif
  double *coord;
  double kx, ky, sx;
  double center[3], kick[2];
  double sigmax, sigmay;
  int flag;

  if ( !USE_MPI || !notSinglePart) {
    if (!np)  
      return;
  }
#if USE_MPI
  else {
    if(isMaster)
      np = 0;
    MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
   if (!np_total)
     return;
  }
#endif
  /* compute bunch center */
  for(i=center[0]=center[1]=center[2]=0; i<np; i++) {
    coord = part[i];
    center[0] += coord[0];
    center[1] += coord[2];
    center[2] += coord[4];
  }
#if USE_MPI
  if (USE_MPI) {
    double center_sum[3];
    MPI_Allreduce (center, center_sum, 3, MPI_DOUBLE, MPI_SUM, workers);
    center[0] = center_sum[0]/np_total;
    center[1] = center_sum[1]/np_total;
    center[2] = center_sum[2]/np_total;
  } 
#else
  center[0] /= np;
  center[1] /= np;
  center[2] /= np;
#endif
 
  /* apply kick to particles */
  if (!nonlinear) {
    for(i=0; i<np; i++) {
      coord = part[i];
      linearSCKick(coord, eptr, center);
    }
  }
  else {
    sigmax = computeRmsCoordinate(part, 0, np);
    sigmay = computeRmsCoordinate(part, 2, np);
    for(i=0; i<np; i++) {
      coord = part[i];

      if ((fabs(coord[0]-center[0])<sigmax) && (fabs(coord[2]-center[1]) < sigmay)) {
        linearSCKick(coord, eptr, center);
        continue;
      }

      if (sigmax/sigmay>0.99 && sigmax/sigmay < 1.01) {
        sx = 0.99 * sigmay;
        flag = nonlinearSCKick(coord, eptr, center, sx, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        kx = kick[0];
        ky = kick[1];
        sx = 1.01 * sigmay;
        flag = nonlinearSCKick(coord, eptr, center, sx, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        kx += kick[0];
        ky += kick[1];
        coord[1] += kx / 2.0;
        coord[3] += ky / 2.0;
      }
      else {
        flag = nonlinearSCKick(coord, eptr, center, sigmax, sigmay, kick); 
        if(!flag) {
          linearSCKick(coord, eptr, center);
          continue;
        }
        coord[1] += kick[0];
        coord[3] += kick[1];
      }
    }
  }

  sc->dmux=sc->dmuy=0.0;       /* reset space charge strength */
}

void linearSCKick(double *coord, ELEMENT_LIST *eptr, double *center)
{
  double k0, kx, ky;
  k0 = sc->c1 * exp(-sqr(coord[4]-center[2])/sqr(sc->sigmaz)/2.0);
  if (sc->horizontal) {
    kx = k0 * sc->dmux / eptr->twiss->betax;	/* From dmux to KL */
    coord[1] += kx*(coord[0]-center[0]);
  }
  if (sc->vertical) {
    ky = k0 * sc->dmuy / eptr->twiss->betay;	/* From dmuy to KL */
    coord[3] += ky*(coord[2]-center[1]);
  }
}

int nonlinearSCKick(double *coord, ELEMENT_LIST *eptr, double *center, 
                     double sigmax, double sigmay, double *kick)
{
  double k0, kx, ky, sqs;
  std::complex <double> wa, wb, w1, w2, w;
  double temp;
  long flag;

  k0 = sc->c1 * exp(-sqr(coord[4]-center[2])/sqr(sc->sigmaz)/2.0) * sqrt(PI/2.0);

  sqs = sqrt(fabs(sqr(sigmax)-sqr(sigmay))*2.0);
  kx = k0 * sc->dmux * sigmax * sqrt(sigmax+sigmay) / sqrt(fabs(sigmax-sigmay)) / eptr->twiss->betax;
  ky = k0 * sc->dmuy * sigmay * sqrt(sigmax+sigmay) / sqrt(fabs(sigmax-sigmay)) / eptr->twiss->betay;

  w1 = std::complex <double> ((coord[0]-center[0])/sqs , (coord[2]-center[1])/sqs);
  w2 = std::complex <double> ((coord[0]-center[0])*sigmay/sigmax/sqs , (coord[2]-center[1])*sigmax/sigmay/sqs);

  temp = exp((-sqr(coord[0]-center[0])/sqr(sigmax)-sqr(coord[2]-center[1])/sqr(sigmay))/2.0);

  
  wa = complexErf(w1, &flag);
  if(!flag) return(0);
  wb = complexErf(w2, &flag);
  if(!flag) return(0);

  w = wa + temp*wb;
  kick[0] = kx * w.imag();
  kick[1] = ky * w.real();
  return(1);
}

void initializeSCMULT(ELEMENT_LIST *eptr, double **part, long np, double Po, long i_pass )
{
  static CHARGE *charge;
	
  if (!eptr->twiss)
    bombElegant((char*)"Twiss parameters must be calculated before SC tracking.", NULL);
		
  if (i_pass==0) {
    while(eptr) {
      if (eptr->type==T_CHARGE) {
        charge = (CHARGE*)eptr->p_elem;
        break;
      }
      eptr =eptr->succ;
    }
    if (charge==NULL) 
      bombElegant((char*)"No charge element is given.", NULL);
	
  }
#if USE_MPI
  /* We set it as single particle case, as the particles have not been distributed 
     when the function is called, all the processors will do the same */
  notSinglePart = 0;  
#endif
  sc->sigmax = computeRmsCoordinate(part, 0, np);
  sc->sigmay = computeRmsCoordinate(part, 2, np);
  sc->sigmaz = computeRmsCoordinate(part, 4, np);
#if USE_MPI
  /* set it back to parallel execution */
  notSinglePart = 1;
#endif
  sc->c0 = sqrt(2.0/PI) * particleRadius * charge->charge / particleCharge;
  sc->c1 = sc->c0/sqr(Po)/sqrt(sqr(Po)+1.0)/sc->sigmaz;
  /*       printf("c0=%.6g, c1=%.6g, sz=%.6g\n\n", sc->c0, sc->c1, sc->sigmaz); */

  sc->dmux=sc->dmuy=0.0;
  sc->length=0.0;
}

void accumulateSCMULT(double **part, long np, ELEMENT_LIST *eptr)
{
  TWISS *twiss0;
  double dmux, dmuy, temp;
  double length;
	
  twiss0 = (eptr->pred)->twiss;
  temp = sc->sigmax + sc->sigmay;
  dmux = twiss0->betax / sc->sigmax / temp;
  dmuy = twiss0->betay / sc->sigmay / temp;
#if USE_MPI
  sc->sigmax = computeRmsCoordinate_p(part, 0, np, eptr);
  sc->sigmay = computeRmsCoordinate_p(part, 2, np, eptr);
#else
  sc->sigmax = computeRmsCoordinate(part, 0, np);
  sc->sigmay = computeRmsCoordinate(part, 2, np);
#endif
  twiss0 = eptr->twiss;
  temp = sc->sigmax + sc->sigmay;
  dmux += twiss0->betax / sc->sigmax / temp;
  dmuy += twiss0->betay / sc->sigmay / temp;
	
  length = ((DRIFT*)eptr->p_elem)->length;
  sc->dmux += dmux * length /2.0;
  sc->dmuy += dmuy * length /2.0;
}

double computeRmsCoordinate(double **coord, long i1, long np)
{
  double vrms=0.0, xc=0.0;
  long i;
#if USE_MPI
  double xc_sum=0.0, vrms_sum=0.0;
  long np_total;
#endif

  if ( !USE_MPI || !notSinglePart) {
    if (!np)
      return(0.0);
  }
#if USE_MPI
  else {
    if(isMaster)
      np = 0;
    MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, workers);
   if (!np_total)
     return(0.0);
  }
#endif

  /* compute centroids */
  for (i=xc=0; i<np; i++) {
    xc  += coord[i][i1];
  }
  /* Compute the sum of xc across all the processors */
  if ( !USE_MPI || !notSinglePart) 
    xc /= np;
#if USE_MPI
  else {
    MPI_Allreduce (&xc, &xc_sum, 1, MPI_DOUBLE, MPI_SUM, workers);
    xc = xc_sum/np_total;
  }
#endif
  for (i=vrms=0; i<np; i++) {
    vrms += sqr(coord[i][i1]-xc );
  }
  if ( !USE_MPI || !notSinglePart)   
    vrms /= np;
#if USE_MPI
  else {
    MPI_Allreduce (&vrms, &vrms_sum, 1, MPI_DOUBLE, MPI_SUM, workers);
    vrms = vrms_sum/np_total;
  }
#endif 
  return(sqrt(vrms));
}

#if USE_MPI
/* We have this new function as we need treat the parallel an serial element separately */
double computeRmsCoordinate_p(double **coord, long i1, long np, ELEMENT_LIST *eptr)
{
  double vrms=0.0, xc=0.0;
  long i, np_total;
  unsigned long classFlags = 0;

  classFlags = entity_description[eptr->type].flags;

  if (classFlags&UNIPROCESSOR) { /* serial element, only master works */
    MPI_Bcast(&np, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (!np)
      return(0.0);

    /* compute centroids */
    if (isMaster) {
      for (i=xc=0; i<np; i++) {
	xc  += coord[i][i1];
      }
      xc  /= np;
    }
    /* Broadcast the xc from master to all the slaves */
    MPI_Bcast(&xc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (isMaster) {
      for (i=vrms=0; i<np; i++) {
	vrms += sqr(coord[i][i1]-xc);
      }
      vrms /= np;    
    }
    /* Broadcast the vrms from master to all the slaves */
    MPI_Bcast(&vrms, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else { /* parallel element, only slaves works */
    double xc_sum=0.0, vrms_sum = 0.0;

    if(isMaster)
      np = 0;
    MPI_Allreduce (&np, &np_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    /* compute centroids */
    if (isSlave) {
      for (i=xc=0; i<np; i++) {
	xc  += coord[i][i1];
      }
    }
    /* Compute the sum of xc across all the processors */
    MPI_Allreduce (&xc, &xc_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    xc = xc_sum/np_total;

    if (isSlave) {
      for (i=vrms=0; i<np; i++) {
	vrms += sqr(coord[i][i1]-xc);
      }
    }
    /* Compute the sum of vrms across all the processors */
    MPI_Allreduce (&vrms, &vrms_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    vrms = vrms_sum/np_total;    
  }
  return(sqrt(vrms));
}
#endif
