/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "mdb.h"
#include "track.h"
#include "match_string.h"
#include "touschekScatter.h"
#include "SDDS.h"
#include "constants.h"
#include "chbook.h"

#define DEBUG 0
#define NDIV 10000;

/* Initialization of the simulation */
static TSCATTER_SPEC *tsSpec = NULL;
void init_TSPEC ();
TSCATTER *initTSCATTER (ELEMENT_LIST *eptr);
/* Calculate local and integrated Touschek scattering rate */
int TouschekRate(LINE_LIST *beamline);
void FIntegral(double tm, double B1, double B2, double *F, char *name);
double Fvalue (double t, double tm, double b1, double b2);
/* Monte Carlo simulation of Touschek scattering */
void TouschekDistribution(RUN *run, LINE_LIST *beamline);

void selectPart(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1);
void bunch2cm(double *p1, double *p2, double *q, double *beta, double *gamma);
void eulertrans(double *v0, double theta, double phi, double *v1, double *v);
void cm2bunch(double *p1, double *p2, double *q, double *beta, double *gamma);
double moeller(double beta0, double theta);
void pickPart(double *weight, long *index, long start, long end, 
              long *iTotal, double *wTotal, double weight_limit, double weight_ave);

char *compose_filename1(char *template, char *root_name, long index);
int str_inp(const char *string1, const char *string2);
void fill_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, 
                 book1 *dp, double *p, double weight);
void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity);

void setupTouschekEffect(NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline) 
{
  ELEMENT_LIST *eptr;
  int flag;

  /* Check if there is TScatter element along beamline. */
  eptr = &(beamline->elem);
  flag = 0;
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      flag = 1;
      break;
    }
    eptr = eptr->succ; 
  }
  if(!flag) 
    bomb("No TSCATTER element along beamline", NULL);                

   /* process the namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  process_namelist(&touschek_scatter, nltext);
  if (echoNamelists) print_namelist(stdout, &touschek_scatter);

  if (!charge)    
    bomb("charge has to be given", NULL);
  if (frequency<1)
    bomb("frequency has to >=1", NULL);
  if (!delta)
    bomb("energy aperture has to be given", NULL);

  if (!p_central_mev)
    bomb("beam energy has to be given", NULL);
  if (!emittance[0] || !emittance[1])
    bomb("bunch emittance has to be given", NULL);
  if (!sigma_dp)
    bomb("energy spread has to be given", NULL);
  if (!sigma_s)
    bomb("bunch length has to be given", NULL);

  /* Initial Touschek calculation then calculate the Twiss function. */
  init_TSPEC ();
  run_twiss_output(run, beamline, NULL, -1);
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  /* Calculate Piwinski's scattering rate. */
  if (TouschekRate(beamline))
    bomb("Touschek scattering rate calculation error", NULL);
  /* Generate scattered particles at TScatter. 
     And track the scattered particles down to beamline. */
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  TouschekDistribution(run, beamline);
#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif
  return;
}

int TouschekRate(LINE_LIST *beamline)
{
  double NP;
  double tm, B1, B2, F, rate, IntR;
  ELEMENT_LIST *eptr;

  double betagamma, gamma; 
  double sp2, sp4, beta2, betagamma2, sh2;
  double sx2, sxb2, dx2, dx_2;
  double sy2, syb2, dy2, dy_2;
  double a0, c0, c1, c2, c3;
  TWISS *twiss0;

  NP = charge/e_mks;
  gamma = p_central_mev/me_mev;
  betagamma = sqrt(sqr(gamma)-1);

  sp2 = sqr(sigma_dp);
  sp4 = sqr(sp2);
  beta2 = sqr(betagamma)/sqr(gamma);
  betagamma2 = sqr(betagamma);
  a0 = sqr(re_mks)*NP/(4*sqrt(PI)*sqr(gamma)*sigma_s);
  tm = beta2*delta*delta;

  IntR = 0.;
  rate = 0.;
  eptr = &(beamline->elem);
  while (eptr) {
    if (eptr->type == T_TSCATTER) {
      ((TSCATTER*)eptr->p_elem)->IntR = IntR;
      ((TSCATTER*)eptr->p_elem)->p_rate = rate*c_mks;
      IntR = 0.;
    }
    if(!(entity_description[eptr->type].flags&HAS_LENGTH) ||
       (eptr->pred && !(((DRIFT*)eptr->p_elem)->length)) ||
       !eptr->succ) {
      eptr = eptr->succ;
      continue;
    }

    twiss0 = eptr->twiss;
    sxb2 = twiss0->betax*emittance[0];
    dx2  = sqr(twiss0->etax);
    sx2  = sxb2 + dx2*sp2;
    dx_2 = sqr(twiss0->alphax * twiss0->etax + twiss0->betax * twiss0->etapx);

    syb2 = twiss0->betay*emittance[1];
    dy2  = sqr(twiss0->etay);
    sy2  = syb2 + dy2*sp2;
    dy_2 = sqr(twiss0->alphay * twiss0->etay + twiss0->betay * twiss0->etapy);

    sh2 = 1/(1/sp2+(dx2+dx_2)/sxb2+(dy2+dy_2)/syb2);
    c0 = sqrt(sh2)/(sigma_dp*betagamma2*emittance[0]*emittance[1]);
    c1 = sqr(twiss0->betax)/(2*betagamma2*sxb2);
    c2 = sqr(twiss0->betay)/(2*betagamma2*syb2);
    c3 = sx2*sy2-sp4*dx2*dy2;
    B1 = c1*(1-sh2*dx_2/sxb2)+c2*(1-sh2*dy_2/syb2);
    B2 = sqr(B1)-sqr(c0)*c3;   	  
    if (B2<0) {
      fprintf(stdout, "B2^2<0 at \"%s\" occurence %ld", eptr->name, eptr->occurence);
      exit(1);
    }
    B2=sqrt(B2);   	  

    FIntegral(tm, B1, B2, &F, eptr->name);

    rate = a0*c0*F*NP;
    IntR += rate * frequency * ((DRIFT*)eptr->p_elem)->length;

    eptr = eptr->succ; 
  }

  return(0);
}

void FIntegral(double tm, double b1, double b2, double *F, char *name) 
{
  double f0, f1, sum;
  double HPI, step;
  long converge=0;
  double t1, k0, k1, km;
  
  HPI = PI/2;
  step = HPI/NDIV;
  k0 = km = atan(sqrt(tm));
  f0 = Fvalue(tm, tm, b1, b2);
  sum = 0;
  f1 = 0;
  
  while (!converge) {
    k1 = k0 + step;
    t1 = sqr(tan(k1));
    if (isnan(t1)) {
      fprintf(stdout, "Integration failed at Element %s", name);
      break;
    }
    
    f1 = Fvalue(t1, tm, b1, b2);

    if (f1>0 && f1*(HPI-k1)<1e-3*sum)
      converge =1;
    
    sum +=(f1+f0)/2*step;
    k0 = k1;
    f0 = f1;    
  }

  *F = sum;

  return;
}

double Fvalue (double t, double tm, double b1, double b2)
{
  double c0, c1, c2, result;

  c0 = (sqr(2*t+1)*(t/tm/(1+t)-1)/t+t-sqrt(t*tm*(1+t))-(2+1/(2*t))*log(t/tm/(1+t)))*sqrt(1+t);
  c1 = exp(-b1*t);
  c2 = dbesi0(b2*t);
  result = c0 * c1 * c2;
  /* If overflow/underflow use approximate equation for modified bessel function. */
  if (isnan(result) || result>FLT_MAX) {
    result=c0*exp(b2*t-b1*t)/sqrt(PIx2*b2*t);
  } 
  return result;
}

char *compose_filename1(char *template, char *root_name, long index)
{
  char *root, *ext, *filename, *inx;
  long i, divp;

  divp = str_inp(template, ".");
  if (divp<0)
    bomb("touschek_scatter: The filename has to be x.x format. x is 1 to anynumber of characters", NULL);
  if (str_in(template, "%s")) {
    root = root_name;
  }
  else {
    root = tmalloc(sizeof(char)*divp);
    for (i=0; i<divp; i++) {
      *(root+i) = *(template+i);
    }
    *(root+i) = '\0';    
  }

  ext =  tmalloc(sizeof(char)*(strlen(template)-divp));
  for(i=0; i<strlen(template)-divp; i++) {
    *(ext+i) = *(template+i+divp);
  }
  *(ext+i) = '\0';

  inx = tmalloc(sizeof(char)*5);
  sprintf(inx, "%04ld", index);

  filename =  tmalloc(sizeof(char)*(strlen(root)+strlen(inx)+strlen(ext)+1));
  strcat(filename,root);
  strcat(filename,inx);
  strcat(filename,ext);
  return(filename);
}

int str_inp(const char *string1, const char *string2)
{
  int i, j, flag;

  if (strlen(string1)<strlen(string2))
    return (-1);

  i=j=0;
  for (i=0;i<strlen(string1); i++) {
    if (*(string1+i) == *string2) {
      if (strlen(string2)+i>strlen(string1))
        return (-1);
      for (j=0; j<strlen(string2); j++) {
        if (*(string1+i+j) != *(string2+j)) {
          flag = 0;
          break;
        }
        flag=1;
      }
      if(flag)
        return (i);
    }
  }

  return (-1);
}

void fill_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, 
                 book1 *dp, double *p, double weight)
{
  chfill1(x, p[0], weight);
  chfill1(y, p[1], weight);
  chfill1(s, p[2], weight);
  chfill1(xp, p[3], weight);
  chfill1(yp, p[4], weight);
  chfill1(dp, p[5], weight);

  return;
}

void print_hbook (book1 *x, book1 *y, book1 *s, book1 *xp, book1 *yp, book1 *dp,
                  TSCATTER *tsptr, char *filename, char *description, int verbosity)
{
  SDDS_DATASET outPage;
  long i;
  double *v1, *v2, *v3, *v4, *v5, *v6;

  v1 = calloc(sizeof(*v1), x->length);
  v2 = calloc(sizeof(*v2), y->length);
  v3 = calloc(sizeof(*v3), s->length);
  v4 = calloc(sizeof(*v4), xp->length);
  v5 = calloc(sizeof(*v5), yp->length);
  v6 = calloc(sizeof(*v6), dp->length);
  for(i=0; i<x->length; i++) {
    v1[i] = ((double)i+0.5)*x->dx + x->xmin;
    v2[i] = ((double)i+0.5)*y->dx + y->xmin;
    v3[i] = ((double)i+0.5)*s->dx + s->xmin;
    v4[i] = ((double)i+0.5)*xp->dx+ xp->xmin;
    v5[i] = ((double)i+0.5)*yp->dx+ yp->xmin;
    v6[i] = ((double)i+0.5)*dp->dx+ dp->xmin;
    x->value[i] *= tsptr->factor;
    y->value[i] *= tsptr->factor;
    s->value[i] *= tsptr->factor;
    xp->value[i] *= tsptr->factor;
    yp->value[i] *= tsptr->factor;
    dp->value[i] *= tsptr->factor;
  }

  /* Open file for writting */
  if (verbosity)
    fprintf( stdout, "Opening \"%s\" for writing...\n", filename);

  if (!SDDS_InitializeOutput(&outPage, SDDS_BINARY, 1, 
                             description, description, filename))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineParameter(&outPage, "Element_Name", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Element_Location", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Piwinski_IntR", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Piwinski_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Ignored_Rate", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "nbins", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Total_count", NULL, NULL, 
                             NULL, NULL, SDDS_LONG, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_1", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_1", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_1", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_2", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_2", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_2", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_3", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_3", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_3", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_4", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_4", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_4", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_5", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_5", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_5", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "VariableName_6", NULL, NULL, 
                             NULL, NULL, SDDS_STRING, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Interval_6", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL) ||
      0>SDDS_DefineParameter(&outPage, "Minimum_6", NULL, NULL, 
                             NULL, NULL, SDDS_DOUBLE, NULL))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (0>SDDS_DefineColumn(&outPage, "x", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "x_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "y", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "y_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "s", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "s_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "xp", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "xp_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "yp", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "yp_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "dp/p", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0) ||
      0>SDDS_DefineColumn(&outPage, "dp/p_rate", NULL, NULL, 
                          NULL, NULL, SDDS_DOUBLE, 0))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_WriteLayout(&outPage) )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* Write to output file */
  if (0>SDDS_StartPage(&outPage, x->length) ||
      !SDDS_SetParameters(&outPage, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "Element_Name", tsptr->name,
                          "Element_Location", tsptr->s, 
                          "VariableName_1", x->vname,
                          "Interval_1", x->dx,
                          "Minimum_1", x->xmin,
                          "VariableName_2", y->vname,
                          "Interval_2", y->dx,
                          "Minimum_2", y->xmin,
                          "VariableName_3", s->vname,
                          "Interval_3", s->dx,
                          "Minimum_3", s->xmin,
                          "VariableName_4", xp->vname,
                          "Interval_4", xp->dx,
                          "Minimum_4", xp->xmin,
                          "VariableName_5", yp->vname,
                          "Interval_5", yp->dx,
                          "Minimum_5", yp->xmin,
                          "VariableName_6", dp->vname,
                          "Interval_6", dp->dx,
                          "Minimum_6", dp->xmin,
                          "Piwinski_IntR", tsptr->IntR,
                          "Piwinski_Rate", tsptr->p_rate,
                          "Total_Rate", tsptr->s_rate,
                          "Ignored_Rate", tsptr->i_rate,
                          "nbins", x->length,
                          "Total_count", x->count, NULL) ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v1,        x->length,  "x") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, x->value,  x->length,  "x_rate") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v2,        y->length,  "y") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, y->value,  y->length,  "y_rate") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v3,        s->length,  "s") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, s->value,  s->length,  "s_rate") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v4,        xp->length, "xp") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, xp->value, xp->length, "xp_rate") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v5,        yp->length, "yp") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, yp->value, yp->length, "yp_rate") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, v6,        dp->length, "dp/p") ||
      !SDDS_SetColumn(&outPage, SDDS_SET_BY_NAME, dp->value, dp->length, "dp/p_rate") ||
      !SDDS_WritePage(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  if (!SDDS_Terminate(&outPage))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  return;
}

void TouschekDistribution(RUN *run, LINE_LIST *beamline)
{
  long i, j, total_event, n_left, eleCount=0;
  ELEMENT_LIST *eptr;
  TSCATTER *tsptr;
  double p1[6], p2[6], dens1, dens2;
  double theta, phi, qa[3], qb[3], beta[3], qabs, gamma;
  double beta0, cross, temp;
  static SDDS_TABLE SDDS_bunch, SDDS_loss;
  BEAM  Beam0, Beam, *beam0, *beam;
  double xrange[6];
  book1 *x, *y, *s, *xp, *yp, *dp;
  book1 *x0, *y0, *s0, *xp0, *yp0, *dp0;
  book1 *lossDis;
  double *weight;
  double ran1[11];
  long *index, iTotal, sTotal;
  double weight_limit, weight_ave, wTotal;

  eptr = &(beamline->elem);
  beam0 = &Beam0;
  beam = &Beam;
  sTotal = (int)beamline->revolution_length;

  while (eptr) {
    if (eptr->type == T_TSCATTER) {
        eleCount++;
        if(eleCount < i_start) {
          eptr = eptr->succ; 
          continue;
        }
        if(eleCount > i_end)
          break;
#if USE_MPI
      MPI_Barrier(MPI_COMM_WORLD); 
      if (myid==0) {
#endif
        weight = (double*)malloc(sizeof(double)*n_simulated);
        beam0->particle = (double**)czarray_2d(sizeof(double), n_simulated, 7);
        beam0->original = beam0->accepted = NULL;
        beam0->n_original = beam0->n_to_track = beam0->n_particle = n_simulated;
        beam0->n_accepted = beam0->n_saved = 0;
        beam0->p0_original = beam0->p0 = tsSpec->betagamma;
        beam0->bunchFrequency = 0.;
        beam0->lostOnPass = tmalloc(sizeof(*(beam0->lostOnPass))*beam0->n_to_track);
	
        tsptr = initTSCATTER (eptr);
        xrange[0] = 0.5*tsSpec->range[0]*tsptr->sigx;
        xrange[1] = 0.5*tsSpec->range[1]*tsptr->sigy;
        xrange[2] = 0.5*tsSpec->range[2]*tsptr->sigz;
        xrange[3] = 0.5*tsSpec->range[0]*sqrt(tsSpec->emit[0]/tsptr->twiss[0][1])*(1.+fabs(tsptr->twiss[0][0]));
        xrange[4] = 0.5*tsSpec->range[1]*sqrt(tsSpec->emit[1]/tsptr->twiss[1][1])*(1.+fabs(tsptr->twiss[1][0]));
	
        xrange[0] = xrange[0]+fabs(0.5*tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[0][0]);
        xrange[1] = xrange[1]+fabs(0.5*tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[1][0]);
        xrange[3] = xrange[3]+fabs(0.5*tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[0][1]);
        xrange[4] = xrange[4]+fabs(0.5*tsSpec->range[2]*tsSpec->sigma_p*tsptr->disp[1][1]);
        xrange[5] = 0.5*tsSpec->range[2]*tsSpec->sigma_p;
        if (initial) {
          tsptr->iniFile = compose_filename1(initial, run->rootname, eptr->occurence);
          x0 = chbook1("x", -xrange[0], xrange[0], tsSpec->nbins);
          y0 = chbook1("y", -xrange[1], xrange[1], tsSpec->nbins);
          s0 = chbook1("s", -xrange[2], xrange[2], tsSpec->nbins);
          xp0 = chbook1("xp", -xrange[3], xrange[3], tsSpec->nbins);
          yp0 = chbook1("yp", -xrange[4], xrange[4], tsSpec->nbins);
          dp0 = chbook1("dp", -xrange[5], xrange[5], tsSpec->nbins);
        }
        if (distribution) {
          tsptr->disFile = compose_filename1(distribution, run->rootname, eptr->occurence);
          x = chbook1("x", -xrange[0], xrange[0], tsSpec->nbins);
          y = chbook1("y", -xrange[1], xrange[1], tsSpec->nbins);
          s = chbook1("s", -xrange[2], xrange[2], tsSpec->nbins);
          xp = chbook1("xp", -xrange[3], xrange[3], tsSpec->nbins);
          yp = chbook1("yp", -xrange[4], xrange[4], tsSpec->nbins);
          dp = chbook1("dp", -0.1, 0.1, tsSpec->nbins);
        }
        if (output) {
          tsptr->outFile = compose_filename1(output, run->rootname, eptr->occurence);
          lossDis = chbook1("los_distribution", 0, sTotal, sTotal);
        }
        if (loss) {
          tsptr->losFile = compose_filename1(loss, run->rootname, eptr->occurence);
          SDDS_BeamScatterLossSetup(&SDDS_loss, tsptr->losFile, SDDS_BINARY, 1, 
                                    "lost particle coordinates", run->runfile,
                                    run->lattice, "touschek_scatter");
        }
        if (bunch) {
          tsptr->bunFile = compose_filename1(bunch, run->rootname, eptr->occurence);
          SDDS_BeamScatterSetup(&SDDS_bunch, tsptr->bunFile, SDDS_BINARY, 1, 
                                "scattered-beam phase space", run->runfile,
                                run->lattice, "touschek_scatter");
        }
        
        i = 0; j=0;total_event=0;
        while(1) {
          if(i>=n_simulated)
            break;
          /* Select the 11 random number then mix them. Use elegant run_setup seed */
          for (j=0; j<11; j++) {
            ran1[j] = random_1_elegant(1);
          }
          randomizeOrder((char*)ran1, sizeof(ran1[0]), 11, 0, random_4);
          
          total_event++;
          selectPart(tsptr, p1, p2, &dens1, &dens2, ran1);
          if (initial) {
            fill_hbook (x0, y0, s0, xp0, yp0, dp0, p1, dens1);
            fill_hbook (x0, y0, s0, xp0, yp0, dp0, p2, dens2);
          }
         /* This is very important. Change from slop to MeV */
          for(j=3; j<5; j++) {
            p1[j] *= tsSpec->pCentral;
            p2[j] *= tsSpec->pCentral;
          }
          p1[5] = (p1[5]+1)*tsSpec->pCentral;
          p2[5] = (p2[5]+1)*tsSpec->pCentral;
	  
          bunch2cm(p1,p2,qa,beta,&gamma);
          
          theta = (ran1[9]*0.998+0.001)*PI;
          phi = ran1[10]*PI;
	  
          temp = dens1*dens2*sin(theta);
          eulertrans(qa,theta,phi,qb,&qabs);
          cm2bunch(p1,p2,qb,beta,&gamma);
          p1[5] -= tsSpec->pCentral;
          p2[5] -= tsSpec->pCentral;
	  
          if(fabs(p1[5])>tsSpec->dp0 || fabs(p2[5])>tsSpec->dp0) {
            beta0=qabs/sqrt(qabs*qabs+me_mev*me_mev);
            cross = moeller(beta0,theta);
            temp *= cross*beta0/gamma/gamma;
	    
            if(fabs(p1[5])>tsSpec->dp0) {
              tsptr->totalWeight += temp;
              p1[3] /= tsSpec->pCentral;
              p1[4] /= tsSpec->pCentral;
              p1[5] /= tsSpec->pCentral;
              if (distribution) {
                fill_hbook (x, y, s, xp, yp, dp, p1, temp);
              }
              tsptr->simuCount++;
	      
              beam0->particle[i][0] = p1[0];
              beam0->particle[i][1] = p1[3];
              beam0->particle[i][2] = p1[1];
              beam0->particle[i][3] = p1[4];
              beam0->particle[i][4] = p1[2];
              beam0->particle[i][5] = p1[5];
              beam0->particle[i][6] = i+1;
              weight[i] = temp;
              i++;
            }
	    
            if(i>=n_simulated)
              break;
	    
            if(fabs(p2[5])>tsSpec->dp0) {
              tsptr->totalWeight += temp;
              p2[3] /= tsSpec->pCentral;
              p2[4] /= tsSpec->pCentral;
              p2[5] /= tsSpec->pCentral;
              if (distribution) {
                fill_hbook (x, y, s, xp, yp, dp, p2, temp);
              }
              tsptr->simuCount++;
	      
              beam0->particle[i][0] = p2[0];
              beam0->particle[i][1] = p2[3];
              beam0->particle[i][2] = p2[1];
              beam0->particle[i][3] = p2[4];
              beam0->particle[i][4] = p2[2];
              beam0->particle[i][5] = p2[5];
              beam0->particle[i][6] = i+1;
              weight[i] = temp;
              i++;
            }
          }
	  
          if (total_event*11 > (long)2e9)  {
            fprintf(stdout, "warning: The total random number used > 2e9. Use less n_simulated or use small delta");
            fflush(stdout);
            break;
          }
        }
        if (total_event/tsptr->simuCount > 20) {
          if (distribution_cutoff[0]<5 || distribution_cutoff[1]<5 ) 
            fprintf(stdout, "waring: Scattering rate is low, please use 5 sigma beam for better simulation.\n");
          else
            fprintf(stdout, "waring: Scattering rate is very low, please ignore the rate from Monte Carlo simulation. Use Piwinski's rate only\n"); 
        }
        tsptr->factor = tsptr->factor / (double)(total_event);
        tsptr->s_rate = tsptr->totalWeight * tsptr->factor;
        tsptr->i_rate = tsptr->s_rate*tsSpec->ignoredPortion;
	
        /* Pick tracking particles from the simulated scattered particles */
        index = (long*)malloc(sizeof(long)*tsptr->simuCount);
        for (i=0; i<tsptr->simuCount; i++) index[i]=i;
        iTotal = 0;
        wTotal =0.;
        weight_limit = tsptr->totalWeight*(1-tsSpec->ignoredPortion);
        weight_ave = tsptr->totalWeight/tsptr->simuCount;
	
        pickPart(weight, index, 0, tsptr->simuCount,  
                 &iTotal, &wTotal, weight_limit, weight_ave);
	
        beam->particle = (double**)czarray_2d(sizeof(double), iTotal, 7);
        beam->original = (double**)czarray_2d(sizeof(double), iTotal, 7);
        beam->accepted = NULL;
        beam->n_original = beam->n_to_track = beam->n_particle = iTotal;
        beam->n_accepted = beam->n_saved = 0;
        beam->p0_original = beam->p0 = tsSpec->betagamma;
        beam->bunchFrequency = 0.;
        beam->lostOnPass = tmalloc(sizeof(*(beam->lostOnPass))*beam->n_to_track);
	
        for (i=0; i<iTotal; i++) {
          beam->original[i][0] = beam->particle[i][0] = beam0->particle[index[i]][0];
          beam->original[i][1] = beam->particle[i][1] = beam0->particle[index[i]][1];
          beam->original[i][2] = beam->particle[i][2] = beam0->particle[index[i]][2];
          beam->original[i][3] = beam->particle[i][3] = beam0->particle[index[i]][3];
          beam->original[i][4] = beam->particle[i][4] = beam0->particle[index[i]][4];
          beam->original[i][5] = beam->particle[i][5] = beam0->particle[index[i]][5];
          beam->original[i][6] = beam->particle[i][6] = i+1;
          weight[i] *= tsptr->factor;
        }
        if (bunch)
          dump_scattered_particles(&SDDS_bunch, beam->particle, (long)iTotal,
                                   weight, tsptr);
        if (distribution) {
          print_hbook(x, y, s, xp, yp, dp, tsptr, tsptr->disFile, 
                      "Distribution of Scattered particles", 1);
          free_hbook1(x);free_hbook1(y);free_hbook1(s);
          free_hbook1(xp);free_hbook1(yp);free_hbook1(dp);
        }
        if (initial) {
          print_hbook(x0, y0, s0, xp0, yp0, dp0, tsptr, tsptr->iniFile, 
                      "Distribution of simulated particles", 1);
          free_hbook1(x0);free_hbook1(y0);free_hbook1(s0);
          free_hbook1(xp0);free_hbook1(yp0);free_hbook1(dp0);
        }
#if USE_MPI
      }
      MPI_Barrier(MPI_COMM_WORLD); 
#endif
      if (do_track) {
        n_left = do_tracking(beam, NULL, (long)iTotal, NULL, beamline, 
                             &tsSpec->betagamma, NULL, NULL, NULL, NULL, run, 1,
                             0, 1, 0, NULL,
                             NULL, NULL, beam->lostOnPass, eptr);
#if USE_MPI
        MPI_Barrier(MPI_COMM_WORLD); 
        if (myid==0) {
#endif
          if (loss) {
            dump_scattered_loss_particles(&SDDS_loss, beam->particle+n_left, beam->original,  
                                          beam->lostOnPass+n_left, beam->n_to_track-n_left, weight, tsptr);
            if (!SDDS_Terminate(&SDDS_loss)) {
              SDDS_SetError("Problem terminating 'losses' file (finish_output)");
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
          }
          if (output) {
            for (i=0; i< beam->n_to_track-n_left; i++) {
              j = (beam->particle+n_left)[i][6]-1;
              chfill1(lossDis, (beam->particle+n_left)[i][4], weight[j]*tsptr->IntR/tsptr->s_rate);
            }
            chprint1(lossDis, tsptr->outFile, "Beam loss distribution in particles/s/m", verbosity);
	    free_hbook1(lossDis);
          }
#if USE_MPI
        }
#endif
      }
#if USE_MPI
      if (myid==0) {
#endif
        free_beamdata(beam);
        free_beamdata(beam0);
        free(weight);
#if USE_MPI
      }
#endif
    }
    eptr = eptr->succ; 
  }
  return;
}

/* Initialize beam parameter at each Scatter element */
TSCATTER *initTSCATTER (ELEMENT_LIST *eptr)
{
  TSCATTER *tsptr;

  tsptr = ((TSCATTER*)eptr->p_elem);

  tsptr->twiss[0][0] = eptr->twiss->alphax;
  tsptr->twiss[0][1] = eptr->twiss->betax;
  tsptr->twiss[0][2] = (1.+ sqr(eptr->twiss->alphax))/eptr->twiss->betax;
  tsptr->twiss[1][0] = eptr->twiss->alphay;
  tsptr->twiss[1][1] = eptr->twiss->betay;
  tsptr->twiss[1][2] = (1.+ sqr(eptr->twiss->alphay))/eptr->twiss->betay;
  tsptr->twiss[2][0] = 0.0;
  tsptr->twiss[2][1] = tsSpec->sigz/tsSpec->sigma_p;
  tsptr->twiss[2][2] = tsSpec->sigma_p/tsSpec->sigz;
  tsptr->disp[0][0] = eptr->twiss->etax;
  tsptr->disp[0][1] = eptr->twiss->etapx;
  tsptr->disp[1][0] = eptr->twiss->etay;
  tsptr->disp[1][1] = eptr->twiss->etapy;

  tsptr->sigx = sqrt(tsSpec->emit[0]*tsptr->twiss[0][1]);
  tsptr->sigy = sqrt(tsSpec->emit[1]*tsptr->twiss[1][1]);
  tsptr->sigz = tsSpec->sigz;
  tsptr->sigxyz = tsptr->sigx * tsptr->sigy * tsptr->sigz;

  tsptr->factor = pow(tsSpec->charge/e_mks, 2.0)*pow(re_mks, 2.0)/4.*c_mks
    *pow(PI, 2.0)/pow(PI, 6.0)/8./8.
    *pow(tsSpec->range[0], 3.0)*pow(tsSpec->range[1], 3.0)*pow(tsSpec->range[2], 3.0)
    /tsptr->sigxyz;

  tsptr->s_rate = tsptr->i_rate = tsptr->totalWeight = 0;
  tsptr->simuCount = 0;
  tsptr->name = eptr->name;
  tsptr->s = eptr->end_pos;
  tsptr->betagamma = tsSpec->betagamma;
  return (tsptr);
}

/* Initialize Touschek Scattering Calculation */
void init_TSPEC ()
{
  if (!(tsSpec = SDDS_Malloc(sizeof(*tsSpec))))
    bomb("memory allocation failure at setup Touscheck scatter", NULL);                

  tsSpec->ebeam = p_central_mev;
  tsSpec->delta = delta;
  tsSpec->sigma_p = sigma_dp;
  tsSpec->sigz = sigma_s;
  tsSpec->charge = charge;
  tsSpec->emit[0] = emittance[0];
  tsSpec->emit[1] = emittance[1];
  tsSpec->emit[2] = sigma_s * sigma_dp;

  tsSpec->gamma = tsSpec->ebeam/me_mev;
  tsSpec->pCentral = sqrt(sqr(tsSpec->ebeam) - me_mev*me_mev);  /*beam momentum in MeV */
  tsSpec->betagamma = tsSpec->pCentral/me_mev;
  tsSpec->dp0 = delta*tsSpec->pCentral;
  tsSpec->range[0] = 2 * distribution_cutoff[0];
  tsSpec->range[1] = 2 * distribution_cutoff[1];
  tsSpec->range[2] = 2 * distribution_cutoff[2];
  tsSpec->nbins = nbins;
  tsSpec->ignoredPortion = ignored_portion;

  return;
}

/************************************************************************************\
 * Modified from S. Khan's code.                                                    *
 *      select two electrons (labelled a and b) for Touschek scattering             *
 * p1[x,y,s,xp,yp,dp/p]                                                             *
\************************************************************************************/

void selectPart(TSCATTER *tsptr, double *p1, double *p2, 
                double *dens1, double *dens2, double *ran1)
{
  int i;
  double U[3], V1[3], V2[3], densa[3], densb[3];

  /* Select random particle's parameter in normal phase space */
  for (i=0; i<3; i++) {
    U[i]  = (ran1[i]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emit[i]);
    V1[i] = (ran1[i+3]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emit[i]);
    V2[i] = (ran1[i+6]-0.5) * tsSpec->range[i]*sqrt(tsSpec->emit[i]);
    densa[i] = exp(-0.5*(U[i]*U[i]+V1[i]*V1[i])/tsSpec->emit[i]);
    if (i==2)
      densb[i] = exp(-0.5*(U[i]*U[i]+V2[i]*V2[i])/tsSpec->emit[i]);
  }
  /* Change particle's parameter from normal phase space to real phase space */
  for (i=0; i<3; i++) {
    p1[i] = p2[i] = sqrt(tsptr->twiss[i][1])*U[i];
    p1[i+3] = (V1[i] - tsptr->twiss[i][0]*U[i])/sqrt(tsptr->twiss[i][1]);
    p2[i+3] = (V2[i] - tsptr->twiss[i][0]*U[i])/sqrt(tsptr->twiss[i][1]);
  }
  /* Dispersion correction */
  p1[0] = p1[0] + p1[5]*tsptr->disp[0][0];
  p1[1] = p1[1] + p1[5]*tsptr->disp[1][0];
  p1[3] = p1[3] + p1[5]*tsptr->disp[0][1];
  p1[4] = p1[4] + p1[5]*tsptr->disp[1][1];

  p2[0] = p1[0] - p2[5]*tsptr->disp[0][0];
  p2[1] = p1[1] - p2[5]*tsptr->disp[1][0];
  U[0] = p2[0]/sqrt(tsptr->twiss[0][1]);
  U[1] = p2[1]/sqrt(tsptr->twiss[1][1]);
  p2[3] = (V2[0] - tsptr->twiss[0][0]*U[0])/sqrt(tsptr->twiss[0][1]);
  p2[4] = (V2[1] - tsptr->twiss[1][0]*U[1])/sqrt(tsptr->twiss[1][1]);
  densb[0] = exp(-0.5*(U[0]*U[0]+V2[0]*V2[0])/tsSpec->emit[0]);
  densb[1] = exp(-0.5*(U[1]*U[1]+V2[1]*V2[1])/tsSpec->emit[1]);

  p2[0] = p1[0];
  p2[1] = p1[1];
  p2[3] = p2[3] + p2[5]*tsptr->disp[0][1];
  p2[4] = p2[4] + p2[5]*tsptr->disp[1][1];

  *dens1 = densa[0] * densa[1] * densa[2]; 
  *dens2 = densb[0] * densb[1] * densb[2]; 

  return;
}

/* This function transfer particle's momentum from lab system to c.o.m system */
/* p1[3]=p1x*c, p1[4]=p1y*c, p1[5]=p1z*c in MeV */
/* p2[3]=p2x*c, p2[4]=p2y*c, p2[5]=p1z*c in MeV */
void bunch2cm(double *p1, double *p2, double *q, double *beta, double *gamma)
{
  double pp1, pp2, e1, e2, ee;
  int i;
  double bb, betap1, factor;

  pp1=0.0;
  pp2=0.0;
  for(i=3; i<6; i++) {
    pp1 = pp1 + sqr(p1[i]);
    pp2 = pp2 + sqr(p2[i]);
  }
  e1=sqrt(me_mev*me_mev+pp1);
  e2=sqrt(me_mev*me_mev+pp2);
  ee=e1+e2;

  betap1=0.0;
  bb=0.0;
  for(i=0; i<3; i++) {
    beta[i]=(p1[i+3]+p2[i+3])/ee;
    betap1=betap1+beta[i]*p1[i+3];
    bb=bb+beta[i]*beta[i];
  }

  *gamma = 1./sqrt(1.-bb);
  factor = ((*gamma)-1.)*betap1/bb;

  for(i=0; i<3; i++) {
    q[i]=p1[i+3]+factor*beta[i]-(*gamma)*e1*beta[i];
  }

  return;
}

/* Rotate scattered p in c.o.m system */
void eulertrans(double *v0, double theta, double phi, double *v1, double *v)
{
  double th, ph, s1, s2, c1, c2;
  double x0, y0, z0; 

  *v=sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
  th=acos(v0[2]/(*v));
  ph=atan2(v0[1],v0[0]);

  s1=sin(th);
  s2=sin(ph);
  c1=cos(th);
  c2=cos(ph);

  x0=cos(theta);
  y0=sin(theta)*cos(phi);
  z0=sin(theta)*sin(phi);

  v1[0] = (*v) * (s1*c2*x0 - s2*y0 - c1*c2*z0);
  v1[1] = (*v) * (s1*s2*x0 + c2*y0 - c1*s2*z0);
  v1[2] = (*v) * (   c1*x0         +    s1*z0);

  return;
}

void cm2bunch(double *p1, double *p2, double *q, double *beta, double *gamma)
{
  int i;
  double pq, e, betaq, bb, factor;

  pq=0.0;
  for(i=0; i<3; i++) {
    pq = pq + q[i]*q[i];
  }

  e=sqrt(me_mev*me_mev+pq);

  betaq=0.0;
  bb=0.0;
  for(i=0; i<3; i++) {
    betaq = betaq + beta[i]*q[i];
    bb = bb + beta[i]*beta[i];
  }

  factor=((*gamma)-1)*betaq/bb;
  for(i=0; i<3; i++) {
    p1[i+3]= q[i] + (*gamma)*beta[i]*e + factor*beta[i];
    p2[i+3]=-q[i] + (*gamma)*beta[i]*e - factor*beta[i];
  }

  return;
}

double moeller(double beta0, double theta)
{
  double cross; 
  double beta2, st2;

  beta2=beta0*beta0;
  st2=sqr(sin(theta));

  cross = (1.-beta2)*(sqr(1.+1./beta2)*(4./st2/st2-3./st2)+1.+4./st2);
 
  return cross;
}

void pickPart(double *weight, long *index, long start, long end, 
              long *iTotal, double *wTotal, double weight_limit, double weight_ave)
{
  long i, i1, i2, N;
  double w1, w2;
  long *index1, *index2;
  double *weight1, *weight2;

  i1=i2=0;
  w1=w2=0.;
  N = end-start;
  if(N<3) return;  /* scattered particles normally appear in pair */
  index2 = (long*)malloc(sizeof(long)*N);
  weight2 = (double*)malloc(sizeof(double)*N);
  index1 = (long*)malloc(sizeof(long)*N);
  weight1 = (double*)malloc(sizeof(double)*N);
  
  for (i=start; i<end; i++) {
    if (weight[i] > weight_ave) {
      weight2[i2] = weight[i];
      index2[i2++] = index[i];
      w2 += weight[i];
    } else {
      weight1[i1] = weight[i];
      index1[i1++] = index[i];
      w1 += weight[i];
    }
  }
  if ((w2+ (*wTotal)) > weight_limit) {
    weight_ave = w2/(double)i2;
    for (i=0; i<i2; i++) {
      index[start+i]=index2[i];
      weight[start+i]=weight2[i];
    }
    free(weight1);
    free(index1);
    free(weight2);
    free(index2);
    pickPart(weight, index, start, start+i2,
             iTotal, wTotal, weight_limit, weight_ave);
    return;
  }

  *iTotal += i2;
  *wTotal += w2;
  weight_ave = w1/(double)i1;
  for (i=0; i<i2; i++) {
    index[start+i]=index2[i];
    weight[start+i]=weight2[i];
  }
  for (i=0; i<i1; i++) {
    index[start+i2+i] = index1[i];
    weight[start+i2+i] = weight1[i];
  }
  free(weight1);
  free(index1);
  free(weight2);
  free(index2);
  pickPart(weight, index, i2+start, end,
           iTotal, wTotal, weight_limit, weight_ave);
  return;
}
