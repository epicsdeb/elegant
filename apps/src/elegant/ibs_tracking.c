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
#include "zibs.h"
#include "SDDS.h"

static double tmp_safe_sqrt;
#define SAFE_SQRT(x) ((tmp_safe_sqrt=(x))<0?0.0:sqrt(tmp_safe_sqrt))
static TWISS *twiss0;

void init_IBS(ELEMENT_LIST *element);
void free_IBS(IBSCATTER *IBS);
void inflateEmittance(double **coord, double Po, 
                      long offset, long istart, long iend, long *index, double factor);
void inflateEmittanceZ(double **coord, double Po, long isRing, double dt,
                       long istart, long iend, long *index, double zRate[3]);
void SDDS_IBScatterSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                         char *command_file, char *lattice_file, char *caller, long isRing);
void dump_IBScatter(SDDS_TABLE *SDDS_table, IBSCATTER *IBS, long pass);
void reset_IBS_output(ELEMENT_LIST *element);

void slicebeam(double **coord, long np, double Po, long nslice, long *index, long *count, double *dt);
void zeroslice (long islice, IBSCATTER *IBS);
long computeSliceParameters(double C[6], double S[6][6], double **part, long *index, long start, long end, double Po);
void forth_propagate_twiss(IBSCATTER *IBS, long islice, double betax0, double alphax0, 
                           double betay0, double alphay0, RUN *run);
void copy_twiss(TWISS *tp0, TWISS *tp1);

void track_IBS(double **coord, long np, IBSCATTER *IBS, double Po, 
               ELEMENT_LIST *element, CHARGE *charge, long i_pass, long n_passes, RUN *run)
{
  long *index, *count;
  long istart, iend, ipart, icoord, ihcoord, islice;
  double aveCoord[6], S[6][6];
  double betax0, alphax0, betay0, alphay0;
  double randomNumber;
  double beta0, beta1, p;
  double RNSigma[3], RNSigmaCheck[3]={0,0,0};
  double bLength, tLength, zRate[3];
  static SDDS_TABLE outPage;
  static long isInit=0, doOut=0;
  double eta[4];

  if (IBS->nslice<1) 
    bombElegant("NSLICE has to be an integer >= 1", NULL);
  if (charge)
    IBS->charge = charge->macroParticleCharge*np;
  if (!IBS->charge)
    bombElegant("bunch charge is not given", NULL);
  if (IBS->isRing && IBS->nslice>1)
    bombElegant("no slice valid for ring beam. NSLICE has to be 1", NULL);

  if (!IBS->s)
    init_IBS(element);
  if (IBS->dT == 0) return;

  eta[0] = IBS->etax[IBS->elements-1];
  eta[1] = IBS->etaxp[IBS->elements-1];
  eta[2] = IBS->etay[IBS->elements-1];
  eta[3] = IBS->etayp[IBS->elements-1];

  index = (long*)malloc(sizeof(long)*np);
  count = (long*)malloc(sizeof(long)*IBS->nslice);
  slicebeam(coord, np, Po, IBS->nslice, index, count, &tLength);
  bLength = IBS->revolutionLength/IBS->dT*tLength;
  if ((IBS->nslice == 1) && (!IBS->isRing))
    bLength /=sqrt(2*PI);

  iend = 0;
  for (islice=0; islice<IBS->nslice; islice++) {
    istart = iend;
    iend += count[islice];

    if (count[islice]<10) {
      fprintf(stdout, "count=%ld, warning: too few particles inside slice #%ld. No IBS taking into account in this slice.\n", count[islice], islice+1);
      zeroslice (islice, IBS);
      continue;
    }

    computeSliceParameters(aveCoord, S, coord, index, istart, iend, Po);
    IBS->emitx0[islice] = correctedEmittance(S, eta, 0, 1, &betax0, &alphax0);
    IBS->emity0[islice] = correctedEmittance(S, eta, 2, 3, &betay0, &alphay0);
    IBS->emitl0[islice] = SAFE_SQRT(S[4][4]*S[5][5]-sqr(S[4][5]));
    IBS->sigmaDelta0[islice] = sqrt(S[5][5]);
    if (IBS->isRing)
      bLength = IBS->revolutionLength/IBS->dT*sqrt(S[4][4]);
    IBS->sigmaz0[islice] = bLength;
    
    if (!IBS->forceMatchedTwiss)
      forth_propagate_twiss(IBS, islice, betax0, alphax0, betay0, alphay0, run);
    IBS->icharge[islice] = IBS->charge * (double)count[islice] / (double)np;
    IBSRate (fabs(IBS->icharge[islice]/particleCharge), 
             IBS->elements, 1, 0, IBS->isRing,
             IBS->emitx0[islice], IBS->emity0[islice], IBS->sigmaDelta0[islice], bLength,
             IBS->s, IBS->pCentral, IBS->betax[islice], IBS->alphax[islice], IBS->betay[islice], IBS->alphay[islice],
             IBS->etax, IBS->etaxp, IBS->etay, IBS->etayp, 
             IBS->xRateVsS[islice], IBS->yRateVsS[islice], IBS->zRateVsS[islice],
             &(IBS->xGrowthRate[islice]), &(IBS->yGrowthRate[islice]), &(IBS->zGrowthRate[islice]), 1);    
    IBS->xGrowthRate[islice] *= IBS->factor;
    IBS->yGrowthRate[islice] *= IBS->factor;
    IBS->zGrowthRate[islice] *= IBS->factor;
  }

  iend = 0;
  for (islice=0; islice<IBS->nslice; islice++) {
    istart = iend;
    iend += count[islice];

    zRate[1] = 1.+IBS->dT*IBS->zGrowthRate[islice];
    if (islice == 0)
      zRate[0] = zRate[1];
    else 
      zRate[0] = 1.+IBS->dT*IBS->zGrowthRate[islice-1];
    if (islice == IBS->nslice-1)
      zRate[2] = zRate[1];
    else
      zRate[2] = 1.+IBS->dT*IBS->zGrowthRate[islice+1];

    RNSigma[0] = RNSigma[1] = RNSigma[2] = RNSigma[3] = 0;
    if (!IBS->smooth) {
      computeSliceParameters(aveCoord, S, coord, index, istart, iend, Po);
      if (IBS->do_x)
        RNSigma[0] = sqrt(fabs(sqr(1 + IBS->dT * IBS->xGrowthRate[islice])-1))*sqrt(S[1][1]);
      if (IBS->do_y)
        RNSigma[1] = sqrt(fabs(sqr(1 + IBS->dT * IBS->yGrowthRate[islice])-1))*sqrt(S[3][3]);
      if (IBS->do_z) {
        RNSigma[2] = sqrt(fabs(sqr(1 + IBS->dT * IBS->zGrowthRate[islice])-1))*sqrt(S[5][5]);
      }
      for (icoord=1, ihcoord=0; icoord<6; icoord+=2, ihcoord++) {
        if (RNSigma[ihcoord]) {
          RNSigmaCheck[ihcoord] = 0;
          if (icoord!=5) {
            for (ipart=istart; ipart<iend; ipart++) {
              randomNumber = gauss_rn_lim(0.0, RNSigma[ihcoord], 3.0, random_2);
              coord[index[ipart]][icoord] += randomNumber;
              RNSigmaCheck[ihcoord] += sqr(randomNumber);
            }
          } else {
            for (ipart=istart; ipart<iend; ipart++) {
              randomNumber = gauss_rn_lim(0.0, RNSigma[ihcoord], 3.0, random_2);
              p = Po*(1+coord[index[ipart]][5]);
              beta0 = p/sqrt(p*p+1);
              coord[index[ipart]][5] += randomNumber;
              p = Po*(1+coord[index[ipart]][5]);
              beta1 = p/sqrt(p*p+1);
              coord[index[ipart]][4] *= beta1/beta0;
              RNSigmaCheck[ihcoord] += sqr(randomNumber);
            }
          }
          RNSigmaCheck[ihcoord] = sqrt(RNSigmaCheck[ihcoord]/(double)(iend-istart)/S[icoord][icoord]+1.);
        }
      }
      /*
      fprintf(stdout,"s=%g,islice=%ld,istart=%ld,iend=%ld,checkz=%g\n", 
              IBS->s[IBS->elements-1], islice, istart, iend, RNSigmaCheck[2]);
      */
    } else {
      /* inflate each emittance by the prescribed factor */
      inflateEmittance(coord, Po, 0, istart, iend, index, (1.+IBS->dT*IBS->xGrowthRate[islice]));
      inflateEmittance(coord, Po, 2, istart, iend, index, (1.+IBS->dT*IBS->yGrowthRate[islice]));
      inflateEmittanceZ(coord, Po, IBS->isRing, tLength, istart, iend, index, zRate);
    }

    /* update beam emittance information after IBS scatter for IBSCATTER */
    if (!IBS->isRing) {
      computeSliceParameters(aveCoord, S, coord, index, istart, iend, Po);
      IBS->emitx[islice] = correctedEmittance(S, eta, 0, 1, 0, 0);
      IBS->emity[islice] = correctedEmittance(S, eta, 2, 3, 0, 0);
      IBS->emitl[islice] = SAFE_SQRT(S[4][4]*S[5][5]-sqr(S[4][5]));
      IBS->sigmaDelta[islice] = sqrt(S[5][5]);
      IBS->sigmaz[islice] = sqrt(S[4][4]);
    }
  }
  free(index);
  free(count);

  if(IBS->filename) {
    if (!isInit) {
      SDDS_IBScatterSetup(&outPage, IBS->filename, SDDS_BINARY, 1, "IBS scatter growth rate output", 
                          run->runfile, run->lattice, "ibs_tracking", IBS->isRing);
      isInit = 1;
    }

    if ((int)i_pass/IBS->interval > doOut) {
      doOut++;
      reset_IBS_output(element);
    }
    if (IBS->output) {
      dump_IBScatter(&outPage, IBS, i_pass);
      IBS->output = 0;
    }
  }
  if (i_pass==n_passes-1)
    free_IBS(IBS);
  return;
}

void inflateEmittance(double **coord, double Po, 
                      long offset, long istart, long iend, long *index, double factor)
{
  long ipart, np;
  double factorSqrt, c[2]={0,0};

  np = iend - istart;
  if (!np)
    return;
  factorSqrt = sqrt(factor);

  for (ipart=istart; ipart<iend; ipart++) {
    c[0] += coord[index[ipart]][offset+0];
    c[1] += coord[index[ipart]][offset+1];
  }
  c[0] /= np;
  c[1] /= np;
  for (ipart=istart; ipart<iend; ipart++) {
    coord[index[ipart]][offset+0] = (coord[index[ipart]][offset+0]-c[0])*factorSqrt+c[0];
    coord[index[ipart]][offset+1] = (coord[index[ipart]][offset+1]-c[1])*factorSqrt+c[1];
  }
}

void inflateEmittanceZ(double **coord, double Po, long isRing, double dt,
                       long istart, long iend, long *index, double zRate[3])
{
  long i, ipart, np;
  double c0, tc, dpc, *time, p, beta0, beta1;

  np = iend - istart;
  if (!np)
    return;
 
  time = tmalloc(sizeof(*time)*np);
  for (i=dpc=tc=0; i<np; i++) {
    ipart = i + istart;
    dpc += coord[index[ipart]][5];
    p = Po*(1+coord[index[ipart]][5]);
    beta0 = p/sqrt(p*p+1);
    time[i] = coord[index[ipart]][4]/(beta0*c_mks);
    tc += time[i];
  }
  tc /= np;
  dpc /= np;

  for (i=0; i<np; i++) {
    ipart = i + istart;

    if (isRing) {
      c0 = sqrt(zRate[1]);
      time[i] = (time[i]-tc)*c0+tc;
      coord[index[ipart]][5] = (coord[index[ipart]][5]-dpc)*c0+dpc;
      p = Po*(1+coord[index[ipart]][5]);
      beta0 = p/sqrt(p*p+1);
      coord[ipart][4] = time[i]*beta0*c_mks;
    } else {
      if (time[i]>tc)
        c0 = (time[i]-tc)/dt*(zRate[2]-zRate[1])+zRate[1];
      else
        c0 = (time[i]-tc)/dt*(zRate[1]-zRate[0])+zRate[1];
      c0 = sqrt(c0);
      p = Po*(1+coord[index[ipart]][5]);
      beta0 = p/sqrt(p*p+1);
      coord[index[ipart]][5] = (coord[index[ipart]][5]-dpc)*c0+dpc;
      p = Po*(1+coord[index[ipart]][5]);
      beta1 = p/sqrt(p*p+1);
      coord[index[ipart]][4] = coord[index[ipart]][4]*beta1/beta0;
    }
  }
  free(time);
  return;
}

/* Set twiss parameter arrays etc. */
void init_IBS(ELEMENT_LIST *element)
{
  long count, nElements, i, j, isRing = 0, init=1;
  double startRingPos, finalPos;
  double s0, s1, dt, delta_s, p0, gamma;
  ELEMENT_LIST *element0, *elementStartRing, *eptr=NULL;
  IBSCATTER *IBS=NULL;
  TWISS *tp;
  
  if (!element->twiss) 
    bombElegant("Twiss parameters must be calculated befor IBS tracking.", NULL);

  /* Find out start point of ring */
  element0 = elementStartRing = element;
  startRingPos = 0;

  count = 0;
  while (element) {
    if (element->type==T_RECIRC) {
      startRingPos = element->end_pos;
      elementStartRing = element->succ;
      count++;
      break;
    }
    count++;
    element = element->succ;
  }

  if (elementStartRing!=element0) {
    element = elementStartRing;
  } else {
    element = element0;
    count = 0;
  }

  nElements =0;
  s1 = startRingPos;
  dt = delta_s = 0.;
  while (element) {
    s0 = s1;
    s1 = element->end_pos;
    if (s1 > s0) {
      if (element->pred)
        p0 = (element->Pref_output + element->pred->Pref_output)/2.;
      else
        p0 = element->Pref_output;
      gamma = sqrt(p0*p0+1);
      dt += (s1-s0)*gamma/p0/c_mks;
      delta_s += (s1-s0);
    }

    if (element->type==T_IBSCATTER) {
      IBS = (IBSCATTER*)element->p_elem; 
      IBS->revolutionLength = delta_s;
      if (nElements<2) 
        bombElegant("you need at least 2 other elements between IBSCATTERS, check twiss file for clue", NULL);
      IBS->dT = dt;
      dt = delta_s = 0.;
      IBS->elements = nElements;
      IBS->offset = count;
      IBS->output = 1;
      count = count + nElements +1;
      isRing = IBS->isRing;
      if (!(IBS->name = SDDS_Calloc(nElements,  sizeof(*(IBS->name)))) ||
          !(IBS->s = SDDS_Calloc(nElements, sizeof(*(IBS->s)))) ||
          !(IBS->pCentral = SDDS_Calloc(nElements, sizeof(*(IBS->pCentral)))))
        bombElegant("memory allocation failure in init_IBS", NULL);
      if (!(IBS->etax = SDDS_Realloc(IBS->etax, sizeof(*(IBS->etax))*nElements)) ||
          !(IBS->etaxp = SDDS_Realloc(IBS->etaxp, sizeof(*(IBS->etaxp))*nElements)) ||
          !(IBS->etay = SDDS_Realloc(IBS->etay, sizeof(*(IBS->etay))*nElements)) ||
          !(IBS->etayp = SDDS_Realloc(IBS->etayp, sizeof(*(IBS->etayp))*nElements)))
        bombElegant("memory allocation failure in init_IBS", NULL);
      if (!(IBS->icharge = SDDS_Realloc(IBS->icharge, sizeof(*(IBS->icharge))*IBS->nslice)) ||
          !(IBS->emitx0 = SDDS_Realloc(IBS->emitx0, sizeof(*(IBS->emitx0))*IBS->nslice)) ||
          !(IBS->emity0 = SDDS_Realloc(IBS->emity0, sizeof(*(IBS->emity0))*IBS->nslice)) ||
          !(IBS->emitl0 = SDDS_Realloc(IBS->emitl0, sizeof(*(IBS->emitl0))*IBS->nslice)) ||
          !(IBS->sigmaz0 = SDDS_Realloc(IBS->sigmaz0, sizeof(*(IBS->sigmaz0))*IBS->nslice)) ||
          !(IBS->sigmaDelta0 = SDDS_Realloc(IBS->sigmaDelta0, sizeof(*(IBS->sigmaDelta0))*IBS->nslice)) ||
          !(IBS->emitx = SDDS_Realloc(IBS->emitx, sizeof(*(IBS->emitx))*IBS->nslice)) ||
          !(IBS->emity = SDDS_Realloc(IBS->emity, sizeof(*(IBS->emity))*IBS->nslice)) ||
          !(IBS->emitl = SDDS_Realloc(IBS->emitl, sizeof(*(IBS->emitl))*IBS->nslice)) ||
          !(IBS->sigmaz = SDDS_Realloc(IBS->sigmaz, sizeof(*(IBS->sigmaz))*IBS->nslice)) ||
          !(IBS->sigmaDelta = SDDS_Realloc(IBS->sigmaDelta, sizeof(*(IBS->sigmaDelta))*IBS->nslice)) ||
          !(IBS->xGrowthRate = SDDS_Realloc(IBS->xGrowthRate, sizeof(*(IBS->xGrowthRate))*IBS->nslice)) ||
          !(IBS->yGrowthRate = SDDS_Realloc(IBS->yGrowthRate, sizeof(*(IBS->yGrowthRate))*IBS->nslice)) ||
          !(IBS->zGrowthRate = SDDS_Realloc(IBS->zGrowthRate, sizeof(*(IBS->zGrowthRate))*IBS->nslice)))
        bombElegant("memory allocation failure in init_IBS", NULL);
      if (!(IBS->betax = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->alphax = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->betay = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->alphay = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->xRateVsS = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->yRateVsS = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)) ||
          !(IBS->zRateVsS = (double**)czarray_2d(sizeof(double), IBS->nslice, nElements)))
        bombElegant("memory allocation failure in init_IBS", NULL);

      for (i=0; i<IBS->elements; i++) {
        cp_str(&IBS->name[i], elementStartRing->name);
        IBS->s[i] = elementStartRing->end_pos;
        IBS->pCentral[i] = elementStartRing->Pref_output;
        IBS->etax[i] = elementStartRing->twiss->etax;
        IBS->etaxp[i] = elementStartRing->twiss->etapx;
        IBS->etay[i] = elementStartRing->twiss->etay;
        IBS->etayp[i] = elementStartRing->twiss->etapy;
        if (init)
          twiss0 = tmalloc(sizeof(*twiss0)*IBS->nslice);
        /*
          twiss0 = SDDS_Realloc(twiss0, sizeof(*twiss0)*IBS->nslice);
        */
        for (j=0; j<IBS->nslice; j++) {
          if (init) {
            tp = &(twiss0[j]);
            copy_twiss(tp, elementStartRing->twiss);
          }
          IBS->betax[j][i] = elementStartRing->twiss->betax;
          IBS->alphax[j][i] = elementStartRing->twiss->alphax;   
          IBS->betay[j][i] = elementStartRing->twiss->betay;   
          IBS->alphay[j][i] = elementStartRing->twiss->alphay;   
        }
        init = 0;
        elementStartRing = elementStartRing->succ;
      }
      IBS->elem = tmalloc(sizeof(*(IBS->elem)));
      IBS->elem->pred = IBS->elem->succ = NULL;
      elementStartRing = elementStartRing->pred;
      for (i=IBS->elements; i>0; i--) {
        add_element(IBS->elem, elementStartRing);
        eptr = IBS->elem->succ;
        /* copy input energy to newly added element. Necessary if beamline contains RF cavity */
        eptr->Pref_input = elementStartRing->Pref_input;
        eptr->Pref_output = elementStartRing->Pref_output;
        elementStartRing = elementStartRing->pred;
      }
      IBS->elem = eptr;
      eptr->pred = eptr->pred->succ = NULL;
      nElements = -1;
      elementStartRing = element->succ;
    }
    nElements ++;
    finalPos = element->end_pos;
    element = element->succ;
  }
  
  if (isRing)
    if (finalPos != IBS->s[IBS->elements-1])
      bombElegant("You must have IBSCATTER at the end of the RING", NULL);
  element = element0;
  return;
}

void free_IBS(IBSCATTER *IBS)
{
  free(IBS->name); free(IBS->s); free(IBS->pCentral); free(IBS->icharge);
  free(IBS->etax); free(IBS->etaxp); free(IBS->etay); free(IBS->etayp);
  free(IBS->emitx0); free(IBS->emity0); free(IBS->emitl0);
  free(IBS->emitx);  free(IBS->emity);  free(IBS->emitl);
  free(IBS->sigmaz0); free(IBS->sigmaDelta0);
  free(IBS->sigmaz);  free(IBS->sigmaDelta);
  free(IBS->xGrowthRate); free(IBS->yGrowthRate); free(IBS->zGrowthRate); 
  free_czarray_2d((void**)IBS->betax, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->betay, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->alphax, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->alphay, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->xRateVsS, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->yRateVsS, IBS->nslice, IBS->elements);
  free_czarray_2d((void**)IBS->zRateVsS, IBS->nslice, IBS->elements);
  free_elements1(IBS->elem);
  return;
}

void slicebeam(double **coord, long np, double Po, long nslice, long *index, long *count, double *dt)
{
  long i, j, islice, total;
  double tMaxAll, tMinAll, tMin, tMax;
  double *time, P, beta;

  for (i=0; i<np; i++)
    index[i] = -1;
  for (islice=0; islice<nslice; islice++)
    count[islice] = 0;

  time = tmalloc(sizeof(*time)*np);
  for (i=0; i<np; i++) {
    P = Po*(1+coord[i][5]);
    beta = P/sqrt(P*P+1);
    time[i] = coord[i][4]/(beta*c_mks);
  }

  /* find limits of bunch longitudinal coordinates */
  tMin = tMaxAll = -DBL_MAX;
  tMax = tMinAll =  DBL_MAX;

  for (i=0; i<np; i++) {
    if (tMinAll>time[i])
      tMinAll = time[i];
    if (tMaxAll<time[i])
      tMaxAll = time[i];
  }
  *dt = (tMaxAll-tMinAll)/(double)nslice;
  
  tMin = tMinAll;
  j = total= 0;
  for (islice=0; islice<nslice; islice++) {
    tMax = tMin + (*dt);
    if (islice == nslice-1)
      tMax = tMaxAll + 0.1*(*dt);
    for (i=0; i<np; i++) {
      if (time[i]>=tMin && time[i]<tMax) {
        count[islice]++;
        index[j++] = i;
      }
    }
    total += count[islice];
    tMin = tMax;
  }
  if (total !=np)
    bombElegant("ibs: slice-beam, total is not equal to np. Report it to code developer", NULL);
  free(time);
  return;
}

void zeroslice (long islice, IBSCATTER *IBS) {
  long i;
  IBS->emitx0[islice] = IBS->emity0[islice] = IBS->emitl0[islice] =0.;
  IBS->sigmaz0[islice] = IBS->sigmaDelta0[islice] = 0.;
  IBS->emitx[islice] = IBS->emity[islice] = IBS->emitl[islice] =0.;
  IBS->sigmaz[islice] = IBS->sigmaDelta[islice] = 0.;
  IBS->xGrowthRate[islice] = IBS->yGrowthRate[islice] = IBS->zGrowthRate[islice] = 0.;
  for (i=0; i<IBS->elements; i++) {
    IBS->xRateVsS[islice][i] = IBS->yRateVsS[islice][i] = IBS->zRateVsS[islice][i] = 0.;
  }
  return;
}

long computeSliceParameters(double C[6], double S[6][6], double **part, long *index, long start, long end, double Po)
{
  long i, j, k, i1;
  double *time, dt, dp, beta, P;

  time = tmalloc(sizeof(*time)*(end-start));

  for (j=0; j<6; j++) {
    C[j] = 0;
    for (k=0; k<6; k++)
      S[j][k] = 0;
  }

  /* compute centroid slice parameters */
  for (i=start; i<end; i++) {
    for (j=0; j<4; j++) {
      C[j] += part[index[i]][j];
    }
    P = Po*(1+part[index[i]][5]);
    beta = P/sqrt(P*P+1);
    i1 = i-start;
    time[i1] = part[index[i]][4]/(beta*c_mks);
    C[4] += time[i1];
    C[5] += part[index[i]][5];
  }

  for (j=0; j<6; j++)
    C[j] /= (double)(end-start);

  for (i=start; i<end; i++) {
    for (j=0; j<4; j++)
      for (k=0; k<=j; k++)
        S[j][k] += (part[index[i]][j]-C[j])*(part[index[i]][k]-C[k]);
    i1 = i-start;
    S[4][4] += sqr(dt = time[i1]          - C[4]);
    S[5][5] += sqr(dp = part[index[i]][5] - C[5]);
    S[4][5] += dt*dp;
  }
  for (j=0; j<6; j++)
    for (k=0; k<=j; k++) {
      S[j][k] /= (double)(end-start); 
      S[k][j] = S[j][k];
    }
  free(time);
  return 0;
}

void forth_propagate_twiss(IBSCATTER *IBS, long islice, double betax0, double alphax0, 
                           double betay0, double alphay0, RUN *run)
{
  ELEMENT_LIST *elem;
  double tune[2];
  long i, waists[2];
  TWISS *tp;

  tp = &(twiss0[islice]);
  propagate_twiss_parameters(tp, tune, waists, NULL, IBS->elem, run, NULL, NULL);
  elem = IBS->elem;
  for (i=0; i<IBS->elements; i++) {
    IBS->betax[islice][i] = elem->twiss->betax;
    IBS->alphax[islice][i] = elem->twiss->alphax;   
    IBS->betay[islice][i] = elem->twiss->betay;   
    IBS->alphay[islice][i] = elem->twiss->alphay;   
    elem = elem->succ;
  }
 
  tp->betax = betax0;
  tp->alphax = alphax0;
  tp->phix = 0;
  tp->etax = IBS->etax[0];
  tp->etapx = IBS->etaxp[0];
  tp->betay = betay0;
  tp->alphay = alphay0;
  tp->phiy = 0;
  tp->etay = IBS->etay[0];
  tp->etapy = IBS->etayp[0];
  return;
}

void reset_IBS_output(ELEMENT_LIST *element)
{
  ELEMENT_LIST *element0;
  IBSCATTER *IBS;

  element0 = element;
  while (element) {
    if (element->type==T_IBSCATTER) {
      IBS = (IBSCATTER*)element->p_elem; 
      IBS->output = 1;
    }
    element = element->succ;
  }
  element = element0;
}

#define IBSCATTER_RING_PARAMETERS 19
#define IBSCATTER_LINAC_PARAMETERS 25
static SDDS_DEFINITION ibscatter_print_parameter[IBSCATTER_LINAC_PARAMETERS] = {
  {"TotalSlice", "&parameter name=TotalSlice, type=long, description=\"Total number of Slices\" &end"},
  {"NSlice", "&parameter name=NSlice, type=long, description=\"The ith slice of the beam\" &end"},
  {"Charge", "&parameter name=Charge, type=double, units=\"nC\", description=\"Slice charge in nC\" &end"},
  {"Particles", "&parameter name=Particles, type=double, description=\"Number of particles in slice\" &end"},
  {"s", "&parameter name=s, type=double, units=\"m\", description=\"IBScatter element location in beamline\" &end"},
  {"Pass", "&parameter name=Pass, type=long, description=\"Pass number\" &end"},
  {"StartPoint", "&parameter name=StartPoint, type=long, description=\"IStart point in the beamline\" &end"},
  {"Elements", "&parameter name=Elements, type=long, description=\"Number of elements to integrate\" &end"},
  {"ds", "&parameter name=ds, type=double, units=\"m\", description=\"Integrated length for IBS rate\" &end"},  
  {"dt", "&parameter name=dt, type=double, units=\"m\", description=\"Integrated time for IBS scattering\" &end"},  
  {"xGrowthRate", "&parameter name=xGrowthRate, symbol=\"g$bIBS,x$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the horizontal plane\" &end"},
  {"yGrowthRate", "&parameter name=yGrowthRate, symbol=\"g$bIBS,y$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the vertical plane\" &end"},
  {"zGrowthRate", "&parameter name=zGrowthRate, symbol=\"g$bIBS,z$n\", units=\"1/s\", type=double, description=\"Accumulated IBS emittance growth rate in the longitudinal plane\" &end"},
  {"enx0", "&parameter name=enx0, symbol=\"$gge$r$bx$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized initial horizontal emittance\" &end"},
  {"eny0", "&parameter name=eny0, symbol=\"$gge$r$by$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized initial vertical emittance\" &end"},
  {"emitx0", "&parameter name=emitx0, symbol=\"$ge$r$bx,Input$n\", units=\"$gp$rm\", type=double, description=\"Initial horizontal emittance\" &end"},
  {"emity0", "&parameter name=emity0, symbol=\"$ge$r$by,Input$n\", units=\"$gp$rm\", type=double, description=\"Initial vertical emittance\" &end"},
  {"sigmaDelta0", "&parameter name=sigmaDelta0, symbol=\"$gs$r$bd,Input$n\", type=double, description=\"Initial momentum spread\" &end"},
  {"sigmaz0", "&parameter name=sigmaz0, symbol=\"$gs$r$bz,Input$n\", units=m, type=double, description=\"Initial bunch length\" &end"},
  {"enx", "&parameter name=enx, symbol=\"$gge$r$bx$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized horizontal emittance with IBS\" &end"},
  {"eny", "&parameter name=eny, symbol=\"$gge$r$by$n\", units=\"m$be$nc $gp$rm\", type=double, description=\"Normalized vertical emittance with IBS\" &end"},
  {"emitx", "&parameter name=emitx, symbol=\"$ge$r$bx$n\", units=\"$gp$rm\", type=double, description=\"Horizontal emittance with IBS\" &end"},
  {"emity", "&parameter name=emity, symbol=\"$ge$r$by$n\", units=\"$gp$rm\", type=double, description=\"Vertical emittance with IBS\" &end"},
  {"sigmaDelta", "&parameter name=sigmaDelta, symbol=\"$gs$r$bd$n\", type=double, description=\"Momentum spread with IBS\" &end"},
  {"sigmaz", "&parameter name=sigmaz, symbol=\"$gs$r$bz$n\", units=m, type=double, description=\"Bunch length with IBS\" &end"},
};
#define IBSCATTER_COLUMNS 13
static SDDS_DEFINITION ibscatter_print_column[IBSCATTER_COLUMNS] = {
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"dIBSRatex", "&column name=dIBSRatex, units=\"1/(m s)\", type=double, description=\"Horizontal IBS Emittance Growth Rate\" &end"},
    {"dIBSRatey", "&column name=dIBSRatey, units=\"1/(m s)\", type=double, description=\"Vertical IBS Emittance Growth Rate\" &end"},
    {"dIBSRatez", "&column name=dIBSRatez, units=\"1/(m s)\", type=double, description=\"Longitudinal IBS Emittance Growth Rate\" &end"},
    {"betax", "&column name=betax, type=double, units=m, symbol=\"$gb$r$bx$n\", description=\"Horizontal beta-function\" &end"},
    {"alphax", "&column name=alphax, type=double, symbol=\"$ga$r$bx$n\", description=\"Horizontal alpha-function\" &end"},
    {"etax", "&column name=etax, type=double, units=m, symbol=\"$gc$r$bx$n\", description=\"Horizontal dispersion\" &end"},
    {"etaxp", "&column name=etaxp, type=double, symbol=\"$gc$r$bx$n$a'$n\", description=\"Slope of horizontal dispersion\" &end"},
    {"betay", "&column name=betay, type=double, units=m, symbol=\"$gb$r$by$n\", description=\"Vertical beta-function\" &end"},
    {"alphay", "&column name=alphay, type=double, symbol=\"$ga$r$by$n\", description=\"Vertical alpha-function\" &end"},
    {"etay", "&column name=etay, type=double, units=m, symbol=\"$gc$r$by$n\", description=\"Vertical dispersion\" &end"},
    {"etayp", "&column name=etayp, type=double, symbol=\"$gc$r$by$n$a'$n\", description=\"Slope of vertical dispersion\" &end"},
};

void SDDS_IBScatterSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                         char *command_file, char *lattice_file, char *caller, long isRing)
{
    log_entry("SDDS_IBScatterSetup");

    if (isRing)
      SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            ibscatter_print_parameter, IBSCATTER_RING_PARAMETERS, ibscatter_print_column, IBSCATTER_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    else
      SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            ibscatter_print_parameter, IBSCATTER_LINAC_PARAMETERS, ibscatter_print_column, IBSCATTER_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    log_exit("SDDS_IBScatterSetup");
}

void dump_IBScatter(SDDS_TABLE *SDDS_table, IBSCATTER *IBS, long pass)
{
  long i, islice;
  double gamma;

  log_entry("dump_IBScatter");

  if (!IBS->elements)
    return;

  gamma = sqrt(ipow(IBS->pCentral[IBS->elements-1], 2)+1.);
  for (islice=0; islice<IBS->nslice; islice++) {
    if (!SDDS_StartTable(SDDS_table, IBS->elements)) {
      SDDS_SetError("Problem starting SDDS table (dump_IBScatter)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    for (i=0; i<IBS->elements; i++) {
      if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                             0, IBS->name[i], 1, IBS->s[i], 
                             2, IBS->xRateVsS[islice][i], 3, IBS->yRateVsS[islice][i], 4, IBS->zRateVsS[islice][i], 
                             5, IBS->betax[islice][i], 6, IBS->alphax[islice][i], 7, IBS->etax[i], 8, IBS->etaxp[i], 
                             9, IBS->betay[islice][i], 10, IBS->alphay[islice][i], 11, IBS->etay[i], 12, IBS->etayp[i], -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_IBScatter)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Particles", fabs(IBS->icharge[islice]/particleCharge), NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Charge", IBS->icharge[islice]*1e9, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "TotalSlice", IBS->nslice, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "NSlice", islice+1, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "s", IBS->s[IBS->elements-1], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Pass", pass, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "StartPoint", IBS->offset, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Elements", IBS->elements, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "ds", IBS->revolutionLength, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "dt", IBS->dT, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "xGrowthRate", IBS->xGrowthRate[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "yGrowthRate", IBS->yGrowthRate[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "zGrowthRate", IBS->zGrowthRate[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "enx0", IBS->emitx0[islice]*gamma, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "eny0", IBS->emity0[islice]*gamma, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emitx0", IBS->emitx0[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emity0", IBS->emity0[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaDelta0", IBS->sigmaDelta0[islice], NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaz0", IBS->sigmaz0[islice], NULL))){
      SDDS_SetError("Problem setting SDDS parameters (dump_IBScatter)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    if (!IBS->isRing) {
      if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emitx", IBS->emitx[islice], NULL))||
          (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "emity", IBS->emity[islice], NULL))||
          (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "enx", IBS->emitx[islice]*gamma, NULL))||
          (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "eny", IBS->emity[islice]*gamma, NULL))||
          (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaDelta", IBS->sigmaDelta[islice], NULL))||
          (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "sigmaz", IBS->sigmaz[islice], NULL))){
        SDDS_SetError("Problem setting SDDS parameters (dump_IBScatter)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }

    if (!SDDS_WriteTable(SDDS_table)) {
      SDDS_SetError("Problem writing SDDS table (dump_IBScatter)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } 
    
    SDDS_UpdatePage(SDDS_table, 0);
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);
  }
  log_exit("dump_IBScatter");    
}

void copy_twiss(TWISS *tp0, TWISS *tp1)
{
  long i;
  tp0->betax = tp1->betax;
  tp0->alphax = tp1->alphax;
  tp0->phix = tp1->phix;
  tp0->etax = tp1->etax;
  tp0->etapx = tp1->etapx;
  tp0->apx = tp1->apx;
  tp0->betay = tp1->betay;
  tp0->alphay = tp1->alphay;
  tp0->phiy = tp1->phiy;
  tp0->etay = tp1->etay;
  tp0->etapy = tp1->etapy;
  tp0->apy = tp1->apy;
  tp0->Cx = tp1->Cx;
  tp0->Cy = tp1->Cy;
  tp0->periodic = tp1->periodic;
  for (i=0; i<6; i++)
    tp0->dI[i] = tp1->dI[i]; 
  return;
}

