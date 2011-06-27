/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/


#include "mdb.h"
#include "madto.h"
#include <ctype.h>

void emitCosyDrift(FILE *fp, char *name, double length)
{
  fprintf(fp, "  PROCEDURE %s ;\n", name);
  fprintf(fp, "    DL %.15g ;\n", length);
  fprintf(fp, "  ENDPROCEDURE ;\n");
}

void emitCosyDipole(FILE *fp, char *name,
                    double length, double angle, 
                    double e1, double e2,
                    double h1, double h2,
                    double K1, double K2, double K3,
                    double hgap, double fint)
{
  double rho, engeCoef[3];
  long i, j;
  double sign = 1;
  
  if (angle<0)
    sign = -1;
  rho = length/(fabs(angle)+1e-300),
  fprintf(fp, "  PROCEDURE %s ;\n", name);
  fprintf(fp, "    VARIABLE W 1 4 4 ;\n");
  fprintf(fp, "    W(1,1) := 1 ;\n");
  fprintf(fp, "    W(2,1) := %.15g ;\n", K1*rho);
  fprintf(fp, "    W(3,1) := %.15g ;\n", K2*rho/2);
  fprintf(fp, "    W(4,1) := %.15g ;\n", K3*rho/6);
  for (i=1; i<=4; i++) 
    for (j=2; j<=4; j++)
      fprintf(fp, "    W(%ld, %ld) := 0 ;\n", i, j);

  if (angle<0) 
    fprintf(fp, "    CB ;\n");
  if (hgap!=0 && fint!=0) {
    computeEngeCoefficients(engeCoef, rho, length, 2*hgap, fint);
    fprintf(fp, "    FR 3 ;\n");
    fprintf(fp, "    FC 1 1 1 %.15g %.15g %.15g\n       0.0 0.0 0.0 ; \n",
            engeCoef[0], engeCoef[1], engeCoef[2]);
    fprintf(fp, "    FC 1 2 1 %.15g %.15g %.15g\n       0.0 0.0 0.0 ; \n",
            engeCoef[0], engeCoef[1], engeCoef[2]);
  } else {
    fprintf(fp, "    FR 1 ;\n");
  }
  fprintf(fp, "    MSS %.15g %.15g %.15g %.15g\n       %.15g %.15g %.15g W; \n",
          rho, 180/PI*fabs(angle), 
          hgap==0?1e-3:hgap, 
          sign*180/PI*e1, h1, 
          sign*180/PI*e2, h2);
  if (angle<0) 
    fprintf(fp, "    CB ;\n");
  fprintf(fp, "    FD ; \n");
  fprintf(fp, "  ENDPROCEDURE; \n");
}


void convert_to_cosy(char *outputfile, LINE_LIST *beamline, 
                     long cosyOrder, double pCentral,
                     double quadBoreRadius, double sextBoreRadius)
{
    ELEMENT_LIST *eptr;
    QUAD  *quad; KQUAD *kquad;
    SEXT  *sext; KSEXT *ksext;
    BEND  *bend;
    HCOR  *hcor;
    VCOR  *vcor;
    DRIFT *drift;
    CSBEND *csbend;
    CSRCSBEND *csrbend;
    CSRDRIFT *csrdrift;
    NIBEND *nibend;
    RFCA *rfca;
    FILE *fp;
    double BRho;
    long i;
    char *ptr;
#define NVAR 6
    char *varName[NVAR] = {
      "A", "B", "G", "R", "MU", "F"
      };
    long varLength[NVAR] = {
      2, 2, 2, 2, 2, 6
      };
    char **nameUsed = NULL;
    long namesUsed = 0;

    fp = fopen_e(outputfile, "w", 0);

    fprintf(fp, "INCLUDE 'COSY';\n");
    fprintf(fp, "PROCEDURE RUN ;\n");
    for (i=0; i<NVAR; i++)
      fprintf(fp, "    VARIABLE %s 100 %ld ; \n", varName[i], varLength[i]);
      
    BRho = pCentral*particleMass*c_mks/particleCharge;

    /* emit procedures describing elements */
    eptr = &(beamline->elem);
    while (eptr) {
      while ((ptr=strchr(eptr->name, ':')))
        *ptr = '_';
      for (i=0; i<namesUsed; i++) {
        if (strcmp(eptr->name, nameUsed[i])==0)
          break;
      }
      if (i!=namesUsed) {
        eptr = eptr->succ;
        continue;
      }
      nameUsed = trealloc(nameUsed, sizeof(*nameUsed)*(namesUsed+1));
      nameUsed[namesUsed++] = eptr->name;
      switch (eptr->type) {
      case T_QUAD:
        quad = (QUAD*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K1 = %e, R = %e, BRho = %e }\n",
                quad->k1, quadBoreRadius, BRho);
        fprintf(fp, "    MQ %.15g %.15g %.15g ; \n",
                quad->length, -quad->k1*BRho*quadBoreRadius, quadBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_KQUAD:
        kquad = (KQUAD*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K1 = %e, R = %e, BRho = %e }\n",
                kquad->k1, quadBoreRadius, BRho);
        fprintf(fp, "    MQ %.15g %.15g %.15g ; \n",
                kquad->length, -kquad->k1*BRho*quadBoreRadius, quadBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_SBEN:
        bend = (BEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, bend->length, bend->angle,
                       bend->e1, bend->e2, bend->h1, bend->h2,
                       bend->k1, bend->k2, 0.0,
                       bend->hgap, bend->fint);
        break;
      case T_NIBEND:
        nibend = (NIBEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, nibend->length, nibend->angle,
                       nibend->e1, nibend->e2, 0.0, 0.0,
                       0.0, 0.0, 0.0,
                       nibend->hgap, nibend->fint);
        break;
      case T_CSBEND:
        csbend = (CSBEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, csbend->length, csbend->angle,
                       csbend->e1, csbend->e2, csbend->h1, csbend->h2,
                       csbend->k1, csbend->k2, csbend->k3,
                       csbend->hgap, csbend->fint);
        break;
      case T_CSRCSBEND:
        csrbend = (CSRCSBEND*)eptr->p_elem;
        emitCosyDipole(fp, eptr->name, csrbend->length, csrbend->angle,
                       csrbend->e1, csrbend->e2, csrbend->h1, csrbend->h2,
                       csrbend->k1, csrbend->k2, csrbend->k3,
                       csrbend->hgap, csrbend->fint);
        break;
      case T_DRIF:
        drift = (DRIFT*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, drift->length);
        break;
      case T_CSRDRIFT:
        csrdrift = (CSRDRIFT*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, csrdrift->length);
        break;
      case T_SEXT:
        sext = (SEXT*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K2 = %e, R = %e, BRho = %e }\n",
                sext->k2, sextBoreRadius, BRho);
        fprintf(fp, "    MSK %.15g %.15g %.15g ;\n",
                sext->length, sext->k2*BRho*sextBoreRadius, sextBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_KSEXT:
        ksext = (KSEXT*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    { K2 = %e, R = %e, BRho = %e }\n",
                ksext->k2, sextBoreRadius, BRho);
        fprintf(fp, "    MSK %.15g %.15g %.15g ;\n",
                ksext->length, ksext->k2*BRho*sextBoreRadius, sextBoreRadius);
        fprintf(fp, "  ENDPROCEDURE; \n");
        break;
      case T_HCOR:
        hcor = (HCOR*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, hcor->length);
        break;
      case T_VCOR:
        vcor = (VCOR*)eptr->p_elem;
        emitCosyDrift(fp, eptr->name, vcor->length);
        break;
      case T_RFCA:
        rfca = (RFCA*)eptr->p_elem;
        fprintf(fp, "  PROCEDURE %s ;\n", eptr->name);
        fprintf(fp, "    VARIABLE V 1 2 2 ;\n");
        fprintf(fp, "    V(1,1) := %.15g ;\n", rfca->volt/1e3);
        fprintf(fp, "    V(1,2) := 0.0 ;\n");
        fprintf(fp, "    V(2,1) := 0.0 ;\n");
        fprintf(fp, "    V(2,2) := 0.0 ;\n");
        if (rfca->length)
          fprintf(fp, "    DL %.15g ;\n", rfca->length/2);
        fprintf(fp, "    RF V 1 %.15g %.15g 0.1 ;\n", rfca->freq, rfca->phase);
        if (rfca->length)
          fprintf(fp, "    DL %.15g ;\n", rfca->length/2);
        fprintf(fp, "  ENDPROCEDURE ;\n");
        break;
      default:
        emitCosyDrift(fp, eptr->name, 0.0);
        break;
      }
      eptr = eptr->succ;
    }

    /* emit beamline definition */
    fprintf(fp, "  PROCEDURE MACH ;");
    eptr = &(beamline->elem);
    i = 0;
    while (eptr) {
      if (i%4==0)
        fprintf(fp, "\n    %s ;", eptr->name);
      else
        fprintf(fp, " %s ;", eptr->name);
      i++;
      eptr = eptr->succ;
    }
    fprintf(fp, "\n  ENDPROCEDURE ;\n");
    fprintf(fp, "    OPENF 7 '%s.map' 'new' ;\n", outputfile);
    fprintf(fp, "    OPENF 8 '%s.tunes' 'new' ;\n", outputfile);
    fprintf(fp, "    OV %ld 2 1 ;\n", cosyOrder);
    fprintf(fp, "    RPE %.15g*PARA(1) ;\n", (sqrt(pCentral*pCentral+1)-1)*particleMass*sqr(c_mks)/particleCharge/1e6);
    fprintf(fp, "    UM ; MACH; \n");
    fprintf(fp, "    PT 7; CLOSEF 7;\n");
    fprintf(fp, "    TP MU; WRITE 8 ' DELTA-DEPENDENT TUNES: '   MU(1) MU(2) ;\n");
    fprintf(fp, "    GT MAP F MU A B G R;\n");
    fprintf(fp, "    WRITE 8 ' DELTA-DEPENDENT FIXED POINT ' F(1) F(2) F(3) F(4) ;\n");
    fprintf(fp, "    WRITE 8 ' DELTA-DEPENDENT ALPHAS ' A(1) A(2) ;\n");
    fprintf(fp, "    WRITE 8 ' DELTA-DEPENDENT BETAS  ' B(1) B(2) ;\n");
    fprintf(fp, "    WRITE 8 ' DELTA-DEPENDENT GAMMAS ' G(1) G(2) ;\n");
    fprintf(fp, "    CLOSEF 8 ;\n");
    fprintf(fp, "ENDPROCEDURE ;\n");
    fprintf(fp, "RUN ;\n");
    fprintf(fp, "END ;\n");
    }

