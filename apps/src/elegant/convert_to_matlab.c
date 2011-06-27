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

void convert_to_EmmaMatlab(char *outputfile, LINE_LIST *beamline, char *header_file, char *ender_file)
{
    ELEMENT_LIST *eptr;
    QUAD  *quad;
    SEXT  *sext;
    BEND  *bend;
    HCOR  *hcor;
    VCOR  *vcor;
    DRIFT *drift;
    CSBEND *csbend;
    CSRCSBEND *csrbend;
    CSRDRIFT *csrdrift;
    FILE *fpi, *fp;
    double length, angle, k1, E1, E2;
    long slices;
    char s[256];
    

    fp = fopen_e(outputfile, "w", 0);

    if (header_file) {
        fpi = fopen_e(header_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    eptr = &(beamline->elem);
    
    while (eptr) {
      length = angle = k1 = E1 = E2 = slices = 0;
      switch (eptr->type) {
      case T_QUAD:
        quad = (QUAD*)eptr->p_elem;
        k1 = quad->k1;
        length = quad->length;
        break;
      case T_SBEN:
        bend = (BEND*)eptr->p_elem;
        angle = bend->angle;
        E1 = bend->e1;
        E2 = bend->e2;
        length = bend->length;
        k1 = bend->k1;
        break;
      case T_CSBEND:
        csbend = (CSBEND*)eptr->p_elem;
        angle = csbend->angle;
        E1 = csbend->e1;
        E2 = csbend->e2;
        length = csbend->length;
        k1 = csbend->k1;
        slices = csbend->n_kicks;
        break;
      case T_CSRCSBEND:
        csrbend = (CSRCSBEND*)eptr->p_elem;
        angle = csrbend->angle;
        E1 = csrbend->e1;
        E2 = csrbend->e2;
        length = csrbend->length;
        k1 = csrbend->k1;
        slices = csrbend->n_kicks;
        break;
      case T_DRIF:
        drift = (DRIFT*)eptr->p_elem;
        length = drift->length;
        break;
      case T_CSRDRIFT:
        csrdrift = (CSRDRIFT*)eptr->p_elem;
        length = csrdrift->length;
        break;
      case T_SEXT:
        sext = (SEXT*)eptr->p_elem;
        length = sext->length;
        break;
      case T_HCOR:
        hcor = (HCOR*)eptr->p_elem;
        length = hcor->length;
        break;
      case T_VCOR:
        vcor = (VCOR*)eptr->p_elem;
        length = vcor->length;
        break;
      default:
        fprintf(stderr, "warning: entity type %s not implemented\n",
               entity_name[eptr->type]);
        break;
      }
      fprintf(fp, "%21.15g %21.15g %21.15g %21.15g %21.15g %ld 0\n",
              length, angle, k1, E1, E2, slices);
      eptr = eptr->succ;
    }

    if (ender_file) {
        fpi = fopen_e(ender_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }
    }

