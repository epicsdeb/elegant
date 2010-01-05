/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: convert_to_xorbit
 * purpose: convert parsed mad lattice to xorbit format 
 * 
 * Michael Borland, 1994
 */
#include "mdb.h"
#include "madto.h"
#include <ctype.h>

#define XORBIT_MIN_LENGTH 2e-6

char *xorbit_name(char *name, char extra_character);
void make_xorbit_friendly(ELEMENT_LIST *elem);
void bracket_with_drifts(ELEMENT_LIST *elem, double length);

void convert_to_xorbit(char *outputfile, LINE_LIST *beamline, long flip_k, 
                    char *header_file, char *ender_file)
{
    ELEMENT_LIST *eptr;
    QUAD  *quad;
    SEXT  *sext;
    BEND  *bend;
    HCOR  *hcor;
    VCOR  *vcor;
    DRIFT *drift;
    MULT *mult;
    char s[100];
    int quad_sign;
    FILE *fpi, *fp;
    double length;

    fp = fopen_e(outputfile, "w", 0);

    quad_sign = flip_k?-1:1;

    if (header_file) {
        fpi = fopen_e(header_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    fputs("name     tn      a         b       slen       xbar      xsig      ybar      ysig\n", fp);
    fputs("a------- ai f---------f---------f---------f---------f---------f---------f---------\n", fp);
    fputs("INJ      i\n", fp);


    eptr = &(beamline->elem);
    make_xorbit_friendly(eptr);
    length = 0;
    while (eptr) {
        switch (eptr->type) {
          case T_QUAD:
            quad = (QUAD*)eptr->p_elem;
            fprintf(fp, "%-8s q %10.6f%10.6f%10.6f\n",
                    eptr->name, quad->k1*quad_sign, 0., quad->length);
            length += quad->length;
            break;
          case T_RBEN: case T_SBEN:
            bend = (BEND*)eptr->p_elem;
            fprintf(fp, "%-8s b %10.6f          %10.6f\n",
                    eptr->name, bend->angle/PI, bend->length);
            fprintf(fp, " &         %10.6f%10.6f%10.6f%10.6f\n",
                    bend->e1/PI, bend->e2/PI, bend->fint, bend->hgap);
            length += bend->length;
            break;
          case T_DRIF:
            drift = (DRIFT*)eptr->p_elem;
            fprintf(fp, "%-8s s %10.6f%10.6f%10.6f\n",
                    eptr->name, 0.0, 0.0, drift->length);
            length += drift->length;
            break;
          case T_SEXT:
            sext = (SEXT*)eptr->p_elem;
            fprintf(fp, "%-8s n3%10.6f%10.6f\n",
                    eptr->name, sext->k2*sext->length/2, 0.0);
            break;
          case T_MULT:
            mult = (MULT*)eptr->p_elem;
            fprintf(fp, "%-8s n%10.6f%10.6f\n",
                    eptr->name, mult->KnL*cos(mult->tilt*(mult->order+1))/factorial(mult->order),
                    -mult->KnL*sin(mult->tilt*(mult->order+1))/factorial(mult->order));
            break;
          case T_HCOR:
            hcor = (HCOR*)eptr->p_elem;
            fprintf(fp, "%-8s n1%10.6f%10.6f\n",
                    eptr->name, -hcor->kick, 0.0);
            break;
          case T_VCOR:
            vcor = (VCOR*)eptr->p_elem;
            fprintf(fp, "%-8s n1%10.6f%10.6f\n",
                    eptr->name, -vcor->kick, 0.0);
            break;
          case T_HMON: case T_VMON: case T_MONI:
            fprintf(fp, "%-8s m\n", eptr->name);
            break;
          case T_MARK:
            fprintf(fp, "%-8s l\n", eptr->name);
            break;
          default:
            printf("warning: entity type %s not implemented\n",
                   entity_name[eptr->type]);
            break;
            }
        eptr = eptr->succ;
        }
    fprintf(fp, " $\n");
    printf("total length = %.15e\n", length);

    if (ender_file) {
        fpi = fopen_e(ender_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }
    }

char *xorbit_name(char *name, char extra_character)
{
    static char buffer[9];
    long length, maxlen;
    buffer[0] = 0;
    
    length = strlen(name);
    if (extra_character)
        maxlen = 7;
    else
        maxlen = 8;
    strncpy(buffer, name, maxlen);
    if (length>maxlen)
        length = maxlen;
    if (extra_character)
        buffer[length] = extra_character;
    buffer[length+1] = 0;
    return(buffer);
    }

void make_xorbit_friendly(ELEMENT_LIST *elem)
{
    ELEMENT_LIST *eptr, *eptr1;
    double total_length;
    MULT *mult;
    SEXT *sext;
    HCOR *hcor;
    VCOR *vcor;
    char buffer[100];

    /* convert non-zero-length multipoles to zero-length multipoles bracketed with drifts */
    eptr = elem;
    while (eptr) {
        switch (eptr->type) {
          case T_MULT:
            mult = (MULT*)eptr->p_elem;
            if (mult->length>XORBIT_MIN_LENGTH)
                bracket_with_drifts(eptr, mult->length/2);
            eptr = eptr->succ;
            break;
          case T_SEXT:
            sext = (SEXT*)eptr->p_elem;
            if (sext->length>XORBIT_MIN_LENGTH)
                bracket_with_drifts(eptr, sext->length/2);
            eptr = eptr->succ;
            break;
          case T_HCOR:
            hcor = (HCOR*)eptr->p_elem;
            if (hcor->length>XORBIT_MIN_LENGTH)
                bracket_with_drifts(eptr, hcor->length/2);
            eptr = eptr->succ;
            break;
          case T_VCOR:
            vcor = (VCOR*)eptr->p_elem;
            if (vcor->length>XORBIT_MIN_LENGTH)
                bracket_with_drifts(eptr, vcor->length/2);
            eptr = eptr->succ;
            break;
          default:
            break;
            }
        eptr = eptr->succ;
        }

    /* concatenate successive drift spaces and give new names */
    eptr = elem;
    while (eptr) {
        switch (eptr->type) {
          case T_DRIF:
            total_length = ((DRIFT*)eptr->p_elem)->length;
            eptr1 = eptr->succ;
            while (eptr1 && eptr1->type==T_DRIF) {
                total_length += ((DRIFT*)eptr1->p_elem)->length;
                eptr1 = eptr1->succ;
                }
            ((DRIFT*)eptr->p_elem)->length = total_length;
            if (!eptr1)
                eptr->succ = NULL;
            else {
                eptr->succ = eptr1;
                eptr1->pred = eptr;
                }
            if (total_length==0) {
                if (eptr->pred)
                    eptr->pred->succ = eptr->succ;
                if (eptr->succ)
                    eptr->succ->pred = eptr->pred;
                }
            if (eptr->pred) {
                sprintf(buffer, "%sl", eptr->pred->name);
                cp_str(&eptr->name, buffer);
                }
            else 
                cp_str(&eptr->name, "STARTl");
            break;
          default:
            break;
            }
        eptr = eptr->succ;
        }
    }

void bracket_with_drifts(ELEMENT_LIST *elem, double length)
{
    ELEMENT_LIST *eptr;
    static char buffer[100];

    sprintf(buffer, "%sl", elem->name);

    eptr = elem->pred;
    elem->pred = tmalloc(sizeof(*elem->pred));
    elem->pred->type = T_DRIF;
    elem->pred->p_elem = tmalloc(sizeof(DRIFT));
    ((DRIFT*)elem->pred->p_elem)->length = length;
    cp_str(&elem->pred->name, buffer);
    elem->pred->succ = elem;
    if (eptr) {
        elem->pred->pred = eptr;
        eptr->succ = elem->pred;
        }
    
    eptr = elem->succ;
    elem->succ = tmalloc(sizeof(*elem->succ));
    elem->succ->type = T_DRIF;
    elem->succ->p_elem = tmalloc(sizeof(DRIFT));
    ((DRIFT*)elem->succ->p_elem)->length = length;
    cp_str(&elem->succ->name, buffer);
    elem->succ->pred = elem;
    if (eptr) {
        eptr->pred = elem->succ;
        elem->succ->succ = eptr;
        }
    }


        
