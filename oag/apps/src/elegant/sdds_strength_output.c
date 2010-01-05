/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: sdds_strength_output.c
 * purpose: SDDS format output of (Element, Parameter, KnL, n, BnL)
 *          for a lattice
 * Michael Borland, 1994
 */
#include "mdb.h"
#include "madto.h"
#include "SDDS.h"

void sdds_strength_output(char *output, LINE_LIST *beamline, char *input)
{
    SDDS_TABLE SDDS_table;
    ELEMENT_LIST *eptr;
    char s[SDDS_MAXLINE], t[SDDS_MAXLINE];
    long row, n;
    double KnL=0.0, L=0.0, KnL2PF=0.0, Kn;
    char *param_name=NULL;

    sprintf(s, "Magnet strengths for beamline %s of file %s",
            beamline->name, input);
    if (!SDDS_InitializeOutput(&SDDS_table, SDDS_ASCII, 1, s, 
                               "magnet strengths", output) ||
        SDDS_DefineParameter(&SDDS_table, "Beamline", NULL, NULL, "Beamline name", NULL, SDDS_STRING,
                             beamline->name)<0 ||
        SDDS_DefineParameter(&SDDS_table, "TimeStamp", NULL, NULL, "TimeStamp", NULL, SDDS_STRING,
                             mtime())<0 ||
        SDDS_DefineColumn(&SDDS_table, "s", NULL, NULL, "s", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "ElementName", NULL, NULL, "Element name", "%10s", SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "ElementOccurence", NULL, NULL, "Occurence number of name", "%4ld", SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "ElementType", NULL, NULL, "Element type", "%10s", SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "ElementParameter", NULL, NULL, "Parameter name", "%10s", SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "KnLToParameterFactor", NULL, "SI", "Factor to convert KnL to parameter", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "KnLName", NULL, NULL, "Name of integrated geometric strength", "%10s", SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "BnLName", NULL, NULL, "Name of integrated field strength", "%10s", SDDS_STRING, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "MultipoleOrder", NULL, NULL, "Multipole order", "%6ld ", SDDS_LONG, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "IdealLength", "L$bideal$n", "m", "Ideal effective length", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "KnL", "K$bn$nL", "1/m$an$n", "Integrated geometric strength", NULL, SDDS_DOUBLE, 0)<0 ||
        SDDS_DefineColumn(&SDDS_table, "Kn", "K$bn$n", "1/m$an$n", "Geometric strength", NULL, SDDS_DOUBLE, 0)<0 ||
        !SDDS_WriteLayout(&SDDS_table) || !SDDS_StartTable(&SDDS_table, beamline->n_elems)) {
        SDDS_SetError("Problem setting up SDDS file for magnet strength output");
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }

    row = 0;
    eptr = &(beamline->elem);
    while (eptr) {
        n = -1;
        switch (eptr->type) {
          case T_QUAD:
            L = ((QUAD*)eptr->p_elem)->length;
            KnL = ((QUAD*)eptr->p_elem)->k1*L;
            n = 1;
            param_name = "K1";
            KnL2PF = 1/L;
            break;
          case T_KQUAD:
            L = ((KQUAD*)eptr->p_elem)->length;
            KnL = ((KQUAD*)eptr->p_elem)->k1*L;
            n = 1;
            param_name = "K1";
            KnL2PF = 1/L;
            break;
          case T_SBEN: case T_RBEN:
            L = ((BEND*)eptr->p_elem)->length;
            KnL = ((BEND*)eptr->p_elem)->angle;
            n = 0;
            param_name = "ANGLE";
            KnL2PF = 1;
            break;
          case T_HCOR:
            L = ((HCOR*)eptr->p_elem)->length;
            KnL = -((HCOR*)eptr->p_elem)->kick;
            n = 0;
            param_name = "KICK";
            KnL2PF = -1;
            break;
          case T_VCOR:
            L = ((VCOR*)eptr->p_elem)->length;
            KnL = ((VCOR*)eptr->p_elem)->kick;
            n = 0;
            param_name = "KICK";
            KnL2PF = 1;
            break;
          case T_SEXT:
            L = ((SEXT*)eptr->p_elem)->length;
            KnL = ((SEXT*)eptr->p_elem)->k2*L;
            n = 2;
            param_name = "K2";
            KnL2PF = 1/L;
            break;
          case T_KSEXT:
            L = ((KSEXT*)eptr->p_elem)->length;
            KnL = ((KSEXT*)eptr->p_elem)->k2*L;
            n = 2;
            param_name = "K2";
            KnL2PF = 1/L;
            break;
          case T_KICKER:
            L = ((KICKER*)eptr->p_elem)->length;
            KnL = ((KICKER*)eptr->p_elem)->angle;
            n = 0;
            param_name = "ANGLE";
            KnL2PF = 1;
            break;
          default:
            break;
            }
        if (L>0)
          Kn = KnL/L;
        else if (KnL!=0)
          Kn = DBL_MAX;
        else 
          Kn = 0;
        if (n>=0) {
          sprintf(s, "B%ldL", n);
          sprintf(t, "K%ldL", n);
          if (!SDDS_SetRowValues(&SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row++,
                                   "s", eptr->end_pos,
                                   "ElementName", eptr->name, 
                                   "ElementOccurence", eptr->occurence, 
                                   "ElementType", entity_name[eptr->type], 
                                   "ElementParameter", param_name,
                                   "KnLToParameterFactor", KnL2PF,
                                   "BnLName", s, 
                                   "KnLName", t,
                                   "MultipoleOrder", n, 
                                   "IdealLength", L, 
                                   "KnL", KnL, "Kn", Kn, NULL)) 
                SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
            }
        eptr = eptr->succ;
        }
        
    if (!SDDS_WriteTable(&SDDS_table) || !SDDS_Terminate(&SDDS_table)) {
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        }

    }


