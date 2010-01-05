/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: output_magnets
 * purpose: write mpl format data to give plot of magnets.
 *
 * Michael Borland, 1988, 1991
 */
#include "mdb.h"
#include "track.h"

void output_magnets(char *filename, char *line_name, LINE_LIST *beamline)
{
    ELEMENT_LIST *eptr;
    QUAD  *qptr; BEND  *bptr; KQUSE *kqsptr;
    KQUAD *kqptr; KSBEND *kbptr; CSBEND *cbptr; CSRCSBEND *csrbptr;
    long iPhase;
    double start, end, total_length, dz, value;
    FILE *fpm;

    log_entry("output_magnets");

    total_length = 0;
    eptr = &(beamline->elem);

    fpm = fopen_e(filename, "w", 0);

    start = end = 0;
    fprintf(fpm, "SDDS1\n&description text=\"magnet layout for beamline %s\" &end\n", line_name);
    fprintf(fpm, "&column name=ElementName, type=string &end\n");
    fprintf(fpm, "&column name=ElementType, type=string &end\n");
    fprintf(fpm, "&column name=s, units=m, type=double &end\n&column name=Profile, type=double &end\n");
    fprintf(fpm, "&data mode=ascii, no_row_counts=1 &end\n");

    eptr = &(beamline->elem);
    fprintf(fpm, "_BEGIN_ MARK 0 0\n");
    while (eptr) {
         switch (eptr->type) {
            case T_QUAD:
                qptr = (QUAD*)eptr->p_elem;
                fprintf(fpm, "%s %s %e  %d\n", eptr->name, entity_name[eptr->type], start, SIGN(qptr->k1));
                end = start+qptr->length;
                fprintf(fpm, "%s %s %e  %d\n%s %s %e  0\n%s %s %e 0 %\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], end, SIGN(qptr->k1), 
                        eptr->name, entity_name[eptr->type], end, 
                        eptr->name, entity_name[eptr->type], start,
                        eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_KQUAD:
                kqptr = (KQUAD*)eptr->p_elem;
                if (kqptr->bore)
                    value = kqptr->B;
                else
                    value = kqptr->k1;
                fprintf(fpm, "%s %s %e  %d\n", eptr->name, entity_name[eptr->type], start, SIGN(value));
                end = start+kqptr->length;
                fprintf(fpm, "%s %s %e  %d\n%s %s %e  0\n%s %s %e 0 %\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], end, SIGN(kqptr->k1), 
                        eptr->name, entity_name[eptr->type], end, 
                        eptr->name, entity_name[eptr->type], start,
                        eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
             case T_KQUSE:
                kqsptr = (KQUSE*)eptr->p_elem;
                value = kqsptr->k1;
                fprintf(fpm, "%s %s %e  %d\n", eptr->name, entity_name[eptr->type], start, SIGN(value));
                end = start+kqsptr->length;
                fprintf(fpm, "%s %s %e  %d\n%s %s %e  0\n%s %s %e 0 %\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], end, SIGN(value), 
                        eptr->name, entity_name[eptr->type], end, 
                        eptr->name, entity_name[eptr->type], start,
                        eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_RBEN:  case T_SBEN:
                bptr = (BEND*)eptr->p_elem;
                end  = start+bptr->length;
                if (bptr->angle>0)
                    fprintf(fpm, 
                            "%s %s %e .33333333\n%s %s %e .33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                else if (bptr->angle<0)
                    fprintf(fpm, 
                            "%s %s %e -.33333333\n%s %s %e -.33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_KSBEND:
                kbptr = (KSBEND*)eptr->p_elem;
                end  = start+kbptr->length;
                if (kbptr->angle>0)
                    fprintf(fpm, 
                            "%s %s %e .33333333\n%s %s %e .33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                else if (kbptr->angle<0)
                    fprintf(fpm, 
                            "%s %s %e -.33333333\n%s %s %e -.33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_SEXT:
                end = start+((SEXT*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e  .5\n%s %s %e .5\n%s %s %e 0\n%s %s %e 0\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                        eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_KSEXT:
                end = start+((KSEXT*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e  .5\n%s %s %e .5\n%s %s %e 0\n%s %s %e 0\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], start, 
                        eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_HCOR:
                end    = start+((HCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e .25\n%s %s %e .25\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_VCOR:
                end    = start+((VCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e -.25\n%s %s %e -.25\n%s %s %e 0\n", 
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_HVCOR:
                end    = start+((HVCOR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e .25\n%s %s %e -.25\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_DRIF:
                start = (end = start+((DRIFT*)eptr->p_elem)->length);
                fprintf(fpm, "%s %s %e  0\n", eptr->name, entity_name[eptr->type], end);
                break;
            case T_HMON:
                dz = ((HMON*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %s %e 0.125\n%s %s %e 0\n%s %s %e 0\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz, 
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz);
                start += 2*dz;
                break;
            case T_VMON:
                dz = ((VMON*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %s %e -0.125\n%s %s %e 0\n%s %s %e 0\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz, 
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz);
                start += 2*dz;
                break;
            case T_MONI:
                dz = ((MONI*)eptr->p_elem)->length/2;
                fprintf(fpm, "%s %s %e 0.125\n%s %s %e 0\n%s %s %e -0.125\n%s %s %e 0\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz, eptr->name, entity_name[eptr->type], start+dz, 
                        eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start+2*dz);
                start += 2*dz;
                break;
            case T_MULT:
                dz = ((MULT*)eptr->p_elem)->length/3;
                fprintf(fpm, "%s %s %e 0.6666\n%s %s %e 0.6666\n%s %s %e 0\n%s %s %e -0.6666\n%s %s %e -0.6666\n%s %s %e 0\n%s %s %e 0\n",
                        eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start+2*dz, eptr->name, entity_name[eptr->type], start+3*dz,
                        eptr->name, entity_name[eptr->type], start+2*dz, eptr->name, entity_name[eptr->type], start+dz, eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start+3*dz);
                start += 3*dz;
                break;
            case T_MARK:    /* zero-length drift */
                break;
            case T_CSBEND:
                cbptr = (CSBEND*)eptr->p_elem;
                end  = start+cbptr->length;
                if (cbptr->angle>0)
                    fprintf(fpm, 
                        "%s %s %e .33333333\n%s %s %e .33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                else if (cbptr->angle<0) 
                    fprintf(fpm, 
                        "%s %s %e -.33333333\n%s %s %e -.33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_CSRCSBEND:
                csrbptr = (CSRCSBEND*)eptr->p_elem;
                end  = start+csrbptr->length;
                if (csrbptr->angle>0)
                    fprintf(fpm, 
                        "%s %s %e .33333333\n%s %s %e .33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                else if (csrbptr->angle<0) 
                    fprintf(fpm, 
                        "%s %s %e -.33333333\n%s %s %e -.33333333\n%s %s %e  0\n%s %s %e  0\n%s %s %e  0\n%s %s %e 0\n",
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end, eptr->name, entity_name[eptr->type], end, 
                            eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], start, eptr->name, entity_name[eptr->type], end);
                start = end;
                break;
            case T_RFCA: case T_TWLA: case T_RAMPRF: case T_RFCW:
            case T_MODRF:
                dz = ((DRIFT*)eptr->p_elem)->length;
                dz /= 8;
                for (iPhase=0; iPhase<9; iPhase++) {
                  fprintf(fpm, "%s %s %e %e\n",
                          eptr->name, entity_name[eptr->type], start+dz*iPhase,
                          0.5*sin((iPhase/8.0)*PIx2));
                }
                start += dz*8;
                break;
            case T_MATR:
                dz = ((MATR*)eptr->p_elem)->length;
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start, 0.5);
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start+dz, -0.5);
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start+dz, 0.5);
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start, -0.5);
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start, 0.0);
                fprintf(fpm, "%s %s %e %e\n", eptr->name, entity_name[eptr->type], start+dz, 0.0);
                start += dz;
                break;
            default:
                if (entity_description[eptr->type].flags&HAS_LENGTH) {
                    dz = ((DRIFT*)eptr->p_elem)->length;
                    fprintf(fpm, "%s %s %e 0\n", eptr->name, entity_name[eptr->type], start+=dz);
                    }
                break;
            }
        eptr = eptr->succ;
        }
    log_exit("output_magnets");
    fclose(fpm);
    }

