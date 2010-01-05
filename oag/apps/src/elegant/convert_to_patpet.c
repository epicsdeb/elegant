/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: convert_to_patpet
 * purpose: convert parsed mad lattice to patpet format (writes data to
 *          file)
 *
 * Michael Borland, 1988
 */
#include "mdb.h"
#include "madto.h"
#include <ctype.h>

#define RAD_TO_DEG(x) ( (int) ((x)*180./PI) )

void convert_to_patpet(char *outputfile, LINE_LIST *beamline, long flip_k, double angle_tolerance,
                    char *header_file, char *ender_file)
{
    FILE *fp;
    ELEMENT_LIST *eptr;
    QUAD  *quad;
    SEXT  *sext;
    BEND  *bend;
    HCOR  *hcor;
    VCOR  *vcor;
    HMON  *hmon;
    VMON  *vmon;
    MONI  *moni;
    DRIFT *drift;
    char s[100], output[300];
    int count, quad_sign, bend_sign;
    FILE *fpi;

    fp = fopen_e(outputfile, "w", 0);

    if (angle_tolerance<=0)
        angle_tolerance = 1e-5;

    if (header_file) {
        fpi = fopen_e(header_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    /* define a zero-length drift space */
    sprintf(output, "%-4s%6d%5d     %20.7e\n", "D000", 1, 0,  0.0);
    do_output(fp, output);

    quad_sign = (flip_k )?(-1):( 1);
    bend_sign = (-1);

    eptr = &(beamline->elem);
    while (eptr) {
        switch (eptr->type) {
            case T_QUAD:
                quad = (QUAD*)eptr->p_elem;
                sprintf(output, 
                    "%-4s%6d%5d     %20.7e%20.7e%5.1f%5.1f\n",
                    eptr->name, 3, RAD_TO_DEG(quad->tilt),  quad->length, 
                    quad->k1*quad_sign, 0., 0.); 
                do_output(fp, output);
                break;
            case T_RBEN:  case T_SBEN:
                bend = (BEND*)eptr->p_elem;
                if (!(bend->e1==0 && bend->e2==0)) {
                    /* deal with non-zero edge angles to the extent of
                     * allowing edge angles to change sector magnet into
                     * rectangular bend
                     */
                    if (eptr->type==T_RBEN) {
                        printf("error: rectangular bends must have zero edge angles\n");
                        exit(1);
                        }
                    if (!nearly_equal(bend->e1, bend->e2, angle_tolerance)) {
                        printf("error: bend %s has unequal edge angles\n",
                            eptr->name);
                        exit(1);
                        }
                    if (!nearly_equal(bend->e1*2, bend->angle, angle_tolerance)) {
                        printf("error: sector bend %s has edge angles of %.4e and  bending angle %.4e.\n",
                            eptr->name, bend->e1, bend->angle);
                        puts("Sector bend may only have edge angles such that they are equivalent to");
                        printf("rectangular magnets (to within a tolerance of %e).\n", 
                                angle_tolerance);
                        puts("You may change this edge-angle tolerance with the -tolerance switch.");
                        exit(1);
                        }
                    else {
                        bend->e1 = bend->e2 = 0;
                        eptr->type = T_RBEN;
                        }
                    }
                sprintf(output,  
                    "%-4s%6d%5d %c   %20.7e%20.7e%5.1f%5.1f%10.7f\n",
                    eptr->name, 2, RAD_TO_DEG(bend->tilt),  
                    ((eptr->type==T_SBEN)?'1':' '), bend->length, 
                    bend_sign*((bend->angle==0)?0.:(bend->length/bend->angle)),
                    0., bend->hgap*2000, bend->k1*quad_sign); 
                do_output(fp, output);
                break;
            case T_DRIF:
                drift = (DRIFT*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e\n",
                    eptr->name, 1, 0,  drift->length);
                do_output(fp, output);
                break;
            case T_SEXT:
                sext = (SEXT*)eptr->p_elem;
                sprintf(output, 
                    "%-4s%6d%5d     %20.7e%20.7e%5.1f%5.1f\n",
                    eptr->name, 4, RAD_TO_DEG(sext->tilt),  sext->length, 
                    sext->length*sext->k2*quad_sign, 0., 0.); 
                do_output(fp, output);
                break;
            case T_HCOR:
                hcor = (HCOR*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e%20.7e\n",
                    eptr->name, 1, RAD_TO_DEG(hcor->tilt), hcor->length, 
                    (hcor->length==0)?0.0:(hcor->kick/hcor->length));
                do_output(fp, output);
                break;
            case T_VCOR:
                vcor = (VCOR*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e%20.7e\n",
                    eptr->name, 1, 90+RAD_TO_DEG(vcor->tilt), vcor->length, 
                    (vcor->length==0)?0.0:(vcor->kick/vcor->length));
                do_output(fp, output);
                break;
            case T_HMON:
                hmon = (HMON*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e\n",
                    eptr->name, 1, 0, hmon->length);
                do_output(fp, output);
                break;
            case T_VMON:
                vmon = (VMON*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e\n",
                    eptr->name, 1, 90, vmon->length);
                do_output(fp, output);
                break;
            case T_MONI:
                moni = (MONI*)eptr->p_elem;
                sprintf(output, "%-4s%6d%5d     %20.7e\n",
                    eptr->name, 1, 90, moni->length);
                do_output(fp, output);
                break;
            case T_MARK:    /* zero-length drift */
                sprintf(output, "%-4s%6d%5d     %20.7e\n",
                    eptr->name, 1, 0,  0.0);
                do_output(fp, output);
                break;
            default:
                printf("warning: entity type %s not implemented\n",
                    entity_name[eptr->type]);
                break;
            }
        eptr = eptr->succ;
        }

    eptr = &(beamline->elem);
    fprintf(fp, "END\n");
    count = 0;

    /* put zero-length drift space at start of lattice */
    fputs("D000", fp);
    count++;

    while (eptr) {
        fprintf(fp, "%-4s", eptr->name);
        if (count==19) {
            fputc('\n', fp);
            count = 0;
            }
        else 
            count++;
        eptr = eptr->succ;
        }

    /* put zero-length drift space at end of lattice */
    fputs("D000", fp);
    if (count==19) {
        fputc('\n', fp);
        count = 0;
        }
    else count++;

    if (count!=0)
        fputc('\n', fp);
    fprintf(fp, "END\n");

    if (ender_file) {
        fpi = fopen_e(ender_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }
    }

