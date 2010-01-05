/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: convert_to_transport
 * purpose: convert parsed mad lattice to transport format (writes output
 *          to a file)
 * 
 * Michael Borland, 1988.
 */
#include "mdb.h"
#include "madto.h"
#include <ctype.h>

#define PC_CONST (((double)1e10)/c_mks)

void convert_to_transport(char *outputfile, LINE_LIST *beamline, long flip_k, double angle_tolerance, 
                          char *header_file, char *ender_file,
                          double quad_aperture, double sext_aperture, double pc)
{
    FILE *fp;
    ELEMENT_LIST *eptr;
    double rho;
    QUAD  *quad;
    SEXT  *sext;
    BEND  *bend;
    DRIFT *drift;
    VCOR *vcor;
    char *quoted_label(), s[100];
    long n_bend_k2=0, n_bend_gap=0, n_bend_fint=0;
    long n_bend_rotated_faces=0, n_bend_curved_faces=0;
    long n_bend_tilts=0, n_quad_tilts=0, n_sext_tilts=0;
    long quad_sign;
    FILE *fpi;
    char buffer[512];

    fp = fopen_e(outputfile, "w", 0);

    quad_sign = flip_k?-1:1;
    
    quad_aperture *= 1e-3;
    sext_aperture *= 1e-3;

    fprintf(fp, "(* pc = %g Gev/c *)\n", pc);

    if (header_file) {
        fpi = fopen_e(header_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    /* unit-change cards */
    fprintf(fp, "15. 1. 'M' ;\n15. 2. 'R' ;\n15. 5. 'M' ;\n15. 6. 'N' ;\n");
    
    /* default beam card */
    sprintf(s, "1. 1E-3 1E-3 1E-3 1E-3 1E-3 1E-3 %.16g;\n", pc);
    str_toupper(s);
    fputs(s, fp);

    /* definitions for use with tplot */
    fprintf(fp, "13 5 'MAPR';\n13 2 'BMPR';\n");

    eptr = &(beamline->elem);
    pc *= PC_CONST;
    while (eptr) {
        switch (eptr->type) {
            case T_QUAD:
                quad = (QUAD*)eptr->p_elem;
                if (quad->tilt!=0) {
                    sprintf(s, "%%QT%02ld", n_quad_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, quad->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                sprintf(buffer, "%8.8s %-5d %21.16e %21.16e %21.16e;\n",
                    quoted_label(eptr->name), 5, quad->length,
                    -quad->k1*quad_aperture*pc*quad_sign, 
                    quad_aperture);
                do_output_transport(fp, buffer);
                if (quad->tilt!=0) {
                    sprintf(s, "%%QT%02ld", n_quad_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, -quad->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                break;
            case T_SBEN: case T_RBEN:
                bend  = (BEND*)eptr->p_elem;
                if (eptr->type==T_RBEN) {
                    rho = bend->length/(2*sin(bend->angle/2));
                    bend->length = rho*bend->angle;
                    }
                else 
                    rho = bend->length/bend->angle;
                /* set special parameters */
                if (bend->hgap!=0) {
                    sprintf(s, "GP%02ld", n_bend_gap++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 5, bend->hgap);
                    do_output_transport(fp, buffer);
                    }
                if (bend->h1!=0) {
                    sprintf(s, "CF%02ld", n_bend_curved_faces++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n", 
                        quoted_label(s), 16, 12, bend->h1);
                    do_output_transport(fp, buffer);
                    }
                if (bend->h2!=0) {
                    sprintf(s, "CF%02ld", n_bend_curved_faces++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n", 
                        quoted_label(s), 16, 13, bend->h2);
                    do_output_transport(fp, buffer);
                    }
                if (bend->k2!=0) {
                    sprintf(s, "BS%02ld", n_bend_k2++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 1, .5*bend->k2*rho);
                    do_output_transport(fp, buffer);
                    }
                if (bend->fint!=.5) {
                    sprintf(s, "FI%02ld", n_bend_fint++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 7, bend->fint);
                    do_output_transport(fp, buffer);
                    }
                /* use 20 cards to tilt and un-tilt the magnet */
                if (bend->tilt!=0) {
                    sprintf(s, "%%BT%02ld", n_bend_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, bend->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                /* print the pole faces and the magnet */
                sprintf(s, "FA%02ld", n_bend_rotated_faces++);
                sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 2, bend->e1*180/PI);
                do_output_transport(fp, buffer);
                sprintf(buffer, "%8.8s %-5d %21.16e %21.16e %21.16e;\n",
                        quoted_label(eptr->name), 4, bend->length, pc/rho, 
                        -bend->k1*sqr(rho));
                do_output_transport(fp, buffer);
                sprintf(s, "FA%02ld", n_bend_rotated_faces++);
                sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 2, bend->e2*180/PI);
                do_output_transport(fp, buffer);
                if (bend->tilt!=0) {
                    sprintf(s, "%%BT%02ld", n_bend_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, -bend->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                /* reset special parameters to default values */
                if (bend->hgap!=0) {
                    sprintf(s, "GP%02ld", n_bend_gap++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 5, 0.);
                    do_output_transport(fp, buffer);
                    }
                if (bend->h1!=0) {
                    sprintf(s, "CF%02ld", n_bend_curved_faces++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n", 
                        quoted_label(s), 16, 12, 0.);
                    do_output_transport(fp, buffer);
                    }
                if (bend->h2!=0) {
                    sprintf(s, "CF%02ld", n_bend_curved_faces++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n", 
                        quoted_label(s), 16, 13, 0.);
                    do_output_transport(fp, buffer);
                    }
                if (bend->k2!=0) {
                    sprintf(s, "BS%02ld", n_bend_k2++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 1, 0.);
                    do_output_transport(fp, buffer);
                    }
                if (bend->fint!=.5) {
                    sprintf(s, "FI%02ld", n_bend_fint++);
                    sprintf(buffer, "%8.8s %-5d %12d %21.16e;\n",
                        quoted_label(s), 16, 7, 0.);
                    do_output_transport(fp, buffer);
                    }
                break;
            case T_DRIF:
                drift = (DRIFT*)eptr->p_elem;
                sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(eptr->name),
                        3, drift->length);
                do_output_transport(fp, buffer);
                break;
            case T_VCOR:
            case T_HCOR:
                vcor = (VCOR*)eptr->p_elem;
                sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(eptr->name),
                        3, vcor->length);
                do_output_transport(fp, buffer);
                break;
            case T_SEXT:
                sext = (SEXT*)eptr->p_elem;
                if (sext->tilt!=0) {
                    sprintf(s, "%%ST%02ld", n_sext_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, sext->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                sprintf(buffer, "%8.8s %-5d %21.16e %21.16e %21.16e;\n",
                        quoted_label(eptr->name), 18, sext->length, 
                        sext->k2*sqr(sext_aperture)*pc/2,
                        sext_aperture);
                do_output_transport(fp, buffer);
                if (sext->tilt!=0) {
                    sprintf(s, "%%ST%02ld", n_sext_tilts++);
                    sprintf(buffer, "%8.8s %-5d %21.16e;\n", 
                        quoted_label(s), 20, -sext->tilt*180/PI);
                    do_output_transport(fp, buffer);
                    }
                break;
            default:
                printf("unknown element type--element '%s' ignored\n",
                        eptr->name);
                break;
            }
        eptr = eptr->succ;
        }
    fputs("SENTINEL\n'MATRICES'\n1\n13 6 'MAPR';\n13 3 'BMPR';\nSENTINEL\n",
        fp);

    if (ender_file) {
        fpi = fopen_e(ender_file, "r", 1);
        while (fgets(s, 256, fpi))
            fputs(s, fp);
        fclose(fpi);
        }

    }

void do_output_transport(FILE *fp, char *s)
{
    str_toupper(s);
    section(s, 75);
    fputs(s, fp);
    }


/* routine: section()
 * purpose: break a string into sections separated by newline characters
 *          string is broken at white-space characters only, with resulting
 *          sections being no longer than n.
 */

void section(char *s, int n)
{
    register char *ptr;

    if (strlen(s)<n) 
        return;

    ptr = s + n;
    while (ptr!=s && !isspace(*ptr))
        ptr--;
    if (ptr==s) 
        return;
    
    while (ptr!=s && isspace(*ptr))
        ptr--;
    if (ptr==s) 
        return;

    *++ptr = '\n';
    section(++ptr, n);
    }

