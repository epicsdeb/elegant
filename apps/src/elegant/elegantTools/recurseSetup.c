/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: recurse.c
 * purpose: work out recursion relations for SLAC 75 for input into mathematica
 * 
 * Michael Borland, 1991, 1996
 $Log: not supported by cvs2svn $
 Revision 1.4  2002/08/14 20:23:46  soliday
 Added Open License

 Revision 1.3  1999/08/05 15:39:53  soliday
 Added WIN32 and Linux support

 Revision 1.2  1999/03/10 17:21:37  borland
 Changed the naming of variables in the mathematica input file so that the
 meaning is clear for higher orders.

 Revision 1.1  1999/03/02 03:48:11  borland
 First version in repository.
 
 */
#include "mdb.h"
#include "scan.h"

void test_and_print(FILE *fpo, int i1, int i2, int factor, char *sfact);


int main(int argc, char **argv)
{
    int nx, ny, n, m, maxord, ix, iy;
    FILE  *fpo;
    char s1[10], s2[10];
    char *filename;

    if (argc>=2)
        fpo = fopen_e(filename=*++argv, "w", 0);
    else {
        fpo = stdout;
        cp_str(&filename, "recurse_results.");
        }

    nx = ny = maxord = query_int("maximum order", 6);

    /* recursion relations for Anm coefficients */
    for (n=0; n<=nx; n++) {
        for (m=0; m<=ny/2; m++) {
            if ((2*m+n)>maxord)
                continue;
            if (2*m+3>9)
                sprintf(s1, "%c", 'A'+2*m+3-9);
            else
                sprintf(s1, "%d", 2*m+3);
            if (n>9)
                sprintf(s2, "%c", 'A'+n-9);
            else
                sprintf(s2, "%d", n);
            fprintf(fpo, "F%s%s := -(", s1, s2);
            test_and_print(fpo, 2*m+1, n+2, 1, "");
            test_and_print(fpo, 2*m+1, n+1, 3*n+1, "h*");
            test_and_print(fpo, 2*m+1,   n, (3*n-1)*n, "h^2*");
            test_and_print(fpo, 2*m+1, n-1, n*(n-1)*(n-1), "h^3*");
            test_and_print(fpo, 2*m+3, n-1, 3*n, "h*");
            test_and_print(fpo, 2*m+3, n-2, 3*n*(n-1), "h^2*");
            test_and_print(fpo, 2*m+3, n-3, n*(n-1)*(n-2), "h^3*");
            fprintf(fpo, ")\n");
            }
        }

    /* expression for A */
    fprintf(fpo, "A := (");
    for (n=0; n<=nx; n++) {
        fprintf(fpo, "%sx^%d/%.0f*(", (n==0?" ":"\n      + "), n, dfactorial(n));
        for (m=0; m<=ny/2; m++) {
            if ((2*m+1+n)>maxord+1)
                break;
            if (2*m+1>9)
                sprintf(s1, "%c", 'A'+(2*m+1-9));
            else
                sprintf(s1, "%d", 2*m+1);
            if (n>9)
                sprintf(s2, "%c", 'A'+n-9);
            else
                sprintf(s2, "%d", n);
            fprintf(fpo, "%sF%s%s*y^%d/%.0f", (m==0?" ":" + "), s1, s2, 2*m+1, dfactorial(2*m+1));
            }
        fprintf(fpo, ")");
        }
    fprintf(fpo, "    )\n\n");

    fprintf(fpo, "Fx := D[A, x]\nFy := D[A, y]\n");
    fprintf(fpo, "F10 := 1\n");
    fprintf(fpo, "F11 := -nh\n");
    for (n=2; n<20; n++) {
      if (n>9)
        fprintf(fpo, "F1%c := %.0f*%ch%d\n", 
                'A'+n-9, dfactorial(n), 'A'+n-1, n);
      else 
        fprintf(fpo, "F1%d := %.0f*%ch%d\n", 
                n, dfactorial(n), 'A'+n-1, n);
    }
    
    fprintf(fpo, "rfile = OpenWrite[\"%sout\"]\n", filename);

    /* expressions for Fx coefficients of powers of x and y */
    for (iy=0; iy<=maxord; iy++)
        for (ix=0; ix<=maxord-iy; ix++) {
            fprintf(fpo, "Fx%02d%02d := Simplify[D[Fx, {x, %d}, {y, %d}]/%.0f /. {x->0, y->0}]\n",
                    ix, iy, ix, iy, dfactorial(ix)*dfactorial(iy));
            fprintf(fpo, "WriteString[rfile, \"Fx%02d%02d = \"]\n", ix, iy);
            fprintf(fpo, "Write[rfile, CForm[Fx%02d%02d]]\n", ix, iy);
            }

    /* expressions for Fy coefficients of powers of x and y */
    for (iy=0; iy<=maxord; iy++)
        for (ix=0; ix<=maxord-iy; ix++) {
            fprintf(fpo, "Fy%02d%02d := Simplify[D[Fy, {x, %d}, {y, %d}]/%.0f /. {x->0, y->0}]\n",
                    ix, iy, ix, iy, dfactorial(ix)*dfactorial(iy));
            fprintf(fpo, "WriteString[rfile, \"Fy%02d%02d = \"]\n", ix, iy);
            fprintf(fpo, "Write[rfile, CForm[Fy%02d%02d]]\n", ix, iy);
            }
    return(0);

}

void test_and_print(FILE *fpo, int i1, int i2, int factor, char *sfact)
{
    char s1[5], s2[5], sf[20];

    if (i1<0 || i2<0 || factor==0)
        return;

    if (i1>9)
        sprintf(s1, "%c", 'A'+(i1-9));
    else
        sprintf(s1, "%d", i1);

    if (i2>9)
        sprintf(s2, "%c", 'A'+(i2-9));
    else
        sprintf(s2, "%d", i2);

    if (factor==1)
        sf[0] = 0;
    else
        sprintf(sf, "%d*", factor);
    
    fprintf(fpo, "+ %s%sF%s%s", sf, sfact, s1, s2);
    }

