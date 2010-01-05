/* Copyright 1994 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: resdiag.c
 * purpose: make mpl format file for resonance diagram
 *
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "scan.h"
#include "table.h"

#define SET_ORDER 0
#define SET_SUPERPERIODICITY 1
#define SET_INTEGER_TUNES 2
#define N_OPTIONS 3

char *option[N_OPTIONS] = {
    "order", "superperiodicity", "integertunes"
    } ;

#define USAGE "sddsresdiag <outputfile> [-order=<upper>[,<lower>]] [-superperiodicity=<integer>] [-integerTunes=<xvalue>,<yvalue>]\n\
resdiag makes an mpl format file that can be plotted to make a resonance diagram.\n\
The resonances are shown up to the order specified (3 is the default).  \n\
Program by Michael Borland. (This is version 1, October 2002)."

#define FL_LOC_BOTTOM 1
#define FL_LOC_TOP    2
#define FL_LOC_RIGHT  4
#define FL_LOC_LEFT   8

int nearly_equal(double x1, double y1, double x2, double y2);
int common_factor(long m, long k, long j);

int main(int argc, char **argv)
{
    char *output = NULL;
    FILE *fp_out;
    long n_super, nsol, m, first;
    long max_order, min_order, order;
    double nux, nuy;
    double nuxi[2], nuyi[2];
    long location[2], j, k;
    SCANNED_ARG *scanned;
    long i_arg, type;
    long nux_int, nuy_int;

    argc = scanargs(&scanned, argc, argv); 
    if (argc<2 || argc>(2+N_OPTIONS)) 
        bomb(NULL, USAGE);

    max_order = 3;
    min_order = 1;
    n_super = 1;
    nux_int = nuy_int = 0;

    for (i_arg=1; i_arg<argc; i_arg++) {
        if (scanned[i_arg].arg_type==OPTION) {
            /* process options here */
            switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
              case SET_SUPERPERIODICITY:
                if (scanned[i_arg].n_items!=2 ||
                    sscanf(scanned[i_arg].list[1], "%ld", &n_super)!=1 ||
                    n_super<=0)
                    bomb("invalid -superperiodicity syntax",  USAGE);
                break;    
              case SET_ORDER:
                if (scanned[i_arg].n_items==2) {
                    if (sscanf(scanned[i_arg].list[1], "%ld", &max_order)!=1 || max_order<0)
                        bomb("invalid -order syntax",  USAGE);
                    }
                else if (scanned[i_arg].n_items==3) {
                    if (sscanf(scanned[i_arg].list[1], "%ld", &max_order)!=1 || max_order<0 ||
                        sscanf(scanned[i_arg].list[2], "%ld", &min_order)!=1 || min_order<0 ||
                        min_order>max_order)
                        bomb("invalid -order syntax",  USAGE);
                    }
                else
                    bomb("invalid -order syntax",  USAGE);
                break;    
              case SET_INTEGER_TUNES:
                if (scanned[i_arg].n_items!=3 ||
                    sscanf(scanned[i_arg].list[1], "%ld", &nux_int)!=1 || nux_int<0 ||
                    sscanf(scanned[i_arg].list[2], "%ld", &nuy_int)!=1 || nuy_int<0 )
                    bomb("invalid -integer_tunes syntax", USAGE);
                break;
              default:
                bomb("unknown option given", USAGE);
                break;
                }
            }
        else {
            if (!output)
                output = scanned[i_arg].list[0];
            else
                bomb("too many files listed", USAGE);
            }
        }

    if (!output)
        bomb("no output file listed", USAGE);
    
    fp_out = fopen_e(output, "w", 0);
    fprintf(fp_out, "SDDS1\n&column name=nux, symbol=\"$gn$r$bx$n\", type=double &end\n");
    fprintf(fp_out, "&column name=nuy, symbol=\"$gn$r$by$n\", type=double &end\n");
    fprintf(fp_out, "&parameter name=Order, type=long &end\n");
    fprintf(fp_out, 
            "&parameter name=Description, type=string, fixed_value=\"%ld-superperiod resonance diagram for order %ld through %ld\" &end\n",
            n_super, min_order, max_order);
    fprintf(fp_out, 
            "&parameter name=ResonanceLabel type=string &end\n");
    fprintf(fp_out, "&data mode=ascii no_row_counts=1 &end\n");
    m = max_order*((nux_int>nuy_int?nux_int:nuy_int)+1);

    first = 1;
    do {
        for (k=-max_order; k<=max_order; k++) {
            for (j=-max_order; j<=max_order; j++) {
/*
                if (m==0 && j>k)
                    continue;
 */
                if (abs(k)+abs(j)>max_order || abs(k)+abs(j)<min_order)
                    continue;
/*                if (common_factor(m, k, j))
                    continue;
 */
                order = abs(k)+abs(j);
                if (k==0) {
                    /* vertical line */
                    if (j==0)
                        continue;
                    if ((nux=(1.0*m*n_super-j*nux_int)/j)<0 || nux>1)
                        continue;
                    nuxi[0] = nux;
                    nuyi[0] = 1;  
                    if (nuxi[0]==1)
                        location[0] = FL_LOC_RIGHT;
                    else
                        location[0] = FL_LOC_TOP;
                    nuxi[1] = nux; 
                    nuyi[1] = 0;
                    location[1] = FL_LOC_BOTTOM;
                    type = 1;
                    nsol = 2;
                    }
                else if (j==0) {
                    /* vertical line */
                    if ((nuy=(1.0*m*n_super-k*nuy_int)/k)<0 || nuy>1)
                        continue;
                    nuxi[0] = 0;
                    nuyi[0] = nuy;
                    location[0] = FL_LOC_TOP;
                    nuxi[1] = 1;
                    nuyi[1] = nuy;
                    location[1] = FL_LOC_RIGHT;
                    type = 2;
                    nsol = 2;
                    }
                else {
                    /* find intersection points with four boundaries of diagram */
                    location[0] = location[1] = 0;
                    nsol = 0;
                    nuyi[0] = 0;
                    if ((nuxi[0] = (m*n_super-j*nux_int-k*nuy_int)/(1.0*j))>=0 && nuxi[0]<=1) {
                        location[0] = FL_LOC_BOTTOM;
                        nsol++;
                        }
                    nuyi[nsol] = 1;
                    if ((nuxi[nsol] = (m*n_super-j*nux_int-k*(1+nuy_int))/(1.0*j))>=0 && nuxi[nsol]<=1 &&
                            !nearly_equal(nuxi[0], nuyi[0], nuxi[1], nuyi[1])) {
                        location[nsol] = FL_LOC_TOP;
                        nsol++;
                        }
                    type = 3;
                    if (nsol!=2) {
                        nuxi[nsol] = 0;
                        if ((nuyi[nsol] = (m*n_super-j*nux_int-k*nuy_int)/(1.0*k))>=0 && nuyi[nsol]<=1 &&
                                !nearly_equal(nuxi[0], nuyi[0], nuxi[1], nuyi[1])) {                          
                            location[nsol] = FL_LOC_LEFT;
                            type = 4;
                            nsol++;
                            }
                        if (nsol!=2) {
                            nuxi[nsol] = 1;
                            if ((nuyi[nsol] = (m*n_super-j*(1+nux_int)-k*nuy_int)/(1.0*k))>=0 && nuyi[nsol]<=1) {
                                location[nsol] = FL_LOC_RIGHT;
                                type = 5;
                                nsol++;
                                }
                            }
                        }
                   }
                if (nsol==2) {
                    if (nuxi[0]>nuxi[1]) {
                        SWAP_DOUBLE(nuxi[0], nuxi[1]);
                        SWAP_DOUBLE(nuyi[0], nuyi[1]);
                        SWAP_LONG(location[0], location[1]);
                        }
                    if (!first) 
                      fprintf(fp_out, "\n");
                    first = 0;
                    fprintf(fp_out, "! %ld NUx +  %ld NUy = %ld        type = %ld\n", j, k, m*n_super, type);
                    fprintf(fp_out, "%ld\n", order);
                    fprintf(fp_out, "%ld$gn$r$bx$n+%ld$gn$r$by$n=%ld\n", j, k, m*n_super);
                    if (location[0]&FL_LOC_BOTTOM || location[0]&FL_LOC_LEFT)
                        fprintf(fp_out, "%lf %lf\n", nux_int+nuxi[0], nuy_int+nuyi[0]);
                    else if (location[0]&FL_LOC_RIGHT)
                        fprintf(fp_out, "%ld %ld\n%lf %lf\n", nux_int+1, nuy_int, nux_int+nuxi[0], nuy_int+nuyi[0]);
                    else if (location[0]&FL_LOC_TOP)
                        fprintf(fp_out, "%ld %ld\n%lf %lf\n", nux_int, nuy_int+1, nux_int+nuxi[0], nuy_int+nuyi[0]);
                    else
                        bomb("unknown location flag for boundary intersection", NULL);

                    if (location[1]&FL_LOC_BOTTOM || location[1]&FL_LOC_LEFT)
                        fprintf(fp_out, "%lf %lf\n%ld %ld\n", nux_int+nuxi[1], nuy_int+nuyi[1], nux_int, nuy_int);
                    else if (location[1]&FL_LOC_RIGHT || location[1]&FL_LOC_TOP) {
                        fprintf(fp_out, "%lf %lf\n%lf %lf\n",
                            nux_int+nuxi[1], nuy_int+nuyi[1], nux_int+nuxi[0], nuy_int+nuyi[0]);
                        if (location[0]&FL_LOC_BOTTOM || location[0]&FL_LOC_LEFT)
                            fprintf(fp_out, "%ld %ld\n", nux_int, nuy_int);
                        else if (location[0]&FL_LOC_RIGHT)
                            fprintf(fp_out, "%ld %ld\n%ld %ld\n", nux_int+1, nuy_int, nux_int, nuy_int);
                        else if (location[0]&FL_LOC_TOP)
                            fprintf(fp_out, "%ld %ld\n%ld %ld\n", nux_int, nuy_int+1, nux_int, nuy_int);
                        }
                    else
                        bomb("unknown location flag for boundary intersection", NULL);
                    }
                }
            }
        } while (--m>=0);
    fclose(fp_out);
    return 0;
    }

int nearly_equal(double x1, double y1, double x2, double y2)
{
    if ((fabs(x1-x2)+fabs(y1-y2))<1e-6)
        return(1);
    return(0);
    }

int common_factor(long m, long k, long j)
{
    long i, n, nd;

    n = (abs(m)>abs(k)?abs(m):abs(k));
    n = (abs(n)>abs(j)?abs(n):abs(j));

    for (i=2; i<=n; i++) {
        if (i*(k/i)-k)
            continue;
        if (i*(j/i)-j)
            continue;
        if (i*(m/i)-m)
            continue;
        return(i);
        }
    return(0);
    }
            

