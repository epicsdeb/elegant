/*************************************************************************\
* Copyright (c) 2003 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2003 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: randpert.c
 * purpose: calculate multipole perturbation effects in magnets by random
 *          perturbation of individual pole positions.
 * reference: K. Halbach, NIM Vol. 74, No. 1, 1969
 * 
 * Michael Borland, 1992
 */
#include "mdb.h"
#include "match_string.h"
#include "namelist.h"

void do_perturbations(NAMELIST_TEXT *nltext);

typedef struct {
    char *name;    /* name of the fundamental */
    long N;        /* fundamental multipole.  N=2 is quadrupole, N=3 is sextupole etc */
    double *am;    /* array of Halbach's n*an/N coefficients */
    double *bm;    /* array of Halbach's n*bn/N coefficients */
    double *rhom;  /* array of Halbach's n*rhon/N coefficients */
    long n;        /* number of entries in the tables */
    } MULTIPOLE_COEFS;

#define N_QUAD_ENTRIES 10
double quad_am[N_QUAD_ENTRIES] = {7.46e-2, 2.14e-1, 2.88e-1, 2.31e-1, 1.08e-1, 2.87e-2, 1.04e-2, 1.56e-2, 
    1.25e-2, 5.81e-3};
double quad_bm[N_QUAD_ENTRIES] = {-4.25e-1, -5.15e-1, -2.88e-1, 6.76e-2, 1.08e-1, 4.45e-2, -1.04e-2, 1.28e-2,
    1.25e-2, 6.37e-3};
double quad_rhom[N_QUAD_ENTRIES] = {1.76e-1, 5.00e-1, 6.60e-1, 5.00e-1, 1.91e-1, 0.0, -3.06e-2, 0.0,
    7.53e-3, 0.0};

#define N_SEXT_ENTRIES 15
double sext_am[N_SEXT_ENTRIES] = {5.09e-2, 1.71e-1, 3.03e-1, 3.90e-1, 3.97e-1, 3.18e-1, 1.95e-1, 9.03e-2,
    2.51e-2, 1.90e-3, 5.49e-3, 1.31e-2, 1.36e-2, 9.85e-3, 4.56e-3};
double sext_bm[N_SEXT_ENTRIES] = {-3.14e-1, -4.95e-1, -5.15e-1, -3.90e-1, -1.73e-1, 6.55e-2, 1.08e-1, 9.03e-2, 4.16e-2,
    -1.90e-3, -1.45e-2, 1.05e-2, 1.07e-2, 9.85e-3, 5.06e-3};
double sext_rhom[N_SEXT_ENTRIES] = {8.47e-2, 2.84e-1, 5.00e-1, 6.39e-1, 6.43e-1, 5.00e-1, 2.88e-1, 1.08e-1, 0.0, -3.38e-2,
    -2.05e-2, 0.0, 7.34e-3, 5.82e-3, 0.0};

#define INDEX_QUAD 0
#define INDEX_SEXT 1
#define N_FUNDAMENTAL_TYPES 2
char *type_name[N_FUNDAMENTAL_TYPES] = {
    "quadrupole", "sextupole"
    } ;
MULTIPOLE_COEFS coefs[N_FUNDAMENTAL_TYPES] = {
    {"quadrupole", 2, quad_am, quad_bm, quad_rhom, N_QUAD_ENTRIES},
    {"sextupole", 3, sext_am, sext_bm, sext_rhom, N_SEXT_ENTRIES},
    } ;


#define MULTI_PERT 0
#define N_COMMANDS 1
char *command[N_COMMANDS] = { "perturbations" };

int main(int argc, char **argv)
{
  FILE *fpin;
  char s[1024];
  NAMELIST_TEXT namelist_text;
  
  if (argc!=2)
    bomb(NULL, "randpert inputfile");
  fpin = fopen_e(argv[1], "r", 0);

  while (get_namelist(s, 1024, fpin)) {
    scan_namelist(&namelist_text, s);
    switch (match_string(namelist_text.group_name, command, N_COMMANDS, EXACT_MATCH)) {
    case MULTI_PERT:
      do_perturbations(&namelist_text);
      break;
    default:
      bomb("unknown namelist/command given", NULL);
      break;
    }
  }
  return 0;
}

void do_perturbations(NAMELIST_TEXT *nltext)
{
#include "sddsrandmult.h"
  long N, iN, m, i, i_case;
  FILE *fpSDDS=NULL, *fpElegant=NULL, *fpKMULT=NULL;
  double *f, *g;        /* normal and skew fractional field errors at r=bore_radius */
  double *frms, *grms;
  double dxh, dyh;      /* offsets for the upper half (lower half is negative of this) */
  double fdx, fdy, alpha, gamma, eps;
  double fdX, fdY, srot, crot;
  double dr, dphi;
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
  complex double exp1, exp2, dHn;
#else
  doublecomplex_sdds exp1, exp2, dHn, tmp3, tmp4, tmp5, tmp6;
  double tmp, tmp2;
#endif

  /* initialize variables */
  SDDS_output = elegant_output = kmult_output = name = NULL;
  type = "quadrupole";
  bore_radius = 0.066;
  effective_length = 0.23;
  n_cases = 1000;
  random_number_seed -= 2;
  n_harm = 0;
  dx_pole = dy_pole = dradius = dx_split = dy_split = dphi_halves = 0;
  
  /* process namelist input */
  process_namelist(&perturbations, nltext);
  print_namelist(stdout, &perturbations);

  if (n_cases<=0)
    bomb("n_cases <= 0", NULL);
  if (bore_radius<=0)
    bomb("bore_radius <= 0", NULL);
  if (!SDDS_output)
    bomb("no output filename given", NULL);
  if ((iN=match_string(type, type_name, N_FUNDAMENTAL_TYPES, 0))<0)
    bomb("unknown magnet type", NULL);
  if (n_harm<=0 || n_harm>coefs[iN].n)
    n_harm = coefs[iN].n;
  random_1(abs(random_number_seed));
  if (!name) {
    if (iN==INDEX_QUAD)
      name = "qh";
    else
      name = "sh";
  }

  /* set up output file */    
  fpSDDS = fopen_e(SDDS_output, "w", 0);
  fprintf(fpSDDS, "SDDS1\n");
  for (i=1; i<=n_harm; i++) 
    fprintf(fpSDDS, "&column name=f%ld, description=\"normal %ld-pole fractional error\", type=double &end\n", i-1, 2*i);
  fprintf(fpSDDS, "&column name=F, description=\"normal total fractional error\", type=double &end\n");
  for (i=1; i<=n_harm; i++) 
    fprintf(fpSDDS, "&column name=g%ld, description=\"skew %ld-pole fractional error\", type=double &end\n", i-1, 2*i);
  fprintf(fpSDDS, "&column name=G, description=\"skew total fractional error\", type=double &end\n");
  fprintf(fpSDDS, "&parameter name=MagnetType, type=string, fixed_value=%s &end\n", type);
  fprintf(fpSDDS, "&parameter name=BoreRadius, type=double, fixed_value=%e, units=m &end\n", bore_radius);
  fprintf(fpSDDS, "&parameter name=EffectiveLength, type=double, fixed_value=%e, units=m &end\n", effective_length);
  fprintf(fpSDDS, "&parameter name=SxPole, type=double, units=m, fixed_value=%e &end\n", dx_pole);
  fprintf(fpSDDS, "&parameter name=SyPole, type=double, units=m, fixed_value=%e &end\n", dy_pole);
  fprintf(fpSDDS, "&parameter name=SxSplit, type=double, units=m, fixed_value=%e &end\n", dx_split);
  fprintf(fpSDDS, "&parameter name=SySplit, type=double, units=m, fixed_value=%e &end\n", dy_split);
  fprintf(fpSDDS, "&parameter name=SphiHalves, type=double, fixed_value=%e &end\n", dphi_halves);
  fprintf(fpSDDS, "&parameter name=SrPole, type=double, units=m, fixed_value=%e &end\n", dradius);
  fprintf(fpSDDS, "&data mode=ascii &end\n");
  fprintf(fpSDDS, "%ld\n", n_cases);

  if (kmult_output) {
    fpKMULT = fopen_e(kmult_output, "w", 0);
    fprintf(fpKMULT, "SDDS1\n");
    fprintf(fpKMULT, "&parameter name=MagnetType, type=string, fixed_value=%s &end\n", type);
    fprintf(fpKMULT, "&parameter name=BoreRadius, type=double, fixed_value=%e, units=m &end\n", bore_radius);
    fprintf(fpKMULT, "&parameter name=EffectiveLength, type=double, fixed_value=%e, units=m &end\n", effective_length);
    fprintf(fpKMULT, "&parameter name=SxPole, type=double, units=m, fixed_value=%e &end\n", dx_pole);
    fprintf(fpKMULT, "&parameter name=SyPole, type=double, units=m, fixed_value=%e &end\n", dy_pole);
    fprintf(fpKMULT, "&parameter name=SxSplit, type=double, units=m, fixed_value=%e &end\n", dx_split);
    fprintf(fpKMULT, "&parameter name=SySplit, type=double, units=m, fixed_value=%e &end\n", dy_split);
    fprintf(fpKMULT, "&parameter name=SphiHalves, type=double, fixed_value=%e &end\n", dphi_halves);
    fprintf(fpKMULT, "&parameter name=SrPole, type=double, units=m, fixed_value=%e &end\n", dradius);
    fprintf(fpKMULT, "&parameter name=referenceRadius, type=double, fixed_value=%e, units=m &end\n", bore_radius);
    fprintf(fpKMULT, "&column name=order, type=long &end\n");
    fprintf(fpKMULT, "&column name=an, type=double &end\n");
    fprintf(fpKMULT, "&column name=bn, type=double &end\n");
    fprintf(fpKMULT, "&data mode=ascii no_row_counts=1 &end\n");
  }
  
  if (elegant_output) {
    fpElegant = fopen_e(elegant_output, "w", 0);
    fprintf(fpElegant, "! elegant lattice and namelist input prepared by sddsrandmult\n");
    fprintf(fpElegant, "! for result of mechanical errors in %s magnets with R=%lgm and L=%lgm\n",
            coefs[iN].name, bore_radius, effective_length);
    fprintf(fpElegant, "! %ld cases where calculated with the following mechanical errors present:\n",
            n_cases);
    if (dx_pole || dy_pole)
      fprintf(fpElegant, "!     individual pole placement:   Sx=%lgmm   Sy=%lgmm\n", dx_pole*1e3, dy_pole*1e3);
    if (dradius)
      fprintf(fpElegant, "!     radius error:                Sr=%lgmm\n", dradius*1e3);
    if (dy_split)
      fprintf(fpElegant, "!     top-bottom splitting:        Sy=%lgmm\n", dy_split*1e3);
    if (dy_split)
      fprintf(fpElegant, "!     left-right spitting:         Sx=%lgmm\n", dx_split*1e3);
    if (dphi_halves)
      fprintf(fpElegant, "!     top-bottom counter-rotation: Sdphi=%lgmrad\n", dphi_halves*1e3);
  } 

  dx_pole /= bore_radius;
  dy_pole /= bore_radius;
  dradius /= bore_radius;
  dx_split /= bore_radius;
  dy_split /= bore_radius;

  f = tmalloc(sizeof(*f)*(2*n_harm+2));
  g = f+n_harm+1;
  frms = tmalloc(sizeof(*frms)*(2*n_harm+2));
  grms = frms+n_harm+1;
  for (i=0; i<n_harm+1; i++)
    frms[i] = grms[i] = 0;

  N = coefs[iN].N;
  crot = cos(PI/2./N);
  srot = sin(PI/2./N);
  for (i_case=0; i_case<n_cases; i_case++) {
    /* calculate the parameters for the halves */
    dxh = gauss_rn_lim(0.0, dx_split, 2.0, random_1);
    dyh = gauss_rn_lim(0.0, dy_split, 2.0, random_1);
    dr = gauss_rn_lim(0.0, dradius, 2.0, random_1);
    dphi = gauss_rn_lim(0.0, dphi_halves, 2.0, random_1);
    for (i=0; i<n_harm+1; i++)
      f[i] = g[i] = 0;
    for (m=0; m<2*N; m++) {
      /* sum over all poles */
      /* -- fractional dx and dy in coordinate system N, for which magnet looks normally oriented */
      fdx = gauss_rn_lim(0.0, dx_pole, 2.0, random_1) + (m<N?dxh:-dxh);
      fdy = gauss_rn_lim(0.0, dy_pole, 2.0, random_1) + (m<N?dyh:-dyh);
      /* -- angle of pole in system S, where magnet has a pole on the x axis */
      alpha = (PI/N)*m;
      /* -- (fdx, fdy) rotated into system S, with radial errors */
      fdX =  crot*fdx + srot*fdy + dr*cos(alpha) ;
      fdY = -srot*fdx + crot*fdy + dr*sin(alpha) ;
      /* -- angle and magnitude of displacement vector in S system */
      if (fdX==0 && fdY==0)
        gamma = eps = 0;
      else {
        gamma = atan2(fdY, fdX);
        eps = hypot(fdY, fdX);
      }
      exp1 = cexpi(gamma-alpha);
#if defined(__USE_ISOC99) || defined(__USE_ISOC94)
      exp2 = conj(exp1);
      for (i=1; i<=n_harm; i++) {
        /* loop over all harmonics */
	dHn = (exp1*(coefs[iN].bm[i-1]-coefs[iN].am[i-1])*eps/2 +
	       exp2*(coefs[iN].bm[i-1]+coefs[iN].am[i-1])*eps/2)
	  *cexpi(-(N+i)*alpha)*cexpi(PI/2-(PI*i)/(2.*N));
        f[i-1] += creal(dHn);
        g[i-1] += cimag(dHn);
      }
#else
      exp2.r = exp1.r;
      exp2.i = -exp1.i;
      for (i=1; i<=n_harm; i++) {
        /* loop over all harmonics */
        tmp = (coefs[iN].bm[i-1]-coefs[iN].am[i-1])*eps/2;
        tmp2 = (coefs[iN].bm[i-1]+coefs[iN].am[i-1])*eps/2;
	tmp3.r = exp1.r * tmp + exp2.r * tmp2;
	tmp3.i = exp1.i * tmp + exp2.i * tmp2;
	tmp4 = cexpi(-(N+i)*alpha);
	tmp5 = cexpi(PI/2-(PI*i)/(2.*N));
	tmp6.r = (tmp3.r * tmp4.r - tmp3.i * tmp4.i);
	tmp6.i = (tmp3.i * tmp4.r + tmp3.r * tmp4.i);
	dHn.r = (tmp5.r * tmp6.r - tmp5.i * tmp6.i);
	dHn.i = (tmp5.i * tmp6.r + tmp5.r * tmp6.i);
        f[i-1] += dHn.r;
        g[i-1] += dHn.i;
      }
#endif
    }
    /* add rotation error contribution */
    if (dphi_halves) {
      for (i=1; i<=n_harm; i++) {
        /* loop over all harmonics */
        if ((i+N)%2 && dphi_halves && coefs[iN].rhom[i-1]) 
          f[i-1] += dphi*coefs[iN].rhom[i-1]/cos(PI*i/(2*N));
      }
    }
    /* find maximum total error for all harmonics greater in harmonic number to the fundamental */
    for (i=N; i<n_harm; i++) {
      f[n_harm] += fabs(f[i]);
      g[n_harm] += fabs(g[i]);
    }
    /* keep statistics on rms for each harmonic */
    for (i=0; i<n_harm; i++) {
      frms[i] += sqr(f[i]);
      grms[i] += sqr(g[i]);
    }
    for (i=1; i<=2*n_harm+2; i++)
      fprintf(fpSDDS, "%13.6le ", f[i-1]);
    fputc('\n', fpSDDS);
  }

  for (i=0; i<n_harm; i++) {
    frms[i] = sqrt(frms[i]/n_cases);
    grms[i] = sqrt(grms[i]/n_cases);
  }
  
  if (kmult_output) {
    for (i=1; i<n_harm; i++)
      fprintf(fpKMULT, "%ld %e %e\n", i, frms[i], grms[i]);
  }
  
  if (elegant_output) {
    switch (iN) {
    case INDEX_QUAD:
      fprintf(fpElegant, "! element definitions:\n%sh: quad,l=%lg\n", name, effective_length/2);
      for (i=1; i<n_harm; i++) {
        fprintf(fpElegant, "%sn%ld: mult,n_kicks=1,order=%ld,knl=%lg,factor=1\n",
                name, i, i, frms[i]*dfactorial(i)*effective_length/ipow(bore_radius, i-1));
        fprintf(fpElegant, "%ss%ld: mult,n_kicks=1,order=%ld,tilt=%.12lf,knl=%lg,factor=1\n",
                name, i, i, PI/4, grms[i]*dfactorial(i)*effective_length/ipow(bore_radius, i-1));
      }
      fprintf(fpElegant, "! beamline definition for element with multipoles:\n%s: line=(%sh,", name, name);
      for (i=1; i<n_harm-1; i++) 
        fprintf(fpElegant, "%sn%ld,%ss%ld%s", name, i, name, i, (i%6==0?",&\n    ":","));
      fprintf(fpElegant, "%sn%ld,%ss%ld,%sh)\n", name, i, name, i, name);
      fprintf(fpElegant, "\n! link definitions between pertubations and fundamental:\n");
      for (i=1; i<n_harm; i++) {
        fprintf(fpElegant, "&link_elements target=\"%sn%ld\", item=\"factor\", source=\"%s\", equation=\"K1\" &end\n",
                name, i, name);
        fprintf(fpElegant, "&link_elements target=\"%ss%ld\", item=\"factor\", source=\"%s\", equation=\"K1\" &end\n",
                name, i, name);
      }
      fprintf(fpElegant, "! error definitions for pertubations:\n");
      fprintf(fpElegant, "&error_element name=\"%sn*\", item=\"knl\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      fprintf(fpElegant, "&error_element name=\"%ss*\", item=\"knl\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      break;
    case INDEX_SEXT:
      fprintf(fpElegant, "! element definitions:\n%sh: sext,l=%lg\n", name, effective_length/2);
      for (i=2; i<n_harm; i++) {
        fprintf(fpElegant, "%sm%ld: mult,n_kicks=1,order=%ld,tilt=%e,knl=%g,factor=1\n",
                name, i, i, (frms[i]?atan2(grms[i], frms[i])/(i+1):0.0),
                sqrt(sqr(frms[i])+sqr(grms[i]))*dfactorial(i)*effective_length/ipow(bore_radius, i-2));
        fprintf(fpElegant, "%sn%ld: mult,n_kicks=1,order=%ld,knl=%g,factor=1\n",
                name, i, i, frms[i]*dfactorial(i)*effective_length/ipow(bore_radius, i-2));
        fprintf(fpElegant, "%ss%ld: mult,n_kicks=1,order=%ld,tilt=%.12f,knl=%lg,factor=1\n",
                name, i, i, PI/6, grms[i]*dfactorial(i)*effective_length/ipow(bore_radius, i-2));
      }
      fprintf(fpElegant, "! beamline definition for element with skew and normal multipoles:\n%s_1: line=(%sh,", name, name);
      for (i=2; i<n_harm-1; i++) 
        fprintf(fpElegant, "%sn%ld,%ss%ld%s", name, i, name, i, (i%6==0?",&\n    ":","));
      fprintf(fpElegant, "%sn%ld,%ss%ld,%sh)\n", name, i, name, i, name);
      fprintf(fpElegant, "! beamline definition for element with tilted multipoles:\n%s_2: line=(%sh,", name, name);
      for (i=2; i<n_harm-1; i++) 
        fprintf(fpElegant, "%sm%ld%s", name, i, (i%12==0?",&\n    ":","));
      fprintf(fpElegant, "%sn%ld,%ss%ld,%sh)\n", name, i, name, i, name);
      fprintf(fpElegant, "\n! link definitions between pertubations and fundamental:\n");
      for (i=2; i<n_harm; i++) {
        fprintf(fpElegant, "&link_elements target=\"%sm%ld\", item=\"factor\", source=\"%s\", equation=\"K2\" &end\n",
                name, i, name);
        fprintf(fpElegant, "&link_elements target=\"%sn%ld\", item=\"factor\", source=\"%s\", equation=\"K2\" &end\n",
                name, i, name);
        fprintf(fpElegant, "&link_elements target=\"%ss%ld\", item=\"factor\", source=\"%s\", equation=\"K2\" &end\n",
                name, i, name);
      }
      fprintf(fpElegant, "! error definitions for pertubations:\n");
      fprintf(fpElegant, "&error_element name=\"%sm*\", item=\"tilt\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      fprintf(fpElegant, "&error_element name=\"%sm*\", item=\"knl\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      fprintf(fpElegant, "&error_element name=\"%sn*\", item=\"knl\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      fprintf(fpElegant, "&error_element name=\"%ss*\", item=\"knl\", amplitude=1, type=\"plus_or_minus\" &end\n", name);
      break;
    default:
      bomb("invalid value of iN encountered--programming error", NULL);
      break;
    }
  }

  if (fpKMULT)
    fclose(fpKMULT);
  if (fpElegant)
    fclose(fpElegant);
  if (fpSDDS)
    fclose(fpSDDS);
}
