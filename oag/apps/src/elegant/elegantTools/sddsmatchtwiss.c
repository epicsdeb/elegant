/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsmatchtwiss.c
 * purpose: take particle coordinate input and convert the Twiss parameters
 *          to user-supplied values.
 *
 * Michael Borland, 2000
 *
 $Log: sddsmatchtwiss.c,v $
 Revision 1.2  2007/03/30 21:20:35  soliday
 Modified to work on WIN32

 Revision 1.1  2007/03/30 16:50:29  soliday
 Moved from directory above.

 Revision 1.21  2007/03/16 02:07:27  borland
 Added option to specify geometric emittance.

 Revision 1.20  2006/10/23 19:49:43  soliday
 Updated to fix an issue with linux-x86_64

 Revision 1.19  2005/11/10 15:38:49  soliday
 Added changes to get it to compile properly with 64 bit compilers.

 Revision 1.18  2003/05/26 17:54:23  borland
 Supplied a missing variable initializatin (x in main()).
 Fixed bug from improper use of spec->element to test for user input of
 an element name (should use spec->flags).

 Revision 1.17  2002/08/14 20:23:49  soliday
 Added Open License

 Revision 1.16  2002/07/01 15:58:27  borland
 Optionally takes specification of twiss parameters from an elegant
 twiss output file.

 Revision 1.15  2002/05/29 20:06:03  borland
 Allow zero emittances.

 Revision 1.14  2002/05/20 23:08:39  borland
 Fixed memory leak.

 Revision 1.13  2002/05/20 21:33:16  borland
 No longer attempts to perform transformation when it isn't requested.
 Fixed bug that resulted when energy spread was exactly zero.

 Revision 1.12  2002/01/30 23:46:44  borland
 Fixed bug in beta/alpha matching when beta and alpha aren't given on the
 command line.

 Revision 1.11  2002/01/22 18:42:46  borland
 Fixed problem with usage message (missing escapes at end of line).

 Revision 1.10  2002/01/21 01:27:02  borland
 Fixed bug with -xplane and -yplane: previously, didn't allow changing
 beta and alpha when changing the emittance.
 Added chirp capability to -zPlane.

 Revision 1.9  2001/10/15 20:37:07  soliday
 Cleaned up for Linux.

 Revision 1.8  2001/10/15 15:42:33  soliday
 Cleaned up for WIN32.

 Revision 1.7  2001/05/30 15:53:27  borland
 -zPlane option now allows specifying alphaz or correlation.

 Revision 1.6  2001/05/15 16:52:54  borland
 Added to -zPlane option the ability to change the central momentum.

 Revision 1.5  2001/05/08 14:07:22  borland
 Added -ztransform option, which permits changing the longitudinal parameters.

 Revision 1.4  2001/04/06 14:44:55  borland
 Added ability to change the emittance of the beam in both planes.

 Revision 1.3  2001/01/08 20:22:16  borland
 Added -oneTransform option, which allows computing the transform only for
 the first page, and using it for all pages.
 Also, added hidden -verbose option.

 Revision 1.2  2000/08/17 13:45:56  borland
 Fixed errors in usage message.

 Revision 1.1  2000/08/14 21:44:56  borland
 First version.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_XPLANE 1
#define SET_YPLANE 2
#define SET_NOWARNINGS 3
#define SET_ONETRANSFORM 4
#define SET_VERBOSE 5
#define SET_ZPLANE 6
#define N_OPTIONS 7

char *option[N_OPTIONS] = {
  "pipe", "xplane", "yplane", "nowarnings", "oneTransform", "verbose", 
  "zplane",
} ;

char *USAGE1="sddsmatchtwiss [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-xPlane=[beta=<meters>,alpha=<value>][,{nemittance=<meters>|emittance=<meters>}]\n\
           [,etaValue=<meters>][,etaSlope=<value>]\n\
           [,filename=<filename>[,element=<name>[,occurrence=<number>]]]]\n\
  [-yPlane=[beta=<meters>,alpha=<value>][,{nemittance=<meters>|emittance=<meters>}]\n\
           [,etaValue=<meters>][,etaSlope=<value>]\n\
           [,filename=<filename>[,element=<name>[,occurrence=<number>]]]]\n\
  [-zPlane=[deltaStDev=<value>][,tStDev=<seconds>]\n\
           [,{correlation=<seconds>|alpha=<value>}][,chirp=<1/seconds>]\n\
           [,betaGamma=<central-value>]]\n\
  [-nowarnings] [-oneTransform]\n";
char *USAGE2="The input file must have columns x, xp, y, yp, and p; for example, an\n\
elegant beam output file is acceptable.  If filename is not given, then\n\
beta and alpha must be given together, or omitted together.  etaValue \n\
and etaSlope may be given individually or together.  If etaValue is not \n\
given, the coordinates have no dispersion adjustment.  If etaSlope is not \n\
given, then slopes have no dispersion adjustment.\n\
If filename is given, the file is expected to be an elegant twiss parameter\n\
output file.  By default, all twiss parameters are taken from the final row\n\
of the first page of this file.  If 'element' is given, however, the twiss\n\
parameters are taken from the last row for which the element name matches the\n\
given string.  However, if 'occurrence=<n>' is given, then the data is taken from\n\
the nth matching row.\n\
For -zPlane, operations are performed as follows: matching of the\n\
fraction momentum spread, bunch length, and correlation/alpha;\n\
application of the chirp; adjustment of central beta-gamma.\n\
If central beta-gamma is adjusted, the normalized emittance is preserved\n\
in x and y planes while the absolute momentum spread is preserved in\n\
the longitudinal plane.\n\
If -oneTransform is given, then the transformation is computed for the\n\
first page only, then reused for all subsequent pages.\n\n\
Program by Michael Borland.  ("__DATE__")\n";

typedef struct {
  double beta, alpha, eta, etap, normEmittance, emittance;
  unsigned long flags;
  double R11, R12, R21, R22;
  double etaBeam, etapBeam;
  char *filename, *element;
  int32_t occurrence;
#define BETA_GIVEN       0x0001U
#define ALPHA_GIVEN      0x0002U
#define ETA_GIVEN        0x0004U
#define ETAP_GIVEN       0x0008U
#define NEMIT_GIVEN      0x0010U
#define FILENAME_GIVEN   0x0020U
#define ELEMENT_GIVEN    0x0040U
#define OCCURRENCE_GIVEN 0x0080U
#define EMIT_GIVEN       0x0100U
} PLANE_SPEC;

typedef struct {
  double deltaStDev, tStDev, correlation, alpha, betaGamma, chirp;
  unsigned long flags;
  double R11, R12, R21, R22;
#define DELTASTDEV_GIVEN  0x0001U
#define TSTDEV_GIVEN      0x0002U
#define CORRELATION_GIVEN 0x0004U
#define BETAGAMMA_GIVEN   0x0008U
#define ALPHAZ_GIVEN      0x0010U
#define CHIRP_GIVEN       0x0020U
} ZPLANE_SPEC;

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC *match,
                           long compute);
long PerformZTransformation(double *t, double *p, double *x, double *xp,
                            double *y, double *yp, long rows, ZPLANE_SPEC *match, long compute);

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);
long LoadTwissFromFile(PLANE_SPEC *spec, long yPlane);

int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout;
  char *inputfile, *outputfile;
  long i_arg, rows, readCode, noWarnings;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  PLANE_SPEC xSpec, ySpec;
  ZPLANE_SPEC zSpec;
  double *x=NULL, *xp=NULL, *y=NULL, *yp=NULL, *p=NULL, *t=NULL;
  long oneTransform, verbose;
  
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) {
    fprintf(stderr, "%s%s\n", USAGE1, USAGE2);
    return(1);
  }

  inputfile = outputfile = NULL;
  pipeFlags = noWarnings = 0;
  xSpec.flags = ySpec.flags = zSpec.flags = 0;
  verbose = oneTransform = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_XPLANE:
        xSpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&xSpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &xSpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &xSpec.alpha, 1, ALPHA_GIVEN,
                          "nemittance", SDDS_DOUBLE, &xSpec.normEmittance, 1, NEMIT_GIVEN,
                          "emittance", SDDS_DOUBLE, &xSpec.emittance, 1, EMIT_GIVEN,
                          "etavalue", SDDS_DOUBLE, &xSpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &xSpec.etap, 1, ETAP_GIVEN,
			  "filename", SDDS_STRING, &xSpec.filename, 1, FILENAME_GIVEN,
			  "element", SDDS_STRING, &xSpec.element, 1, ELEMENT_GIVEN,
			  "occurrence", SDDS_LONG, &xSpec.occurrence, 1, 
			  OCCURRENCE_GIVEN,
                          NULL) ||
            (xSpec.flags&BETA_GIVEN && !(xSpec.flags&ALPHA_GIVEN)) ||
            (!(xSpec.flags&BETA_GIVEN) && xSpec.flags&ALPHA_GIVEN) ||
            (xSpec.flags&BETA_GIVEN && xSpec.beta<=0) || 
            (xSpec.flags&NEMIT_GIVEN && xSpec.normEmittance<0) ||
            (xSpec.flags&EMIT_GIVEN && xSpec.emittance<0) ||
	    ((xSpec.flags&ELEMENT_GIVEN || xSpec.flags&OCCURRENCE_GIVEN) &&
	     !(xSpec.flags&FILENAME_GIVEN)))
          SDDS_Bomb("invalid -xPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
	if (xSpec.flags&FILENAME_GIVEN && !LoadTwissFromFile(&xSpec, 0)) 
	  SDDS_Bomb("invalid -xPlane syntax/values---problem loading data from twiss file");
        break;
      case SET_YPLANE:
        ySpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&ySpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "beta", SDDS_DOUBLE, &ySpec.beta, 1, BETA_GIVEN,
                          "alpha", SDDS_DOUBLE, &ySpec.alpha, 1, ALPHA_GIVEN,
                          "nemittance", SDDS_DOUBLE, &ySpec.normEmittance, 1, NEMIT_GIVEN,
                          "emittance", SDDS_DOUBLE, &ySpec.emittance, 1, EMIT_GIVEN,
                          "etavalue", SDDS_DOUBLE, &ySpec.eta, 1, ETA_GIVEN,
                          "etaslope", SDDS_DOUBLE, &ySpec.etap, 1, ETAP_GIVEN,
			  "filename", SDDS_STRING, &ySpec.filename, 1, FILENAME_GIVEN,
			  "element", SDDS_STRING, &ySpec.element, 1, ELEMENT_GIVEN,
			  "occurrence", SDDS_LONG, &ySpec.occurrence, 1, 
			  OCCURRENCE_GIVEN,
                          NULL) ||
            (ySpec.flags&BETA_GIVEN && !(ySpec.flags&ALPHA_GIVEN)) ||
            (!(ySpec.flags&BETA_GIVEN) && ySpec.flags&ALPHA_GIVEN) ||
            (ySpec.flags&BETA_GIVEN && ySpec.beta<=0) ||
            (ySpec.flags&NEMIT_GIVEN && ySpec.normEmittance<0) ||
            (ySpec.flags&EMIT_GIVEN && ySpec.emittance<0) ||
	    ((ySpec.flags&ELEMENT_GIVEN || ySpec.flags&OCCURRENCE_GIVEN) &&
	     !(ySpec.flags&FILENAME_GIVEN)))
          SDDS_Bomb("invalid -yPlane syntax/values---watch out for abbreviations of etaValue and etaSlope");
	if (ySpec.flags&FILENAME_GIVEN && !LoadTwissFromFile(&ySpec, 1)) 
	  SDDS_Bomb("invalid -yPlane syntax/values---problem loading data from twiss file");
        break;
      case SET_ZPLANE:
        zSpec.flags = 0;
        s_arg[i_arg].n_items--;
        if (!scanItemList(&zSpec.flags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                          "tStDev", SDDS_DOUBLE, &zSpec.tStDev, 1, TSTDEV_GIVEN,
                          "deltaStDev", SDDS_DOUBLE, &zSpec.deltaStDev, 1, DELTASTDEV_GIVEN,
                          "correlation", SDDS_DOUBLE, &zSpec.correlation, 1, CORRELATION_GIVEN,
                          "alpha", SDDS_DOUBLE, &zSpec.alpha, 1, ALPHAZ_GIVEN,
                          "chirp", SDDS_DOUBLE, &zSpec.chirp, 1, CHIRP_GIVEN,
                          "betagamma", SDDS_DOUBLE, &zSpec.betaGamma, 1, BETAGAMMA_GIVEN,
                          NULL) ||
            bitsSet(zSpec.flags)<1 ||
	    bitsSet(zSpec.flags&(ALPHAZ_GIVEN+CORRELATION_GIVEN))>1 ||
            (zSpec.flags&TSTDEV_GIVEN &&  zSpec.tStDev<0) ||
            (zSpec.flags&DELTASTDEV_GIVEN && zSpec.deltaStDev<0) ||
            (zSpec.flags&BETAGAMMA_GIVEN && zSpec.betaGamma<=0))
          SDDS_Bomb("invalid -zPlane syntax/values");
        break;
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_NOWARNINGS:
        noWarnings = 1;
        break;
      case SET_ONETRANSFORM:
        oneTransform = 1;
        break;
      case SET_VERBOSE:
        verbose = 1;
        break;
      default:
        fprintf(stdout, "error: unknown switch: %s\n", s_arg[i_arg].list[0]);
        fflush(stdout);
        exit(1);
        break;
      }
    }
    else {
      if (inputfile==NULL)
        inputfile = s_arg[i_arg].list[0];
      else if (outputfile==NULL)
        outputfile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }

  processFilenames("sddsmatchbeam", &inputfile, &outputfile, pipeFlags, noWarnings, NULL);

  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "x", "m") || !check_sdds_beam_column(&SDDSin, "y", "m") ||
      !check_sdds_beam_column(&SDDSin, "xp", NULL) || !check_sdds_beam_column(&SDDSin, "yp", NULL) ||
      (!check_sdds_beam_column(&SDDSin, "p", "m$be$nc") && !check_sdds_beam_column(&SDDSin, "p", NULL))) {
    fprintf(stderr, 
            "sddsmatchtwiss: one or more data quantities (x, xp, y, yp, p) have the wrong units or are not present in %s", 
            inputfile);
    exit(1);
  }

  if (!SDDS_InitializeCopy(&SDDSout, &SDDSin, outputfile, "w") ||
      !SDDS_SaveLayout(&SDDSout) || !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if ((rows=SDDS_RowCount(&SDDSin))==0) {
      if (!SDDS_CopyPage(&SDDSin, &SDDSout) || !SDDS_WritePage(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      continue;
    }
    if (!(x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
        !(xp = SDDS_GetColumnInDoubles(&SDDSin, "xp")) ||
        !(y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
        !(yp = SDDS_GetColumnInDoubles(&SDDSin, "yp")) ||
        !(t = SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
        !(p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (xSpec.flags)
      PerformTransformation(x, xp, p, rows, &xSpec, readCode==1?1:!oneTransform);
    if (ySpec.flags)
      PerformTransformation(y, yp, p, rows, &ySpec, readCode==1?1:!oneTransform);
    if (zSpec.flags)
      PerformZTransformation(t, p, x, xp, y, yp, rows, &zSpec, readCode==1?1:!oneTransform);
    if (verbose) {
      if (xSpec.flags)
        fprintf(stderr, "x transformation: %le, %le, %le, %le, %le, %le\n",
                xSpec.R11, xSpec.R12, xSpec.R21, xSpec.R22, xSpec.etaBeam,
                xSpec.etapBeam);
      if (ySpec.flags)
        fprintf(stderr, "y transformation: %le, %le, %le, %le, %le, %le\n",
                ySpec.R11, ySpec.R12, ySpec.R21, ySpec.R22, ySpec.etaBeam,
                ySpec.etapBeam);
      if (zSpec.flags)
        fprintf(stderr, "z transformation: %le, %le, %le, %le\n",
                zSpec.R11, zSpec.R12, zSpec.R21, zSpec.R22);
    }
    if (!SDDS_CopyPage(&SDDSout, &SDDSin) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, x, rows, "x")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, y, rows, "y")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, xp, rows, "xp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, yp, rows, "yp")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, t, rows, "t")) ||
        !(SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME, p, rows, "p")) ||
        !SDDS_WritePage(&SDDSout))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    free(x);
    free(xp);
    free(y);
    free(yp);
    free(t);
    free(p);
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 0;
}

long PerformTransformation(double *x, double *xp, double *p, long rows, PLANE_SPEC *match,
                           long computeTransform)
{
  long i;
  double pAve;
  double S11, S22, S12, S16, S26, S66;
  double R11, R12, R21, R22;
  double x0, xp0;
  double emit, beta1, beta2, alpha1, alpha2;
  double eta1, etap1, eta2, etap2;

  pAve = arithmeticAverage(p, rows);
  for (i=0; i<rows; i++)
    p[i] = p[i]/pAve - 1;

  if (match->flags&EMIT_GIVEN) {
    match->flags |= NEMIT_GIVEN;
    match->normEmittance = match->emittance*pAve;
  }

  if (computeTransform) {
    computeCorrelations(&S11, &S16, &S66, x, p, rows);
    if (S66)
      match->etaBeam = eta1 = S16/S66;
    else
      match->etaBeam = eta1 = 0;
    computeCorrelations(&S22, &S26, &S66, xp, p, rows);
    if (S66)
      match->etapBeam = etap1 = S26/S66;
    else
      match->etapBeam = etap1 = 0;
  } else {
    eta1 = match->etaBeam;
    etap1 = match->etapBeam;
  }
  
  for (i=0; i<rows; i++)
    x[i] -= p[i]*eta1;
  for (i=0; i<rows; i++)
    xp[i] -= p[i]*etap1;
  
  if (match->flags&BETA_GIVEN || match->flags&NEMIT_GIVEN) {
    if (computeTransform || match->flags&NEMIT_GIVEN) {
      computeCorrelations(&S11, &S12, &S22, x, xp,rows);
      if ((emit = S11*S22-sqr(S12))<=0) {
        SDDS_Bomb("emittance is zero");
      }
      else
        emit = sqrt(emit);
      beta1 = S11/emit;
      alpha1 = -S12/emit;
      if (match->flags&BETA_GIVEN) {
        beta2 = match->beta;
        alpha2 = match->alpha;
      }
      else {
        beta2 = beta1;
        alpha2 = alpha1;
      }
      match->R11 = R11 = beta2/sqrt(beta1*beta2);
      match->R12 = R12 = 0;
      match->R21 = R21 = (alpha1-alpha2)/sqrt(beta1*beta2);
      match->R22 = R22 = beta1/sqrt(beta1*beta2);
      if (match->flags&NEMIT_GIVEN) {
	double factor;
	factor = sqrt(match->normEmittance/(emit*pAve));
	match->R11 = (R11 *= factor);
	match->R12 = (R12 *= factor);
	match->R22 = (R22 *= factor);
	match->R21 = (R21 *= factor);
      }
    }
    else {
      R11 = match->R11;
      R12 = match->R12;
      R21 = match->R21;
      R22 = match->R22;
    }
    for (i=0; i<rows; i++) {
      x0 = x[i];
      xp0 = xp[i];
      x[i] = R11*x0 + R12*xp0;
      xp[i] = R21*x0 + R22*xp0;
    }
  }

  if (match->flags&ETA_GIVEN)
    eta2 = match->eta;
  else
    eta2 = eta1;
  if (match->flags&ETAP_GIVEN)
    etap2 = match->etap;
  else
    etap2 = etap1;
  for (i=0; i<rows; i++) {
    x[i] += eta2*p[i];
    xp[i] += etap2*p[i];
    p[i] = (p[i]+1)*pAve;
  }
  return 1;
}

long PerformZTransformation(double *t, double *p, 
                            double *x, double *xp,
                            double *y, double *yp,
                            long rows, ZPLANE_SPEC *match,
                            long computeTransform)
{
  long i;
  double pAve, tAve;
  double S11, S22, S12;
  double R11, R12, R21, R22;
  double delta0, t0;
  double emit1, beta1, beta2, alpha1, alpha2;
  double emit2, ratio;

  pAve = arithmeticAverage(p, rows);
  for (i=0; i<rows; i++)
    p[i] = p[i]/pAve - 1;
  tAve = arithmeticAverage(t, rows);
  for (i=0; i<rows; i++)
    t[i] = t[i] - tAve;

  if (computeTransform) {
    computeCorrelations(&S11, &S12, &S22, t, p, rows);
    if ((emit1 = S11*S22-sqr(S12))<=0)
      SDDS_Bomb("longitudinal emittance is zero");
    else
      emit1 = sqrt(emit1);
    beta1 = S11/emit1;
    alpha1 = -S12/emit1;

    if (match->flags&TSTDEV_GIVEN)
      S11 = sqr(match->tStDev);
    if (match->flags&DELTASTDEV_GIVEN)
      S22 = sqr(match->deltaStDev);
    if (match->flags&CORRELATION_GIVEN) 
      S12 = match->correlation*sqrt(S11*S22);
    else if (match->flags&ALPHAZ_GIVEN)
      S12 = -match->alpha*sqrt(S11*S22)/sqrt(1+sqr(match->alpha));
    if ((emit2 = S11*S22-sqr(S12))<=0)
      SDDS_Bomb("longitudinal emittance is zero");
    else
      emit2 = sqrt(emit2);
    beta2 = S11/emit2;
    if (!(match->flags&ALPHAZ_GIVEN))
      alpha2 = -S12/emit2;
    else
      alpha2 = match->alpha;
    ratio = sqrt(emit2/emit1);
    
    match->R11 = R11 = ratio*beta2/sqrt(beta1*beta2);
    match->R12 = R12 = 0;
    match->R21 = R21 = ratio*(alpha1-alpha2)/sqrt(beta1*beta2);
    match->R22 = R22 = ratio*beta1/sqrt(beta1*beta2);
  }
  else {
    R11 = match->R11;
    R12 = match->R12;
    R21 = match->R21;
    R22 = match->R22;
  }

  for (i=0; i<rows; i++) {
    t0 = t[i];
    delta0 = p[i];
    t[i] = R11*t0 + R12*delta0;
    p[i] = R21*t0 + R22*delta0;
  }

  if (match->flags&CHIRP_GIVEN) {
    for (i=0; i<rows; i++)
      p[i] += t[i]*match->chirp;
  }

  if (match->flags&BETAGAMMA_GIVEN) {
    double ratio;
    ratio = sqrt(pAve/match->betaGamma);
    for (i=0; i<rows; i++) {
      x[i] *= ratio;
      xp[i] *= ratio;
      y[i] *= ratio;
      yp[i] *= ratio;
      t[i] += tAve;
      p[i] = p[i]*pAve + match->betaGamma;
    }
  } else {
    for (i=0; i<rows; i++) {
      p[i] = (p[i]+1)*pAve;
      t[i] += tAve;
    }
  }
  
  return 1;
}

long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units)
{
  char *units1;
  if (SDDS_GetColumnIndex(SDDS_table, name)<0)
    return(0);
  if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name)!=SDDS_STRING) {
    SDDS_SetError("units field of column has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!units) {
    if (!units1)
      return(1);
    if (SDDS_StringIsBlank(units1)) {
      free(units1);
      return(1);
    }
    return(0);
  }
  if (!units1)
    return(0);
  if (strcmp(units, units1)==0) {
    free(units1);
    return(1);
  }
  free(units1);
  return(0);
}

long LoadTwissFromFile(PLANE_SPEC *spec, long yPlane)
{
  SDDS_DATASET SDDSin;
  long rows=0, rowOfInterest;
  double *betaData=NULL, *alphaData=NULL;
  double *etaData=NULL, *etapData=NULL;
  char *name[8] = {"betax", "alphax", "etax", "etaxp",
		   "betay", "alphay", "etay", "etayp"};

  if (!SDDS_InitializeInput(&SDDSin, spec->filename) || 
      SDDS_ReadPage(&SDDSin)!=1) {
    SDDS_SetError("problem reading Twiss reference file");
    return 0;
  }
  if (SDDS_CheckColumn(&SDDSin, name[4*yPlane+0], "m", 
		       SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, name[4*yPlane+1], NULL, 
		       SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, name[4*yPlane+2], "m", 
		       SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, name[4*yPlane+3], NULL, 
		       SDDS_ANY_FLOATING_TYPE, stdout)!=SDDS_CHECK_OK ||
      SDDS_CheckColumn(&SDDSin, "ElementName", NULL, 
		       SDDS_STRING, stdout)!=SDDS_CHECK_OK) {
    SDDS_SetError("invalid/missing columns in Twiss reference file");
    return 0;
  }
  if (spec->flags&ELEMENT_GIVEN) {
    if (!SDDS_SetRowFlags(&SDDSin, 1) ||
        (rows=SDDS_MatchRowsOfInterest(&SDDSin, "ElementName", spec->element,
				       SDDS_AND))<=0) {
      SDDS_SetError("Problem finding data for beta function reference.  Check for existence of element.");
      return 0;
    }
    if (spec->occurrence>0 && spec->occurrence>rows) {
      SDDS_SetError("Too few occurrences of reference element in beta function reference file.");
      return 0;
    }
  } 
  if ((rows=SDDS_CountRowsOfInterest(&SDDSin))<1) {
    SDDS_SetError("No data in beta function reference file.");
    return 0;
  }
  if (!(betaData=SDDS_GetColumnInDoubles(&SDDSin, name[4*yPlane+0])) ||
      !(alphaData=SDDS_GetColumnInDoubles(&SDDSin, name[4*yPlane+1])) ||
      !(etaData=SDDS_GetColumnInDoubles(&SDDSin, name[4*yPlane+2])) ||
      !(etapData=SDDS_GetColumnInDoubles(&SDDSin, name[4*yPlane+3])) ) {
    SDDS_SetError("Problem getting data for beta function reference.");
    return 0;
  }
  if (spec->flags&ELEMENT_GIVEN && spec->element && spec->occurrence>0)
    rowOfInterest = spec->occurrence-1;
  else
    rowOfInterest = rows-1;
  if (!(spec->flags&BETA_GIVEN))
    spec->beta = betaData[rowOfInterest];
  if (!(spec->flags&ALPHA_GIVEN))
    spec->alpha = alphaData[rowOfInterest];
  if (!(spec->flags&ETA_GIVEN))
    spec->eta = etaData[rowOfInterest];
  if (!(spec->flags&ETAP_GIVEN))
    spec->etap = etapData[rowOfInterest];
  spec->flags |= BETA_GIVEN+ALPHA_GIVEN+ETA_GIVEN+ETAP_GIVEN;

  free(betaData);
  free(alphaData);
  free(etaData);
  free(etapData);
  return 1;
}

