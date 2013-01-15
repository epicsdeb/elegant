/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsemitmeas
 * purpose: process elegant runs to determine emittance of beam.
 *          Also processes experimental data using elegant-computed
 *          matrices.
 * Michael Borland, 1997.
 */
/*
 * $Log: not supported by cvs2svn $
 * Revision 1.1  2007/03/30 16:50:29  soliday
 * Moved from directory above.
 *
 * Revision 1.6  2005/12/15 00:36:04  borland
 * Fixed bug in parsing error settings.  Fixed bug in setting output parameters
 * (S33 not set).
 *
 * Revision 1.5  2005/11/22 23:21:20  borland
 * Added momentum aperture search, which necessitated adding an argument to
 * do_tracking, resulting in changes in many files.
 * Also improved convergence of orbit finder, adding a second iteration using
 * tracking if the matrix-based method fails.
 *
 * Revision 1.4  2005/11/10 15:38:49  soliday
 * Added changes to get it to compile properly with 64 bit compilers.
 *
 * Revision 1.3  2004/12/30 04:39:04  borland
 * Standard deviations in the output are no longer labeled "Sigma".
 *
 * Revision 1.2  2004/12/30 03:47:07  borland
 * Units of emittance are just m, not pi*m.
 *
 * Revision 1.1  2004/12/30 03:29:45  borland
 * Simplified version of sddsemitmeas, with only the Monte-Carlo-based
 * error analysis.  More reliable and preferred over sddsemitmeas.
 *
 * Revision 1.19  2004/12/29 23:33:56  borland
 * Fixed bug in previous change.
 *
 * Revision 1.18  2004/12/29 23:31:19  borland
 * Now outputs the sigma data and fit even with errors.  In this case, any
 * first perturbation is provided.
 *
 * Revision 1.17  2003/02/10 20:18:55  borland
 * Will now use Step as the independent variable if no quad-like names are
 * found.
 *
 * Revision 1.16  2002/08/14 20:23:48  soliday
 * Added Open License
 *
 * Revision 1.15  2002/04/25 21:27:52  borland
 * Fixed an error message. SHould now fit y plane even if x plane fit fails.
 *
 * Revision 1.14  2001/10/15 20:37:06  soliday
 * Cleaned up for Linux.
 *
 * Revision 1.13  2000/08/17 13:45:41  borland
 * Fixed error in usage message for -errorLevels .
 *
 * Revision 1.12  2000/08/02 01:36:17  borland
 * Fixed an error message that was sent to stdout.
 *
 * Revision 1.11  2000/05/12 20:59:30  borland
 * Made some changes to try to accomodate cause when there is a matrix error.
 * Not successful.
 *
 * Revision 1.10  1999/10/12 21:50:01  borland
 * All printouts now go to the stdout rather than stderr.  fflush statements,
 * some unnecessary, were added in a mostly automated fashion.
 *
 * Revision 1.9  1999/08/05 15:41:23  soliday
 * Added WIN32 and Linux support
 *
 * Revision 1.8  1999/07/02 15:57:55  borland
 * Replaced all instances of the "HUGE" macro with "DBL_MAX".  The former
 * is only defined on Solaris.
 *
 * Revision 1.7  1998/12/23 19:47:51  borland
 * Improved robustness of error handling when, for example, the orbit or
 * tunes can't be found.  By default, this is turned off, but can be
 * activated by using the soft_failure=0 flag in the track namelist.
 * For the find_aperture and analyze_map features, this flag is not
 * yet implemented, so the behavior has changed.
 *
 * Revision 1.6  1998/11/16 22:18:27  borland
 * Added string parameter with emittance label to output file.
 *
 * Revision 1.5  1998/10/14 17:12:57  borland
 * Modified the output file.  If there are no error sets requested, then the
 * file contains the results (emittances, etc.) as parameters and the
 * data+fit as columns.  If there are error sets requested, the file contains
 * just results as parameters.
 *
 * Revision 1.4  1998/08/11 19:47:39  borland
 * Fixed bug with variable name (obsolete option that is still needed
 * internally).  Should remove this.
 *
 * Revision 1.3  1998/08/11 18:53:43  borland
 * Initialized pipeFlags in main() to prevent error in filename parsing.
 *
 * Revision 1.2  1998/02/26 16:16:29  borland
 * This version works.
 *
 * Revision 1.1  1997/09/25 19:32:25  borland
 * First version in repository.  Doesn't work.  Probably doesn't even compile.
 *
 *
 */
#include "mdb.h"
#include "matlib.h"
#include "match_string.h"
#include "SDDS.h"
#include "scan.h"
#include "table.h"
#include <time.h>
#include <memory.h>

#define SET_ERROR_LEVEL 0
#define SET_N_ERROR_SETS 1
#define SET_SEED 2
#define SET_VERBOSITY 3
#define SET_DEVIATION_LIMIT 4
#define SET_RESOLUTION 5
#define SET_LIMIT_MODE 6
#define SET_IGNORE_PLANE 7
#define SET_SIGMA_DATA 8
#define SET_PIPE 9
#define SET_INCLUDE_DISPERSION 10
#define SET_ERROR_DATA 11
#define N_OPTIONS 12

char *option[N_OPTIONS] = {
    "errorlevel", "nerrorsets", "seed", "verbosity",
    "deviationlimit", "resolution", "limitmode", "ignoreplane",
    "sigmadata", "pipe", "includedispersion", "errordata",
    } ;

#define USAGE "sddsemitproc\n\
 [<inputfile>] [<outputfile>] [-pipe=[input][,output]]\n\
 -sigmaData=<xName>,<yName> [-includeDispersion[=vertical]] \n\
 [-errorData=<xName>,<yName> | \n\
 [-errorLevel=<x_valueInm>,<y_valueInm>,[{gaussian,<nSigmas> | uniform}]]\n\
 [-nErrorSets=<number>]\n\
 [-limitMode={resolution | zero}[{,reject}]\n\
 [-deviationLimit=<xLevelm>,<yLevelm>]\n\
 [-resolution=<xResolutionm>,<yResolutionm>]\n\
 [-seed=integer] [-ignorePlane={x | y}]\n\
 [-verbosity=level]\n\n\
Program by Michael Borland. (This is version 2, November 2005.)"

static char *additional_help[] = {
USAGE,
"\nThis program computes the emittance from simulated or experimental",
"beamsize data, using simulated data for response matrices.",
"inputfile is an elegant \"final\" parameters output file (SDDS format),",
"which contains the R matrix and simulated beam sigmas as a function",
"of some variable or variables.  It may also contain experimental data",
"that has been obtained separately and put into the elegant output file",
"(e.g., using sddsxref). \n",
"-sigmaData is used to name the columns in <inputfile> that contain the",
"    beam sizes to be used in fitting.",
"-errorData and -errorLevel are two ways to specify the rms errors in the",
"    measurements.  -nErrorSets specifies the number of ensembles of errors",
"    to add to the sigmas. -deviationLimit allows exclusion from the ",
"    fit of bad data.  -limitMode allows you to specify what is done",
"    with data at or below the resolution limit.",
"-resolution allows specification of the measurement resolution,",
"    which is subtracted in quadrature from the sigma or width values.",
"-includeDispersion is used if there is dispersion in the system. For this",
"    to work, you must have a dipole in the beamline that you are simulating.",
"    Otherwise, the program cannot separate dispersive effects from emittance.",
NULL
    } ; 

#define N_INCREMENT 10

double solve_normal_form(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double *s2_fit);
double solve_normal_form_opt(MATRIX *F, MATRIX *sF, MATRIX *P, MATRIX *M, MATRIX *C, double dev_limit,
    long *n_used, double *s2_fit);
void print_fit(char *filename, double *variable_data, char *variable_name,
    double *sigma2_fit, double resol, char *sigma_name, long n_pts);
double propagate_errors_for_emittance(double **Sigma, double **Covar);
void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long n_configs, long equal_weights);
double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, MATRIX *K, 
    double dev_limit, long n_configs, double uncert_min, double *fit_sig2_return);
long make_tweeked_data_set(MATRIX *s2, double *sigma, double *errorData, double error_level, double error_sigmas, long error_type_code, 
    long n_configs, double resol, long reject_at_limit, double limit, long *n_at_resol);
long SetSigmaData(SDDS_DATASET *SDDSout, char *dataName, MATRIX *s2, char *fitName, double *fitSigSqr, 
                  long configs);

#define GAUSSIAN_ERRORS 0
#define UNIFORM_ERRORS  1
#define N_ERROR_TYPES   2
char *error_type[N_ERROR_TYPES] = {
    "gaussian", "uniform"
    } ;

#define LIMIT_AT_RESOLUTION 0
#define LIMIT_AT_ZERO       1
#define N_LIMIT_OPTIONS     2
char *limit_option[N_LIMIT_OPTIONS] = {
    "resolution", "zero"
    } ;

#define MAX_N_TRIES 100

int main(
     int argc,
     char **argv
     )
{
  SDDS_TABLE SDDSin, SDDSout;
  double *R11;        /* R11 matrix element for ith configuration */
  double *R12=NULL;   /* R12 matrix element for ith configuration */
  double *R16=NULL;   /* R16 matrix element for ith configuration */
  double *sigmax, *uncertx;               /* sigma in x plane for ith configuration */
  double *R33=NULL, *R34=NULL, *R36=NULL, *sigmay, *uncerty;   /* similar data for y plane */
  long n_configs, i_config;
  MATRIX *Rx, *Ry, *s2x, *s2y;
  MATRIX *Sx, *Sy, *sSx, *sSy, *Kx, *Ky;
  double S11_sum, S11_sum2, S12_sum, S12_sum2, S22_sum, S22_sum2;
  double S66_sum, S66_sum2, S16_sum, S16_sum2, S26_sum, S26_sum2;
  double S33_sum, S33_sum2, S34_sum, S34_sum2, S44_sum, S44_sum2;
  double S36_sum, S46_sum, S36_sum2, S46_sum2;
  double betax, alphax, betax_sum, betax_sum2, alphax_sum, alphax_sum2;
  double betay, alphay, betay_sum, betay_sum2, alphay_sum, alphay_sum2;
  long i_variable;
  double emitx, emity;
  SCANNED_ARG *scanned;
  long i_arg, i;
  char *input, *output;
  char *variable_name;
  double *variable_data, *x_fit_sig2, *y_fit_sig2;
  long n_error_sets, seed, i_error, error_type_code;
  double error_sigmas;
  double emitx_max, emitx_min;
  double emity_max, emity_min;
  double emitx_sum, emity_sum;
  double emitx2_sum, emity2_sum;
  double x_error_level, y_error_level;
  double md_x, md_y, *dev_limit, contrib;
  long i_dev, n_dev_limits;
  long nx_used, nx_used_sum;
  long ny_used, ny_used_sum;
  long n_good_fits_x, n_good_fits_y;
  double x_resol, y_resol;
  double x_limit=0.0, y_limit=0.0;
  long limit_code, reject_at_limit;
  long n_xresol, n_yresol, ignore_x, ignore_y, error_output;
  double x_uncert_min, y_uncert_min;
  double includeDispersion;
  long equal_weights_x_fit=0, equal_weights_y_fit=0;
  long verbosity;
  char *x_width_name, *y_width_name;
  char *x_error_name, *y_error_name;
  unsigned long pipeFlags;
  char emitLabel[256], xEmitLabel[256], yEmitLabel[256];
  
  argc = scanargs(&scanned, argc, argv);
  if (argc<2 || argc>(2+N_OPTIONS)) {
    for (i=0; ; i++) {
      if (!additional_help[i])
        break;
      puts(additional_help[i]);
    }
    exit(1);
  }

  input = output = NULL;
  dev_limit = NULL;
  n_error_sets = 1;
  n_dev_limits = 0;
  x_error_level = y_error_level = 0;
  seed = -1;
  variable_name = NULL;
  x_resol = y_resol = 0;
  error_type_code = GAUSSIAN_ERRORS;
  error_sigmas = 2;
  limit_code = LIMIT_AT_ZERO;
  reject_at_limit = 0;
  ignore_x = ignore_y = 0;
  x_uncert_min = y_uncert_min = 0;
  verbosity = 0;
  x_width_name = "Sx";
  y_width_name = "Sy";
  x_error_name = y_error_name = NULL;
  includeDispersion = 0;
  pipeFlags = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (scanned[i_arg].arg_type==OPTION) {
      switch (match_string(scanned[i_arg].list[0], option,
                           N_OPTIONS, 0)) {
        /* process options here */
      case SET_ERROR_LEVEL:
        if (scanned[i_arg].n_items<3)
          bomb("invalid -errorLevel syntax (too few qualifiers)", USAGE);
        if (!sscanf(scanned[i_arg].list[1], "%lf", &x_error_level) ||
            !sscanf(scanned[i_arg].list[2], "%lf", &y_error_level)) 
          bomb("invalid -errorLevel syntax (invalid x or y level)", USAGE);
        if (scanned[i_arg].n_items>=4 && scanned[i_arg].n_items<=5) {
          if ((error_type_code=match_string(scanned[i_arg].list[3], error_type, N_ERROR_TYPES, 0))<0)
            bomb("unknown error type (invalid error type)", USAGE);
          if (scanned[i_arg].n_items==5) {
            if (error_type_code!=GAUSSIAN_ERRORS || 
                !sscanf(scanned[i_arg].list[4], "%lf", &error_sigmas) ||
                error_sigmas<0)
              bomb("invalid -errorLevel syntax (invalid error sigmas level)", USAGE);
          }
        }
        break;
      case SET_N_ERROR_SETS:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &n_error_sets) ||
            n_error_sets<0)
          bomb("invalid -n_error_sets syntax", USAGE);
        break;
      case SET_SEED:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &seed) ||
            seed<0)
          bomb("invalid -seed syntax", USAGE);
        break;
      case SET_DEVIATION_LIMIT:
        if ((n_dev_limits=scanned[i_arg].n_items-1)<1)
          bomb("invalid -deviation_limit syntax", USAGE);
        dev_limit = tmalloc(sizeof(double)*n_dev_limits);
        for (i=0; i<n_dev_limits; i++) {
          if (!sscanf(scanned[i_arg].list[i+1], "%lf", dev_limit+i) ||
              dev_limit[i]==0) 
            bomb("invalid -deviation_limit syntax", USAGE);
        }
        break;
      case SET_RESOLUTION:
        if (scanned[i_arg].n_items!=3 ||
            !sscanf(scanned[i_arg].list[1], "%lf", &x_resol) ||
            !sscanf(scanned[i_arg].list[2], "%lf", &y_resol) ||
            x_resol<0 || y_resol<0) 
          bomb("invalid -resolution syntax", USAGE);
        break;
      case SET_LIMIT_MODE:
        if (scanned[i_arg].n_items<2 || scanned[i_arg].n_items>4 ||
            (limit_code=match_string(scanned[i_arg].list[1], limit_option,
                                     N_LIMIT_OPTIONS, 0))<0)
          bomb("invalid -limit_mode syntax", USAGE);
        if (scanned[i_arg].n_items==3) {
          if (scanned[i_arg].list[2][0]=='r')
            reject_at_limit = 1;
          else
            bomb("invalid -limit_mode syntax", USAGE);
        }
        break;
      case SET_IGNORE_PLANE:
        if (scanned[i_arg].n_items!=2 ||
            !((ignore_x=(scanned[i_arg].list[1][0]=='x')) ||
              (ignore_y=(scanned[i_arg].list[1][0]=='y')) ) )
          bomb("invalid -ignore_plane syntax", USAGE);
        break;
      case SET_VERBOSITY:
        if (scanned[i_arg].n_items!=2 ||
            !sscanf(scanned[i_arg].list[1], "%ld", &verbosity) ||
            verbosity<0)
          bomb("invalid -verbosity syntax", USAGE);
        break;
      case SET_SIGMA_DATA:
        if (scanned[i_arg].n_items!=3 ||
            !strlen(x_width_name = scanned[i_arg].list[1]) ||
            !strlen(y_width_name = scanned[i_arg].list[2]))
          bomb("invalid -sigmaData syntax", USAGE);
        break;
      case SET_ERROR_DATA:
        if (scanned[i_arg].n_items!=3 ||
            !strlen(x_error_name = scanned[i_arg].list[1]) ||
            !strlen(y_error_name = scanned[i_arg].list[2]))
          bomb("invalid -errorData syntax", USAGE);
        break;
      case SET_PIPE:
	if (!processPipeOption(scanned[i_arg].list+1, scanned[i_arg].n_items-1, &pipeFlags))
	  SDDS_Bomb("invalid -pipe syntax");
	break;
      case SET_INCLUDE_DISPERSION:
        includeDispersion = 1;
        if (scanned[i_arg].n_items==2) {
          if (strncmp(scanned[i_arg].list[1], "vertical", strlen(scanned[i_arg].list[1]))==0)
            includeDispersion = 2;
          else
            SDDS_Bomb("invalid -includeDispersion syntax");
         } else if (scanned[i_arg].n_items>2)  
           SDDS_Bomb("invalid -includeDispersion syntax");
        SDDS_Bomb("The -includeDispersion option isn't reliable.  Suggest subtracting off energy spread contribution to beam size manually.");
        break;
      default:
        bomb("unknown option given", USAGE);
        break;
      }
    }
    else {
      if (!input)
        input = scanned[i_arg].list[0];
      else if (!output)
        output = scanned[i_arg].list[0];
      else
        bomb("too many filenames given", USAGE);
    }
  }

  processFilenames("sddsemitmeas", &input, &output, pipeFlags, 0, NULL);

  if (ignore_x && ignore_y)
    bomb("can't ignore both x and y planes!", USAGE);

  x_limit = (limit_code==LIMIT_AT_ZERO?0.0:x_resol);
  y_limit = (limit_code==LIMIT_AT_ZERO?0.0:y_resol);

  if (seed<0) {
    /* generate seed from system clock */
    seed = (int)time((time_t)0);
  }
  random_1(-seed);

  if (n_error_sets<=2)  
    SDDS_Bomb("number of error sets must be >2");

  if (n_dev_limits<=0) {
    dev_limit = tmalloc(sizeof(double));
    n_dev_limits = 1;
    dev_limit[0] = 0;
  }

  error_output = 1;
  
  if (!SDDS_InitializeInput(&SDDSin, input))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 1, NULL, NULL, output))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  if (!ignore_x) {
    if (!SDDS_DefineSimpleParameter(&SDDSout, "ex", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S11", "m$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S12", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S22", "", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "betax", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "alphax", "", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "xRMSDeviation", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "xGoodFits", "", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "xAverageFitPoints", "", SDDS_DOUBLE) ||
        (error_output && 
         (!SDDS_DefineSimpleParameter(&SDDSout, "exStDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S11StDev", "m$a2$n", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S12StDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S22StDev", "", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "betaxStDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "alphaxStDev", "", SDDS_DOUBLE)))) 
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (!SDDS_DefineSimpleColumn(&SDDSout, "xSigmaData", "m", SDDS_DOUBLE) ||
         !SDDS_DefineSimpleColumn(&SDDSout, "xSigmaFit", "m", SDDS_DOUBLE)) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  if (!ignore_y) {
    if (!SDDS_DefineSimpleParameter(&SDDSout, "ey", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S33", "m$a2$n", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S34", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "S44", "", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "betay", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "alphay", "", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "yRMSDeviation", "m", SDDS_DOUBLE) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "yGoodFits", "", SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(&SDDSout, "yAverageFitPoints", "", SDDS_DOUBLE) ||
        (error_output &&
         (!SDDS_DefineSimpleParameter(&SDDSout, "eyStDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S33StDev", "m$a2$n", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S34StDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "S44StDev", "", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "betayStDev", "m", SDDS_DOUBLE) ||
          !SDDS_DefineSimpleParameter(&SDDSout, "alphayStDev", "", SDDS_DOUBLE))))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (!SDDS_DefineSimpleColumn(&SDDSout, "ySigmaData", "m", SDDS_DOUBLE) ||
         !SDDS_DefineSimpleColumn(&SDDSout, "ySigmaFit", "m", SDDS_DOUBLE)) {
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
  }
  if (includeDispersion && 
      (!SDDS_DefineSimpleParameter(&SDDSout, "Sdelta", NULL, SDDS_DOUBLE) ||
       (error_output && !SDDS_DefineSimpleParameter(&SDDSout, "SdeltaSigma", NULL, SDDS_DOUBLE))))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (SDDS_GetColumnIndex(&SDDSin, x_width_name)<0 ||
      SDDS_GetColumnIndex(&SDDSin, y_width_name)<0 ||
      SDDS_GetColumnIndex(&SDDSin, "R11")<0 ||
      SDDS_GetColumnIndex(&SDDSin, "R12")<0 ||
      SDDS_GetColumnIndex(&SDDSin, "R33")<0 ||
      SDDS_GetColumnIndex(&SDDSin, "R34")<0)
    SDDS_Bomb("input file missing required quantities");
  if (x_error_name) {
    if (SDDS_GetColumnIndex(&SDDSin, x_error_name)<0 ||
        SDDS_GetColumnIndex(&SDDSin, y_error_name)<0)
      SDDS_Bomb("input file missing error quantities");
  }
  if (includeDispersion &&
      (SDDS_GetColumnIndex(&SDDSin, "R16")<0 ||
       SDDS_GetColumnIndex(&SDDSin, "R36")<0))
      SDDS_Bomb("input file missing R16 or R36");
  
  if (variable_name) {
    if ((i_variable=SDDS_GetColumnIndex(&SDDSin, variable_name))<0)
      bomb("no match for variable for fit output/filtering", NULL);
  }
  else {
    char **column_name;
    int32_t column_names;
    if (!(column_name = SDDS_GetColumnNames(&SDDSin, &column_names)))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    for (i_variable=0; i_variable<column_names; i_variable++)
      if (column_name[i_variable][0]=='Q')
        break;
    if (i_variable==column_names) {
      for (i_variable=0; i_variable<column_names; i_variable++)
      if (strchr(column_name[i_variable], 'Q'))
        break;
    }
    if (i_variable==column_names) {
      for (i_variable=0; i_variable<column_names; i_variable++)
      if (strcmp(column_name[i_variable], "Step")==0)
        break;
    }
    if (i_variable==column_names)
      bomb("you did not specify -variable_name, and there is no obvious choice in the input data", NULL);
    variable_name = column_name[i_variable];
  }
  if (!SDDS_TransferColumnDefinition(&SDDSout, &SDDSin, variable_name, NULL) ||
      !SDDS_DefineSimpleParameter(&SDDSout, "EmittanceLabel", "", SDDS_STRING) ||
      !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  uncertx = uncerty = x_fit_sig2 = y_fit_sig2 = NULL;
  while (SDDS_ReadTable(&SDDSin)>0) {
    n_configs = SDDS_CountRowsOfInterest(&SDDSin);
    if (!(R11 = SDDS_GetColumn(&SDDSin, "R11")) ||
        !(R12 = SDDS_GetColumn(&SDDSin, "R12")) ||
        !(R33 = SDDS_GetColumn(&SDDSin, "R33")) ||
        !(R34 = SDDS_GetColumn(&SDDSin, "R34")) ) 
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (includeDispersion &&
        (!(R16 = SDDS_GetColumn(&SDDSin, "R16")) ||
         !(R36 = SDDS_GetColumn(&SDDSin, "R36"))))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    if (!(sigmax = SDDS_GetColumn(&SDDSin, x_width_name)))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    if (!(sigmay = SDDS_GetColumn(&SDDSin, y_width_name)))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    if (!(variable_data = SDDS_GetColumn(&SDDSin, variable_name)))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    if (!x_error_name) {
      uncertx = SDDS_Realloc(uncertx, sizeof(double)*n_configs);
      uncerty = SDDS_Realloc(uncerty, sizeof(double)*n_configs);
      for (i_config=0 ; i_config<n_configs; i_config++) {
        uncertx[i_config] = uncerty[i_config] = 0;
      }
    } else {
      if (uncertx)
        free(uncertx);
      if (uncerty)
        free(uncerty);
      if (!(uncertx = SDDS_GetColumn(&SDDSin, x_error_name)))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      if (!(uncerty = SDDS_GetColumn(&SDDSin, y_error_name)))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    } 

    x_fit_sig2 = SDDS_Realloc(x_fit_sig2, sizeof(double)*n_configs);
    y_fit_sig2 = SDDS_Realloc(y_fit_sig2, sizeof(double)*n_configs);

    if (n_configs<4)
      continue;

    if (!SDDS_StartPage(&SDDSout, n_configs) || !SDDS_CopyColumns(&SDDSout, &SDDSin)) 
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    if (!ignore_x) {
      if (includeDispersion==1) {
        for (i_config=0; i_config<n_configs; i_config++) 
          if (R16[i_config]) 
            break;
        if (i_config==n_configs)
          SDDS_Bomb("you asked to include dispersion, but R16=0 for all configurations");
      }
      if (!x_error_name) {
        equal_weights_x_fit = 0;
        for (i_config=0; i_config<n_configs; i_config++)
          uncertx[i_config] = 1;
      }
    }
    if (!ignore_y) {
      if (includeDispersion==2) {
        for (i_config=0; i_config<n_configs; i_config++) 
          if (R36[i_config]) 
            break;
        if (i_config==n_configs)
          SDDS_Bomb("you asked to include vertical dispersion, but R36=0 for all configurations");
      }
      if (!y_error_name) {
        equal_weights_y_fit = 0;
        for (i_config=0; i_config<n_configs; i_config++)
          uncerty[i_config] = 1;
      }
    }

    m_alloc(&Rx,  n_configs, includeDispersion==1?6:3);
    m_alloc(&Ry,  n_configs, includeDispersion==2?6:3);
    m_alloc(&s2x, n_configs, 1);
    m_alloc(&s2y, n_configs, 1);
    m_alloc(&Sx,  includeDispersion==1?6:3, 1);
    m_alloc(&Sy,  includeDispersion==2?6:3, 1);
    m_alloc(&sSx,  includeDispersion==1?6:3, includeDispersion==1?6:3); 
    m_alloc(&sSy,  includeDispersion==2?6:3, includeDispersion==2?6:3); 
    m_alloc(&Kx,  n_configs, n_configs);  m_zero(Kx);
    m_alloc(&Ky,  n_configs, n_configs);  m_zero(Ky);

    for (i_config=0; i_config<n_configs; i_config++) {
      Rx->a[i_config][0]  =  sqr(R11[i_config]);
      Rx->a[i_config][1]  =  2*R11[i_config]*R12[i_config];
      Rx->a[i_config][2]  =  sqr(R12[i_config]);
      if (includeDispersion==1) {
        Rx->a[i_config][3] = sqr(R16[i_config]);
        Rx->a[i_config][4] = 2*R11[i_config]*R16[i_config];
        Rx->a[i_config][5] = 2*R12[i_config]*R16[i_config];
      }
      Ry->a[i_config][0]  =  sqr(R33[i_config]);
      Ry->a[i_config][1]  =  2*R33[i_config]*R34[i_config];
      Ry->a[i_config][2]  =  sqr(R34[i_config]);
      if (includeDispersion==2) {
        Ry->a[i_config][3] = sqr(R36[i_config]);
        Ry->a[i_config][4] = 2*R33[i_config]*R36[i_config];
        Ry->a[i_config][5] = 2*R34[i_config]*R36[i_config];
      }
    }

    for (i_dev=0; i_dev<n_dev_limits; i_dev++) {
      emitx_sum = emitx2_sum = 0;
      emity_sum = emity2_sum = 0;
      emitx_max = -(emitx_min = DBL_MAX);
      emity_max = -(emity_min = DBL_MAX);
      S11_sum = S11_sum2 = S12_sum = S12_sum2 = S22_sum = S22_sum2 = 0;
      S33_sum = S33_sum2 = S34_sum = S34_sum2 = S44_sum = S44_sum2 = 0;
      S16_sum = S16_sum2 = S26_sum = S26_sum2 = 0;
      S36_sum = S36_sum2 = S46_sum = S46_sum2 = 0;
      S66_sum = S66_sum2 = 0;
      betax_sum = betax_sum2 = alphax_sum = alphax_sum2 = 0;
      betay_sum = betay_sum2 = alphay_sum = alphay_sum2 = 0;
      md_x = md_y = nx_used_sum = ny_used_sum = 0;
      n_good_fits_x = n_good_fits_y = 0;
      n_xresol = n_yresol = 0; 
      for (i_error=0; i_error<n_error_sets; i_error++) {
        if (verbosity>1) {
          fprintf(stderr, "Error set %ld...\n", i_error);
        }
        for (i_config=0; i_config<n_configs; i_config++) {
          s2x->a[i_config][0] =  sqr( sigmax[i_config] ) - sqr(x_resol);
          s2y->a[i_config][0] =  sqr( sigmay[i_config] ) - sqr(y_resol);
	  if (verbosity>3) 
	    fprintf(stderr, "Data point %ld: sigmax^2=%e, s2x=%e, sigmay^2=%e, s2y=%e\n",
		    i_config, sqr(sigmax[i_config]),  s2x->a[i_config][0],
		    sqr(sigmay[i_config]), s2y->a[i_config][0]);
        }

        if (!ignore_x) {
          if (x_error_level==0 && i_error==0 && !x_error_name) {
            /* must do initial fit to find error level */
            set_up_covariance_matrix(Kx, sigmax, uncertx, n_configs, equal_weights_x_fit);
            x_error_level = estimate_uncertainty(uncertx, Sx, sSx, Rx, s2x, Kx, dev_limit[i_dev], 
                                                 n_configs, x_uncert_min, x_fit_sig2);
	    if (verbosity>2) 
	      fprintf(stderr, "x error level estimate: %le\n", x_error_level);
	  }
          if (!make_tweeked_data_set(s2x, sigmax, x_error_name?uncertx:NULL, x_error_level, error_sigmas, error_type_code, n_configs, 
                                     x_resol, reject_at_limit, x_limit, &n_xresol))  {
            bomb("fatal error: failed to get acceptable error set\n", NULL);
          }
          set_up_covariance_matrix(Kx, sigmax, uncertx, n_configs, equal_weights_x_fit);
        }

        if (!ignore_y) {
          if (y_error_level==0 && i_error==0 && !y_error_name) {
            /* must do initial fit to find error level */
            set_up_covariance_matrix(Ky, sigmay, uncerty, n_configs, equal_weights_y_fit);
            y_error_level = estimate_uncertainty(uncerty, Sy, sSy, Ry, s2y, Ky, dev_limit[i_dev], 
                                                 n_configs,  y_uncert_min, y_fit_sig2);
	    if (verbosity>2) 
	      fprintf(stderr, "y error level estimate: %le\n", y_error_level);
          }
          if (!make_tweeked_data_set(s2y, sigmay, y_error_name?uncerty:NULL, y_error_level, error_sigmas, error_type_code, n_configs, 
                                     y_resol, reject_at_limit, y_limit, &n_yresol))  {
            bomb("fatal error: failed to get acceptable error set\n", NULL);
          }
          set_up_covariance_matrix(Ky, sigmay, uncerty, n_configs, equal_weights_y_fit);
        }

        if (!ignore_x) {
          if ((contrib = solve_normal_form_opt(Sx, sSx, Rx, s2x, Kx, dev_limit[i_dev], &nx_used, x_fit_sig2))<0) {
            fprintf(stderr, "Problem fitting for x\n");
          } else {
            md_x += contrib;
            if (nx_used) {
	      if (verbosity>3) {
		long it;
		for (it=0; it<(includeDispersion==1?6:3); it++)
		  fprintf(stderr, "Sx->a[%ld][0] = %e\n", it, Sx->a[it][0]);
	      }
              if (includeDispersion==1) {
                Sx->a[0][0] -= sqr(Sx->a[4][0])/Sx->a[3][0];
                Sx->a[1][0] -= Sx->a[4][0]*Sx->a[5][0]/Sx->a[3][0];
                Sx->a[2][0] -= sqr(Sx->a[5][0])/Sx->a[3][0];
              }
              S11_sum += Sx->a[0][0]; S11_sum2 += sqr(Sx->a[0][0]);
              S12_sum += Sx->a[1][0]; S12_sum2 += sqr(Sx->a[1][0]);
              S22_sum += Sx->a[2][0]; S22_sum2 += sqr(Sx->a[2][0]);
              if ((emitx = Sx->a[0][0]*Sx->a[2][0]-sqr(Sx->a[1][0]))>0) {
                emitx = sqrt(emitx);
                if (verbosity>2) {
                  fprintf(stderr, "Horizontal emittance: %e\n", emitx);
                }
                if (includeDispersion==1) {
                  S66_sum += Sx->a[3][0]; S66_sum2 += sqr(Sx->a[3][0]);
                  S16_sum += Sx->a[4][0]; S16_sum2 += sqr(Sx->a[4][0]);
                  S26_sum += Sx->a[5][0]; S26_sum2 += sqr(Sx->a[5][0]);
                }
                betax = Sx->a[0][0]/emitx; 
                betax_sum += betax;  betax_sum2  += sqr(betax);
                alphax = -Sx->a[1][0]/emitx;
                alphax_sum += alphax; alphax_sum2 += sqr(alphax);
                emitx_sum  += emitx;
                emitx2_sum += sqr(emitx);
                if (emitx_max<emitx)
                  emitx_max = emitx;
                if (emitx_min>emitx)
                  emitx_min = emitx;
                nx_used_sum += nx_used;
                n_good_fits_x++;
	      } else {
		fprintf(stderr, "Horizontal emittance couldn't be determined\n");
              }
	    }
          }
        }
        

        if (!ignore_y) {       
          if ((contrib = solve_normal_form_opt(Sy, sSy, Ry, s2y, Ky, dev_limit[i_dev], &ny_used, y_fit_sig2))<0) {
            fprintf(stderr, "Problem fitting for y\n");
          } else {
            md_y += contrib;
            if (ny_used) {
              if (includeDispersion==2) {
                Sy->a[0][0] -= sqr(Sy->a[4][0])/Sy->a[3][0];
                Sy->a[1][0] -= Sy->a[4][0]*Sy->a[5][0]/Sy->a[3][0];
                Sy->a[2][0] -= sqr(Sy->a[5][0])/Sy->a[3][0];
              }
              S33_sum += Sy->a[0][0]; S33_sum2 += sqr(Sy->a[0][0]);
              S34_sum += Sy->a[1][0]; S34_sum2 += sqr(Sy->a[1][0]);
              S44_sum += Sy->a[2][0]; S44_sum2 += sqr(Sy->a[2][0]);
              if ((emity = Sy->a[0][0]*Sy->a[2][0]-sqr(Sy->a[1][0]))>0) {
                emity = sqrt(emity);
                if (verbosity>2) {
                  fprintf(stderr, "Vertical emittance: %e\n", emity);
                }
                if (includeDispersion==2) {
                  S66_sum += Sy->a[3][0]; S66_sum2 += sqr(Sy->a[3][0]);
                  S36_sum += Sy->a[4][0]; S36_sum2 += sqr(Sy->a[4][0]);
                  S46_sum += Sy->a[5][0]; S46_sum2 += sqr(Sy->a[5][0]);
                }
                betay = Sy->a[0][0]/emity; 
                betay_sum += betay;  betay_sum2  += sqr(betay);
                alphay = -Sy->a[1][0]/emity;
                alphay_sum += alphay; alphay_sum2 += sqr(alphay);
                emity_sum  += emity;
                emity2_sum += sqr(emity);
                if (emity_max<emity)
                  emity_max = emity;
                if (emity_min>emity)
                  emity_min = emity;
                ny_used_sum += ny_used;
                n_good_fits_y++;
	      } else {
		fprintf(stderr, "Vertical emittance couldn't be determined\n");
              }
            }
          }
        }
      }
      
      xEmitLabel[0] = yEmitLabel[1] = 0;
      if (!ignore_x && n_good_fits_x) {
        if (!error_output) 
          sprintf(xEmitLabel, "$ge$r$bx$n = %.3g um", emitx_sum/n_good_fits_x*1e6);
        else 
          sprintf(xEmitLabel, "$ge$r$bx$n = %.3g $sa$e %.3g um",
                  1e6*emitx_sum/n_good_fits_x, 
                  1e6*sqrt(emitx2_sum/n_good_fits_x-sqr(emitx_sum/n_good_fits_x)));
        if (!SDDS_SetParameters
            (&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
             "ex", emitx_sum/n_good_fits_x, "S11", S11_sum/n_good_fits_x,
             "S12", S12_sum/n_good_fits_x, "S22", S22_sum/n_good_fits_x,
             "betax", betax_sum/n_good_fits_x, "alphax", alphax_sum/n_good_fits_x,
             "xRMSDeviation", md_x/n_good_fits_x, "xGoodFits", n_good_fits_x,
             "xAverageFitPoints", (double)nx_used_sum/n_good_fits_x,
             NULL) ||
            (error_output &&
             !SDDS_SetParameters
             (&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
              "exStDev", sqrt(emitx2_sum/n_good_fits_x-sqr(emitx_sum/n_good_fits_x)),
              "S11StDev", sqrt(S11_sum2/n_good_fits_x-sqr(S11_sum/n_good_fits_x)),
              "S12StDev", sqrt(S12_sum2/n_good_fits_x-sqr(S12_sum/n_good_fits_x)),
              "S22StDev", sqrt(S22_sum2/n_good_fits_x-sqr(S22_sum/n_good_fits_x)),
              "betaxStDev", sqrt(betax_sum2/n_good_fits_x-sqr(betax_sum/n_good_fits_x)),
              "alphaxStDev", sqrt(alphax_sum2/n_good_fits_x-sqr(alphax_sum/n_good_fits_x)),
              NULL)))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (SetSigmaData(&SDDSout, "xSigmaData", s2x, "xSigmaFit", x_fit_sig2, n_configs)) 
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (!ignore_y && n_good_fits_y) {
        if (!error_output) 
          sprintf(yEmitLabel, "$ge$r$by$n = %.3g um", emity_sum/n_good_fits_y*1e6);
        else 
          sprintf(yEmitLabel, "$ge$r$by$n = %.3g $sa$e %.3g um",
                  1e6*emity_sum/n_good_fits_y, 
                  1e6*sqrt(emity2_sum/n_good_fits_y-sqr(emity_sum/n_good_fits_y)));
        if (!SDDS_SetParameters
            (&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
             "ey", emity_sum/n_good_fits_y, "S33", S33_sum/n_good_fits_y,
             "S34", S34_sum/n_good_fits_y, "S44", S44_sum/n_good_fits_y,
             "betay", betay_sum/n_good_fits_y, "alphay", alphay_sum/n_good_fits_y,
             "yRMSDeviation", md_y/n_good_fits_y, "yGoodFits", n_good_fits_y,
             "yAverageFitPoints", (double)ny_used_sum/n_good_fits_y,
             NULL) ||
            (error_output && 
             !SDDS_SetParameters
             (&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
              "eyStDev", sqrt(emity2_sum/n_good_fits_y-sqr(emity_sum/n_good_fits_y)),
              "S33StDev", sqrt(S33_sum2/n_good_fits_y-sqr(S33_sum/n_good_fits_y)),
              "S34StDev", sqrt(S34_sum2/n_good_fits_y-sqr(S34_sum/n_good_fits_y)),
              "S44StDev", sqrt(S44_sum2/n_good_fits_y-sqr(S44_sum/n_good_fits_y)),
              "betayStDev", sqrt(betay_sum2/n_good_fits_y-sqr(betay_sum/n_good_fits_y)),
              "alphayStDev", sqrt(alphay_sum2/n_good_fits_y-sqr(alphay_sum/n_good_fits_y)),
              NULL)))
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
        if (!SetSigmaData(&SDDSout, "ySigmaData", s2y, "ySigmaFit", y_fit_sig2, n_configs)) 
          SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
      if (includeDispersion &&
          (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               "Sdelta", sqrt(S66_sum/(includeDispersion==1?n_good_fits_x:n_good_fits_y)), NULL) ||
           (error_output &&
            !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                                "SdeltaSigma", 
                                sqrt(S66_sum2/(includeDispersion==1?n_good_fits_x:n_good_fits_y)
                                     -sqr(S66_sum/(includeDispersion==1?n_good_fits_x:n_good_fits_y))),
                                NULL))))SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      sprintf(emitLabel, "%s %s", xEmitLabel, yEmitLabel);
      if (!SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                              "EmittanceLabel", emitLabel, NULL) ||
          !SDDS_WritePage(&SDDSout))
        SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }
    
    m_free(&Rx);
    m_free(&Ry);
    m_free(&s2x);
    m_free(&s2y);
    m_free(&Sx);
    m_free(&Sy);
    m_free(&sSx);
    m_free(&sSy);
    m_free(&Kx);
    m_free(&Ky);
  }

  if (!SDDS_Terminate(&SDDSin) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  return(0);
}

void print_fit(
               char *filename, double *variable_data, char *variable_name,
               double *sigma2_fit, double resol, char *sigma_name,
               long n_pts
               )
{
  long i, n_good;
  FILE *fp;

  for (i=n_good=0; i<n_pts; i++)
    if (sigma2_fit[i]>=0)
      n_good++;

  fp = fopen_e(filename, "w", 0);
  fprintf(fp, "%s\n%s\nEMITMEAS fit output\n\n%ld\n", variable_name, sigma_name,
          n_good);

  for (i=0; i<n_pts; i++) 
    if (sigma2_fit[i]>=0)
      fprintf(fp, "%e %e\n", variable_data[i], sqrt(sigma2_fit[i]+sqr(resol)));

  fclose(fp);
}

double solve_normal_form_opt(
                             MATRIX *F,       /* Mx1 matrix of Fit coefficients (returned) */
                             MATRIX *sF,      /* MxM matrix of errors in Fit coefficients (returned) */
                             MATRIX *P,       /* NxM matrix of Parameters of fit (provided) */
                             MATRIX *M,       /* Nx1 column matrix of Measured quantities (provided) */
                             /* M = P.F is the equation being solved */
                             MATRIX *K,      /* NxN inverse covariance matrix for Measured quantities.
                                                K[i][j] = delta(i,j)/uncert[i]^2 */
                             double dev_limit,/* limit on deviation for any point used in final fit */
                             long *n_used,    /* number of points used in final fit */
                             double *s2_fit  /* sigma^2 from fit (returned) */
                             )
{
  long i, j, i_good;
  double rms_error, error;
  double *s2_fit2;
  long *index, n_good, *good;
  MATRIX *Pc, *Mc, *Kc;

  rms_error = solve_normal_form(F, sF, P, M, K, s2_fit);

  if (dev_limit==0) {
    *n_used = M->n;
    return(rms_error);
  }

  good = tmalloc(sizeof(*good)*M->n);
  if (dev_limit<0)
    dev_limit = -dev_limit*rms_error;
  for (i=n_good=0; i<M->n; i++) {
    error = fabs(sqrt(s2_fit[i]) - sqrt(M->a[i][0]));
    if (dev_limit>error) {
      n_good++;
      good[i] = 1;
    }
    else {
      s2_fit[i] = -1;  /* used to mark excluded points for caller */
      good[i] = 0;
    }
  }
  if (n_good==0) {
    *n_used = 0;
    free(good);
    return(0.0);
  }

  m_alloc(&Pc, n_good, P->m);
  m_alloc(&Mc, n_good, 1);
  m_alloc(&Kc, n_good, n_good);
  m_zero(Kc);
  index = tmalloc(sizeof(int)*n_good);
  s2_fit2 = tmalloc(sizeof(*s2_fit2)*n_good);

  for (i=i_good=0; i<M->n; i++) {
    if (good[i]) {
      index[i_good] = i;
      for (j=0; j<P->m; j++)
        Pc->a[i_good][j] = P->a[i][j];
      for (j=0; j<Mc->m; j++)
        Mc->a[i_good][j] = M->a[i][j];
      Kc->a[i_good][i_good] = K->a[i][i];
      i_good++;
    }
  }            

  *n_used = n_good;
  rms_error = solve_normal_form(F, sF, Pc, Mc, Kc, s2_fit2);
  for (i=0; i<n_good; i++)
    s2_fit[index[i]] = s2_fit2[i];

  free(good);
  free(s2_fit2);
  free(index);
  m_free(&Pc);
  m_free(&Mc);
  m_free(&Kc);

  return(rms_error);
}

double solve_normal_form(
                         MATRIX *F,       /* Mx1 matrix of Fit coefficients (returned) */
                         MATRIX *sF,      /* MxM covariance matrix for Fit coefficients (returned) */
                         MATRIX *P,       /* NxM matrix of Parameters of fit (provided) */
                         MATRIX *M,       /* Nx1 column matrix of Measured sigma^2 (provided) */
                         MATRIX *K,       /* NxN inverse covariance matrix for Measured quantities (provided) .
                                             K[i][j] = delta(i,j)/uncert[i]^2 */
                         /* M = P.F is the equation being solved                          */
                         /* the solution is F = inverse(transpose(P) K P) tranpose(P) K M = T.M */
                         /*           and  sF = T.K.transpose(T) */
                         double *s2_fit      /* sigma^2 from fit (returned) */
                         )
{
  long n, m, i;
  MATRIX *Pt, *Pt_K, *Pt_K_P, *Inv_Pt_K_P, *Inv_PtKP_PtK;
  MATRIX *Mp, *Tt, *TC, *C;
  double error, rms_error;

  n = M->n;        /* n is the number of experimental measurements */
  m = P->m;        /* m is 3 or 6, the number of unknowns to determine */
  m_alloc(&Pt , m, n);
  m_alloc(&Pt_K, m, n);
  m_alloc(&Pt_K_P, m, m);
  m_alloc(&Inv_Pt_K_P, m, m);
  m_alloc(&Inv_PtKP_PtK, m, n);
  m_alloc(&Mp, M->n, M->m);
  m_alloc(&Tt, n, m);
  m_alloc(&TC, m, n);
  m_alloc(&C, K->n, K->m);

  /* find the fit */
  if (!m_trans(Pt, P)) {
    fprintf(stderr, "matrix error--call was: m_trans(Pt, P)\n");
    return -1;
  }
  if (!m_mult(Pt_K, Pt, K)) {
    fprintf(stderr, "matrix error--call was: m_mult(Pt_K, Pt, K)\n");
    return -1;
  }
  if (!m_mult(Pt_K_P, Pt_K, P)) {
    fprintf(stderr, "matrix error--call was: m_mult(Pt_K_P, Pt_K, P)\n");
    return -1;
  }
  if (!m_invert(Inv_Pt_K_P, Pt_K_P)) {
    fprintf(stderr, "matrix error--call was: m_invert(Inv_Pt_K_P, Pt_K_P)\n");
    return -1;
  }

  if (!m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K)) {
    fprintf(stderr, "matrix error--call was: m_mult(Inv_PtKP_PtK, Inv_Pt_K_P, Pt_K)\n");
    return -1;
  }
  if (!m_mult(F, Inv_PtKP_PtK, M)) {
    fprintf(stderr, "matrix error--call was: m_mult(F, Inv_PtKP_PtK, M)\n");
    return -1;
  }
  m_zero(sF);
  if (m_invert(C, K)) {
    if (!m_trans(Tt, Inv_PtKP_PtK)) {
      fprintf(stderr, "matrix error--call was: m_trans(Tt, Inv_PtKP_PtK)\n");
      return -1;
    }
    if (!m_mult(TC, Inv_PtKP_PtK, C)) {
      fprintf(stderr, "matrix error--call was: m_mult(TC, Inv_PtKP_PtK, C)\n");
      return -1;
    }
    if (!m_mult(sF, TC, Tt)) {
      fprintf(stderr, "matrix error--call was: m_mult(sF, TC, Tt)\n");
      fprintf(stderr, "sF: %d x %d\n", sF->n, sF->m);
      fprintf(stderr, "TC: %d x %d\n", TC->n, TC->m);
      fprintf(stderr, "Tt: %d x %d\n", Tt->n, Tt->m);
      return -1;
    }
  }

  /* evaluate the fit */
  if (!m_mult(Mp, P, F)) {
    fprintf(stderr, "matrix error--call was: m_mult(Mp, P, F)\n");
    return -1;
  }
  for (i=rms_error=0; i<Mp->n; i++) {
    if ((s2_fit[i] = Mp->a[i][0])<0) {
      fprintf(stderr, "bad fit--negative sigma^2!: fit is %e, data is %e\n",
	      s2_fit[i], M->a[i][0]);
      s2_fit[i] = Mp->a[i][0] = 0;
    }
    error = sqrt(Mp->a[i][0]) - sqrt(M->a[i][0]);
    rms_error += sqr(error);
  }
  if (Mp->n>3)
    rms_error = sqrt(rms_error/(Mp->n-3));

  m_free(&Pt);
  m_free(&Pt_K);
  m_free(&Pt_K_P);
  m_free(&Inv_Pt_K_P);
  m_free(&Inv_PtKP_PtK);
  m_free(&Mp);
  m_free(&Tt);
  m_free(&TC);
  m_free(&C);

  return(rms_error);
}

long make_tweeked_data_set(MATRIX *s2, double *sigma, double *uncert, double error_level, double error_sigmas, long error_type_code, 
                          long n_configs, double resol, long reject_at_limit, double limit, long *n_at_resol)
{
  long i_config, n_tries;
  double tweek;

  for (i_config=n_tries=0; i_config<n_configs && n_tries<MAX_N_TRIES; i_config++) {
    if (uncert)
      error_level = uncert[i_config];
    tweek = (error_type_code==UNIFORM_ERRORS?
             2*error_level*(random_1(1)-.5):
             gauss_rn_lim(0.0, error_level, error_sigmas, random_1));
    s2->a[i_config][0] = sqr( sigma[i_config] + tweek ) - sqr(resol);
    if (s2->a[i_config][0]<sqr(limit)) {
      if (reject_at_limit) {
        n_tries++;
        i_config--;
        continue;
      }
      else {
        *n_at_resol += 1;
        s2->a[i_config][0] = sqr(limit);
      }
    }
  }
  if (n_tries==MAX_N_TRIES)
    return(0);
  return(1);
}

void set_up_covariance_matrix(MATRIX *K, double *sigma, double *uncert, long n_configs, long equal_weights)
/* actually, K is the inverse of the covariance matrix */
{
  long i_config;

  for (i_config=0; i_config<n_configs; i_config++) {
    if (sigma[i_config]==0 || uncert[i_config]==0 || equal_weights)
      K->a[i_config][i_config] = 1;
    else
      K->a[i_config][i_config] = 1./sqr(2*sigma[i_config]*uncert[i_config]);
  }
}

double estimate_uncertainty(double *uncert, MATRIX *S, MATRIX *sS, MATRIX *R, MATRIX *s2, 
                            MATRIX *K, double dev_limit, long n_configs, double uncert_min, double *fit_sig2_return)
{
  double md, *fit_sig2;
  long n_used, i_config;

  if (fit_sig2_return)
    fit_sig2 = fit_sig2_return;
  else
    fit_sig2 = tmalloc(sizeof(*fit_sig2)*n_configs);

  /* find initial fit with supplied covariance matrix */
  md = solve_normal_form_opt(S, sS, R, s2, K, 0.0, &n_used, fit_sig2);
  if (!md || !n_used)
    bomb("unable to find initial fit (1)", NULL);

  /* calculate new covariance matrix */
  for (i_config=0; i_config<n_configs; i_config++)
    K->a[i_config][i_config] = 1./sqr(2*md)/s2->a[i_config][0];

  /* do second fit, excluding points that lie to far out */
  md = solve_normal_form_opt(S, sS, R, s2, K, dev_limit, &n_used, fit_sig2);
  if (!md || !n_used)
    bomb("unable to find initial fit (2)", NULL);

  /* calculate new covariance matrix */
  if (uncert_min && md<uncert_min)
    md = uncert_min;
  for (i_config=0; i_config<n_configs; i_config++) {
    K->a[i_config][i_config] = 1./sqr(2*md)/s2->a[i_config][0];
    uncert[i_config] = md;
  }

  if (!fit_sig2_return)
    free(fit_sig2);

  return(md);
}

long SetSigmaData(SDDS_DATASET *SDDSout, char *dataName, MATRIX *s2, 
                  char *fitName, double *fitSigSqr, long configs)
{
  static double *buffer = NULL;
  long i;

  if (!(buffer=SDDS_Realloc(buffer, configs*sizeof(*buffer))))
    return 0;

  for (i=0; i<configs; i++) {
    if (s2->a[i][0]<0)
      buffer[i] = 0;
    else
      buffer[i] = sqrt(s2->a[i][0]);
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, buffer, configs, dataName))
    return 0;

  for (i=0; i<configs; i++) {
    if (fitSigSqr[i]<0) 
      buffer[i] = 0;
    else
      buffer[i] = sqrt(fitSigSqr[i]);
  }
  if (!SDDS_SetColumn(SDDSout, SDDS_SET_BY_NAME, buffer, configs, fitName))
    return 0;
  
  return 1;
}
