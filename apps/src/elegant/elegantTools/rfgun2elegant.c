/* Copyright 1995 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: rfgun2elegant
 * purpose: convert RFGUN output to a form that can be used
 *          with elegant.  The RFGUN output must be in SDDS protocol.
 * issues:  elegant wants a beam of equally-charged particles, which RFGUN 
 *          does not supply.  Also, RFGUN is often run with y=yp=0, to save 
 *          cpu when the cylindrical symmetry is unbroken.
 * 
 * The program assumes the charge of any macro-particle can be adequately
 * expressed as a multiple of the charge in the least-charged macro-particle.
 * Hence, the number of new macro-particles is 
 *             N_new = (int)(multiplier*Charge_old/Charge_min)
 * where multiplier is specified with the -multiplier switch, and defaults to
 * 1.
 *
 * If -symmetrize is given, the new particles are redistributed according to
 *             x_new(i) =  x_old*cos(i*dtheta+theta0)
 *            xp_new(i) = xp_old*cos(i*dtheta+theta0)
 *             y_new(i) =  y_old*sin(i*dtheta+theta0)
 *            yp_new(i) = yp_old*sin(i*dtheta+theta0)
 * where theta is a random number between 0 and 2*pi, i runs from 0 to N_new-1,
 * and dtheta = 2*pi/N_new.
 * 
 * 
 * Michael Borland, 1991. 
 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"
#include "match_string.h"

#define SET_MULTIPLIER 0
#define SET_SYMMETRIZE 1
#define SET_SEED 2
#define SET_SAMPLE_FRACTION 3
#define SET_VERBOSE 4
#define SET_PMIN 5
#define SET_DRIFT 6
#define N_OPTIONS 7

char *option[N_OPTIONS] = {
    "multiplier", "symmetrize", "seed", "sample_fraction", 
    "verbose", "pmin", "drift,"
    } ;

char *USAGE="rfgun2elegant inputfile outputfile [-multiplier=number] [-symmetrize]\n\
[-seed=random_number_seed] [-sample_fraction=value] [-pmin=value]  [-verbose] [-drift=<meters>]\n\n\
Program by Michael Borland.  (This is version 2, January 2010)";

main(int argc, char **argv)
{
  SDDS_TABLE inTable, outTable;
  char *input, *output;
  long i, i_arg;
  SCANNED_ARG *scanned;
  long multiplier, symmetrize, seed;
  double sample_fraction;
  double charge_min, theta0, theta, dtheta;
  double cos_theta, sin_theta, p;
  char s[200];
  long n_new, n_mp, i_new, j, verbose, inPoints, outPoints;
  double *pxIn, *tIn, *xIn, *yIn, *pyIn, *pzIn, *charge;
  double pmin = 0, xp, yp;
  double drift = 0;

  argc = scanargs(&scanned, argc, argv); 
  if (argc<3 || argc>(3+N_OPTIONS)) 
    bomb(NULL, USAGE);

  input = output = NULL;
  multiplier = 1;
  symmetrize = 0;
  seed = 123456789;
  sample_fraction = 1;
  verbose = 0;

  for (i_arg=1; i_arg<argc; i_arg++) {
    if (scanned[i_arg].arg_type==OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_MULTIPLIER:
        if (scanned[i_arg].n_items!=2 ||
            1!=sscanf(scanned[i_arg].list[1], "%ld", &multiplier) ||
            multiplier<=0)
          bomb("invalid -multiplier syntax", NULL);
        break;
      case SET_SYMMETRIZE:
        symmetrize = 1;
        break;
      case SET_SEED:
        if (scanned[i_arg].n_items!=2 ||
            1!=sscanf(scanned[i_arg].list[1], "%ld", &seed) ||
            seed<=0)
          bomb("invalid -seed syntax", NULL);
        break;
      case SET_SAMPLE_FRACTION:
        if (scanned[i_arg].n_items!=2 ||
            1!=sscanf(scanned[i_arg].list[1], "%lf", &sample_fraction) ||
            sample_fraction<=0)
          bomb("invalid -sample_fraction syntax", NULL);
        break;
      case SET_VERBOSE:
        verbose = 1;
        break;
      case SET_PMIN:
        if (scanned[i_arg].n_items!=2 ||
            1!=sscanf(scanned[i_arg].list[1], "%lf", &pmin) ||
            pmin<0)
          bomb("invalid -pmin syntax", NULL);
        break;
      case SET_DRIFT:
        if (scanned[i_arg].n_items!=2 ||
            1!=sscanf(scanned[i_arg].list[1], "%lf", &drift))
          bomb("invalid -drift syntax", NULL);
        break;
      default:
        bomb("unknown option in command line", NULL);
        break;
      }
    }
    else {
      if (!input)
        input = scanned[i_arg].list[0];
      else if (!output)
        output = scanned[i_arg].list[0];
      else
        bomb("too many filenames listed", NULL);
    }
  }
  if (!input || !output) 
    bomb("required arguments not supplied", NULL);

  if (!SDDS_InitializeInput(&inTable, input))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_CheckColumn(&inTable, "x", "m", 0, stderr)!=SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid x data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "y", "m", 0, stderr)!=SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid y data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "px", "m$be$nc", 0, stderr)!=SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid xp data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "py", "m$be$nc", 0, stderr)!=SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid yp data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "t", "s", 0, stderr)!=SDDS_CHECK_OKAY) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid t data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "pz", "m$be$nc", 0, stderr)!=SDDS_CHECK_OKAY ) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid p data in rfgun input file", NULL);
  }
  if (SDDS_CheckColumn(&inTable, "charge", "C", 0, stderr)!=SDDS_CHECK_OKAY ) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    bomb("missing or invalid charge data in rfgun input file", NULL);
  }


  /* The transfer column statements must appear in the order x, xp, y, yp, t, p if the
     rest of the program is to work!
     */
  if (!SDDS_InitializeOutput(&outTable, SDDS_BINARY, 0, NULL, NULL, output) ||
      !SDDS_DefineSimpleColumn(&outTable, "x", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&outTable, "xp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&outTable, "y", "m", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&outTable, "yp", NULL, SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&outTable, "t", "s", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(&outTable, "p", "m$be$nc", SDDS_DOUBLE) ||
      !SDDS_SaveLayout(&outTable) || !SDDS_WriteLayout(&outTable))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  charge_min = 0;
  n_mp = 0;
  if (seed<0)
    seed = -seed;
  seed = 2*((long)(seed/2))+1;
  random_1(-seed);
  while (SDDS_ReadTable(&inTable)>0) {
    if (!(inPoints = SDDS_CountRowsOfInterest(&inTable)))
      continue;
    if (!(charge = SDDS_GetColumnInDoubles(&inTable, "charge")) ||
        !(xIn = SDDS_GetColumnInDoubles(&inTable, "x")) ||
        !(yIn = SDDS_GetColumnInDoubles(&inTable, "y")) ||
        !(tIn =  SDDS_GetColumnInDoubles(&inTable, "t")) ||
        !(pxIn = SDDS_GetColumnInDoubles(&inTable, "px")) ||
        !(pyIn = SDDS_GetColumnInDoubles(&inTable, "py")) ||
        !(pzIn =  SDDS_GetColumnInDoubles(&inTable, "pz")) )
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (charge_min==0) {
      charge_min = DBL_MAX;
      for (i=0; i<inPoints; i++) {
        if (charge_min>charge[i])
          charge_min = charge[i];
      }
      if (charge_min==0)
        bomb("charge==0 for particle", NULL);
    }
    /* calculate the number of particles that will be generated */
    outPoints = 0;
    for (i=0; i<inPoints; i++)
      outPoints += (1.0*multiplier*charge[i])/charge_min+0.5;
    if (verbose) {
      sprintf(s, "creating %ld macro-particles for a total of %ld\n", 
              (long)(sample_fraction*outPoints), n_mp+(long)(sample_fraction*outPoints));
      fputs(s, stdout);
    }
    if (!outPoints)
      continue;
    if (!SDDS_StartTable(&outTable, ((long)(outPoints+100))))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    i_new = 0;
    for (i=0; i<inPoints; i++) {
      if (pzIn[i]<=0)
        continue;
      p = sqrt(sqr(pxIn[i])+sqr(pyIn[i])+sqr(pzIn[i]));
      if (p<pmin)
	continue;
      n_new = (1.0*multiplier*charge[i])/charge_min+0.5;
      xp = pxIn[i]/pzIn[i];
      yp = pyIn[i]/pzIn[i];
      if (symmetrize) {
        dtheta = PIx2/n_new;
        theta0 = random_1(0)*PIx2;
        for (j=0; j<n_new; j++) {
          if (sample_fraction==1 || random_1(0)<sample_fraction) {
            cos_theta = cos(theta = theta0 + j*dtheta);
            sin_theta = sin(theta);
            if (!SDDS_SetRowValues(&outTable, SDDS_BY_INDEX|SDDS_PASS_BY_VALUE, i_new,
                                   0, (xIn[i]+drift*xp)*cos_theta, 1, xp*cos_theta,
                                   2, (xIn[i]+drift*xp)*sin_theta, 3, xp*sin_theta,
                                   4, tIn[i],  5, p, -1))
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            i_new++;
          }
        }                
      }
      else {
        for (j=0; j<n_new; j++) {
          if (sample_fraction==1 || random_1(0)<sample_fraction) {
            if (!SDDS_SetRowValues(&outTable, SDDS_BY_INDEX|SDDS_PASS_BY_VALUE, i_new,
                                   0, xIn[i]+xp*drift,  1, xp, 
                                   2, yIn[i]+yp*drift,  3, yp,
                                   4, tIn[i],  5, p, -1))
              SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            i_new++;
          }
        }                
      }
    }
    if (i_new) {
      n_mp += i_new;
      if (!SDDS_WriteTable(&outTable))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  free(charge);
  free(xIn);
  free(yIn);
  free(tIn);
  free(pxIn);
  free(pyIn);
  free(pzIn);
  if (verbose)
    printf("\n%ld macro-particles created\n", n_mp);
}


