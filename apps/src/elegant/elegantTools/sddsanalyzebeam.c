/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsanalyzebeam.c
 * purpose: take particle coordinate input and analyze to give
 *          moments and twiss parameters
 *
 * Michael Borland, 2000
 *
 $Log: not supported by cvs2svn $
 Revision 1.5  2009/07/21 14:29:49  borland
 Added symbols to the columns.

 Revision 1.4  2008/08/21 16:58:16  borland
 Incremented version number and date.

 Revision 1.3  2008/08/21 14:54:06  borland
 Added option to use canonical coordinates.

 Revision 1.2  2007/07/17 21:02:26  soliday
 Changed strcpy_s to strcpy_ss because of a conflict with VC2005.

 Revision 1.1  2007/03/30 16:50:28  soliday
 Moved from directory above.

 Revision 1.16  2007/02/08 16:51:46  ywang25
 Parallelized several new elements.

 Revision 1.15  2006/08/23 13:20:33  borland
 Now outputs the longitudinal emittance even with only doing corrected quantities.

 Revision 1.14  2005/11/10 15:38:49  soliday
 Added changes to get it to compile properly with 64 bit compilers.

 Revision 1.13  2005/02/24 19:27:53  borland
 Added -generate option to sddsanalyzebeam.  This option can generate a new
 phase-space distribution with the same 6D rms properties as the input beam.
 Uses meschach library.

 Revision 1.12  2004/09/08 20:04:36  borland
 Fixed problem with units for sij output.

 Revision 1.11  2004/08/10 19:02:24  borland
 Fixed uninitialized variable problem that resulted in wrong sigma values
 in output file for single-particle beam.

 Revision 1.10  2003/11/24 17:20:07  borland
 Added computation of longitudinal emittance.

 Revision 1.9  2003/02/26 18:00:12  borland
 Copies input parameters through to output columns.
 Think I fixed a memory leak.

 Revision 1.8  2003/01/09 15:58:51  borland
 Added -correctedOnly option.

 Revision 1.7  2002/08/14 20:23:48  soliday
 Added Open License

 Revision 1.6  2001/10/15 20:37:06  soliday
 Cleaned up for Linux.

 Revision 1.5  2001/10/15 15:42:33  soliday
 Cleaned up for WIN32.

 Revision 1.4  2001/08/06 13:56:40  borland
 Fixed memory leak.

 Revision 1.3  2001/07/02 18:09:46  borland
 Usage and error messages now mention the requirement for 't' column.

 Revision 1.2  2000/12/06 02:14:47  borland
 Now provides "corrected" and "uncorrected" values for beta and alpha.
 The corrected values are the best ones if one really knows the dispersion.

 Revision 1.1  2000/08/14 21:44:57  borland
 First version.

 */
#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_NOWARNINGS 1
#define SET_CORRECTED_ONLY 2
#define SET_GENERATE 3
#define SET_CANONICAL 4
#define N_OPTIONS 5

char *option[N_OPTIONS] = {
  "pipe", "nowarnings", "correctedonly", "generate", "canonical"
} ;

char *USAGE="sddsanalyzebeam [-pipe=[input][,output]] [<SDDSinputfile>] [<SDDSoutputfile>]\n\
  [-nowarnings] [-correctedOnly] [-canonical] [-generate=<outputFile>,<particles>,<cutoff>]\n\
Computes Twiss parameters and other properties of a particle beam.\n\
The input file must have columns x, xp, y, yp, t, and p; for example, an elegant\n\
beam output file is acceptable.\nIf -correctedOnly is given, then only the\n\
\"corrected values\" (with dispersion terms correctly subtracted) are given.\n\
Use this if you plan to use the output as input to the twiss_output command in\n\
elegant.\n\
If -canonical is given, all computations are performed for canonical variables.\n\
If -generate is given, a new gaussian distribution with the specified number of\n\
particles is generated based on the analysis.\n\n\
Program by Michael Borland.  (This is version 4, August 2008.)\n";


long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile, long correctedOnly, long canonical, SDDS_DATASET *SDDSin);
long check_sdds_beam_column(SDDS_TABLE *SDDS_table, char *name, char *units);
long SetUpGenerateFile(SDDS_DATASET *SDDSgen, char *generateFile);
void GenerateAndDumpParticles(SDDS_DATASET *SDDSgen, double C[6], double S[6][6], double pAve,
			      double cutoff, long np);
void TransformToCanonicalMomenta(double **data, long np, double pAve);

char *CenName[6];
char *CorName[6][6];
char *psName[6] = {"x", "xp", "y", "yp", "t", "delta" };
char *qsName[6] = {"x", "qx", "y", "qy", "t", "delta" };
/* columns in output file derived from parameters in input file */
int32_t pColumns;
char **pColumn;
  
int main(int argc, char **argv)
{
  SDDS_DATASET SDDSin, SDDSout, SDDSgen;
  char *inputfile, *outputfile, *generateFile;
  long iPart, particles, i_arg, readCode, noWarnings, row, nToGenerate;
  long i, j;
  SCANNED_ARG *s_arg;
  unsigned long pipeFlags;
  double *x, *xp, *y, *yp, *p=NULL, *t;
  double pAve=0.0, sum;
  double S[6][6], C[6], beta[3], alpha[3], eta[4], emit[3], emitcor[3], beamsize[6];
  double betacor[3], alphacor[3];
  double *data[6], Sbeta[6][6];
  long tmpFileUsed, correctedOnly, canonical=0;
  double cutoff = 3;

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);

  inputfile = outputfile = generateFile = NULL;
  pipeFlags = noWarnings = correctedOnly = 0;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
      case SET_PIPE:
        if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
          SDDS_Bomb("invalid -pipe syntax");
        break;
      case SET_NOWARNINGS:
        noWarnings = 1;
        break;
      case SET_CORRECTED_ONLY:
        correctedOnly = 1;
        break;
      case SET_GENERATE:
	if (s_arg[i_arg].n_items!=4 ||
	    !strlen(generateFile=s_arg[i_arg].list[1]) ||
	    sscanf(s_arg[i_arg].list[2], "%ld", &nToGenerate)!=1 ||
	    sscanf(s_arg[i_arg].list[3], "%lf", &cutoff)!=1 ||
	    nToGenerate<1)
	  SDDS_Bomb("invalid -generate syntax");
	break;
      case SET_CANONICAL:
        canonical = 1;
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

  processFilenames("sddsanalyzebeam", &inputfile, &outputfile, pipeFlags, noWarnings, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  if (!SDDS_InitializeInput(&SDDSin, inputfile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

  /* check the input file for valid data */
  if (!check_sdds_beam_column(&SDDSin, "x", "m") || !check_sdds_beam_column(&SDDSin, "y", "m") ||
      !check_sdds_beam_column(&SDDSin, "xp", NULL) || !check_sdds_beam_column(&SDDSin, "yp", NULL) ||
      !check_sdds_beam_column(&SDDSin, "t", "s") ||
      (!check_sdds_beam_column(&SDDSin, "p", "m$be$nc") && !check_sdds_beam_column(&SDDSin, "p", NULL))) {
    fprintf(stderr, 
            "sddsanalyzebeam: one or more of (x, xp, y, yp, t, p) have the wrong units or are not present in %s", 
            inputfile);
    exit(1);
  }

  if (!SetUpOutputFile(&SDDSout, outputfile, correctedOnly, canonical, &SDDSin) ||
      !SDDS_SetParameters(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                          "canonical", canonical, NULL))
    SDDS_Bomb("problem setting up output file");
  
  if (generateFile) {
    if (!SetUpGenerateFile(&SDDSgen, generateFile))
      SDDS_Bomb("problem setting up -generate output file");
    random_1(-987654321);
  }

  for (i=0; i<6; i++)
    data[i] = NULL;

  row = 0;
  while ((readCode=SDDS_ReadPage(&SDDSin))>0) {
    if (readCode!=1 && !SDDS_LengthenTable(&SDDSout, 1))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (i=0; i<2; i++)
      beta[i] = alpha[i] = emit[i] = emitcor[i] = eta[i] = eta[i+2] = betacor[i] = alphacor[i] = 0;
    for (i=0; i<6; i++) {
      for (j=C[i]=beamsize[i]=0; j<=i; j++)
        S[i][j] = 0;
      if (data[i]) {
        free(data[i]);
        data[i] = NULL;
      }
    }
    x = xp = y = yp = t = p = NULL;
    if ((particles=SDDS_RowCount(&SDDSin))>2) {
      if (!(data[0] = x = SDDS_GetColumnInDoubles(&SDDSin, "x")) ||
          !(data[1] = xp = SDDS_GetColumnInDoubles(&SDDSin, "xp")) ||
          !(data[2] = y = SDDS_GetColumnInDoubles(&SDDSin, "y")) ||
          !(data[3] = yp = SDDS_GetColumnInDoubles(&SDDSin, "yp")) ||
          !(data[4] = t =  SDDS_GetColumnInDoubles(&SDDSin, "t")) ||
          !(data[5] = p = SDDS_GetColumnInDoubles(&SDDSin, "p")))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);

      /* convert momentum to (p-po)/po */
      pAve = arithmeticAverage(p, particles);
      for (i=0; i<particles; i++)
        p[i] = p[i]/pAve-1;

      if (canonical) 
        TransformToCanonicalMomenta(data, particles, pAve);

      /* compute and subtract off average values */
      for (i=0; i<6; i++) {
        C[i] = arithmeticAverage(data[i], particles);
        for (j=0; j<particles; j++)
          data[i][j] -= C[i];
      }

      /* compute correlations */
      for (i=0; i<6; i++) {
        for (j=0; j<=i; j++) {
          for (iPart=sum=0; iPart<particles; iPart++)
            sum += data[i][iPart]*data[j][iPart];
          S[j][i] = S[i][j] = sum/particles;
        }
      }
      if (generateFile)
	GenerateAndDumpParticles(&SDDSgen, C, S, pAve, cutoff, nToGenerate);

      for (i=0; i<6; i++)
        beamsize[i] = sqrt(S[i][i]);

      /* compute correlations with energy correlations removed */
      if (S[5][5])
        for (i=0; i<4; i++) 
          eta[i] = S[i][5]/S[5][5];
      else 
        for (i=0; i<4; i++) 
          eta[i] = 0;
      for (i=0; i<4; i++) {
        for (iPart=0; iPart<particles; iPart++)
          data[i][iPart] -= eta[i]*data[5][iPart];
      }
      for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
          for (iPart=sum=0; iPart<particles; iPart++)
            sum += data[i][iPart]*data[j][iPart];
          Sbeta[j][i] = Sbeta[i][j] = sum/particles;
        }
      }
      /* compute beta functions etc */
      for (i=0; i<3; i++) {
        emitcor[i] = emit[i] = beta[i] = alpha[i] = 0;
        if ((emit[i] = S[2*i+0][2*i+0]*S[2*i+1][2*i+1]-sqr(S[2*i+0][2*i+1]))>0) {
          emit[i] = sqrt(emit[i]);
          beta[i] = S[2*i+0][2*i+0]/emit[i];
          alpha[i] = -S[2*i+0][2*i+1]/emit[i];
        } else
          emit[i] = 0;
        if ((emitcor[i] = Sbeta[2*i+0][2*i+0]*Sbeta[2*i+1][2*i+1]-sqr(Sbeta[2*i+0][2*i+1]))>0) {
          emitcor[i] = sqrt(emitcor[i]);
          betacor[i] = Sbeta[2*i+0][2*i+0]/emitcor[i];
          alphacor[i] = -Sbeta[2*i+0][2*i+1]/emitcor[i];
        } else
          emitcor[i] = 0;
      }
    }
    /* set centroids and sigmas */
    for (i=0; i<6; i++) {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                             row, CenName[i], C[i], NULL)) {
        SDDS_SetError("Problem setting centroid value");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      for (j=0; j<=i; j++) {
        if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                               row, CorName[i][j], S[i][j], NULL)) {
          SDDS_SetError("Problem setting sigma value");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row,
                           "Sx", beamsize[0], "Sxp", beamsize[1],
                           "Sy", beamsize[2], "Syp", beamsize[3],
                           "St", beamsize[4], "Sdelta", beamsize[5], NULL)) {
      SDDS_SetError("Problem setting sigma value");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    /* set twiss parameter, etc */
    if (!correctedOnly) {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row, 
                             "betax", beta[0], "alphax", alpha[0], "etax", eta[0],
                             "etaxp", eta[1], "ex", emit[0], 
                             "enx", emit[0]*pAve, 
                             "betay", beta[1], "alphay", alpha[1], "etay", eta[2],
                             "etayp", eta[3], "ey", emit[1], 
                             "eny", emit[1]*pAve, 
                             "pAverage", pAve, 
                             "ecx", emitcor[0], "ecnx", emitcor[0]*pAve,
                             "betacx", betacor[0], "alphacx", alphacor[0],
                             "ecy", emitcor[1], "ecny", emitcor[1]*pAve,
                             "betacy", betacor[1], "alphacy", alphacor[1],
                             "el", emit[2], "pAverage", pAve, "ElementName", "None", NULL)) {
        SDDS_SetError("Problem setting Twiss values");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    else {
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, row, 
                             "etax", eta[0], "etaxp", eta[1], 
                             "etay", eta[2], "etayp", eta[3], 
                             "ex", emitcor[0], "enx", emitcor[0]*pAve,
                             "betax", betacor[0], "alphax", alphacor[0],
                             "ey", emitcor[1], "eny", emitcor[1]*pAve,
                             "betay", betacor[1], "alphay", alphacor[1], 
                             "el", emit[2], "pAverage", pAve, NULL)) {
        SDDS_SetError("Problem setting Twiss values");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }

    /* copy parameters from the input file to columns in the output file */
    for (i=0; i<pColumns; i++) {
      long buffer[16];
      if (!pColumn[i])
        continue;
      if (!SDDS_GetParameter(&SDDSin, pColumn[i], buffer)) {
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }
      if (!SDDS_SetRowValues(&SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_REFERENCE, row,
                             pColumn[i], buffer, NULL)) {
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }
    }
    row++;
  }
  if (readCode==0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSin) ||
      !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return 0;
}

long SetUpOutputFile(SDDS_DATASET *SDDSout, char *outputfile, long correctedOnly, long canonical, SDDS_DATASET *SDDSin)
{
  char units[128], name[128];
  char *ppUnits[6] = {"m", 0, "m", 0, "s", 0} ;
  long i, j;
  
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, NULL, NULL, outputfile) ||
      !SDDS_DefineSimpleParameter(SDDSout, "canonical", NULL, SDDS_LONG) ||
      !SDDS_DefineSimpleColumn(SDDSout, "ElementName", NULL, SDDS_STRING) ||
      SDDS_DefineColumn(SDDSout, "ex", "$ge$r$bx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "enx", "$ge$r$bnx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "betax", "$gb$r$bx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "alphax", "$ga$r$bx$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "etax", "$gc$r$bx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "etaxp", "$gc$r$bx$n$a'$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "ey", "$ge$r$by$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "eny", "$ge$r$bny$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "betay", "$gb$r$by$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "alphay", "$ga$r$by$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "etay", "$gc$r$by$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "etayp", "$gc$r$by$n$a'$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Sx", "$gs$r$bx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Sxp", "$gs$r$bx'$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Sy", "$gs$r$by$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Syp", "$gs$r$by'$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "St", "$gs$r$bt$n", "s", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "Sdelta", "$gs$bd$n$r", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "pAverage", "<p>", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(SDDSout, "el", "$ge$r$bl$n", "s", NULL, NULL, SDDS_DOUBLE, 0)<0) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!correctedOnly &&
      (SDDS_DefineColumn(SDDSout, "ecx", "$ge$r$bcx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "ecnx", "$ge$r$bcnx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "betacx", "$gb$r$bcx$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "alphacx", "$ga$r$bcx$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "ecy", "$ge$r$bcy$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "ecny", "$ge$r$bcny$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "betacy", "$gb$r$bcy$n", "m", NULL, NULL, SDDS_DOUBLE, 0)<0 ||
       SDDS_DefineColumn(SDDSout, "alphacy", "$ga$r$bcy$n", NULL, NULL, NULL, SDDS_DOUBLE, 0)<0)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (i=0; i<6; i++) {
    sprintf(name, "C%s", canonical==0?psName[i]:qsName[i]);
    if (!SDDS_CopyString(&CenName[i], name))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    if (!SDDS_DefineSimpleColumn(SDDSout, name, ppUnits[i], SDDS_DOUBLE))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    for (j=0; j<=i; j++) {
      sprintf(name, "s%ld%ld", j+1, i+1);
      units[0] = 0;
      if (ppUnits[i]) {
        if (ppUnits[j]) {
          if (strcmp(ppUnits[i], ppUnits[j])==0 && strlen(ppUnits[i]))
            sprintf(units, "%s$a2$n", ppUnits[i]);
          else if (strlen(ppUnits[i]))
            sprintf(units, "%s %s", ppUnits[i], ppUnits[j]);
        }
        else {
          strcpy_ss(units, ppUnits[i]);
        }
      } else if (ppUnits[j])
        strcpy_ss(units, ppUnits[j]);
      if (!SDDS_DefineSimpleColumn(SDDSout, name, units, SDDS_DOUBLE))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      if (!SDDS_CopyString(&CorName[i][j], name))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!(pColumn = SDDS_GetParameterNames(SDDSin, &pColumns))) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  for (i=0; i<pColumns; i++) {
    if (SDDS_GetColumnIndex(SDDSout, pColumn[i])>=0) {
      free(pColumn[i]);
      pColumn[i] = NULL;
    }
    if (!SDDS_DefineColumnLikeParameter(SDDSout, SDDSin, pColumn[i], NULL)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!SDDS_SaveLayout(SDDSout) || !SDDS_WriteLayout(SDDSout) ||
      !SDDS_StartPage(SDDSout, 1)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  return(1);
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

#define PHASE_SPACE_COLUMNS 7
static SDDS_DEFINITION column_definition[PHASE_SPACE_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"t", "&column name=t, units=s, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    } ;

#define PHASE_SPACE_PARAMETERS 2
static SDDS_DEFINITION parameter_definition[PHASE_SPACE_PARAMETERS] = {
  {"pCentral", "&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double &end"},
  {"Particles", "&parameter name=Particles, type=long &end"},
  };

long SetUpGenerateFile(SDDS_DATASET *SDDSgen, char *generateFile)
{
  long i;

  if (!SDDS_InitializeOutput(SDDSgen, SDDS_BINARY, 0, NULL, NULL, generateFile)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  
  for (i=0; i<PHASE_SPACE_PARAMETERS; i++) 
    if (!SDDS_ProcessParameterString(SDDSgen, parameter_definition[i].text, 0) ||
	SDDS_GetParameterIndex(SDDSgen, parameter_definition[i].name)<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  
  /* define SDDS columns */
  for (i=0; i<PHASE_SPACE_COLUMNS; i++) 
    if (!SDDS_ProcessColumnString(SDDSgen, column_definition[i].text, 0) ||
	SDDS_GetColumnIndex(SDDSgen, column_definition[i].name)<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }

  if (!SDDS_WriteLayout(SDDSgen)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  return 1;
}

#include       "matrix.h"
#include       "matrix2.h"   
#include       <setjmp.h>
#include        "err.h"


void GenerateAndDumpParticles(SDDS_DATASET *SDDSgen, double C[6], double S[6][6], double pAve,
			      double cutoff, long np)
{
  long i, j, ip;
  MAT *Sigma, *U, *V, *Ut;
  VEC *SValue, *Coord1, *Coord0;
  double sigma[6];

  Sigma = m_get(6, 6);
  for (i=0; i<6; i++)
    for (j=0; j<6; j++)
      Sigma->me[i][j] = S[i][j]/sqrt(S[i][i]*S[j][j]);

  U = m_get(6, 6);
  V = m_get(6, 6);
  SValue = v_get(6);
  svd(Sigma, U, V, SValue);
  Ut = m_transp(U, MNULL);

#ifdef DEBUG
  for (i=0; i<6; i++) {
    printf("U%ld: ", i);
    for (j=0; j<6; j++)
      printf("%10.3e ", U->me[i][j]);
    printf("\n");
  }
  for (i=0; i<6; i++) {
    printf("V%ld: ", i);
    for (j=0; j<6; j++)
      printf("%10.3e ", V->me[i][j]);
    printf("\n");
  }
  for (i=0; i<6; i++)
    printf("S%ld: %10.3e\n", i, SValue->ve[i]);
#endif
  
  for (i=0; i<6; i++)
    sigma[i] = sqrt(SValue->ve[i]);

  if (!SDDS_StartPage(SDDSgen, np)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }

  Coord0 = v_get(6);
  Coord1 = v_get(6);
  for (ip=0; ip<np; ip++) {
    for (i=0; i<6; i++)
      Coord0->ve[i] = gauss_rn_lim(0.0, 1.0, cutoff, random_1)*sigma[i];
    vm_mlt(V, Coord0, Coord1);
    for (i=0; i<5; i++)
      if (!SDDS_SetRowValues(SDDSgen, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip,
			     i, Coord1->ve[i]*sqrt(S[i][i])+C[i], -1)) {
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }

    if (!SDDS_SetRowValues(SDDSgen, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ip,
			   i, pAve*(1+Coord1->ve[i]*sqrt(S[i][i])+C[i]), -1)) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
  if (!SDDS_WritePage(SDDSgen)) {
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exit(1);
  }
  
}

void TransformToCanonicalMomenta(double **data, long np, double pAve)
{
  double qz, xp, yp, delta;
  long ip;

  /* Compute canonical momenta, save slopes */
  for (ip=0; ip<np; ip++) {
    xp    = data[1][ip];
    yp    = data[3][ip];
    delta = data[5][ip];
    qz = (1+delta)/sqrt(1+sqr(xp)+sqr(yp));
    data[1][ip] *= qz;  /* qx = px/p0 = xp*(1+delta)/sqrt(1+xp^2+yp^2) = xp*qz*/
    data[3][ip] *= qz;  /* qy = yp*qz */
    data[5][ip]  = qz;
  }
}

