/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: sddsfindresonances.c
   Hairong Shang, Oct 2005
   *
   $Log: not supported by cvs2svn $
   Revision 1.2  2010/07/13 21:37:07  borland
   Added ability to increase the number of multipoles using all=<number> syntax.

   Revision 1.1  2007/03/30 16:50:29  soliday
   Moved from directory above.

   Revision 1.4  2006/10/23 19:49:43  soliday
   Updated to fix an issue with linux-x86_64

   Revision 1.3  2005/11/28 16:01:01  borland
   Added -variables option.

   Revision 1.2  2005/11/18 21:05:19  shang
   updated the comments.

   Revision 1.1  2005/11/18 20:58:51  shang
   first version.

   */


#include "mdb.h"
#include "scan.h"
#include "SDDS.h"

#define SET_PIPE 0
#define SET_MULTIPOLES 1
#define SET_TYPE 2
#define SET_VARIABLES 3
#define N_OPTIONS 4
char *option[N_OPTIONS] = {
  "pipe",  "multipoles", "type", "variables",
} ;

#define USE_DIPOLE 0x0001U
#define USE_QUADRUPOLE 0x0002U
#define USE_SEXTUPOLE 0x0004U
#define USE_OCTUPOLE 0x0008U
#define USE_ALL_MULTIPOLES 0x0010U

#define USE_NORMAL_TYPE 0x0001U
#define USE_SKEW_TYPE 0x0002U

char *USAGE="sddsfindresonances [-pipe=[input][,output]] [<inputFile>] [<outputfile>]\n\
 -multipoles=[all=<integer>]|[dipole,][quadrupole,][sextupole,][octupole,] \n\
 [-type=[normal,][skew]] [-variables=<firstColumn>,<secondColumn>]\n\
-multipoles    for user to choose which multipoles (dipole, quadrupole, \n\
               sextupole, and/or octupole) to use or compute all multipoles by choosing all.\n\
               By default, all multipoles will be computed up to octupole order (4).\n\
               For multipoles up to 2*n-pole, use all=<n>\n\
-type          for user to choose to use normal and/or skew for each \n\
               multipoles. default is normal and skew. \n\
-variables     Give names for the variable columns.  The defaults are x and y.\n\n\
Program by H. Shang, November 2005 (This is Version 2, November 2005.) \n";

static double *sortData[2];
int SDDS_CompareDoubleRows(const void *vrow1, const void *vrow2);
/*comment on algorithms 
  This program (resonance finder) interpolate solutions of the equation 
  M * nux + N * nuy = P and M * nux - N * nuy = P (M, N and P are positive integers)
  for various values of M, N and P, nux and nuy are between 0 and 1.
  Q = 2 * (M +N) is the multipole (2=dipole, 4=quadrupole 6=sextupole, 8=octupole) that drives the resonance.
  That is, the combinations of M, and N should satisfy (M + N = Q/2)
         for dipole:  M + N = 1
         for quadrupole:  M + N =2
         for sextupole:  M + N = 3
         for octupole:   M + N = 4

   For each multipole, there are two types -- normal and/or skew; for each Q, we have
   M      N         multipole type
   Q/2    0         normal
   Q/2-1  1         skew
   Q/2-2  2         normal
   ect.
   for example, for quadrupole, Q=4, there are following possible values of M and N
           M    N    type
           2    0    normal
           1    1    skew
           0    2    normal
   the type of normal and/or skew is choosed by user. User also may chose any one mutipole or
   any combinations of multipoles.
   where P varies from min(M,N) (or 0 if min(M,N)>0) to max(M, N) for each pair of M and N.
*/
/* comments on coding:
   inputFile should have x, y, nux, and nuy columns, which compose the indepData.
   i.e., indepData[0]=x, indepData[1]=y, indepData[2]=nux, indepData[3]=nuy
   two sets of data are interpolated,
   first set data is sorted  by x and then by y
   second set data is sorted by y and then by x
   multx in the code is M, and multy is N, and offset is P
   both multx*nux + multy*nuy = offset and multx*nux - multy*nuy = offset are interpolated.
   the outputData (0 for x, 1 for y, 2 for nux, and 3 for nuy) are the interpolated values where
   two consecutives values of (multx*nux +- multy*nuy - offset) have opsitive signs (i.e, 
   value of "multx*nux +- multy*nuy - offset" crosses zero).
   
   The coding algorithm is
   foreach selected multipole do {
        multx = Q/2 (which is the value of pole[i_pole])
        count=0
        while (multx>=0) {
             multy = Q/2 - multx 
             if (both normal and skew type are requested do following
             else if (count%2==0) do following only when normal type is requested
             else if (count%2!=0) do following only when skew type is requested
             #if multy is 0, no need to change sign, so set signs=1
             if (multy==0) signs=1 
             else signs=2
             for (i=0; i<signs; i++) {
                 if (sign) multy = -multy (change sign)
                 start = min(multx, multy)
                 end = max(multx, multy)
                 if (start>0) start=0
                 for (offset=start; offset<=end; offset++) {
                     for (plane=0; plane<2; plane++) {
                         #plane=0 -- first set of sorted data, sort by x and then y
                         #plane=1 -- second set of sorted data, sort by y and then x
                         #find the resonace of following equation
                         #multx *nux + multy * nuy = offset
                         for (i_data=0; i_data<datapoints-1; i++) {
                             #skip if the two consecutive values of x (first set data)
                               or y (second set data) are different.
                             delta = multx *nux[i_data] + multy * nuy[i_data] - offset
                             delta1 = multx *nux[i_data+1] + multy * nuy[i_data+1] - offset
                             if delta and delta1 have different signs -- interpolate the values of
                                 x, y, nux and nuy (resonance) through delta and delta1 and write to
                              output data array
                         }
                     } 
                 }
             }
             multx --;
             count ++;
        }
    }   
*/

int main(int argc, char **argv)
{
  SCANNED_ARG *s_arg;
  char *inputFile, *outputFile;
  SDDS_DATASET inData, outData;
  long i_arg, tmpFileUsed, rows, i, k, sign, signs, multx, multy, start, end, offset;
  long count, datapoints=0, outputpoints=0, **sortIndex, plane, index, index1, *pole=NULL, i_pole, poles;
  long maximumOrder;
  char *columnName[4] = {NULL, NULL, "nux", "nuy"};  
  char *multipole0[4]={"dipole", "quadrupole", "sextupole", "octupole"};
  char **multipole = NULL;
  double **outputData, *tmpData, **indepData, slope, delta, delta1;  
  unsigned long multipoleFlags, typeFlags, pipeFlags;
  char label[1024];
  long inputPage;

  tmpFileUsed=rows=0;
  inputFile=outputFile=NULL;
  multipoleFlags = typeFlags = pipeFlags = 0;
  outputData = (double**)malloc(sizeof(*outputData)*4);
  indepData = (double**)malloc(sizeof(*indepData)*4);
  /* indepData[i] i=0,1,2,3 is the input data of x, y, nux, and nuy respectively */
  /* outputData[i] i=0,1,2,3 is the interpolated values of x, y, nux, and nuy respectively */
  outputData[0] = outputData[1] = outputData[2] = outputData[3]= NULL;
  indepData[0] = indepData[1] = indepData[2] = indepData[3] = NULL;
  tmpData=NULL;
  /*sortIndex is for sorting x (and then y) or y (and then y) data*/
  /* 0 (sort x and then y), 1 (sort y and then x) */
  sortIndex = malloc(sizeof(*sortIndex)*2);
  sortIndex[0] = sortIndex[1] = NULL;
  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv);
  if (argc<2) 
    bomb(NULL, USAGE);
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
       switch (match_string(s_arg[i_arg].list[0], option, N_OPTIONS, 0)) {
       case SET_PIPE:
         if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags))
           SDDS_Bomb("invalid -pipe syntax");
         break;
       case SET_MULTIPOLES:
         s_arg[i_arg].n_items--;
         maximumOrder = 4;
         if (!scanItemList(&multipoleFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                           "dipole", -1, NULL, 0, USE_DIPOLE,
                           "quadrupole", -1, NULL, 0, USE_QUADRUPOLE,
                           "sextupole", -1, NULL, 0, USE_SEXTUPOLE,
                           "octupole", -1, NULL, 0, USE_OCTUPOLE,
                           "all", SDDS_LONG, &maximumOrder, 1, USE_ALL_MULTIPOLES,
                           NULL))
           SDDS_Bomb("Invalid -multipoles syntax");
         s_arg[i_arg].n_items++;
         break;
       case SET_TYPE:
         s_arg[i_arg].n_items--;
         if (!scanItemList(&typeFlags, s_arg[i_arg].list+1, &s_arg[i_arg].n_items, 0,
                           "normal", -1, NULL, 0, USE_NORMAL_TYPE,
                           "skew", -1, NULL, 0, USE_SKEW_TYPE,
                           NULL))
           SDDS_Bomb("Invalid -multipoles type");
         s_arg[i_arg].n_items++;
         break;
       case SET_VARIABLES:
         if (s_arg[i_arg].n_items!=3)
           SDDS_Bomb("Invalid -variables syntax");
         columnName[0] = s_arg[i_arg].list[1];
         columnName[1] = s_arg[i_arg].list[2];
         break;
       default:
         fprintf(stderr, "unknown option - %s provided.\n%s", s_arg[i_arg].list[0], USAGE);
         exit(1);
       }  
     } else {
      if (inputFile==NULL)
        inputFile = s_arg[i_arg].list[0];
      else if (outputFile==NULL)
        outputFile = s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many filenames");
    }
  }
  if (!typeFlags) 
    typeFlags |= USE_NORMAL_TYPE | USE_SKEW_TYPE;
  if (!multipoleFlags || multipoleFlags&USE_ALL_MULTIPOLES) 
    multipoleFlags |= USE_DIPOLE | USE_QUADRUPOLE | USE_SEXTUPOLE | USE_OCTUPOLE;
  if (columnName[0]==NULL && (!SDDS_CopyString(columnName+0, "x") || !SDDS_CopyString(columnName+1, "y")))
    SDDS_Bomb("Problem copying column name strings");

  poles = maximumOrder < 4 ? 4 : maximumOrder ;
  pole = calloc(sizeof(*pole), poles);
  multipole = calloc(sizeof(*multipole), poles);
  for (i=0; i<4; i++)
    multipole[i] = multipole0[i];
  /* pole[i] i=0,1,2,3 represents dipole, quadrupole, sextupole and octupole respectively,
     if it is 0, means the corresponding multipole is not chosen and will not be processed.
     if it is not zero, the resonance for the corresponding multipole will be found and write
     to output file.
     the value of pole[i] = Q/2 */
  if (multipoleFlags&USE_DIPOLE) 
    pole[0]=1;
  if (multipoleFlags&USE_QUADRUPOLE)
    pole[1]=2;
  if (multipoleFlags&USE_SEXTUPOLE)
    pole[2]=3;
  if (multipoleFlags&USE_OCTUPOLE)
    pole[3]=4;
  if (maximumOrder>4) {
    char s[1024];
    for (i=4; i<maximumOrder; i++) {
      pole[i] = i+1;
      sprintf(s, "%ld-pole", 2*(i+1));
      cp_str(multipole+i, s);
    }
  }
  
  processFilenames("sddsfindresonances", &inputFile, &outputFile, pipeFlags, 0, &tmpFileUsed);
  if (tmpFileUsed)
    SDDS_Bomb("can't overwrite input file");
  
  
  if (!SDDS_InitializeInput(&inData, inputFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  /* check the input file for valid data */
  for (i=0;i<4;i++) {
    if (SDDS_GetColumnIndex(&inData, columnName[i])<0) {
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      SDDS_Bomb("something wrong with input columns, should have x, y, nux and nuy columns");
    }
  }
  
  if (!SDDS_InitializeOutput(&outData, SDDS_BINARY, 1, NULL, NULL, outputFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"multipole",NULL, NULL, NULL,NULL,SDDS_STRING, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"type",NULL, NULL, NULL,NULL,SDDS_STRING, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"resonanceLabel",NULL, NULL, NULL,NULL,SDDS_STRING, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"multx",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"multy",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"offset",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&outData,"inputPage",NULL, NULL, NULL,NULL,SDDS_LONG, 0)<0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (i=0; i<4; i++) {
    if (!SDDS_TransferColumnDefinition(&outData, &inData, columnName[i], NULL))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(&outData))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  
  /* read indepData (x, y, nux, and nuy) */
  inputPage = 0;
  while (SDDS_ReadPage(&inData)>0) {
    inputPage ++;
    rows=0;
    if (!(rows=SDDS_RowCount(&inData)))
      continue;
    for (i=0; i<4; i++) {
      if (indepData[i])
	free(indepData[i]);
      if (!(indepData[i] = SDDS_GetColumnInDoubles(&inData, columnName[i])))
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    datapoints = rows;

    for (i=0; i<2; i++) {
      if (!i) {
	/* sort by x first and then sort by y; the sorted index is represented by sortIndex[0]*/
	sortData[0] = indepData[0];
	sortData[1] = indepData[1];
      } else {
	/* sort by y first and then sort by x; the sorted index is represented by sortIndex[1]*/
	sortData[0] = indepData[1];
	sortData[1] = indepData[0];
      }
      sortIndex[i] = malloc(sizeof(**sortIndex)*datapoints);
      for (k=0;k<datapoints;k++)
	sortIndex[i][k]=k;
      qsort((void*)sortIndex[i], datapoints, sizeof(**sortIndex), SDDS_CompareDoubleRows);
    }
    
    for (i_pole=0; i_pole<poles; i_pole++) {
      /* process each selected multipole */
      if (!pole[i_pole])
	continue;
      /* multx (M) starts from pole[i_pole] (i.e., Q/2) and 
	 loop with -1 interval until multx reaches 0 */
      multx=pole[i_pole];
      count=0;
      outputpoints=0;
      while (multx>=0) {
	/* multy (N) = Q/2 - multx */
	multy = pole[i_pole] - multx;
	if ((count%2==0 && typeFlags&USE_NORMAL_TYPE) ||
	    (count%2 && typeFlags&USE_SKEW_TYPE)) {
	  /* if multy==0, no need to change the sign, just find resonances for multx*nux = offset
	     otherwise, find resonances for both multx*nux + multy *nuy = offset and
	     multx*nux - multy *nuy = offset */
	  if (!multy) {
	    /* multy is zero */
	    if (multx<0)
	      signs = 0; /* skip multx<0 as it is the same as |multx| */
	    else
	      signs = 1;
	  } else {
	    /* multy is nonzero */
	    if (multx)
	      signs = 2;
	    else {
	      signs = 1;
	      multy = abs(multy);
	    }
	  }
	  for (sign=0; sign<signs; sign++) {
	    if (sign)
	      /* change sign of multy */
	      multy = -multy;
	    start=multx<multy?multx:multy;
	    if (start>0)
	      start=0;
	    end=multx>multy?multx:multy;
	    rows = 0;
	    /* find the resonances for each offset*/
	    for (offset=start; offset<=end; offset++) {
	      /* process two sets of data, first set is the data sorted by x and then y (plane=0),
		 second set is the data sorted by y and then x (plane=1)*/
	      for (plane=0; plane<2; plane++) {
		for (i=0; i<datapoints-1; i++) {
		  index = sortIndex[plane][i];
		  index1 = sortIndex[plane][i+1];
		  if (indepData[plane][index]!=indepData[plane][index1])
		    continue;
		  delta = multx * indepData[2][index] + multy * indepData[3][index] - offset * 1.0;
		  delta1 = multx * indepData[2][index1] + multy * indepData[3][index1] - offset * 1.0;
		  if ((delta<=0 && delta1>=0) || (delta>=0 && delta1<=0)) {
		    if (!outputpoints) {
		      outputpoints=20;
		      for (k=0; k<4; k++) 
			outputData[k]=SDDS_Realloc(outputData[k], sizeof(**outputData)*outputpoints);
		    }
		    if (rows>=outputpoints) {
		      outputpoints +=rows+10;
		      for (k=0; k<4; k++) 
			outputData[k]=SDDS_Realloc(outputData[k], sizeof(**outputData)*outputpoints);
		    }
		    for (k=0; k<4; k++) {
		      if (indepData[k][index] == indepData[k][index1])
			outputData[k][rows] = indepData[k][index];
		      else {
			slope = (delta1 - delta)/(indepData[k][index1] - indepData[k][index]);
			if (slope) 
			  outputData[k][rows] = indepData[k][index] - delta/slope;
			else
			  outputData[k][rows] = (indepData[k][index] + indepData[k][index1])/2;
		      }
		    }
		    rows++;
		  } /*end of if deltaTune[i] [i+1] cross zero */
		} /*end of for i=o to datapoints loop */
	      } /*end of for plane loop */
	    } /* end of for offset loop */
	    if (rows) {
	      sprintf(label, "%ld*$gn$r$bx$n + %ld*$gn$r$by$n", multx, multy);
	      if (!SDDS_StartPage(&outData, rows))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "resonanceLabel", label, NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "multipole", multipole[i_pole], NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "multx", multx, NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "multy", multy, NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (count%2==0) {
		if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
					"type", "normal", NULL))
		  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      } else {
		if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
					"type", "skew", NULL))
		  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      }
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "offset", offset, NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      if (!SDDS_SetParameters(&outData, SDDS_BY_NAME|SDDS_PASS_BY_VALUE,
				      "inputPage", inputPage, NULL))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      for (k=0; k<4; k++) {
		if (!SDDS_SetColumn(&outData,  SDDS_SET_BY_NAME, outputData[k], rows, columnName[k]))
		  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	      }
	      if (!SDDS_WritePage(&outData))
		SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	    }
	  } /*end of for sign  loop */
	} /*end of if */
	count++;
	multx --;
      }
    } /*end of for i_pole */
  }


  for (k=0;k<4;k++)
    if (outputData[k])
      free(outputData[k]);

  if (!SDDS_Terminate(&outData))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  for (i=0; i<2; i++)
    free(sortIndex[i]);
  free(sortIndex);
  for (k=0;k<4;k++)
    free(indepData[k]);
  free(indepData);
  free(outputData);
  free(pole);
  free_scanargs(&s_arg,argc);
  return 0;
}

int SDDS_CompareDoubleRows(const void *vrow1, const void *vrow2)
{
  static long row1, row2;
  double diff;
  long i;
  
  row1 = *(long*)vrow1;
  row2 = *(long*)vrow2;
  for (i=0; i<2; i++) {
    diff = sortData[i][row1] - sortData[i][row2];
    if (diff>0)
      return 1;
    else if (diff<0)
      return -1;
  }
  return 0;
}


