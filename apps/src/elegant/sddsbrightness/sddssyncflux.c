/*
 *
 *  $Log: not supported by cvs2svn $
 *  Revision 1.4  2012/03/01 21:22:12  soliday
 *  Fixed a typo in the usage message.
 *
 *  Revision 1.3  2011/03/16 19:11:47  borland
 *  Fixed units for on-axis flux.
 *
 *  Revision 1.2  2011/03/16 19:09:10  borland
 *  Fixed calculation of total flux for dipole.  Cf. Attwood, 5.8 (multiply by
 *  2*pi since we want the flux in a full circle, not per mrad).
 *
 *  Revision 1.1  2011/03/09 20:44:37  shang
 *  added description
 *
 *  Revision 1.7  2011/03/08 16:25:48  shang
 *  fixed -pipe option bug which did not work.
 *
 *  Revision 1.6  2009/04/27 19:00:59  shang
 *  now uses the GSL bessel functions instead of NR.
 *
 *  Revision 1.5  2009/04/07 22:17:37  shang
 *  modified to print out message for undulator since it is not implemented yet.
 *
 *  Revision 1.4  2007/01/04 16:41:01  shang
 *  modified to the CVS version of bessel G1 value file.
 *
 *  Revision 1.3  2006/10/10 18:19:01  soliday
 *  Updated for WIN32.
 *
 *  Revision 1.2  2005/11/04 22:46:19  soliday
 *  Updated code to be compiled by a 64 bit processor.
 *
 *  Revision 1.1  2004/01/12 20:21:13  shang
 *  first version for calculating photon flux of bend, wigger and undulator magnet. The
 *  calculation for undulator has not been implemented yet.
 *
    first version H. SHANG
*/

#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include <ctype.h>
#ifndef USE_GSL
#error The GSL library must be available to build sddssyncflux.
#endif
#include "gsl/gsl_sf_bessel.h"

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define ALPHA 1/137.03604
#define ECHARGE 1.6021892e-19
#define ESPA 1.0e-10

#define ENERGY_MODE 0
#define WAVELENGTH_MODE 1
#define COLUMN_MODES 2
char *column_mode[COLUMN_MODES]={
  "energy","wavelength",
};


#define LINEAR_MODE 0
#define LOGARITHMIC_MODE 1
#define CHANGE_MODES 2
char *change_mode[CHANGE_MODES] = {
  "linear","logarithmic",
};


#define BENDMAGNET_SOURCE 0
#define WIGGLER_SOURCE 1
#define UNDULATOR_SOURCE 2
#define NUM_SOURCES 3
char *device_source[NUM_SOURCES]= {
  "bendMagnet","wiggler","undulator",
};

#define CLO_VERBOSE 0
#define CLO_FILEVALUES 1
#define CLO_PIPE 2
#define CLO_MODE 3
#define CLO_EBEAMENERGY 4
#define CLO_EBEAMGAMMA 5
#define CLO_EBEAMCURRENT 6
#define CLO_SOURCE 7
#define CLO_G1VALUEFILE 8
#define COMMANDLINE_OPTIONS 9
char *commandline_option[COMMANDLINE_OPTIONS] = {
  "verbose","fileValues","pipe","mode","eBeamEnergy","eBeamGamma","eBeamCurrent","source","g1ValueFile",
};

char *USAGE1 = "sddssyncflux <outputFile> -verbose \n\
   [-fileValues=<filename>[,energy=<columnName or wavelength=columnName>] \n\
   [-mode=energy(wavelength),linear(logarithmic),start=<value>,end=<value>,step(factor)=<value>] \n\
   [-source=bendMagnet[,field=xx[,radius=yy][,criticalEnergy=ZZ]] \n\
     [-source=wiggler(undulator),period=xx[,field=yy][,K=zz],numberOfPeriods=<n>] \n\
   [-eBeamEnergy=<value> [-eBeamCurrent=<value>] [-eBeamGamma=<value>] \n\
<outputFile>       the results are written into the output file. \n\
-verbose           flag for printing messages. \n\
-pipe              output result to the pipe. \n\
-fileValues=<filename>[,energy=<columnName or wavelength=columnName>] \n\
                   get the energy or waveformlength from provided file instead of by -mode option. \n\
-mode=energy,linear,start=<value>,end=<value>,step=<value> \n\
                   Generate photon energy column in eV linearly, from start to end in steps. \n\
-mode=lenergy,ogarithmic,start=<value>,end=<value>,factor=<value> \n\
                   Generate photon energy column logarithmically, from start to end by multiplying\n\
                   factor from point to point.\n";
char *USAGE2 = "-mode=wavelength,linear,start=<value>,end=<value>,step=<value> \n\
                   Generate photon wavelength column in nm linearly, from start to end in steps. \n\
-mode=wavelength,logarithmic,start=<value>,end=<value>,factor=<value> \n\
                   Generate photon wavelength column in nm logarithmically, from start \n\
                   to end by multiplying factor from point to point.\n\
-source=bendMagnet[,field=xx][,radius=yy][,critialEnergy=zz] \n\
                   bend magnet source, magnetic filed=xx Tesla (default=0.6 Teslas). \n\
                   bend radius= yy meter (no default value), criticalEnergy=zz eV (no default value).\n\
                   only one of field, radius and K is needed to be provided. \n\
-source=wiggler,period=xx[,field=yy][,K=zz],numberOfPeriods=n \n\
                   Wiggler source, period=xx cm (default=5 cm). \n\
                   Peak magnetic field=yy Tesla (default=1 Teslas).\n\
                   Undulator parameter, K=zz (no default values) \n\
                   only two of period, field and K is needed to be provided. \n\
                   numberOfPeriods needs to be provided. \n";
char *USAGE3 = "-source=undulator,period=xx[,field=yy][,K=zz],numberOfPerios=n \n\
                   undulator source, period=xx cm (default=5 cm). \n\
                   Peak magnetic field=yy Tesla (default=1 Teslas).\n\
                   Undulator parameter, K=zz (no default values) \n\
                   only two of period, field and K is needed to be provided. \n\
                   numberOfPeriods needs to be provided. \n\
note that only one source is accepted at one time. \n\
-eBeamEnergy       Electron beam energy in Gev, default 7Gev \n\
-eBeamGamma        Electron beam gamma. \n\
-eBeamCurrent      Electron beam current in A. \n\
-g1ValueFile       give the file which contains the values of y and yGy, where \n\
                   yGy=y*intergration of K5/3 from y to infinity. \n\n\
sddssyncflux calculates the photo flux of synchrotron radiation of bending magnet and wiggler, undulator is not being implemented. \n\n";

long ObtainPhotonEnergyAndWavelengthFromFile(char *fileValuesFile, char *ColumnName, long mode, 
                                             double **Energy, double **Wavelength);
long CalculatePhotonEnergyAndWavelengthFromRange(double start, double end, double step, double factor,
                                                 long columnmode, long changemode,
                                                 double **Energy, double **Wavelength);
long CalculateTotalPhotonFluxOfBendMagnet(char *g1ValueFile, double *energy, double **totalFlux,
                                          double criticalEnergy,long n,double eGamma, double eCurrent);
long CalculateOnAxisPhotonFluxOfBendMagnet(double *energy, double **onAxisFlux, double criticalEnergy, long n,double eGamma, double eCurrent);
void GetG1ValuesFromFile (char *g1ValueFile, double **y, double **yGy,long *rows);
long CalculatePhotonFluxOfWiggler(double *energy, double **onAxisFlux, double **totalFlux,
                                 double **PowerPerEV, double **xWidth, double **yWidth,long n,
                                 double eCurrent, double eGamma, double eEnergy, double K, long nPeriods,
                                 double field, char *g1ValueFile);

int main(int argc, char **argv)
{
  SDDS_TABLE SDDSout;
  SCANNED_ARG *s_arg;
  char *outputFile,*fileValuesFile,*energyColumn,*wavelengthColumn, *columnMode,*changeMode,*device,*columnName,*g1ValueFile;
  long i_arg, i,eGammaGiven,eCurrentGiven,eEnergyGiven,fieldGiven;
  long radiusGiven,criticalGiven,periodGiven,sourceGiven,mode,rows,changemode,columnmode,verbose,tmpfile_used;
  unsigned long pipeFlags,flags;
  double start,end,step,factor,eEnergy,eCurrent,eGamma,field,radius,criticalEnergy,period,K,totalPower;
  double *energy,*wavelength,*xWidth,*yWidth,*totalFlux,*onAxisFlux, *PowerPerEV;
  char *input="obset";
  int32_t nPeriods=0; 
  long modeGiven=0;
  

  SDDS_RegisterProgramName(argv[0]);
  argc = scanargs(&s_arg, argc, argv); 
  if (argc<3 || argc>(3+COMMANDLINE_OPTIONS)) {
    fprintf(stderr, "%s%s%s\n", USAGE1, USAGE2, USAGE3);
    exit(1);
  }
  start=end=0;
  outputFile=fileValuesFile=g1ValueFile=NULL;
  energyColumn=wavelengthColumn=columnName=NULL;
  eEnergy=7.0;
  eCurrent=1.0;
  pipeFlags=flags=0;
  eGammaGiven=eCurrentGiven=eEnergyGiven=0;
  columnMode=changeMode=device=NULL;
  mode=changemode=columnmode=0;
  step=factor=0;
  field=radius=criticalEnergy=period=0;
  fieldGiven=radiusGiven=periodGiven=criticalGiven=sourceGiven=0;
  energy=wavelength=xWidth=yWidth=totalFlux=onAxisFlux=NULL;
  rows=tmpfile_used=0;
  K=0;
  PowerPerEV=NULL;

  for (i_arg=1; i_arg<argc; i_arg++) {
    if (s_arg[i_arg].arg_type==OPTION) {
      switch (match_string(s_arg[i_arg].list[0],commandline_option,COMMANDLINE_OPTIONS,0)) {
      case CLO_VERBOSE:
        verbose=1;
        break;
      case CLO_PIPE:
	if (!processPipeOption(s_arg[i_arg].list+1, s_arg[i_arg].n_items-1, &pipeFlags)) {
          fprintf(stderr, "Error (%s): invalid -pipe syntax\n", argv[0]);
          return(1);
        }
	if (pipeFlags!=USE_STDOUT)
          SDDS_Bomb("only -pipe=out syntax is valid!");
        break;
      case CLO_G1VALUEFILE:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -fileValues syntax");
        g1ValueFile=s_arg[i_arg].list[1];
        break;
      case CLO_FILEVALUES:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -fileValues syntax");
        fileValuesFile = s_arg[i_arg].list[1];
        s_arg[i_arg].n_items -= 2;
        if (!scanItemList(&flags, s_arg[i_arg].list+2, &s_arg[i_arg].n_items, 0,
                          "energy", SDDS_STRING, &energyColumn, 1, 0,
                          "wavelength", SDDS_STRING, &wavelengthColumn, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -fileValues syntax");
        if (!energyColumn && !wavelengthColumn)
          SDDS_Bomb("invalid -fileValues syntax, neither energyColumn nor wavelengthColumn are provided");
        if (energyColumn) {
          columnName=energyColumn;
          mode=ENERGY_MODE;
        } else {
          columnName=wavelengthColumn;
          mode=WAVELENGTH_MODE;
        }
        break;
      case CLO_MODE:
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -mode syntax");
        columnMode = s_arg[i_arg].list[1];
        changeMode = s_arg[i_arg].list[2];
        columnmode=match_string(columnMode,column_mode,COLUMN_MODES,0);
        changemode=match_string(changeMode,change_mode,CHANGE_MODES,0);
        if (columnmode<0 || changemode<0) 
          SDDS_Bomb("invalid -mode syntax, it should be -mode=energy(wavelength),linear(logarithmic)...");
        s_arg[i_arg].n_items -= 3;
        if (!scanItemList(&flags, s_arg[i_arg].list+3, &s_arg[i_arg].n_items, 0,
                          "start", SDDS_DOUBLE, &start, 1, 0,
                          "end", SDDS_DOUBLE, &end, 1, 0,
                          "step", SDDS_DOUBLE, &step, 1, 0,
                          "factor",SDDS_DOUBLE, &factor, 1, 0,
                          NULL))
          SDDS_Bomb("invalid -mode syntax");
        if (changemode==LINEAR_MODE && !step)
          SDDS_Bomb("invalid -mode syntax, step value is not given for linear changing mode!");
        if (columnmode==LOGARITHMIC_MODE && !factor)
          SDDS_Bomb("invalid -mode syntax, factor value is not given for logarithmic changing mode!");
        if (!start || !end)
          SDDS_Bomb("invalid -mode syntax,  the start and end value should be provided!");
        modeGiven=1;
        break;
      case CLO_SOURCE:
        if (sourceGiven) {
          fprintf(stderr,"Warning: multiply sources provided, only the first one will be taken.\n");
          break;
        } 
        if (s_arg[i_arg].n_items<2)
          SDDS_Bomb("invalid -source syntax");
        device=s_arg[i_arg].list[1];
        if (match_string(device,device_source,NUM_SOURCES,0)<0)
          SDDS_Bomb("Invalid -source syntax1!");
        s_arg[i_arg].n_items -= 2;
        if (!scanItemList(&flags, s_arg[i_arg].list+2, &s_arg[i_arg].n_items, 0,
                          "period", SDDS_DOUBLE, &period, 1, 0,
                          "field", SDDS_DOUBLE, &field, 1, 0,
                          "radius", SDDS_DOUBLE, &radius, 1, 0,
                          "K",SDDS_DOUBLE, &K, 1, 0,
                          "numberOfPeriods",SDDS_LONG,&nPeriods,1,0,
                          "criticalEnergy",SDDS_DOUBLE,&criticalEnergy,1,0,
                          NULL))
          SDDS_Bomb("invalid -source syntax2");
        sourceGiven=1;
        break;
      case CLO_EBEAMENERGY:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -eBeamEnergy syntax");
        if (sscanf(s_arg[i_arg].list[1], "%lf", &eEnergy)!=1)
          SDDS_Bomb("invalid -eBeamEnergy value");
        eEnergyGiven=1;
        break;
      case CLO_EBEAMGAMMA:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -eBeamGamma syntax");
        if (sscanf(s_arg[i_arg].list[1], "%lf", &eGamma)!=1)
          SDDS_Bomb("invalid -eBeamGamma value");
        eGammaGiven=1;
        break;
      case CLO_EBEAMCURRENT:
        if (s_arg[i_arg].n_items!=2)
          SDDS_Bomb("invalid -eBeamCurrent syntax");
        if (sscanf(s_arg[i_arg].list[1], "%lf", &eCurrent)!=1)
          SDDS_Bomb("invalid -eBeamCurrent value");
        eCurrentGiven=1;
        break;
      default:
        SDDS_Bomb("Unknown argmuments");
        break;
      }
    } else {
      if (!outputFile)
        outputFile=s_arg[i_arg].list[0];
      else
        SDDS_Bomb("too many files given!");
    }
  }
  if (!modeGiven && !fileValuesFile) 
    SDDS_Bomb("Either fileValues or mode option has to be provided!");
  if (!outputFile && !pipeFlags)
    SDDS_Bomb("Either output file or -pipe[=out] should be given!");
  if (outputFile && pipeFlags)
    SDDS_Bomb("Only one of output file and -pipe should be given!");
  processFilenames("sddssyncflux",&input, &outputFile, pipeFlags,1,&tmpfile_used);
  if (eEnergyGiven && eGammaGiven && verbose)  {
    fprintf(stderr,"Warning: electron gamma will be recalculated when electron energy is provided!");
    eGammaGiven=0;
  }
  if (!eEnergyGiven && !eGammaGiven) {
    eEnergyGiven=1;
    eEnergy=7.0;
  }  
  if (eEnergyGiven) 
    eGamma=eEnergy/0.000511003;
  if (eGammaGiven)
    eEnergy=eGamma*0.000511003;
  switch (match_string(device,device_source,NUM_SOURCES,0)) {
  case BENDMAGNET_SOURCE:
    if (field) fieldGiven=1;
    else if (radius) radiusGiven=1;
    else if (criticalEnergy) criticalGiven=1;
    else {
      field=0.6; /*default value of bend magnet field in Tesla */
      fieldGiven=1;
    }
    if (fieldGiven) {
      radius=3.3*eEnergy/field;
      criticalEnergy=665*eEnergy*eEnergy*field;  
    } else if (radiusGiven) {
      field=3.3*eEnergy/radius;
      criticalEnergy=665*eEnergy*eEnergy*field;  
    } else if (criticalGiven) {
      field=criticalEnergy/665.0/pow(eEnergy,2);
      radius=3.3*eEnergy/field;
    }
    break;
  case WIGGLER_SOURCE:
  case UNDULATOR_SOURCE:
    if (!nPeriods) 
      SDDS_Bomb("The number of periods is not provided for wiggler/undulator source!");
    if (!period && !field && !K) {
      period=5; /*default value  period in cm */
      field=1; /*default value of field in Tesla */
    }
    if (period && field)
      K=0.934*period*field;
    else if (period && K)
      field=K/0.934/period;
    else if (K && field)
      period=K/0.934/field;
    criticalEnergy=665*eEnergy*eEnergy*field;
    break;
  }
  /*setup output file */
  if (!SDDS_InitializeOutput(&SDDSout, SDDS_BINARY, 0, NULL, NULL, outputFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_DefineParameter(&SDDSout, "eBeamEnergy", NULL, "Gev", "Electron Beam Energy",
                                   NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "eBeamCurrent", NULL, "A", "Electron Beam Current",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "eBeamGamma", NULL, NULL, "Electron Beam Gamma",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "Device", NULL, NULL, NULL,
                           NULL, SDDS_STRING, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "Field", NULL, "Tesla", NULL,
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "Radius", NULL, "m", "Bend radius of bend magnet",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "Period", NULL, "cm", "Period length of wiggler or undulator",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "K", NULL, NULL, "Wiggler or Undulator parameter",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "CriticalEnergy", NULL, "eV", "Photon critical energy",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineParameter(&SDDSout, "TotalPower", NULL, "W", "Total Power",
                           NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "Energy", NULL, "eV", "Photon Energy",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "Wavelength", NULL, "nm", "Photon Wavelength",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "TotalFlux", NULL, "photons/second/0.1%BW", 
                        "angle-integrated photon flux",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "OnAxisFlux", NULL, "photons/mrad$a2$n/second/0.1%BW", 
                        "on-axis angular photon flux density",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "XWidth", NULL, "mrad", 
                        "horizontal rms width of the x-ray beam",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "YWidth", NULL, "mrad", 
                        "vertical rms width of the x-ray beam",NULL, SDDS_DOUBLE, 0)<0 ||
      SDDS_DefineColumn(&SDDSout, "PowerPerEV", NULL, NULL, 
                        NULL,NULL, SDDS_DOUBLE, 0)<0 ||
      !SDDS_WriteLayout(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (fileValuesFile) 
    rows=ObtainPhotonEnergyAndWavelengthFromFile(fileValuesFile,columnName,mode,&energy,&wavelength);
  else
    rows=CalculatePhotonEnergyAndWavelengthFromRange(start,end,step,factor,columnmode,changemode,&energy,&wavelength);
  if (!rows)
    SDDS_Bomb("Unable to get the photon energy and wavelength!");
  switch (match_string(device,device_source,NUM_SOURCES,0)) {
  case  BENDMAGNET_SOURCE:
    if (!CalculateTotalPhotonFluxOfBendMagnet(g1ValueFile,energy,&totalFlux,criticalEnergy,rows,eGamma,eCurrent))
      SDDS_Bomb("Error in calculating total flux of bend magnet!");
    if (!CalculateOnAxisPhotonFluxOfBendMagnet(energy,&onAxisFlux,criticalEnergy,rows,eGamma,eCurrent))
      SDDS_Bomb("Error in calculating on-axis flux of bend magnet!");
    xWidth=(double*)malloc(sizeof(*xWidth)*rows);
    yWidth=(double*)malloc(sizeof(*yWidth)*rows);
    PowerPerEV=(double*)malloc(sizeof(*PowerPerEV)*rows);
    for (i=0;i<rows;i++) {
      xWidth[i]=1;
      yWidth[i]=totalFlux[i]/sqrt(2*PI)/onAxisFlux[i];
      PowerPerEV[i]=1000*ECHARGE*totalFlux[i];
    }
    break;
  case WIGGLER_SOURCE:
    CalculatePhotonFluxOfWiggler(energy,&onAxisFlux,&totalFlux,&PowerPerEV,&xWidth,&yWidth,rows,eCurrent,
                                eGamma,eEnergy,K,nPeriods,field,g1ValueFile);
    break;
  case UNDULATOR_SOURCE:
    fprintf(stderr, "Udulator not implemented yet.\n");
    exit(1);
    break;
  }
  
  totalPower=1.27*eEnergy*eEnergy*field*field*eCurrent*radius*0.001;
  if (!SDDS_StartPage(&SDDSout, rows) || 
      !SDDS_SetParameters(&SDDSout,SDDS_BY_NAME|SDDS_PASS_BY_VALUE,"eBeamEnergy",eEnergy,
                          "eBeamCurrent",eCurrent,"eBeamGamma",eGamma,"Device",device,
                          "Field",field,"Radius",radius,"Period",period,"K",K,
                          "CriticalEnergy",criticalEnergy,"TotalPower",totalPower,NULL) ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,energy,rows,"Energy") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,wavelength,rows,"Wavelength") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,totalFlux,rows,"TotalFlux") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,onAxisFlux,rows,"OnAxisFlux") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,xWidth,rows,"XWidth") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,yWidth,rows,"YWidth") ||
      !SDDS_SetColumnFromDoubles(&SDDSout, SDDS_SET_BY_NAME,PowerPerEV,rows,"PowerPerEV") )
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WritePage(&SDDSout) || !SDDS_Terminate(&SDDSout))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  free(energy);
  free(wavelength);
  free(xWidth);
  free(yWidth);
  free(totalFlux);
  free(onAxisFlux);
  if (energyColumn) free(energyColumn);
  if (wavelengthColumn) free(wavelengthColumn);
  
  free_scanargs(&s_arg,argc);
  
  return 0;
}
  
long ObtainPhotonEnergyAndWavelengthFromFile(char *fileValuesFile, char *columnName, long mode,
                                             double **Energy, double **Wavelength)
{
  SDDS_TABLE inTable;
  long rows,i;
  
  *Energy=*Wavelength=NULL;
  if (!SDDS_InitializeInput(&inTable,fileValuesFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_ReadPage(&inTable)<=0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  rows=SDDS_CountRowsOfInterest(&inTable);
  switch (mode) {
    case ENERGY_MODE:
    if (!(*Energy=SDDS_GetColumnInDoubles(&inTable, columnName)))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
   *Wavelength=(double*)malloc(sizeof(**Wavelength)*rows);
    for (i=0;i<rows;i++) 
      (*Wavelength)[i]=1239.85/(*Energy)[i];
    break;
  case WAVELENGTH_MODE:
    if (!(*Wavelength=SDDS_GetColumnInDoubles(&inTable, columnName)))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    *Energy=(double*)malloc(sizeof(**Energy)*rows);
    for (i=0;i<rows;i++) 
      (*Energy)[i]=1239.85/(*Wavelength)[i];
    break;
  default:
    SDDS_Bomb("Unkonw mode provided!");
    break;
  }
  if (!SDDS_Terminate(&inTable))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  return rows;
}

long CalculatePhotonEnergyAndWavelengthFromRange(double start, double end, double step, double factor,
                                              long columnmode, long changemode,
                                              double **Energy, double **Wavelength)
{
  long rows;
  double currentVal;
  
  *Energy=*Wavelength=NULL;
  rows=0;
  currentVal=start;
  do {
    *Energy=SDDS_Realloc(*Energy,sizeof(**Energy)*(rows+1));
    *Wavelength=SDDS_Realloc(*Wavelength,sizeof(**Wavelength)*(rows+1));
    switch (columnmode) {
    case ENERGY_MODE:
      (*Energy)[rows]=currentVal;
      (*Wavelength)[rows]=1239.85/currentVal;
      break;
    case WAVELENGTH_MODE:
      (*Wavelength)[rows]=currentVal;
      (*Energy)[rows]=1239.85/currentVal;
      break;
    }
    switch (changemode) {
    case LINEAR_MODE:
      if (!step)
        SDDS_Bomb("The step value is not given for linear changing mode!");
      if (step<0)
        SDDS_Bomb("The step should be greater than 0, otherwise, infinity loop expected!");
      currentVal +=step;
      break;
    case LOGARITHMIC_MODE:
      if (!factor)
        SDDS_Bomb("The factor value is not given for logarithmic changing mode!");
      if (factor<1)
        SDDS_Bomb("The factor value should be greater than 1, otherwise, infinity loop expected!");
      currentVal *=factor;
      break;
    }
    rows++;
  } while (currentVal<=end);
  
  return (rows);
}

long CalculateTotalPhotonFluxOfBendMagnet(char *g1ValueFile, double *energy, double **totalFlux,
  double criticalEnergy,long n, double eGamma, double eCurrent) 
{
  double *y,*yGy;
  double g1;
  long rows=0,i;
  unsigned long interpCode;
  OUTRANGE_CONTROL belowRange,aboveRange;

  if (!n || !energy)
    return 0;
  GetG1ValuesFromFile(g1ValueFile,&y,&yGy,&rows);
  aboveRange.flags = belowRange.flags = OUTRANGE_VALUE;
  belowRange.value=aboveRange.value=0;
  *totalFlux=(double*)malloc(sizeof(**totalFlux)*n);
  for (i=0;i<n;i++) {
    g1=interpolate(yGy,y,rows,energy[i]/criticalEnergy,&belowRange,&aboveRange,1,&interpCode,1);
    (*totalFlux)[i]=1.0e-6*sqrt(3)*ALPHA/ECHARGE*eGamma*eCurrent*g1;
  }
  return 1;
}

void GetG1ValuesFromFile (char *g1ValueFile, double **y, double **yGy, long *rows)
{
  SDDS_TABLE g1Table;
  
  *y=*yGy=NULL;
  if (!g1ValueFile) 
    SDDS_CopyString(&g1ValueFile,"/usr/local/oag/apps/configData/syncflux/besselG1.result");
  if (!SDDS_InitializeInput(&g1Table,g1ValueFile))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (SDDS_ReadPage(&g1Table)<=0)
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!(*rows=SDDS_CountRowsOfInterest(&g1Table)))
    SDDS_Bomb("No data found in the g1 value file!");
  if (!(*y=SDDS_GetColumnInDoubles(&g1Table, "y")) || !(*yGy=SDDS_GetColumnInDoubles(&g1Table, "yGy")))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_Terminate(&g1Table))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
}



long CalculateOnAxisPhotonFluxOfBendMagnet(double *energy, double **onAxisFlux, double criticalEnergy,
                                           long n,double eGamma, double eCurrent)
{
  double rk;
  double xnu,x,factor;
  long i;
  
  xnu=2.0/3.0;  
  if (!n || !energy)
    return 0;
  *onAxisFlux=(double*)malloc(sizeof(**onAxisFlux)*n);
  factor=3.0e-9*ALPHA/4.0/(PI*PI)/ECHARGE;
  for (i=0;i<n;i++) {
    x=energy[i]/criticalEnergy;
    rk = gsl_sf_bessel_Knu(xnu, x/2);
    (*onAxisFlux)[i]=factor*eGamma*eGamma*eCurrent*x*x*rk*rk;
  }
  return 1;
}


long CalculatePhotonFluxOfWiggler(double *energy, double **onAxisFlux, double **totalFlux,
                                 double **PowerPerEV, double **xWidth, double **yWidth,long n,
                                 double eCurrent, double eGamma, double eEnergy, double K, long nPeriods,
                                 double field, char *g1ValueFile)
{
  long N,i,j,rows=0;
  unsigned long interpCode;
  double deltaPhi, thetaMax, pasin,totalflux,onaxisflux,offaxisflux,phi,theta,B,criticalEnergy;
  double totalflux1, onaxisflux1,g1,xnu,h2,sigma,deltaTheta;
  double *y,*yGy;
  OUTRANGE_CONTROL belowRange,aboveRange;

  N=max(K*15, 30);
  deltaPhi=PI/(4*N);
  thetaMax=1000*K/eGamma;
  pasin=500/eGamma;
  xnu=2.0/3.0;

  *onAxisFlux=(double*)malloc(sizeof(**onAxisFlux)*n);
  *totalFlux=(double*)malloc(sizeof(**totalFlux)*n);
  *PowerPerEV=(double*)malloc(sizeof(**PowerPerEV)*n);
  *xWidth=(double*)malloc(sizeof(**xWidth)*n);
  *yWidth=(double*)malloc(sizeof(**yWidth)*n);
  
  aboveRange.flags = belowRange.flags = OUTRANGE_VALUE;
  belowRange.value=aboveRange.value=0;
  GetG1ValuesFromFile(g1ValueFile,&y,&yGy,&rows);

  for (i=0;i<n;i++) {
    totalflux=0;
    onaxisflux=0;
    offaxisflux=0;
    for (j=0;j<N;j++) {
      phi=j*deltaPhi;
      theta=thetaMax*sin(phi);
      B=field*cos(phi);
      criticalEnergy=665*eEnergy*eEnergy*B;
      g1=interpolate(yGy,y,rows,energy[i]/criticalEnergy,&belowRange,&aboveRange,1,&interpCode,1);
      totalflux1=1.0e-6*sqrt(3)*ALPHA/2/PI/ECHARGE*eGamma*eCurrent*g1;
      h2=gsl_sf_bessel_Knu(xnu,energy[i]/criticalEnergy/2.0)*energy[i]/criticalEnergy;
      onaxisflux1=3.0e-9*ALPHA/4/(PI*PI)/ECHARGE*eGamma*eGamma*eCurrent*h2*h2;
      sigma=totalflux1/sqrt(2*PI)/onaxisflux1;
      deltaTheta=thetaMax*cos(phi)*deltaPhi;
      totalflux +=totalflux1*deltaTheta;
      onaxisflux +=onaxisflux1*deltaTheta*exp(-theta*theta/(2*sigma*sigma));
      offaxisflux += onaxisflux1*deltaTheta*exp(-(theta*theta+pasin*pasin)/(2*sigma*sigma));
    }
    (*totalFlux)[i]=totalflux*4*nPeriods;
    (*onAxisFlux)[i]=onaxisflux*4*nPeriods;
    (*PowerPerEV)[i]=1000*ECHARGE*(*totalFlux)[i];
    (*yWidth)[i]=pasin/sqrt(2*log((*onAxisFlux)[i]/offaxisflux));
    (*xWidth)[i]=(*totalFlux)[i]/2/PI/(*yWidth)[i]/(*onAxisFlux)[i];
  }
  return 0;
}
