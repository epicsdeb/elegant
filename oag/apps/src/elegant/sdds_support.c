/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: sdds_support.c
 * purpose: routines for working with SDDS files
 * 
 * Michael Borland, 1993.
 */
#include "mdb.h"
#include "fftpackC.h"
#include "track.h"
#include "SDDS.h"

void SDDS_ElegantOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                             char *contents, char *command_file, char *lattice_file, SDDS_DEFINITION *parameter_definition,
                             long n_parameters, SDDS_DEFINITION *column_definition, long n_columns,
                             char *caller, long flags)
{
    static char description[SDDS_MAXLINE];
    long i, last_index, index;

    log_entry("SDDS_ElegantOutputSetup");
    last_index = -1;

#if USE_MPI && !SDDS_MPI_IO 
    if (myid!=0)
      return;
#endif

    if (flags&SDDS_EOS_NEWFILE) {
        if (SDDS_IsActive(SDDS_table)==1) {
            if (!SDDS_Terminate(SDDS_table)) {
                fprintf(stdout, "Unable to terminate SDDS output (%s)\n", caller);
                fflush(stdout);
                SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
                exit(1);
                }
            }
        sprintf(description, "%s--input: %s  lattice: %s", contents, command_file, lattice_file);
#if SDDS_MPI_IO
	/* We still need serial IO for some operations */
	if ((SDDS_table->parallel_io && (!SDDS_Parallel_InitializeOutput(SDDS_table, description, contents, filename))) ||
	    (!SDDS_table->parallel_io && (!SDDS_InitializeOutput(SDDS_table, mode, 1, description, contents, filename))))
#else
        if (!SDDS_InitializeOutput(SDDS_table, mode, 1, description, contents, filename)) 
#endif
	{
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
	}
    }

    /* define SDDS parameters */
    for (i=0; i<n_parameters; i++) {
        if (!SDDS_ProcessParameterString(SDDS_table, parameter_definition[i].text, 0) ||
            (index=SDDS_GetParameterIndex(SDDS_table, parameter_definition[i].name))<0) {
            fprintf(stdout, "Unable to define SDDS parameter for %s--string was:\n%s\n", caller,
                    parameter_definition[i].text);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        if (last_index==-1 && flags&SDDS_EOS_NEWFILE && index!=0)
            fprintf(stdout, "\7\7\7WARNING: first defined parameter index for SDDS file %s is not 0--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        if (last_index!=-1 && index!=(last_index+1))
            fprintf(stdout, "\7\7\7WARNING: parameter indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        last_index = index;
        }

    last_index = -1;
    /* define SDDS columns */
    for (i=0; i<n_columns; i++) {
        if (!SDDS_ProcessColumnString(SDDS_table, column_definition[i].text, 0) ||
            (index=SDDS_GetColumnIndex(SDDS_table, column_definition[i].name))<0) {
            fprintf(stdout, "Unable to define SDDS column for %s--string was:\n%s\n", caller,
                    column_definition[i].text);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        if (last_index==-1 && flags&SDDS_EOS_NEWFILE && index!=0)
            fprintf(stdout, "\7\7\7WARNING: first defined column index for SDDS file %s is not 0--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        if (last_index!=-1 && index!=(last_index+1))
            fprintf(stdout, "\7\7\7WARNING: column indices for SDDS file %s are not sequential--this will probably cause unexpected results\n", filename);
            fflush(stdout);
        last_index = index;
        }

#if SDDS_MPI_IO 
  /* In the case of parallel IO, the WiteLayout will be called at the time it dumps data or setups the output,
     as the communicator information is required */
    if (!SDDS_table->parallel_io)
#endif
    if (flags&SDDS_EOS_COMPLETE) {
        if (!SDDS_WriteLayout(SDDS_table)) {
            fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
            fflush(stdout);
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
            exit(1);
            }
        }

    log_exit("SDDS_ElegantOutputSetup");
    }


#define STANDARD_PARAMETERS 1
static SDDS_DEFINITION standard_parameter[STANDARD_PARAMETERS] = {
    {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
    } ;

#define ELEMENT_COLUMNS 4
static SDDS_DEFINITION element_column[ELEMENT_COLUMNS] = {
    {"s", "&column name=s, units=m, type=double, description=\"Distance\" &end"},
    {"ElementName", "&column name=ElementName, type=string, description=\"Element name\", format_string=%10s &end"},
    {"ElementOccurence", 
         "&column name=ElementOccurence, type=long, description=\"Occurence of element\", format_string=%6ld &end"},
    {"ElementType", "&column name=ElementType, type=string, description=\"Element-type name\", format_string=%10s &end"},
    } ;

#define PHASE_SPACE_COLUMNS 7
static SDDS_DEFINITION phase_space_column[PHASE_SPACE_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"t", "&column name=t, units=s, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    } ;

#define PHASE_SPACE_PARAMETERS 4
static SDDS_DEFINITION phase_space_parameter[PHASE_SPACE_PARAMETERS] = {
  {"Step", "&parameter name=Step, type=long, description=\"Simulation step\" &end"},
  {"pCentral", "&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double, description=\"Reference beta*gamma\" &end"},
  {"Charge", "&parameter name=Charge, type=double, units=C, description=\"Beam charge\" &end"},
  {"Particles", "&parameter name=Particles, type=long, description=\"Number of particles\" &end"},
};
  
void SDDS_PhaseSpaceSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_PhaseSpaceSetup");
#if SDDS_MPI_IO
    SDDS_table->parallel_io = 1;
    /* set up parallel IO information */      
    SDDS_MPI_Setup(SDDS_table, 1, n_processors, myid, MPI_COMM_WORLD, 0);
#endif
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            phase_space_parameter, PHASE_SPACE_PARAMETERS, phase_space_column, PHASE_SPACE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    log_exit("SDDS_PhaseSpaceSetup");
    }


#define BEAM_LOSS_COLUMNS 8
static SDDS_DEFINITION beam_loss_column[BEAM_LOSS_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"s", "&column name=s, units=m, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    {"Pass", "&column name=Pass, type=long &end"},
    } ;

void SDDS_BeamLossSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_BeamLossSetup");
#if SDDS_MPI_IO  
    SDDS_table->parallel_io = 1;
    /* set up parallel IO information */      
    SDDS_MPI_Setup(SDDS_table, 1, n_processors, myid, MPI_COMM_WORLD, 0);
#endif

    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, beam_loss_column, BEAM_LOSS_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    log_exit("SDDS_BeamLossSetup");
    }


#define CENTROID_COLUMNS 8
static SDDS_DEFINITION centroid_column[CENTROID_COLUMNS] = {
    {"Cx", "&column name=Cx, symbol=\"<x>\", units=m, type=double, description=\"x centroid\" &end"},
    {"Cxp", "&column name=Cxp, symbol=\"<x'>\", type=double, description=\"x' centroid\" &end"},
    {"Cy", "&column name=Cy, symbol=\"<y>\", units=m, type=double, description=\"y centroid\" &end"},
    {"Cyp", "&column name=Cyp, symbol=\"<y'>\", type=double, description=\"y' centroid\" &end"},
    {"Cs", "&column name=Cs, symbol=\"<s>\", units=m, type=double, description=\"mean distance traveled\" &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"<$gd$r>\", type=double, description=\"delta centroid\" &end"},
    {"Particles", "&column name=Particles, description=\"Number of particles\", type=long &end"},
    {"pCentral", "&column name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double, description=\"Reference beta*gamma\" &end"},
    } ;

void SDDS_CentroidOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_CentroidOutputSetup");
#if USE_MPI  
    if (myid!=0)
      return;
#if SDDS_MPI_IO
     SDDS_table->parallel_io = 0;
#endif
#endif

    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, element_column, ELEMENT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE);
    SDDS_ElegantOutputSetup(SDDS_table, NULL, 0, 0, NULL, NULL, NULL, NULL, 0, 
                            centroid_column, CENTROID_COLUMNS, caller, SDDS_EOS_COMPLETE);
    log_exit("SDDS_CentroidOutputSetup");
    }

#define SIGMA_MATRIX_COLUMNS 57
static SDDS_DEFINITION sigma_matrix_column[SIGMA_MATRIX_COLUMNS] = {
    {"s1",    "&column name=s1, symbol=\"$gs$r$b1$n\", units=m, type=double, description=\"sqrt(<x*x>)\" &end"},
    {"s12",    "&column name=s12, symbol=\"$gs$r$b12$n\", units=m, type=double, description=\"<x*xp'>\" &end"},
    {"s13",    "&column name=s13, symbol=\"$gs$r$b13$n\", units=\"m$a2$n\", type=double, description=\"<x*y>\" &end"},
    {"s14",    "&column name=s14, symbol=\"$gs$r$b14$n\", units=m, type=double, description=\"<x*y'>\" &end"},
    {"s15",    "&column name=s15, symbol=\"$gs$r$b15$n\", units=\"m$a2$n\", type=double, description=\"<x*s>\" &end"},
    {"s16",    "&column name=s16, symbol=\"$gs$r$b16$n\", units=m, type=double, description=\"<x*delta>\" &end"},
    {"s2",    "&column name=s2, symbol=\"$gs$r$b2$n\", type=double, description=\"sqrt(<x'*x'>)\" &end"},
    {"s23",    "&column name=s23, symbol=\"$gs$r$b23$n\", units=m, type=double, description=\"<x'*y>\" &end"},
    {"s24",    "&column name=s24, symbol=\"$gs$r$b24$n\", type=double, description=\"<x'*y'>\" &end"},
    {"s25",    "&column name=s25, symbol=\"$gs$r$b25$n\", units=m, type=double, description=\"<x'*s>\" &end"},
    {"s26",    "&column name=s26, symbol=\"$gs$r$b26$n\", type=double, description=\"<x'*delta>\" &end"},
    {"s3",    "&column name=s3, symbol=\"$gs$r$b3$n\", units=m, type=double, description=\"sqrt(<y*y>)\" &end"},
    {"s34",    "&column name=s34, symbol=\"$gs$r$b34$n\", units=m, type=double, description=\"<y*y'>\" &end"},
    {"s35",    "&column name=s35, symbol=\"$gs$r$b35$n\", units=\"m$a2$n\", type=double, description=\"<y*s>\" &end"},
    {"s36",    "&column name=s36, symbol=\"$gs$r$b36$n\", units=m, type=double, description=\"<y*delta>\" &end"},
    {"s4",    "&column name=s4, symbol=\"$gs$r$b4$n\", type=double, description=\"sqrt(<y'*y')>\" &end"},
    {"s45",    "&column name=s45, symbol=\"$gs$r$b45$n\", units=m, type=double, description=\"<y'*s>\" &end"},
    {"s46",    "&column name=s46, symbol=\"$gs$r$b46$n\", type=double, description=\"<s'*delta>\" &end"},
    {"s5",    "&column name=s5, symbol=\"$gs$r$b5$n\", units=m, type=double, description=\"sqrt(<s*s>)\" &end"},
    {"s56",    "&column name=s56, symbol=\"$gs$r$b56$n\", units=m, type=double, description=\"<s*delta>\" &end"},
    {"s6",    "&column name=s6, symbol=\"$gs$r$b6$n\", type=double, description=\"sqrt(<delta*delta>)\" &end"},
    {"ma1",    "&column name=ma1, symbol=\"max$sb$ex$sb$e\", units=m, type=double, description=\"maximum absolute value of x\" &end"},
    {"ma2",    "&column name=ma2, symbol=\"max$sb$ex'$sb$e\", type=double, description=\"maximum absolute value of x'\" &end"},
    {"ma3",    "&column name=ma3, symbol=\"max$sb$ey$sb$e\", units=m, type=double, description=\"maximum absolute value of y\" &end"},
    {"ma4",    "&column name=ma4, symbol=\"max$sb$ey'$sb$e\", type=double, description=\"maximum absolute value of y'\" &end"},
    {"ma5",    "&column name=ma5, symbol=\"max$sb$e$gD$rs$sb$e\", type=double, units=m, description=\"maximum absolute value of s\" &end"},
    {"ma6",    "&column name=ma6, symbol=\"max$sb$e$gd$r$sb$e\", type=double, description=\"maximum absolute value of delta\" &end"},
    {"minimum1",   "&column name=minimum1, symbol=\"x$bmin$n\", type=double, units=m &end"},
    {"minimum2",   "&column name=minimum2, symbol=\"x'$bmin$n\", type=double, units=m &end"},
    {"minimum3",   "&column name=minimum3, symbol=\"y$bmin$n\", type=double, units=m &end"},
    {"minimum4",   "&column name=minimum4, symbol=\"y'$bmin$n\", type=double, units=m &end"},
    {"minimum5",   "&column name=minimum5, symbol=\"$gD$rs$bmin$n\", type=double, units=m &end"},
    {"minimum6",   "&column name=minimum6, symbol=\"$gd$r$bmin$n\", type=double, units=m &end"},
    {"maximum1",   "&column name=maximum1, symbol=\"x$bmax$n\", type=double, units=m &end"},
    {"maximum2",   "&column name=maximum2, symbol=\"x'$bmax$n\", type=double, units=m &end"},
    {"maximum3",   "&column name=maximum3, symbol=\"y$bmax$n\", type=double, units=m &end"},
    {"maximum4",   "&column name=maximum4, symbol=\"y'$bmax$n\", type=double, units=m &end"},
    {"maximum5",   "&column name=maximum5, symbol=\"$gD$rs$bmax$n\", type=double, units=m &end"},
    {"maximum6",   "&column name=maximum6, symbol=\"$gd$r$bmax$n\", type=double, units=m &end"},
    {"Sx",    "&column name=Sx, symbol=\"$gs$r$bx$n\", units=m, type=double, description=\"sqrt(<(x-<x>)^2>)\" &end"},
    {"Sxp",    "&column name=Sxp, symbol=\"$gs$r$bx'$n\", type=double, description=\"sqrt(<(x'-<x'>)^2>)\" &end"},
    {"Sy",    "&column name=Sy, symbol=\"$gs$r$by$n\", units=m, type=double, description=\"sqrt(<(y-<y>)^2>)\" &end"},
    {"Syp",    "&column name=Syp, symbol=\"$gs$r$by'$n\", type=double, description=\"sqrt(<(y'-<y'>)^2>)\" &end"},
    {"Ss",    "&column name=Ss, units=m, symbol=\"$gs$r$bs$n\", type=double, description=\"sqrt(<(s-<s>)^2>)\" &end"},
    {"Sdelta",    "&column name=Sdelta, symbol=\"$gs$bd$n$r\", type=double, description=\"sqrt(<(delta-<delta>)^2>)\" &end"},
    {"ex", "&column name=ex, symbol=\"$ge$r$bx$n\", units=m, type=double, description=\"geometric horizontal emittance\" &end"},
    {"enx", "&column name=enx, symbol=\"$ge$r$bx,n$n\", type=double, units=m, description=\"normalized horizontal emittance\"  &end"},
    {"ecx", "&column name=ecx, symbol=\"$ge$r$bx,c$n\", units=m, type=double, description=\"geometric horizontal emittance less dispersive contributions\" &end"},
    {"ecnx", "&column name=ecnx, symbol=\"$ge$r$bx,cn$n\", type=double, units=m, description=\"normalized horizontal emittance less dispersive contributions\" &end"},
    {"ey", "&column name=ey, symbol=\"$ge$r$by$n\", units=m, type=double, description=\"geometric vertical emittance\" &end"},
    {"eny", "&column name=eny, symbol=\"$ge$r$by,n$n\", type=double, units=m, description=\"normalized vertical emittance\" &end"},
    {"ecy", "&column name=ecy, symbol=\"$ge$r$by,c$n\", units=m, type=double, description=\"geometric vertical emittance less dispersive contributions\" &end"},
    {"ecny", "&column name=ecny, symbol=\"$ge$r$by,cn$n\", type=double, units=m, description=\"normalized vertical emittance less dispersive contributions\" &end"},
    {"betaxBeam", "&column name=betaxBeam, symbol=\"$gb$r$bx,beam$n\", type=double, units=m, description=\"betax for the beam, excluding dispersive contributions\" &end"},
    {"alphaxBeam", "&column name=alphaxBeam, symbol=\"$ga$r$bx,beam$n\", type=double, description=\"alphax for the beam, excluding dispersive contributions\" &end"},
    {"betayBeam", "&column name=betayBeam, symbol=\"$gb$r$by,beam$n\", type=double, units=m, description=\"betay for the beam, excluding dispersive contributions\" &end"},
    {"alphayBeam", "&column name=alphayBeam, symbol=\"$ga$r$by,beam$n\", type=double, description=\"alphay for the beam, excluding dispersive contributions\" &end"},
    } ;

void SDDS_SigmaOutputSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row,
                           char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_SigmaOutputSetup");
#if USE_MPI  
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "sigma matrix", command_file, lattice_file, 
                            standard_parameter, STANDARD_PARAMETERS, element_column, ELEMENT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE);
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, NULL, NULL, NULL,
                            NULL, 0, sigma_matrix_column, SIGMA_MATRIX_COLUMNS,
                            caller, SDDS_EOS_COMPLETE);
    log_exit("SDDS_SigmaOutputSetup");
    }


#define WATCH_PARAMETER_MODE_COLUMNS 25
#define WATCH_CENTROID_MODE_COLUMNS 15
static SDDS_DEFINITION watch_parameter_mode_column[WATCH_PARAMETER_MODE_COLUMNS] = {
    {"Step", "&column name=Step, type=long &end"},
    {"Pass", "&column name=Pass, type=long &end"},
    {"Cx", "&column name=Cx, symbol=\"<x>\", units=m, type=double, description=\"x centroid\" &end"},
    {"Cxp", "&column name=Cxp, symbol=\"<x'>\", type=double, description=\"x' centroid\" &end"},
    {"Cy", "&column name=Cy, symbol=\"<y>\", units=m, type=double, description=\"y centroid\" &end"},
    {"Cyp", "&column name=Cyp, symbol=\"<y'>\", type=double, description=\"y' centroid\" &end"},
    {"Cs", "&column name=Cs, symbol=\"<s>\", units=m, type=double, description=\"mean distance traveled\" &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"<$gd$r>\", type=double, description=\"delta centroid\" &end"},
    {"Ct", "&column name=Ct, symbol=\"<t>\", units=s, type=double, description=\"mean time of flight\" &end"},
    {"dCt", "&column name=dCt, symbol=\"$gD$r<t>\", units=s, type=double, description=\"change in mean time of flight from last pass\" &end"},
    {"Particles", "&column name=Particles, description=\"Number of particles\", type=long, &end"},
    {"Transmission", "&column name=Transmission, description=Transmission, type=double &end"},
    {"pCentral", "&column name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", type=double, description=\"Reference beta*gamma\" &end"},
    {"pAverage", "&column name=pAverage, symbol=\"p$bave$n\", units=\"m$be$nc\", type=double, description=\"Mean beta*gamma\" &end"},
    {"KAverage",  "&column name=KAverage, symbol=\"K$bave$n\", units=MeV, type=double, description=\"Mean kinetic energy\" &end"},
    {"Sx",    "&column name=Sx, symbol=\"$gs$r$bx$n\", units=m, type=double, description=\"sqrt(<(x-<x>)^2>)\" &end"},
    {"Sxp",    "&column name=Sxp, symbol=\"$gs$r$bx'$n\", type=double, description=\"sqrt(<(x'-<x'>)^2>)\" &end"},
    {"Sy",    "&column name=Sy, symbol=\"$gs$r$by$n\", units=m, type=double, description=\"sqrt(<(y-<y>)^2>)\" &end"},
    {"Syp",    "&column name=Syp, symbol=\"$gs$r$by'$n\", type=double, description=\"sqrt(<(y'-<y'>)^2>)\" &end"},
    {"Ss",    "&column name=Ss, units=m, symbol=\"$gs$r$bs$n\", type=double, description=\"sqrt(<(s-<s>)^2>)\" &end"},
    {"Sdelta",    "&column name=Sdelta, symbol=\"$gs$bd$n$r\", type=double, description=\"sqrt(<(delta-<delta>)^2>)\" &end"},
    {"St", "&column name=St, symbol=\"$gs$r$bt$n\", units=s, type=double, description=\"sqrt(<(t-<t>)^2>)\" &end"},
    {"ex", "&column name=ex, symbol=\"$ge$r$bx$n\", units=m, type=double, description=\"geometric horizontal emittance\" &end"},
    {"ey", "&column name=ey, symbol=\"$ge$r$by$n\", units=m, type=double, description=\"geometric vertical emittance\" &end"},
    {"el", "&column name=el, symbol=\"$ge$r$bl$n\", units=s, type=double, description=\"longitudinal emittance\" &end"},
    } ;

#define WATCH_POINT_FFT_COLUMNS 4
static SDDS_DEFINITION watch_point_fft_column[WATCH_POINT_FFT_COLUMNS] = {
    {"f", "&column name=f, symbol=\"f/f$bsample$n\", description=\"Normalized frequency\", type=double &end"},
    {"Cx", "&column name=Cx, symbol=\"FFT <x>\", units=m, type=double, description=\"FFT of x centroid\" &end"},
    {"Cy", "&column name=Cy, symbol=\"FFT <y>\", units=m, type=double, description=\"FFT of y centroid\" &end"},
    {"Cdelta", "&column name=Cdelta, symbol=\"FFT <$gd$r>\", type=double, description=\"FFT of delta centroid\" &end"},
    } ;

void SDDS_WatchPointSetup(WATCH *watch, long mode, long lines_per_row,
                          char *command_file, char *lattice_file, char *caller, char *qualifier,
                          char *previousElementName)
{
  char s[100];
  SDDS_TABLE *SDDS_table;
  char *filename;
  long watch_mode, columns;
  
#if USE_MPI && !SDDS_MPI_IO 
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

  SDDS_table = &watch->SDDS_table;
  filename = watch->filename;
  watch_mode = watch->mode_code;

#if SDDS_MPI_IO
  SDDS_table->parallel_io = 0; /* By default, I/O will be done in serial */
#endif
  switch (watch_mode) {
    case WATCH_COORDINATES:
#if SDDS_MPI_IO
    SDDS_table->parallel_io = 1; 
    SDDS_MPI_Setup(SDDS_table, 1, n_processors, myid, MPI_COMM_WORLD, 0);
#endif
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point phase space",
                            command_file, lattice_file,
                            phase_space_parameter, PHASE_SPACE_PARAMETERS,
                            NULL, 0, caller, SDDS_EOS_NEWFILE);
#if SDDS_MPI_IO
    /* Open file here for parallel IO */
    if (!SDDS_MPI_File_Open(SDDS_table->MPI_dataset, SDDS_table->layout.filename, SDDS_MPI_WRITE_ONLY)) 
      SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDS_table->MPI_dataset->MPI_file);
#endif
    columns = 0;
    if (watch->xData) {
      if ((watch->xIndex[0]=SDDS_DefineColumn(SDDS_table, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (!watch->excludeSlopes &&
           (watch->xIndex[1]=SDDS_DefineColumn(SDDS_table, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0)) {
        fprintf(stdout, "Unable to define SDDS columns x and/or xp for file %s (%s)\n",
                filename, caller);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if (watch->yData) {
      if ((watch->yIndex[0]=SDDS_DefineColumn(SDDS_table, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (!watch->excludeSlopes &&
           (watch->yIndex[1]=SDDS_DefineColumn(SDDS_table, "yp",  NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0)) {
              fprintf(stdout, "Unable to define SDDS columns y and yp for file %s (%s)\n",
                filename, caller);
              fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if (watch->longitData) {
      if ((watch->longitIndex[0]
           =SDDS_DefineColumn(SDDS_table, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->longitIndex[1]
           =SDDS_DefineColumn(SDDS_table, "p", NULL, "m$be$nc", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
          (watch->longitIndex[2]
           =SDDS_DefineColumn(SDDS_table, "dt", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0) {
        fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
                filename, caller);
        fflush(stdout);
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
        exit(1);
      }
    }
    if ((watch->IDIndex = SDDS_DefineColumn(SDDS_table, "particleID", NULL, NULL, NULL, NULL, SDDS_LONG, 0))<0) {
      fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
              filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    
    if (!SDDS_DefineSimpleParameter(SDDS_table, "Pass", NULL, SDDS_LONG) ||
        !SDDS_DefineSimpleParameter(SDDS_table, "PassLength", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDS_table, "PassCentralTime", "s", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDS_table, "s", "m", SDDS_DOUBLE)) {
      fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (SDDS_DefineParameter(SDDS_table, "Description", NULL, NULL, NULL, "%s", SDDS_STRING, watch->label)<0) {
      fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (SDDS_DefineParameter(SDDS_table, "PreviousElementName", NULL, NULL, NULL, "%s", SDDS_STRING, 
                             previousElementName?previousElementName:"_BEG_")<0) {
      fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
#if SDDS_MPI_IO
    if (!SDDS_MPI_WriteLayout(SDDS_table))  
#else     
    if (!SDDS_WriteLayout(SDDS_table))
#endif
    {
      fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
#if SDDS_MPI_IO
    if ((watch->useDisconnect) && (!SDDS_DisconnectFile(SDDS_table)))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
#endif
    break;
  case WATCH_CENTROIDS:
    if (isMaster)
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point centroids",
                            command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS,
                            watch_parameter_mode_column, WATCH_CENTROID_MODE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    break;
  case WATCH_PARAMETERS:
    if (isMaster)
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "watch-point parameters", 
                            command_file, lattice_file,
                            standard_parameter, STANDARD_PARAMETERS, 
                            watch_parameter_mode_column, WATCH_PARAMETER_MODE_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    break;
  case WATCH_FFT:
    if (isMaster) {
    if (!qualifier)
      watch->window_code = FFT_HANNING;
    else if ((watch->window_code=match_string(qualifier, fft_window_name, N_FFT_WINDOWS, 0))<0) 
      bomb("invalid window mode for WATCH fft", NULL);
    sprintf(s,  "watch-point centroid FFT (%s window)", fft_window_name[watch->window_code]);
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, s,
                            command_file, lattice_file, standard_parameter, STANDARD_PARAMETERS, 
                            watch_point_fft_column, WATCH_POINT_FFT_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    }
    break;
  default:
    break;
  }
}

void SDDS_HistogramSetup(HISTOGRAM *histogram, long mode, long lines_per_row,
                          char *command_file, char *lattice_file, char *caller)
{
  SDDS_TABLE *SDDS_table=NULL;
  char *filename=NULL;
  long column, columns;

  if (isMaster) {
    SDDS_table = &histogram->SDDS_table;
    filename = histogram->filename;
    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, "histograms of phase-space coordinates",
			    command_file, lattice_file,
			    phase_space_parameter, PHASE_SPACE_PARAMETERS,
			    NULL, 0, caller, SDDS_EOS_NEWFILE);
  }
  for (column=0; column<7; column++)
    histogram->columnIndex[column][0] = 
      histogram->columnIndex[column][1] = -1;

  if (isMaster) {
    columns = 0;
    if (histogram->xData) {
      if ((histogram->columnIndex[0][0]
	   =SDDS_DefineColumn(SDDS_table, "x", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[0][1]
	   =SDDS_DefineColumn(SDDS_table, "xFrequency", NULL, 
			      histogram->normalize?"1/m":"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[1][0]
	   =SDDS_DefineColumn(SDDS_table, "xp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[1][1]
	   =SDDS_DefineColumn(SDDS_table, "xpFrequency", NULL, 
			      histogram->normalize?NULL:"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0) {
	fprintf(stdout, "Unable to define SDDS columns x and xp for file %s (%s)\n",
		filename, caller);
	fflush(stdout);
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }
      columns += 4;
    }
    if (histogram->yData) {
      if ((histogram->columnIndex[2][0]
	   =SDDS_DefineColumn(SDDS_table, "y", NULL, "m", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[2][1]
	   =SDDS_DefineColumn(SDDS_table, "yFrequency", NULL, 
			      histogram->normalize?"1/m":"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[3][0]
	   =SDDS_DefineColumn(SDDS_table, "yp", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[3][1]
	   =SDDS_DefineColumn(SDDS_table, "ypFrequency", NULL, 
			      histogram->normalize?NULL:"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0) {
	fprintf(stdout, "Unable to define SDDS columns y and yp for file %s (%s)\n",
		filename, caller);
	fflush(stdout);
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }
      columns += 4;
    }
    if (histogram->longitData) {
      if ((histogram->columnIndex[4][0]
	   =SDDS_DefineColumn(SDDS_table, "t", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[4][1]
	   =SDDS_DefineColumn(SDDS_table, "tFrequency", NULL, 
			      histogram->normalize?"1/s":"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[5][0]
	   =SDDS_DefineColumn(SDDS_table, "delta", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[5][1]
	   =SDDS_DefineColumn(SDDS_table, "deltaFrequency", NULL, 
			      histogram->normalize?NULL:"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[6][0]
	   =SDDS_DefineColumn(SDDS_table, "dt", NULL, "s", NULL, NULL, SDDS_DOUBLE, 0))<0 ||
	  (histogram->columnIndex[6][1]
	   =SDDS_DefineColumn(SDDS_table, "dtFrequency", NULL, 
			      histogram->normalize?"1/s":"Particles/Bin", NULL, NULL, SDDS_DOUBLE, 0))<0 ) {
	fprintf(stdout, "Unable to define SDDS columns t, dt, and p for file %s (%s)\n",
		filename, caller);
	fflush(stdout);
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	exit(1);
      }
      columns += 6;
    }
    if (!columns)
      bomb("no output selected for histogram", NULL);
    if (!SDDS_DefineSimpleParameter(SDDS_table, "Pass", NULL, SDDS_LONG) ||
	!SDDS_DefineSimpleParameter(SDDS_table, "PassLength", "m", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDS_table, "PassCentralTime", "s", SDDS_DOUBLE) ||
	!SDDS_DefineSimpleParameter(SDDS_table, "s", "m", SDDS_DOUBLE)) {
      fprintf(stdout, "Unable define SDDS parameter for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
    if (!SDDS_WriteLayout(SDDS_table)) {
      fprintf(stdout, "Unable to write SDDS layout for file %s (%s)\n", filename, caller);
      fflush(stdout);
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
      exit(1);
    }
  }
#if USE_MPI
  /* Broadcast the column index information to all the slaves */
  MPI_Bcast (&histogram->columnIndex[0][0], 14, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
}

void dump_watch_particles(WATCH *watch, long step, long pass, double **particle, long particles, 
                          double Po, double length, double charge, double z)
{
  long i, row;
  double p, t0, t;
#if SDDS_MPI_IO
  long total_row;
#endif

  log_entry("dump_watch_particles");
#if SDDS_MPI_IO
  if (isMaster && notSinglePart)   /* No particle will be dumped by master */
    particles = 0; 
  else
#endif
  {
  if (!watch->initialized)
    bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_particles)", NULL);
  if (!particle)
    bomb("NULL coordinate pointer passed to dump_watch_particles", NULL);
  if (watch->fraction>1)
    bomb("logic error--fraction>1 in dump_watch_particles", NULL);
  for (i=0; i<particles; i++)
    if (!particle[i]) {
      fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_particles)\n", i);
      fflush(stdout);
      abort();
    }
  }
  if (!SDDS_StartTable(&watch->SDDS_table, particles)) {
    SDDS_SetError("Problem starting SDDS table (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  row = 0;
  t0 = (pass-watch->passLast)*length*sqrt(Po*Po+1)/(c_mks*(Po+1e-32)) + watch->t0Last;
  watch->t0Last = t0;
  watch->passLast = pass;
#if SDDS_MPI_IO
  if ((isSlave&&notSinglePart)||(!notSinglePart&&isMaster))
#endif
  for (i=0; i<particles; i++) {
    if (watch->fraction==1 || random_2(0)<watch->fraction) {
      if (watch->xData && 
          !SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                             watch->xIndex[0], particle[i][0], 
                             watch->excludeSlopes?-1:watch->xIndex[1], particle[i][1],
                             -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (watch->yData &&
          !SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                             watch->yIndex[0], particle[i][2], 
                             watch->excludeSlopes?-1:watch->yIndex[1], particle[i][3],
                             -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (watch->longitData) {
        p = Po*(1+particle[i][5]);
        t = particle[i][4]/(c_mks*p/sqrt(sqr(p)+1));
        if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row,
                               watch->longitIndex[0], t,
                               watch->longitIndex[1], p,
                               watch->longitIndex[2], t-t0,
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
      if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, row++,
                             watch->IDIndex, (long)particle[i][6], -1)) {
        SDDS_SetError("Problem setting SDDS row values (dump_watch_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
  }
#if SDDS_MPI_IO
  if (USE_MPI&&notSinglePart){

    MPI_Allreduce (&row, &total_row, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (isMaster)
      row = total_row;
  }
  if (total_row)
#endif
  if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "Step", step, "Pass", pass, "Particles", row, "pCentral", Po,
                          "PassLength", length, 
                          "Charge", charge,
                          "PassCentralTime", t0, "s", z,
                          NULL)) {
    SDDS_SetError("Problem setting SDDS parameters (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
#if SDDS_MPI_IO
  if ((watch->useDisconnect) && (!SDDS_ReconnectFile(&watch->SDDS_table)))
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (total_row)
  if (!SDDS_MPI_WriteTable(&watch->SDDS_table))
#else
  if (SDDS_IsDisconnected(&(watch->SDDS_table)) &&  (!SDDS_ReconnectFile(&watch->SDDS_table)))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_WriteTable(&watch->SDDS_table))
#endif
  {
    SDDS_SetError("Problem writing SDDS table (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  } 
  if (!inhibitFileSync)
    SDDS_DoFSync(&watch->SDDS_table);
  if ((watch->useDisconnect) && (!SDDS_DisconnectFile(&watch->SDDS_table)))
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  if (!SDDS_ShortenTable(&watch->SDDS_table, 1))
  {
    SDDS_SetError("Problem shortening SDDS table (dump_watch_particles)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  log_exit("dump_watch_particles");
}

double tmp_safe_sqrt;
#define SAFE_SQRT(x) ((tmp_safe_sqrt=(x))<0?(double)0.0:sqrt(tmp_safe_sqrt))

void dump_watch_parameters(WATCH *watch, long step, long pass, long n_passes, double **particle, 
                           long particles, long original_particles,  double Po, 
                           double revolutionLength)
{
    long sample, i, j, watchStartPass=watch->start_pass;
    double tc, tc0, p_sum, gamma_sum, sum, p=0.0;
    
#if USE_MPI  
    long particles_total;
    double emittance_l;
#ifdef  USE_MPE /* use the MPE library */
  int event1a, event1b;
  event1a = MPE_Log_get_event_number();
  event1b = MPE_Log_get_event_number();
  if(isMaster) {
    MPE_Describe_state(event1a, event1b, "End Watch", "pink");
  }
#endif

    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid==0)
      particles = 0;

    MPI_Allreduce(&particles, &particles_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (isMaster) {
      log_entry("dump_watch_parameters");
 
      if (!watch->initialized)
        bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_parameters)", NULL);
    
      if (watch->fraction>1)
        bomb("logic error--fraction>1 in dump_watch_parameters", NULL);
    }
    if (isSlave) {
      if (!particle)
        bomb("NULL coordinate pointer passed to dump_watch_parameters", NULL);
    }

    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_parameters)\n", i);
            fflush(stdout);
            abort();
            }
    if (isMaster) 
    if ((watchStartPass==pass) && !SDDS_StartTable(&watch->SDDS_table, (n_passes-watchStartPass)/watch->interval+1)) {
        SDDS_SetError("Problem starting SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }

    sample = (pass-watchStartPass)/watch->interval;
    tc0 = (pass-watch->passLast)*revolutionLength*sqrt(Po*Po+1)/(c_mks*(Po+1e-32)) + watch->t0Last;
    watch->t0Last = tc0;
    watch->passLast = pass;
    /* compute centroids, sigmas, and emittances for x, y, and s */
    for (i=0; i<6; i+=2 ) {
      double x, xp, cx, cxp, sum_x = 0, sum_xp = 0, sum_x2, sum_xxp, sum_xp2;
      for (j=0; j<particles; j++) {
        sum_x   += particle[j][i];
        sum_xp  += particle[j][i+1];
      }
#if SDDS_MPI_IO
      if (USE_MPI) {
	double sum_x_total, sum_xp_total;
 
        if (particles_total) {
	  MPI_Allreduce(&sum_x, &sum_x_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&sum_xp, &sum_xp_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  sum_x = sum_x_total;
	  sum_xp = sum_xp_total;
	  cx = sum_x/particles_total;
	  cxp = sum_xp/particles_total;
	}
	else
	  cx = cxp = 0;
      }
#else
      if (particles) {
        cx  = sum_x/particles;
        cxp = sum_xp/particles;
      }
      else 
        cx = cxp = 0;
#endif
      if (isMaster) 
	if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
			       i==0?"Cx" :(i==2?"Cy" :"Cs" ), cx,
			       i==0?"Cxp":(i==2?"Cyp":"Cdelta"), cxp, NULL)) {
	  SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      if (watch->mode_code==WATCH_PARAMETERS) {
        for (j=sum_x2=sum_xp2=sum_xxp=0; j<particles; j++) {
          x  = particle[j][i]   - cx;
          xp = particle[j][i+1] - cxp;
          sum_x2  += sqr(x);
          sum_xxp += x*xp;
          sum_xp2 += sqr(xp);
        }
#if SDDS_MPI_IO
	if (USE_MPI) {
          double sum[3], sum_total[3];
          sum[0] = sum_x2; sum[1] = sum_xxp; sum[2] = sum_xp2;
	  MPI_Allreduce(sum, sum_total, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          sum_x2 = sum_total[0];
          sum_xxp = sum_total[1];
          sum_xp2 = sum_total[2];
	}
        if (isMaster) {
	  if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
				 i==0?"Sx" :(i==2?"Sy" :"Ss" ), particles_total?sqrt(sum_x2/particles_total):0.0,
				 i==0?"Sxp":(i==2?"Syp":"Sdelta"), particles_total?sqrt(sum_xp2/particles_total):0.0, NULL)) {
	    SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	  if (i!=4)
	    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
				   i==0?"ex":"ey", SAFE_SQRT(sum_x2*sum_xp2-sqr(sum_xxp))/particles_total, NULL)) {
	      SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	    }
	}
#else 
	if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
			       i==0?"Sx" :(i==2?"Sy" :"Ss" ), particles?sqrt(sum_x2/particles):0.0,
			       i==0?"Sxp":(i==2?"Syp":"Sdelta"), particles?sqrt(sum_xp2/particles):0.0, NULL)) {
	  SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
	if (i!=4)
	  if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
				 i==0?"ex":"ey", SAFE_SQRT(sum_x2*sum_xp2-sqr(sum_xxp))/particles, NULL)) {
	    SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
#endif
      }
    }

#if SDDS_MPI_IO 
    emittance_l = rms_longitudinal_emittance_p(particle, particles, Po);  
    if (isMaster) {
      if (watch->mode_code==WATCH_PARAMETERS)
	if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
			       "el", emittance_l, NULL)) {
	  SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
    }
#else
    if (watch->mode_code==WATCH_PARAMETERS)
      if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
			     "el", rms_longitudinal_emittance(particle, particles, Po), NULL)) {
	SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
#endif
    /* time centroid and sigma */
    for (i=p_sum=gamma_sum=sum=0; i<particles; i++) {
      p = Po*(1+particle[i][5]);
      p_sum     += p;
      gamma_sum += sqrt(sqr(p)+1);
      sum += particle[i][4]/(p/sqrt(sqr(p)+1)*c_mks) ;
    }
#if SDDS_MPI_IO
    if (USE_MPI) {
      double p_sum_total, gamma_sum_total, sum_total;
      MPI_Allreduce (&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce (&p_sum, &p_sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce (&gamma_sum, &gamma_sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
      if (particles_total)
	tc = sum_total/particles_total;
      else
	tc = 0;
      if (isMaster) {
	if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
			       "Ct", tc, "dCt", tc-tc0,
			       "pAverage", p_sum_total/particles_total, "pCentral", Po, 
			       "KAverage", (gamma_sum_total/particles_total-1)*me_mev, NULL)) {
	  SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      }
    }
#else
    if (particles)
      tc = sum/particles;
    else
      tc = 0;
 
    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                           "Ct", tc, "dCt", tc-tc0,
                           "pAverage", p_sum/particles, "pCentral", Po, 
                           "KAverage", (gamma_sum/particles-1)*particleMassMV, NULL)) {
      SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
#endif
    if (watch->mode_code==WATCH_PARAMETERS) {
      for (i=sum=0; i<particles; i++) {
        p = Po*(1+particle[i][5]);
        sum += sqr(particle[i][4]/(p/sqrt(sqr(p)+1)*c_mks) - tc);
      }
#if SDDS_MPI_IO
      if (USE_MPI) {
	double sum_total;
	MPI_Allreduce(&sum, &sum_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (isMaster)
	  if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
				 "St", particles_total?sqrt(sum_total/particles_total):0.0, NULL)) {
	    SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	if (isMaster)
          particles = particles_total; /* make sure the partiles will not be used in a loop */
      }
#else
      if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,
                             "St", particles?sqrt(sum/particles):0.0, NULL)) {
        SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
#endif
    }

    if (isMaster) {
      /* number of particles */
      if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, sample,			     
#if SDDS_MPI_IO
			     "Particles", particles_total,
			     "Transmission", (original_particles?((double)particles_total)/original_particles:(double)0.0),
#else 
			     "Particles", particles,
			     "Transmission", (original_particles?((double)particles)/original_particles:(double)0.0),
#endif
			     "Pass", pass, 
			     "Step", step, NULL)) {
	SDDS_SetError("Problem setting row values for SDDS table (dump_watch_parameters)");
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
			      "Step", step, NULL)) {
	SDDS_SetError("Problem setting parameter values for SDDS table (dump_watch_parameters)");
	SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }


      if (sample==(n_passes-1)/watch->interval) {
	if (watch->flushInterval>0) {
	  if (sample!=watch->flushSample && !SDDS_UpdatePage(&watch->SDDS_table, 0)) {
	    SDDS_SetError("Problem writing data for SDDS table (dump_watch_parameters)");
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	  }
	}
	else if (!SDDS_WriteTable(&watch->SDDS_table)) {
	  SDDS_SetError("Problem writing data for SDDS table (dump_watch_parameters)");
	  SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
	}
      } else {
      if (watch->flushInterval>0 && sample%watch->flushInterval==0 &&
          !SDDS_UpdatePage(&watch->SDDS_table, 0)) {
        SDDS_SetError("Problem flushing data for SDDS table (dump_watch_parameters)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
      watch->flushSample = sample;
    }
#ifdef  USE_MPE
	      MPE_Log_event(event1a, 0, "start watch"); /* record time spent on I/O operations */
#endif
      if (!inhibitFileSync)
	SDDS_DoFSync(&watch->SDDS_table);
#ifdef  USE_MPE
	      MPE_Log_event(event1b, 0, "end watch");
#endif   
      log_exit("dump_watch_parameters");
    }
}


void dump_watch_FFT(WATCH *watch, long step, long pass, long n_passes, double **particle, long particles,
                           long original_particles,  double Po)
{
    long sample_index, i, samples;
    double sum_x, sum_y, sum_dp, *sample[4];

#if USE_MPI  
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

    log_entry("dump_watch_FFT");

    samples = n_passes/watch->interval;
    if (!power_of_2(samples)) 
        bomb("number of samples must be a power of 2 for FFT (dump_watch_FFT)", NULL);
    if (!watch->initialized)
        bomb("uninitialized watch-point (coordinate mode) encountered (dump_watch_FFT)", NULL);
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_watch_FFT", NULL);
    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_watch_FFT)\n", i);
            fflush(stdout);
            abort();
            }

    if (pass==0 && !SDDS_StartTable(&watch->SDDS_table, 2*samples)) {
        SDDS_SetError("Problem starting SDDS table (dump_watch_FFT)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    sample_index = pass/watch->interval;
    for (i=sum_x=sum_y=sum_dp=0; i<particles; i++) {
        sum_x  += particle[i][0];
        sum_y  += particle[i][2];
        sum_dp += particle[i][5];
        }
    if (particles) {
        sum_x  /= particles;
        sum_y  /= particles;
        sum_dp /= particles;
        }
    /* SDDS is used to store the centroid vs pass data prior to FFT */
    if (!SDDS_SetRowValues(&watch->SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, sample_index,
                           0, ((double)sample_index)/samples, 1, sum_x, 2, sum_y, 3, sum_dp, -1)) {
        SDDS_SetError("Problem setting row values for SDDS table (dump_watch_FFT)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (sample_index==samples-1) {
        /* sample/samples stored temporarily in frequency column: */
        if (!(sample[0]=(double*)SDDS_GetColumn(&watch->SDDS_table, "f"))  || 
            !(sample[1]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cx")) ||
            !(sample[2]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cy")) ||
            !(sample[3]=(double*)SDDS_GetColumn(&watch->SDDS_table, "Cdelta"))) {
            SDDS_SetError("Problem retrieving sample values from SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        do_watch_FFT(sample, samples, 1, watch->window_code);
        do_watch_FFT(sample, samples, 2, watch->window_code);
        do_watch_FFT(sample, samples, 3, watch->window_code);
        if (!SDDS_EraseData(&watch->SDDS_table) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[0], samples/2, 0) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[1], samples/2, 1) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[2], samples/2, 2) ||
            !SDDS_SetColumn(&watch->SDDS_table, SDDS_SET_BY_INDEX, sample[3], samples/2, 3) ) {
            SDDS_SetError("Problem setting columns of SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        for (i=0; i<4; i++)
            free(sample[i]);
        if (!SDDS_SetParameters(&watch->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                                "Step", step, NULL)) {
            SDDS_SetError("Problem setting parameter values for SDDS table (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!SDDS_WriteTable(&watch->SDDS_table)) {
            SDDS_SetError("Problem writing data to SDDS file (dump_watch_FFT)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        if (!inhibitFileSync)
          SDDS_DoFSync(&watch->SDDS_table);
        }
    log_exit("dump_watch_FFT");
    }
 
void do_watch_FFT(double **data, long n_data, long slot, long window_code) 
{
    double *real_imag, r1, r2;
    long i;

#if USE_MPI  
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

    log_entry("do_watch_FFT");
    if (!data[slot])
        bomb("data pointer for specified slot is NULL (do_watch_FFT)", NULL);

    /* load into temporary array */
    real_imag = tmalloc(sizeof(*real_imag)*(n_data+2));
    for (i=0; i<n_data; i++)
        real_imag[i] = data[slot][i];

    /* apply the window */
    switch (window_code) {
      case FFT_HANNING:
        r1 = PIx2/(n_data-1);
        for (i=0; i<n_data; i++) 
            real_imag[i] *= (1-cos(i*r1))/2;
        break;
      case FFT_PARZEN:
        r1 = (n_data-1)/2.0;
        for (i=0; i<n_data; i++) 
            real_imag[i] *= 1-FABS((i-r1)/r1);
        break;
      case FFT_WELCH:
        r1 = (n_data-1)/2.0;
        r2 = sqr((n_data+1)/2.0);
        for (i=0; i<n_data; i++) 
            real_imag[i] *= 1 - sqr(i-r1)/r2;
        break;
      case FFT_UNIFORM:
        break;
      default:
        break;
        }

    /* do the FFT */
    if (!realFFT2(real_imag, real_imag, n_data, 0)) {
        fprintf(stdout, "error: FFT fails for slot %ld of watch point\n", slot);
        fflush(stdout);
        abort();
        }

    /* put the amplitudes back into the data array.  Multiply the non-DC
     * terms by 2 to account for amplitude in negative frequencies. 
     */
    data[slot][0] = fabs(real_imag[0]);
    for (i=1; i<n_data/2; i++) 
        data[slot][i] = sqrt(sqr(real_imag[2*i]) + (i==0?0:sqr(real_imag[2*i+1])))*2;
    
    free(real_imag); 
    log_exit("do_watch_FFT");
    }


void dump_particle_histogram(HISTOGRAM *histogram, long step, long pass, double **particle, long particles, 
                          double Po, double length, double charge, double z)
{
  long icoord, ipart, ibin;
  double p, t0;
  static long maxBins = 0, maxParticles = 0;
  static double *frequency=NULL, *coordinate=NULL, *histData=NULL;
  double center, range, lower, upper;
  
#if USE_MPI  
  static double *buffer=NULL;
  long particles_total;
#endif

  log_entry("dump_particle_histogram");
  if (!histogram->initialized)
    bomb("uninitialized histogram encountered (dump_particle_histogram)", NULL);
  if (isSlave)
    if (!particle)
      bomb("NULL coordinate pointer passed to dump_particle_histogram", NULL);
  /*
  for (ipart=0; ipart<particles; ipart++)
    if (!particle[ipart]) {
      fprintf(stdout, "error: coordinate slot %ld is NULL (dump_particle_histogram)\n", ipart);
      fflush(stdout);
      abort();
    }
  */
  if (isMaster) {
    if (!SDDS_StartTable(&histogram->SDDS_table, histogram->bins)) {
      SDDS_SetError("Problem starting SDDS table (dump_particle_histogram)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (histogram->bins>maxBins) {
    maxBins = histogram->bins;
    if (!(frequency=SDDS_Realloc(frequency, sizeof(*frequency)*maxBins)) ||
        !(coordinate=SDDS_Realloc(coordinate, sizeof(*coordinate)*maxBins))) {
      SDDS_Bomb("Memory allocation failure (dump_particle_histogram)");
    }
#if USE_MPI
    if (!(buffer=SDDS_Realloc(buffer, sizeof(*buffer)*maxBins))) {
      SDDS_Bomb("Memory allocation failure (dump_particle_histogram)");
    }
    MPI_Reduce (&particles, &particles_total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  }
  if (isSlave) {
    if (particles>maxParticles) {
      maxParticles = particles;
      if (!(histData = SDDS_Realloc(histData, sizeof(*histData)*maxParticles)))
	SDDS_Bomb("Memory allocation failure (dump_particle_histogram)");
    }
  }
  
  t0 = pass*length*sqrt(Po*Po+1)/(c_mks*(Po+1e-32));  
  for (icoord=0; icoord<7; icoord++) {
    upper = -(lower=DBL_MAX);
    if (histogram->columnIndex[icoord][0]==-1)
      continue;
    if (isSlave) {
      if (icoord==4) {
	/* t */
	for (ipart=0; ipart<particles; ipart++) {
	  p = Po*(1+particle[ipart][5]);
	  histData[ipart] = particle[ipart][4]/(c_mks*p/sqrt(sqr(p)+1));
	}
      } else if (icoord==6) {
	/* dt (could reuse computations for t, but would need to store another array) */
	for (ipart=0; ipart<particles; ipart++) {
	  p = Po*(1+particle[ipart][5]);
	  histData[ipart] = particle[ipart][4]/(c_mks*p/sqrt(sqr(p)+1)) - t0;
	}
      } else {
	for (ipart=0; ipart<particles; ipart++)
	  histData[ipart] = particle[ipart][icoord];
      }
      find_min_max(&lower, &upper, histData, particles);
      if (lower==upper) {
	lower -= 1e-16;
	upper += 1e-16;
      }
    }
#if USE_MPI
    if (isMaster)    /* The particles are on slave */
      particles = 0;    
    find_global_min_max(&lower, &upper, particles, MPI_COMM_WORLD); 
#endif
    if (histogram->count==0 || !histogram->fixedBinSize) {
      range = histogram->binSizeFactor*(upper-lower);
      histogram->binSize[icoord] = range/histogram->bins;
    } else
      range = histogram->binSize[icoord]*histogram->bins;
    center = (lower+upper)/2;
    range *= (1+1e-12); /* Rounding error might cause particle lost in the first and last bins */ 
    lower = center-range/2;
    upper = center+range/2;
    for (ibin=0; ibin<histogram->bins; ibin++) 
      coordinate[ibin] = (ibin+0.5)*histogram->binSize[icoord] + lower;
    make_histogram(frequency, histogram->bins, lower, upper, histData, particles, 1);
#if USE_MPI
    if (USE_MPI) /* This will update the number of particles locally on master if there are lost on slaves */
      MPI_Reduce (&particles, &particles_total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (isMaster)
      particles = particles_total;
    MPI_Allreduce(frequency, buffer, histogram->bins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(frequency, buffer, sizeof(*buffer)*histogram->bins);    
    if (isMaster && particles_total) 
#endif
    {
      if (histogram->normalize)
	for (ibin=0; ibin<histogram->bins; ibin++)
	  frequency[ibin] /= particles*(upper-lower)/histogram->bins;
      if (!SDDS_SetColumn(&histogram->SDDS_table, SDDS_SET_BY_INDEX, 
			  coordinate, histogram->bins, histogram->columnIndex[icoord][0]) ||
	  !SDDS_SetColumn(&histogram->SDDS_table, SDDS_SET_BY_INDEX,
			  frequency, histogram->bins, histogram->columnIndex[icoord][1])) {
	SDDS_Bomb("Problem setting column values (dump_particle_histogram)");
      }
    }
  }
#if SDDS_MPI_IO
  if (isMaster && particles_total)
#endif
  {
    if (!SDDS_SetParameters(&histogram->SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
			    "Step", step, "Pass", pass, "Particles", particles, "pCentral", Po,
			    "PassLength", length, "Charge", charge,
			    "PassCentralTime", t0, "s", z,
			    NULL)) {
      SDDS_SetError("Problem setting SDDS parameters (dump_particle_histogram)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    if (!SDDS_WriteTable(&histogram->SDDS_table)) {
      SDDS_SetError("Problem writing SDDS table (dump_particle_histogram)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } 
    if (!inhibitFileSync)
      SDDS_DoFSync(&histogram->SDDS_table);
  }
  histogram->count++;
}


void dump_phase_space(SDDS_TABLE *SDDS_table, double **particle, long particles, long step, double Po,
        double charge)
{
    long i;
    double p;

#if SDDS_MPI_IO
    long total_particles;

    /* Open file here for parallel IO */
    if (!SDDS_table->layout.layout_written) { /* Check if the file has been opened already */
      if (!SDDS_MPI_File_Open(SDDS_table->MPI_dataset, SDDS_table->layout.filename, SDDS_MPI_WRITE_ONLY)) 
	SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDS_table->MPI_dataset->MPI_file);
      if (!SDDS_MPI_WriteLayout(SDDS_table))  
	SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDS_table->MPI_dataset->MPI_file);
    }
#endif

    log_entry("dump_phase_space");
    if (isSlave) {
      if (!particle)
        bomb("NULL coordinate pointer passed to dump_phase_space", NULL);
    
      for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_phase_space)\n", i);
            fflush(stdout);
            abort();
	}
    }
    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
#if SDDS_MPI_IO
    if (notSinglePart) {
      MPI_Reduce (&particles, &total_particles, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      if (isMaster)
	particles = total_particles;
    }
#endif  
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                            "Step", step, "pCentral", Po, "Particles", particles,
                            "Charge", charge, NULL)) {
        SDDS_SetError("Problem setting parameter values for SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
#if SDDS_MPI_IO
    if ((notSinglePart&&isSlave)||(!notSinglePart&&isMaster))
#endif
    for (i=0; i<particles; i++) {
        p = Po*(1+particle[i][5]);
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particle[i][0], 1, particle[i][1], 2, particle[i][2], 3, particle[i][3],
                               4, particle[i][4]/(c_mks*p/sqrt(sqr(p)+1)), 5, p, 
                               6, (long)particle[i][6], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_phase_space)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step", step, NULL)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
#if SDDS_MPI_IO
    if (!SDDS_MPI_WriteTable(SDDS_table) || !SDDS_ShortenTable(SDDS_table, 1))
#else
    if (!SDDS_WriteTable(SDDS_table) || !SDDS_ShortenTable(SDDS_table, 1))
#endif
      {
        SDDS_SetError("Problem writing SDDS table (dump_phase_space)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } 
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);

    log_exit("dump_phase_space");
    }

#ifdef SORT
static int comp_IDs1(const void **coord1, const void **coord2)
  {
    long v1, v2;
    v1 = (*((long**)coord1))[1];
    v2 = (*((long**)coord2))[1];
    if (v1 < v2) 
      return -1;
    else if (v1 > v2)
      return  1;
    else 
      return 0;
  }
#endif

void dump_lost_particles(SDDS_TABLE *SDDS_table, double **particle, long *lostOnPass, 
			 long particles, long step)
{
    long i;
#ifdef SORT   /* sort for comparing the serial and parallel versions */
    long **tmp, j;
#endif
#if SDDS_MPI_IO

    /* Open file here for parallel IO */
    if (!SDDS_table->layout.layout_written) { /* Check if the file has been opened already */
      if (!SDDS_MPI_File_Open(SDDS_table->MPI_dataset, SDDS_table->layout.filename, SDDS_MPI_WRITE_ONLY)) 
	SDDS_MPI_BOMB("SDDS_MPI_File_Open failed.", &SDDS_table->MPI_dataset->MPI_file);
      if (!SDDS_MPI_WriteLayout(SDDS_table))  
	SDDS_MPI_BOMB("SDDS_MPI_WriteLayout failed.", &SDDS_table->MPI_dataset->MPI_file);
    }
#endif

    log_entry("dump_lost_particles");
    if (isSlave) {
      if (!particle)
        bomb("NULL coordinate pointer passed to dump_lost_particles", NULL);
    }

#ifdef SORT   /* sort for comparing the serial and parallel versions */
    if (SORT && particles) {
        tmp = tmalloc(sizeof(*tmp) * particles);
        for (j=0; j<particles; j++) {
          tmp[j] = tmalloc(sizeof(*(tmp[j])) * 2);
	  tmp[j][0] = lostOnPass[j]; 
          tmp[j][1] = particle[j][6];
	}    
        /* This call looks strange, but works because the particle data is contiguous.  The sort
         * actually moves the data around without changing the pointers.
         */
	qsort(particle[0], particles, COORDINATES_PER_PARTICLE*sizeof(double), comp_IDs); 
        /* This call is more standard.  We sort the pointers and leave the data in place. */
	qsort(tmp, particles, sizeof(long*), comp_IDs1);
	for (j=0; j<particles; j++) {
	  lostOnPass[j] = tmp[j][0];
	}
      }
    if (SORT && particles) {
      for (j=0; j<particles; j++) {
	free(tmp[j]);
      }
      free(tmp);
    }
#endif

    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_lost_particles)\n", i);
            fflush(stdout);
            abort();
            }

    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    for (i=0; i<particles; i++) {
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particle[i][0], 1, particle[i][1], 2, particle[i][2], 3, particle[i][3],
                               4, particle[i][4], 5, particle[i][5],
                               6, (long)particle[i][6], 7, lostOnPass[i], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_lost_particles)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Step", step, NULL)) {
        SDDS_SetError("Problem setting SDDS parameters (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
#if SDDS_MPI_IO
    if (!SDDS_MPI_WriteTable(SDDS_table))
#else
    if (!SDDS_WriteTable(SDDS_table))
#endif
      {
        SDDS_SetError("Problem writing SDDS table (dump_lost_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      } 
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);

    log_exit("dump_lost_particles");
    
}

void dump_centroid(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums, LINE_LIST *beamline, long n_elements, long step,
                          double p_central)
{
    long i, j;
    BEAM_SUMS *beam;
    ELEMENT_LIST *eptr;
    double *cent;
    char *name, *type_name;
    long s_index, Cx_index=0, occurence;

#if USE_MPI  
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

    log_entry("dump_centroid");

    if (!SDDS_table)
        bomb("SDDS_TABLE pointer is NULL (dump_centroid)", NULL);
    if (!sums)
        bomb("BEAM_SUMS pointer is NULL (dump_centroid)", NULL);
    if (!beamline)
        bomb("LINE_LIST pointer is NULL (dump_centroid)", NULL);
    if ((s_index=SDDS_GetColumnIndex(SDDS_table, "s"))<0 ||
        (Cx_index=SDDS_GetColumnIndex(SDDS_table, "Cx"))<0) {
        SDDS_SetError("Problem getting index of SDDS columns (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }

    if (!SDDS_StartTable(SDDS_table, n_elements+1)) {
        SDDS_SetError("Problem starting SDDS table (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                            "Step", step, NULL)) {
        SDDS_SetError("Problem setting parameter values for SDDS table (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
 
    if (!(beamline->flags&BEAMLINE_CONCAT_DONE))
        eptr = &(beamline->elem);
    else 
        eptr = &(beamline->ecat);
    name = "_BEG_";
    type_name = "MARK";
    occurence = 1;
    cent = tmalloc(sizeof(*cent)*6);
    for (i=0; i<n_elements; i++) {
        if (!eptr) {
            fprintf(stdout, "element pointer is NULL, i=%ld (dump_centroid)", i);
            fflush(stdout);
            abort();
            }
        beam = sums+i;
        if (!beam->centroid) {
            fprintf(stdout, "beam->centroid is NULL, i=%ld (dump_centroid)\n", i);
            fflush(stdout);
            abort();
            }
        if (beam->n_part)
            for (j=0; j<6; j++) {
                cent[j] = beam->centroid[j];
                if (isnan(cent[j]) || isinf(cent[j]))
                    cent[j] = DBL_MAX;
                }
        else
            for (j=0; j<6; j++)
                cent[j] = 0;
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               Cx_index  , cent[0], Cx_index+1, cent[1], Cx_index+2, cent[2],
                               Cx_index+3, cent[3], Cx_index+4, cent[4], Cx_index+5, cent[5],
                               Cx_index+6, beam->n_part, Cx_index+7, beam->p0, 
                               s_index, beam->z, s_index+1, name, s_index+2, occurence,
                               s_index+3, type_name, -1)) {
            SDDS_SetError("Problem setting row values for SDDS table (dump_centroid)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        
        name = eptr->name;
        type_name = entity_name[eptr->type];
        occurence = eptr->occurence;
        if (eptr->succ)
            eptr = eptr->succ;
        else if (beamline->elem_recirc)
            eptr = beamline->elem_recirc;
        else
            eptr = &(beamline->elem);
        }

    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Unable to write centroid data (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);
    if (!SDDS_EraseData(SDDS_table)) {
        SDDS_SetError("Unable to erase centroid data (dump_centroid)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    free(cent);
    log_exit("dump_centroid");
    }

void dump_sigma(SDDS_TABLE *SDDS_table, BEAM_SUMS *sums, LINE_LIST *beamline, long n_elements, long step,
                double p_central)
{
  long i, j, ie, offset, plane, index;
  BEAM_SUMS *beam;
  ELEMENT_LIST *eptr;
  double emit, emitNorm, emitc, emitcNorm, beta, alpha;
  char *name, *type_name;
  long s_index=0, ma1_index=0, min1_index=0, max1_index=0, Sx_index=0, occurence, ex_index=0, betax_index=0;
  long sNIndex[6]={0,0,0,0,0,0};

#if USE_MPI  
    if (myid<0)
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid!=0)
      return;
#endif

  if (!SDDS_table)
    bomb("SDDS_TABLE pointer is NULL (dump_sigma)", NULL);
  if (!sums)
    bomb("BEAM_SUMS pointer is NULL (dump_centroid)", NULL);
  if (!beamline)
    bomb("LINE_LIST pointer is NULL (dump_centroid)", NULL);
  if ((s_index=SDDS_GetColumnIndex(SDDS_table, "s"))<0 ||
      (sNIndex[0]=SDDS_GetColumnIndex(SDDS_table, "s1"))<0 ||
      (sNIndex[1]=SDDS_GetColumnIndex(SDDS_table, "s2"))<0 ||
      (sNIndex[2]=SDDS_GetColumnIndex(SDDS_table, "s3"))<0 ||
      (sNIndex[3]=SDDS_GetColumnIndex(SDDS_table, "s4"))<0 ||
      (sNIndex[4]=SDDS_GetColumnIndex(SDDS_table, "s5"))<0 ||
      (sNIndex[5]=SDDS_GetColumnIndex(SDDS_table, "s6"))<0 ||
      (ma1_index=SDDS_GetColumnIndex(SDDS_table, "ma1"))<0 ||
      (min1_index=SDDS_GetColumnIndex(SDDS_table, "minimum1"))<0 ||
      (max1_index=SDDS_GetColumnIndex(SDDS_table, "maximum1"))<0 ||
      (Sx_index=SDDS_GetColumnIndex(SDDS_table, "Sx"))<0 ||
      (ex_index=SDDS_GetColumnIndex(SDDS_table, "ex"))<0 ||
      (betax_index=SDDS_GetColumnIndex(SDDS_table, "betaxBeam"))<0) {
    SDDS_SetError("Problem getting index of SDDS columns (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!SDDS_StartTable(SDDS_table, n_elements+1)) {
    SDDS_SetError("Problem starting SDDS table (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
                          "Step", step, NULL)) {
    SDDS_SetError("Problem setting parameter values for SDDS table (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  if (!(beamline->flags&BEAMLINE_CONCAT_DONE))
    eptr = &(beamline->elem);
  else 
    eptr = &(beamline->ecat);
  name = "_BEG_";
  type_name = "MARK";
  occurence = 1;
  for (ie=0; ie<n_elements; ie++) {
    beam  = sums+ie;
    if (!beam->sigma) {
      fprintf(stdout, "error: sigma element of BEAM_SUMS is undefined for element %ld (%s)\n",
              ie, name);
      fflush(stdout);
      exit(1);
    }
    if (beam->n_part) {
      for (i=0; i<6; i++) {
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                               Sx_index+i, sqrt(beam->sigma[i][i]),
                               sNIndex[i], sqrt(beam->sigma[i][i]),
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        offset = 1;
        for (j=i+1; j<6; j++, offset++) {
          if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                                 sNIndex[i]+offset, beam->sigma[i][j], -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
          }
        }
      }
      
      /* Set values for maximum amplitudes of coordinates */
      for (i=0; i<6; i++) {
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, ma1_index+i, beam->maxabs[i], -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, min1_index+i, beam->min[i], -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, max1_index+i, beam->max[i], -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
      for (plane=0; plane<=2; plane+=2) {
        /* emittance and beam twiss parameters */
        computeEmitTwissFromSigmaMatrix(&emit, &emitc, &beta, &alpha, beam->sigma, plane);
        emitNorm = emit*beam->p0*(1+beam->centroid[5]);
        emitcNorm = emitc*beam->p0*(1+beam->centroid[5]);
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                               ex_index+2*plane, emit, 
                               ex_index+1+2*plane, emitNorm,
                               ex_index+2+2*plane, emitc, 
                               ex_index+3+2*plane, emitcNorm,
                               betax_index+plane+0, beta,
                               betax_index+plane+1, alpha,
                               -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      }
    }
    else {
      for (index=sNIndex[0]; index<sNIndex[0]+21; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=ma1_index; index<ma1_index+4; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=Sx_index; index<Sx_index+6; index++) 
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
      for (index=ex_index; index<ex_index+4; index++)
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, index, (double)0.0, -1)) {
          SDDS_SetError("Problem setting SDDS row values (dump_sigma)");
          SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    }
    
    if (!SDDS_SetRowValues(SDDS_table,  SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, ie, 
                           s_index, beam->z, s_index+1, name, s_index+2, occurence,
                           s_index+3, type_name, -1)) {
      SDDS_SetError("Problem setting row values for SDDS table (dump_sigma)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    
    name = eptr->name;
    type_name = entity_name[eptr->type];
    occurence = eptr->occurence;
    if (eptr->succ)
      eptr = eptr->succ;
    else if (beamline->elem_recirc)
      eptr = beamline->elem_recirc;
    else
      eptr = &(beamline->elem);
  }

  if (!SDDS_WriteTable(SDDS_table)) {
    SDDS_SetError("Unable to write sigma data (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!inhibitFileSync)
    SDDS_DoFSync(SDDS_table);
  if (!SDDS_EraseData(SDDS_table)) {
    SDDS_SetError("Unable to erase sigma data (dump_sigma)");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

long check_sdds_column(SDDS_TABLE *SDDS_table, char *name, char *units)
{
  char *units1;
  if (SDDS_GetColumnIndex(SDDS_table, name)<0)
    return(0);
  if (SDDS_GetColumnInformation(SDDS_table, "units", &units1, SDDS_GET_BY_NAME, name)!=SDDS_STRING) {
    SDDS_SetError("units field of column has wrong data type!");
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (!units || SDDS_StringIsBlank(units)) {
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

#define BEAM_SCATTER_PARAMETERS 5
static SDDS_DEFINITION beam_scatter_parameter[BEAM_SCATTER_PARAMETERS] = {
  {"Particles", "&parameter name=Particles, type=long, description=\"Total simulated scatted particles\" &end"},
  {"IntRate", "&parameter name=IntRate, type=double, units=\"1/s\", description=\"Integrated Scattering Rate\" &end"},
  {"PLocalRate", "&parameter name=PLocalRate, type=double, units=\"1/s/m\", description=\"Piwinski's Local Rate\" &end"},
  {"SLocalRate", "&parameter name=SLocalRate, type=double, units=\"1/s/m\", description=\"Simulated Local Rate\" &end"},
  {"IgnoredRate", "&parameter name=IgnoredRate, type=double, units=\"1/s/m\", description=\"Ignored Scattering Rate\" &end"},
};
#define BEAM_SCATTER_COLUMNS 9
static SDDS_DEFINITION beam_scatter_column[BEAM_SCATTER_COLUMNS] = {
    {"x", "&column name=x, units=m, type=double &end"},
    {"xp", "&column name=xp, symbol=\"x'\", type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"yp", "&column name=yp, symbol=\"y'\", type=double &end"},
    {"t", "&column name=t, units=s, type=double &end"},
    {"p", "&column name=p, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    {"Rate", "&column name=Rate, units=\"1/s/m\", type=double &end"},
    {"IntRate", "&column name=IntRate, units=\"1/s\", type=double &end"},
};

void SDDS_BeamScatterSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_BeamScatterSetup");

    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            beam_scatter_parameter, BEAM_SCATTER_PARAMETERS, beam_scatter_column, BEAM_SCATTER_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);

    log_exit("SDDS_BeamScatterSetup");
}

void dump_scattered_particles(SDDS_TABLE *SDDS_table, double **particle, 
                              long particles, double *weight, TSCATTER *tsptr)
{
    long i;
    double rate;

    log_entry("dump_scattered_particles");
    if (!particle)
        bomb("NULL coordinate pointer passed to dump_scattered_particles", NULL);

    for (i=0; i<particles; i++)
        if (!particle[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_scattered_particles)\n", i);
            fflush(stdout);
            abort();
        }

    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_scattered_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    for (i=0; i<particles; i++) {
      rate = weight[i]*tsptr->IntR/tsptr->s_rate;
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particle[i][0], 1, particle[i][1], 2, particle[i][2], 3, particle[i][3],
                               4, particle[i][4]/c_mks, 5, (particle[i][5]+1.)*tsptr->betagamma,
                               6, (long)particle[i][6], 7, weight[i], 8, rate, -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_scattered_particles)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    }

    if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Particles", particles, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "IntRate", tsptr->IntR, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "PLocalRate", tsptr->p_rate, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "SLocalRate", tsptr->s_rate, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "IgnoredRate", tsptr->i_rate, NULL))) {
        SDDS_SetError("Problem setting SDDS parameters (dump_scattered_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }

    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Problem writing SDDS table (dump_scattered_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    } 
       
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);
    if (!SDDS_Terminate(SDDS_table)) {
      SDDS_SetError("Problem terminating 'distribution' file (TouschekDistribution)");
      SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    log_exit("dump_scattered_particles");    
}
  
#define BEAM_SCATTER_LOSS_COLUMNS 13
static SDDS_DEFINITION beam_scatter_loss_column[BEAM_SCATTER_LOSS_COLUMNS] = {
    {"x0", "&column name=x0, units=m, type=double &end"},
    {"xp0", "&column name=xp0, symbol=\"x'\", type=double &end"},
    {"y0", "&column name=y0, units=m, type=double &end"},
    {"yp0", "&column name=yp0, symbol=\"y'\", type=double &end"},
    {"s0", "&column name=s0, units=m, type=double &end"},
    {"p0", "&column name=p0, units=\"m$be$nc\", type=double &end"},
    {"particleID", "&column name=particleID, type=long &end"},
    {"x", "&column name=x, units=m, type=double &end"},
    {"y", "&column name=y, units=m, type=double &end"},
    {"s", "&column name=s, units=m, type=double &end"},
    {"Pass", "&column name=Pass, type=long &end"},
    {"Rate", "&column name=Rate, units=\"1/s/m\", type=double &end"},
    {"IntRate", "&column name=IntRate, units=\"1/s\", type=double &end"},
} ;
#define BEAM_SCATTER_LOSS_PARAMETERS 6
static SDDS_DEFINITION beam_scatter_loss_parameter[BEAM_SCATTER_LOSS_PARAMETERS] = {
  {"Particles", "&parameter name=Particles, type=long, description=\"Total lost simulated scatted particles\" &end"},
  {"S0", "&parameter name=S0, type=double, units=m, description=\"Scatter location\" &end"},
  {"TotalRate", "&parameter name=TotalRate, type=double, units=\"1/s/m\", description=\"Total simulated scatter rate\" &end"},
  {"LossRate", "&parameter name=LossRate, type=double, units=\"1/s/m\", description=\"Total loss rate including the ignored rate\" &end"},
  {"IntRate", "&parameter name=IntRate, type=double, units=\"1/s\", description=\"Integrated Scattering Rate\" &end"},
  {"IntLossRate", "&parameter name=IntLossRate, type=double, units=\"1/s\", description=\"Total integrated loss rate\" &end"},
};

void SDDS_BeamScatterLossSetup(SDDS_TABLE *SDDS_table, char *filename, long mode, long lines_per_row, char *contents, 
                          char *command_file, char *lattice_file, char *caller)
{
    log_entry("SDDS_BeamScatterLossSetup");

    SDDS_ElegantOutputSetup(SDDS_table, filename, mode, lines_per_row, contents, command_file, lattice_file,
                            beam_scatter_loss_parameter, BEAM_SCATTER_LOSS_PARAMETERS, 
                            beam_scatter_loss_column, BEAM_SCATTER_LOSS_COLUMNS,
                            caller, SDDS_EOS_NEWFILE|SDDS_EOS_COMPLETE);
    log_exit("SDDS_BeamScatterLossSetup");
}

void dump_scattered_loss_particles(SDDS_TABLE *SDDS_table, double **particleLos, double **particleOri, 
                                   long *lostOnPass, long particles, double *weight, TSCATTER *tsptr)
{
    long i, j;
    double rate, lossRate, intLossRate;

    log_entry("dump_scattered_loss_particles");
    if (!particleLos)
        bomb("NULL coordinate pointer passed to dump_scattered_loss_particles", NULL);

    for (i=0; i<particles; i++)
        if (!particleLos[i]) {
            fprintf(stdout, "error: coordinate slot %ld is NULL (dump_scattered_loss_particles)\n", i);
            fflush(stdout);
            abort();
            }

    if (!SDDS_StartTable(SDDS_table, particles)) {
        SDDS_SetError("Problem starting SDDS table (dump_scattered_loss_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    lossRate=intLossRate=0.;
    for (i=0; i<particles; i++) {
        j = (long)particleLos[i][6]-1;
        rate = weight[j]*tsptr->IntR/tsptr->s_rate;
        lossRate += weight[j];
        intLossRate += rate;
        if (!SDDS_SetRowValues(SDDS_table, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                               0, particleOri[j][0], 1, particleOri[j][1], 2, particleOri[j][2], 3, particleOri[j][3],
                               4, tsptr->s, 5, (1.+particleOri[j][5])*tsptr->betagamma, 6, (long)particleOri[j][6],
                               7, particleLos[i][0],8, particleLos[i][2],9, particleLos[i][4],
                               10, lostOnPass[i], 11, weight[j], 12, rate, -1)) {
            SDDS_SetError("Problem setting SDDS row values (dump_scattered_loss_particles)");
            SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
            }
        }
    if ((!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "Particles", particles, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "S0", tsptr->s, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "TotalRate", tsptr->s_rate, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "LossRate", lossRate, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "IntRate", tsptr->IntR, NULL))||
        (!SDDS_SetParameters(SDDS_table, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, "IntLossRate", intLossRate, NULL))) {
        SDDS_SetError("Problem setting SDDS parameters (dump_scattered_loss_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        }
    if (!SDDS_WriteTable(SDDS_table)) {
        SDDS_SetError("Problem writing SDDS table (dump_scattered_loss_particles)");
        SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
        } 
    if (!inhibitFileSync)
      SDDS_DoFSync(SDDS_table);

    log_exit("dump_scattered_loss_particles");
    
}

void computeEmitTwissFromSigmaMatrix(double *emit, double *emitc, double *beta, double *alpha, double sigma[6][6], long plane)
{
  if (plane<0 || plane>4 || plane%2!=0)
    bomb("invalid value for plane in computeEmitTwissFromSigmaMatrix", NULL);
  if (!emit || (!emitc && (alpha || beta)))
    bomb("invalid pointers passed to computeEmitTwissFromSigmaMatrix", NULL);
  
  /* emittance */
  *emit = SAFE_SQRT(sigma[0+plane][0+plane]*sigma[1+plane][1+plane] - sqr(sigma[0+plane][1+plane]));

  if (emitc) {
    /* corrected emittance */
    if (sigma[5][5] && plane!=4)
      *emitc = SAFE_SQRT(sqr(*emit) - 
                         (sqr(sigma[0+plane][5])*sigma[1+plane][1+plane] -
                          2*sigma[0+plane][1+plane]*sigma[0+plane][5]*sigma[1+plane][5] +
                          sqr(sigma[1+plane][5])*sigma[0+plane][0+plane])/sigma[5][5]);
    else
      *emitc = *emit;
    
    /* beta and alpha */
    if (beta && alpha) {
      *beta = *alpha = 0;
      if (*emitc) {
        if (sigma[5][5] && plane!=4) {
          double s11c, s12c;
          s11c = sigma[0+plane][0+plane] - sqr(sigma[0+plane][5])/sigma[5][5];
          s12c = sigma[0+plane][1+plane] - sigma[0+plane][5]*sigma[1+plane][5]/sigma[5][5];
          *beta = s11c/(*emitc);
          *alpha = -s12c/(*emitc);
        } else {
          *beta = sigma[0+plane][0+plane]/(*emitc);
          *alpha = -sigma[0+plane][1+plane]/(*emitc);
        }
      }
    }
  }
}


