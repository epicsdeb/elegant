/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* program: elegant
 * purpose: accelerator simulation
 * Michael Borland, 1989-2004
 */
#include "mdb.h"
#include "mdbsun.h"
#include "track.h"
#include "elegant.h"
#include "SDDS.h"
#include "scan.h"
#include <ctype.h>
#include "match_string.h"
#include <signal.h>
#include <time.h>
#if defined(UNIX) || defined(_WIN32)
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#endif

#if defined(linux)
#include <sched.h>
#include <linux/version.h>
#else
#define KERNEL_VERSION(a,b,c) 0
#endif

#include "chromDefs.h"
#include "correctDefs.h"
#include "tuneDefs.h"

void traceback_handler(int code);
void createSemaphoreFile(char *filename);
void readApertureInput(NAMELIST_TEXT *nltext, RUN *run);
void initialize_structures(RUN *run_conditions, VARY *run_control, ERRORVAL *error_control, CORRECTION *correct, 
                           BEAM *beam, OUTPUT_FILES *output_data, OPTIMIZATION_DATA *optimize,
                           CHROM_CORRECTION *chrom_corr_data, TUNE_CORRECTION *tune_corr_data,
                           ELEMENT_LINKS *links);
void do_semaphore_setup(char **semaphoreFile, NAMELIST_TEXT *nltext);
void printFarewell(FILE *fp);
void closeBeamlineOutputFiles(LINE_LIST *beamline);
void setSigmaIndices();
void process_particle_command(NAMELIST_TEXT *nltext);
void processGlobalSettings(NAMELIST_TEXT *nltext);
void freeInputObjects();

#define DESCRIBE_INPUT 0
#define DEFINE_MACRO 1
#define DEFINE_CPU_LIST 2
#define N_OPTIONS 3
char *option[N_OPTIONS] = {
    "describeinput",
    "macro",
    "cpulist",
        };

#define SHOW_USAGE    0x0001
#define SHOW_GREETING 0x0002

void showUsageOrGreeting (unsigned long mode)
{
#if USE_MPI
  char *USAGE="usage: mpirun -np <number of processes> Pelegant <inputfile> [-macro=<tag>=<value>,[...]]";
  char *GREETING="This is elegant 24.1Beta1, "__DATE__", by M. Borland, W. Guo, V. Sajaev, Y. Wang, Y. Wu, and A. Xiao.\nParallelized by Y. Wang, H. Shang, and M. Borland.";
#else
  char *USAGE="usage: elegant <inputfile> [-macro=<tag>=<value>,[...]]";
  char *GREETING="This is elegant 24.1Beta1, "__DATE__", by M. Borland, W. Guo, V. Sajaev, Y. Wang, Y. Wu, and A. Xiao.";
#endif
  if (mode&SHOW_GREETING)
    puts(GREETING);
  if (mode&SHOW_USAGE)
    puts(USAGE);
}


#define RUN_SETUP        0
#define RUN_CONTROL      1
#define VARY_ELEMENT     2
#define ERROR_CONTROL    3
#define ERROR_ELEMENT    4
#define SET_AWE_BEAM     5
#define SET_BUNCHED_BEAM 6
#define CORRECTION_SETUP 7
#define MATRIX_OUTPUT    8
#define TWISS_OUTPUT     9
#define TRACK           10
#define STOP            11
#define OPTIMIZATION_SETUP  12
#define OPTIMIZE        13
#define OPTIMIZATION_VARIABLE   14
#define OPTIMIZATION_CONSTRAINT 15
#define OPTIMIZATION_COVARIABLE 16
#define SAVE_LATTICE    17
#define RPN_EXPRESSION  18
#define PROGRAM_TRACE   19
#define CHROMATICITY 20
#define CLOSED_ORBIT    21
#define FIND_APERTURE   22
#define ANALYZE_MAP     23
#define CORRECT_TUNES   24
#define LINK_CONTROL    25
#define LINK_ELEMENTS   26
#define STEERING_ELEMENT 27
#define AMPLIF_FACTORS 28
#define PRINT_DICTIONARY 29
#define FLOOR_COORDINATES 30
#define CORRECTION_MATRIX_OUTPUT 31
#define LOAD_PARAMETERS 32
#define SET_SDDS_BEAM   33
#define SUBPROCESS      34
#define FIT_TRACES      35
#define SASEFEL_AT_END  36
#define ALTER_ELEMENTS  37
#define OPTIMIZATION_TERM 38
#define SLICE_ANALYSIS  39
#define DIVIDE_ELEMENTS 40
#define TUNE_SHIFT_WITH_AMPLITUDE 41
#define TRANSMUTE_ELEMENTS 42
#define TWISS_ANALYSIS 43
#define SEMAPHORES     44
#define FREQUENCY_MAP  45
#define INSERT_SCEFFECTS 46 
#define MOMENTUM_APERTURE 47
#define APERTURE_INPUT    48
#define COUPLED_TWISS_OUTPUT 49
#define LINEAR_CHROMATIC_TRACKING_SETUP 50
#define RPN_LOAD 51
#define MOMENTS_OUTPUT 52
#define TOUSCHEK_SCATTER 53
#define INSERT_ELEMENTS 54
#define CHANGE_PARTICLE 55
#define GLOBAL_SETTINGS 56
#define REPLACE_ELEMENTS 57
#define APERTURE_DATAX   58
#define MODULATE_ELEMENTS 59
#define PARALLEL_OPTIMIZATION_SETUP 60
#define RAMP_ELEMENTS 61
#define N_COMMANDS      62

char *command[N_COMMANDS] = {
    "run_setup", "run_control", "vary_element", "error_control", "error_element", "awe_beam", "bunched_beam",
    "correct", "matrix_output", "twiss_output", "track", "stop", 
    "optimization_setup", "optimize", "optimization_variable", "optimization_constraint",
    "optimization_covariable", "save_lattice", "rpn_expression", "trace", "chromaticity", "closed_orbit",
    "find_aperture", "analyze_map", "correct_tunes", "link_control", "link_elements",
    "steering_element", "amplification_factors", "print_dictionary", "floor_coordinates", "correction_matrix_output",
    "load_parameters", "sdds_beam", "subprocess", "fit_traces", "sasefel", "alter_elements",
    "optimization_term", "slice_analysis", "divide_elements", "tune_shift_with_amplitude",
    "transmute_elements", "twiss_analysis", "semaphores", "frequency_map", "insert_sceffects", "momentum_aperture", 
    "aperture_input", "coupled_twiss_output", "linear_chromatic_tracking_setup", "rpn_load",
    "moments_output", "touschek_scatter", "insert_elements", "change_particle", "global_settings","replace_elements",
    "aperture_data", "modulate_elements", "parallel_optimization_setup", "ramp_elements",
  } ;

char *description[N_COMMANDS] = {
    "run_setup                        defines lattice input file, primary output files, tracking order, etc.",
    "run_control                      defines number of steps, number of passes, number of indices, etc.",
    "vary_element                     defines element family and item to vary with an index",
    "error_control                    sets up and control random errors",
    "error_element                    defines an element family and item to add errors to",
    "awe_beam                         defines name of input beam data file, type of data, and some preprocessing",
    "bunched_beam                     defines beam distribution",
    "correct                          requests orbit or trajectory correction and define parameters",
    "matrix_output                    requests SDDS-format output or printed output of the matrix",
    "twiss_output                     requests output of Twiss parameters, chromaticity, and acceptance",
    "track                            command to begin tracking",
    "stop                             command to stop reading input file and end the run",
    "optimization_setup               requests running of optimization mode and sets it up",
    "optimize                         command to begin optimization",
    "optimization_variable            defines an element family and item to vary for optimization",
    "optimization_constraint          defines a constraint on optimization",
    "optimization_covariable          defines an element family and item to compute from optimization variables",
    "optimization_term                specifies an individual term in the optimization equation",
    "save_lattice                     command to save the current lattice",
    "rpn_expression                   command to execute an rpn expression (useful for optimization)",
    "trace                            requests tracing of program calls and defines parameters of trace",
    "chromaticity                     requests correction of the chromaticity",
    "closed_orbit                     requests output of the closed orbit",
    "find_aperture                    command to do aperture searches",
    "analyze_map                      command to do map analysis",
    "correct_tunes                    requests correction of the tunes",
    "link_control                     sets up and control element links",
    "link_elements                    defines a link between two sets of elements",
    "steering_element                 defines an element (group) and item for steering the beam",
    "amplification_factors            computes orbit/trajectory amplification factors",
    "print_dictionary                 prints a list of accelerator element types and allowed parameters",
    "floor_coordinates                computes the floor coordinates for the ends of elements",
    "correction_matrix_output         prints response matrices and their inverses",
    "load_parameters                  sets up loading of parameter values for elements",
    "sdds_beam                        defines name of input beam data file",
    "subprocess                       executes a string in a sub-shell",
    "fit_traces                       obtains a lattice model by fitting to multiple tracks through a beamline",
    "sasefel                          computes parameters of SASE FEL at end of system",
    "alter_elements                   alters a common parameter for one or more elements",
    "slice_analysis                   computes and outputs slice analysis of the beam",
    "divide_elements                  sets up parser to automatically divide specified elements into parts",
    "tune_shift_with_amplitude        sets up twiss module for computation of tune shift with amplitude",
    "transmute_elements               defines transmutation of one element type into another",
    "twiss_analysis                   requests twiss analysis of regions of a beamline (for optimization)",
    "semaphores                       requests use of semaphore files to indicate run start and end",
    "frequency_map                    command to perform frequency map analysis",
    "insert_sceffects                 add space charge element to beamline and set calculation flags", 
    "momentum_aperture                determine momentum aperture from tracking as a function of position in the ring",
    "aperture_data                    provide an SDDS file with the physical aperture vs s", 
    "coupled_twiss_output             compute coupled beamsizes and twiss parameters",
    "linear_chromatic_tracking_setup  set up chromatic derivatives for linear chromatic tracking",
    "rpn_load                         load SDDS data into rpn variables",
    "moments_output                   perform moments computations and output to file",
    "touschek_scatter                 calculate Touschek lifetime, simulate touschek scattering effects, find out momentum aperture through tracking", 
    "insert_elements                  insert elements into already defined beamline", 
    "change_particle                  change the particle type",
    "global_settings                  change various global settigs",
    "replace_elements                 remove or replace elements inside beamline",
    "aperture_input                   provide an SDDS file with the physical aperture vs s (same as aperture_data)", 
    "modulate_elements                modulate values of elements as a function of time",
    "parallel_optimization_setup      requests running of parallel optimization mode and sets it up",
    "ramp_elements                    ramp values of elements as a function of pass",
  } ;

#define NAMELIST_BUFLEN 65536

#define DEBUG 0

static VARY run_control;
static char *semaphoreFile[3];
 
long writePermitted = 1;
/* For serial version, isMaster and isSlave are both set to 1 */
long isMaster = 1;
long isSlave = 1;
long notSinglePart=0; /* All the processors will do the same thing by default */
long runInSinglePartMode = 0; /* This flag will be set as true under some special cases, such as genetic optimization */
double factor = 1.0;    /* In serial version, the memory will be allocted for all the particles */ 
long do_find_aperture = 0;

#if USE_MPI
int n_processors = 1;
int myid;
int dumpAcceptance = 0;
parallelMode parallelStatus = trueParallel; 
int partOnMaster = 1; /* indicate if the particle information is available on master */
long watch_not_allowed = 0;
long lessPartAllowed = 0; /* By default, the number of particles is required to be at least n_processors-1 */
MPI_Comm workers;
int fd; /* save the duplicated file descriptor stdout to use it latter */
long last_optimize_function_call = 0;
int min_value_location = 0;
#endif

#ifdef SET_DOUBLE
void 
set_fpu (unsigned int mode)
{
  __asm__ ("fldcw %0" : : "m" (*&mode));
}
#endif
 
int run_coupled_twiss_output(RUN *run, LINE_LIST *beamline, double *starting_coord);
void finish_coupled_twiss_output();
void run_rpn_load(NAMELIST_TEXT *nltext, RUN *run);
void setupLinearChromaticTracking(NAMELIST_TEXT *nltext, LINE_LIST *beamline);
void setup_coupled_twiss_output(NAMELIST_TEXT *nltext, RUN *run, 
                                LINE_LIST *beamline, long *do_coupled_twiss_output,
                                long default_order);



int main(argc, argv)
int argc;
char **argv;
{
  char **macroTag, **macroValue=NULL;
  long macros;
  LINE_LIST *beamline=NULL;        /* pointer to root of linked list */
  FILE *fp_in=NULL;
  char *inputfile;
  SCANNED_ARG *scanned;
  char s[NAMELIST_BUFLEN], *ptr;
  long i;
  RUN run_conditions;
  CORRECTION correct;
  ERRORVAL error_control;
  BEAM beam;
  OUTPUT_FILES output_data;
  OPTIMIZATION_DATA optimize;
  CHROM_CORRECTION chrom_corr_data;
  TUNE_CORRECTION tune_corr_data;
  ELEMENT_LINKS links;
  char *saved_lattice = NULL;
  long correction_setuped, run_setuped, run_controled, error_controled, beam_type, commandCode;
  long do_chromatic_correction = 0, do_twiss_output = 0, fl_do_tune_correction = 0, do_coupled_twiss_output = 0;
  long do_moments_output = 0;
  long do_closed_orbit = 0, do_matrix_output = 0, do_response_output = 0;
  long last_default_order = 0, new_beam_flags, links_present, twiss_computed = 0, moments_computed = 0;
  long correctionDone, linear_chromatic_tracking_setup_done = 0;
  double *starting_coord, finalCharge;
  long namelists_read = 0, failed, firstPass;
  double apertureReturn;
#if USE_MPI
#ifdef MPI_DEBUG
  FILE *fpError; 
  char fileName[15];
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;  
#endif 
  MPI_Group world_group, worker_group;
  int ranks[1];
#endif

  semaphoreFile[0] = semaphoreFile[1] = semaphoreFile[2] = NULL;

#if USE_MPI
  MPI_Init(&argc,&argv);
  /* get the total number of processors */
  MPI_Comm_size(MPI_COMM_WORLD, &n_processors);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  /* create a new communicator group with the slave processors only */
  ranks[0] = 0;  /* first process is master */
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_excl(world_group, 1, ranks, &worker_group);
  MPI_Comm_create(MPI_COMM_WORLD, worker_group, &workers);

  if (myid!=0) {
#ifdef MPI_DEBUG
    sprintf(fileName, "error.%d", myid);
    fpError = fopen(fileName, "w");
    freopen(fileName, "w", stdout);
    freopen(fileName, "w", stderr);
    MPI_Get_processor_name(processor_name,&namelen);
    fprintf(stdout, "Process %d on %s\n", myid, processor_name);
#else   
    /* duplicate an open file descriptor, tested with gcc */
    fd = dup(fileno(stdout));
    /* redirect output, only the master processor will write on screen or files */
#if defined(_WIN32)
    freopen("NUL","w",stdout);
    /* freopen("NUL","w",stderr); */
#else
    freopen("/dev/null","w",stdout);
    /*    freopen("/dev/null","w",stderr); */
#endif
#endif
    writePermitted = isMaster = 0;
    isSlave = 1;
  }
  else
    isSlave = 0;

  if (sizeof(int)<4) { /* The size of integer is assumed to be 4 bytes to handle a large number of particles */
    printf("Warning!!! The INT_MAX could be too small to record the number of particles.\n"); 
  }
#if !SDDS_MPI_IO
  if (isSlave && (n_processors>3))   /* This will avoid wasting memory on a laptop with a small number of cores */    
    factor = 2.0/(n_processors-1);   /* In parallel version, a portion of memory will be allocated on each slave */
#endif
#endif
 
#ifdef SET_DOUBLE
  set_fpu (0x27F);  /* use double-precision rounding */
#endif

  signal(SIGINT, traceback_handler);
  signal(SIGILL, traceback_handler);
  signal(SIGABRT, traceback_handler);
  signal(SIGFPE, traceback_handler);
  signal(SIGSEGV, traceback_handler);
#ifndef _WIN32
  signal(SIGHUP, traceback_handler);
  signal(SIGQUIT, traceback_handler);
  signal(SIGTRAP, traceback_handler);
  signal(SIGBUS, traceback_handler);
#endif
  
  log_entry("main");
  if (!SDDS_CheckTableStructureSize(sizeof(SDDS_TABLE))) {
    fprintf(stderr, "table structure size is inconsistent\n");
    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
    exitElegant(1);
  }
  load_hash = NULL;     
  compute_offsets();
  set_max_name_length(100);
  setSigmaIndices();
  
  macros = 1;
  if (!(macroTag = malloc(sizeof(*macroTag))) ||
      !(macroValue = malloc(sizeof(*macroValue)))) 
    bombElegant("memory allocation failure setting up default macro tag/value pairs", NULL);
  macroTag[0] = "INPUTFILENAME";
  macroValue[0] = NULL;  /* will fill in later */

#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
  init_stats();
#endif

  /*
    #ifdef SUNOS4
    malloc_debug(1);
    mallopt(M_MXFAST, 1024);
    mallopt(M_NLBLKS, 8);
    #endif
    */
  
  argc = scanargs(&scanned, argc, argv);
  if (argc<2) {
    showUsageOrGreeting(SHOW_USAGE|SHOW_GREETING);
    fflush(stdout);
    link_date();
    exitElegant(1);
  }
  
  showUsageOrGreeting(SHOW_GREETING);
  fflush(stdout);
  link_date();
  if (getenv("RPN_DEFNS")) {
    rpn(getenv("RPN_DEFNS"));
    if (rpn_check_error()) exitElegant(1);
  }

  inputfile = NULL;
  for (i=1; i<argc; i++) {
    if (scanned[i].arg_type==OPTION) {
      switch (match_string(scanned[i].list[0], option, N_OPTIONS, 0)) {
      case DESCRIBE_INPUT:
        show_namelists_fields(stdout, namelist_pointer, namelist_name, n_namelists);
        if (argc==2)
          exitElegant(0);
        break;
      case DEFINE_MACRO:
        if ((scanned[i].n_items-=1)<1) {
          fprintf(stdout, "Invalid -macro syntax\n");
          showUsageOrGreeting(SHOW_USAGE);
          exitElegant(1);
        }
        if (!(macroTag=SDDS_Realloc(macroTag, sizeof(*macroTag)*(macros+scanned[i].n_items))) ||
            !(macroValue=SDDS_Realloc(macroValue, sizeof(*macroValue)*(macros+scanned[i].n_items))))
          bombElegant("memory allocation failure (-macro)", NULL);
        else {
          long j;
          for (j=0; j<scanned[i].n_items; j++) {
            macroTag[macros] = scanned[i].list[j+1];
            if (!(macroValue[macros] = strchr(macroTag[macros], '='))) {
              fprintf(stdout, "Invalid -macro syntax\n");
              showUsageOrGreeting(SHOW_USAGE);
              exitElegant(1);
            }
            macroValue[macros][0] = 0;
            macroValue[macros] += 1;
            macros++;
          }
        }
        break;
      case DEFINE_CPU_LIST:
#if (!USE_MPI && defined(linux))
#if (LINUX_VERSION_CODE >= KERNEL_VERSION(2,6,0))
        if (scanned[i].n_items<2) {
          fprintf(stdout, "Invalid -cpuList syntax\n");
          showUsageOrGreeting(SHOW_USAGE);
          exitElegant(1);
        }
        else {
          long j, k;
          unsigned long processorMask = 0;
          for (j=1; j<scanned[i].n_items; j++) {
            if (sscanf(scanned[i].list[j], "%ld", &k)!=1 && k<0)  {
              fprintf(stdout, "Invalid -cpuList syntax\n");
              showUsageOrGreeting(SHOW_USAGE);
              exitElegant(1);
            }
            processorMask |= (unsigned long)(ipow(2, k)+0.5);
          }
          printf("processorMask = %lx\n", processorMask);
          sched_setaffinity(0, sizeof(processorMask), &processorMask);
        }
#else /* (LINUX_VERSION_CODE >= KERNEL_VERSION(2,6,0)) */
        printf("warning: CPU list ignored\n");
#endif /* (LINUX_VERSION_CODE >= KERNEL_VERSION(2,6,0)) */
#else /*(!USE_MPI && defined(linux)) */
        printf("warning: CPU list ignored\n");
#endif /*(!USE_MPI && defined(linux)) */
        break;
      default:
        fprintf(stdout, "Unknown option given.\n");
        showUsageOrGreeting(SHOW_USAGE);
        break;
      }
    }
    else { 
      /* filenames */
      if (!inputfile)
        fp_in = fopen_e(inputfile = scanned[i].list[0], "r", 0);
      else {
        fprintf(stdout, "Too many file names listed.\n");
        showUsageOrGreeting(SHOW_USAGE);
        exitElegant(1);
      }
    }
  }
  
  if (!inputfile) {
    fprintf(stdout, "No input file was given.\n");
    showUsageOrGreeting(SHOW_USAGE);
    exitElegant(1);
  }
  cp_str(&macroValue[0], inputfile);
  
#if defined(CONDOR_COMPILE)
  sprintf(s, "%s.ckpt", inputfile);
  init_image_with_file_name(s);
#else
#ifndef _WIN32
  signal(SIGUSR2, SIG_IGN);
#endif
#endif
  
  initialize_structures(&run_conditions, &run_control, &error_control, &correct, &beam, &output_data,
                        &optimize, &chrom_corr_data, &tune_corr_data, &links);
  
  run_setuped = run_controled = error_controled = correction_setuped = 0;
  
  starting_coord = tmalloc(sizeof(*starting_coord)*7);
  
  beam_type = -1;

  while (get_namelist(s, NAMELIST_BUFLEN, fp_in)) {
    substituteTagValue(s, NAMELIST_BUFLEN, macroTag, macroValue, macros);
#if DEBUG
    fprintf(stderr, "%s\n", s);
#endif
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
    report_stats(stdout, "statistics: ");
    fflush(stdout);
#endif
    if (namelists_read)
      free_namelist_text(&namelist_text);
    scan_namelist(&namelist_text, s);
    namelists_read = 1;
#if USE_MPI /* synchronize all the processes before execute an input statement */ 
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    switch ((commandCode=match_string(namelist_text.group_name, command, N_COMMANDS, EXACT_MATCH))) {
    case CHANGE_PARTICLE:
      if (run_setuped)
        bombElegant("particle command should precede run_setup", NULL);  /* to ensure nothing is inconsistent */
      process_particle_command(&namelist_text);
      break;
    case RUN_SETUP:
      beam_type = -1;
      
      initialize_structures(NULL, &run_control, &error_control, &correct, &beam, &output_data,
                            &optimize, &chrom_corr_data, &tune_corr_data, &links);
      clearSliceAnalysis();
      finish_load_parameters();
      run_setuped = run_controled = error_controled = correction_setuped = 0;
      
      run_setuped = run_controled = error_controled = correction_setuped = do_closed_orbit = do_chromatic_correction = 
        fl_do_tune_correction = 0;
      do_twiss_output = do_matrix_output = do_response_output = do_coupled_twiss_output = do_moments_output = do_find_aperture = 0;
      linear_chromatic_tracking_setup_done = 0;

      set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
      set_print_namelist_flags(0);
      if (processNamelist(&run_setup, &namelist_text)==NAMELIST_ERROR)
        bombElegant(NULL, NULL);
      if (echoNamelists) print_namelist(stdout, &run_setup);
      setSearchPath(search_path);
      /* check for validity of namelist inputs */
      if (lattice==NULL) {
        if (!saved_lattice)
          bombElegant("no lattice supplied", NULL);
        if (default_order!=last_default_order)
          delete_matrix_data(NULL);
      }
      else {
        /* free previous lattice info */
	if (beamline)
          closeBeamlineOutputFiles(beamline);
        free_elements(NULL);
        free_beamlines(NULL);
        freeInputObjects();
        saved_lattice = lattice;
      }
      free_beamdata(&beam);
      if (default_order<1 || default_order>3)
        bombElegant("default_order is out of range", NULL);
      if (concat_order>3)
        bombElegant("concat_order is out of range", NULL);
      if (p_central && p_central_mev)
        bombElegant("give only one of p_central and p_central_mev", NULL);
      if (p_central_mev!=0 && p_central==0)
        p_central = p_central_mev/particleMassMV;
      if (p_central<=0 && !expand_for)
        bombElegant("p_central<=0 and p_central_mev<=0", NULL);
      if (expand_for)
        p_central = find_beam_p_central(expand_for);
      if (random_number_seed==0) {
        random_number_seed = (long)time((time_t)0);
        random_number_seed = 2*(random_number_seed/2) + 1;
        fprintf(stdout, "clock-generated random_number_seed = %ld\n", random_number_seed);
        fflush(stdout);
      }
      
      seedElegantRandomNumbers(random_number_seed, 0);

      /* copy run data into run_conditions structure */
      run_conditions.ideal_gamma = sqrt(sqr(p_central)+1);
      run_conditions.p_central = p_central;
      run_conditions.default_order = default_order;
      run_conditions.concat_order = concat_order;
      run_conditions.print_statistics = print_statistics;
      run_conditions.combine_bunch_statistics = combine_bunch_statistics;
      run_conditions.wrap_around = wrap_around;
      if ((run_conditions.final_pass = final_pass))
        run_conditions.wrap_around = 1;
      run_conditions.tracking_updates = tracking_updates;
      run_conditions.always_change_p0 = always_change_p0;
      run_conditions.load_balancing_on = load_balancing_on;
      run_conditions.random_sequence_No = random_sequence_No;
      remaining_sequence_No = random_sequence_No; /* For Pelegant regression test */

      /* extract the root filename from the input filename */
      strcpy_ss(s, inputfile);
      if (rootname==NULL) {
        clean_filename(s);
        if ((ptr=strrchr(s, '.')))
          *ptr = 0;
        cp_str(&rootname, s);
        run_conditions.rootname = rootname;
        run_conditions.runfile  = inputfile;
      }
      else {
        run_conditions.rootname = rootname;
        run_conditions.runfile  = compose_filename(inputfile, rootname);
      }
      /* In the version with parallel I/O, these file names need to be known by all the processors, otherwise there
	 will be a synchronization issue when calculated accumulated sum in do_tracking */
      run_conditions.acceptance = compose_filename(acceptance, rootname);
#if USE_MPI
     if (run_conditions.acceptance)
        dumpAcceptance = 1;
#endif
      run_conditions.centroid   = compose_filename(centroid, rootname);
      run_conditions.sigma      = compose_filename(sigma, rootname);

      if (run_conditions.apertureData.initialized && !run_conditions.apertureData.persistent)
        resetApertureData(&(run_conditions.apertureData)); 

#if !SDDS_MPI_IO
      if (isMaster)
#endif
      {
	run_conditions.output     = compose_filename(output, rootname);
	run_conditions.losses     = compose_filename(losses, rootname);
	run_conditions.final      = compose_filename(final, rootname);
      }
#if !SDDS_MPI_IO
      else
	run_conditions.final = run_conditions.output = run_conditions.losses = NULL;
#endif     

      if (isMaster) { /* These files will be written by master */
	magnets                   = compose_filename(magnets, rootname);
	semaphore_file            = compose_filename(semaphore_file, rootname);
	parameters                = compose_filename(parameters, rootname);        
      }
      else {
        magnets = semaphore_file = parameters = NULL;
      }
      
      if (semaphore_file && fexists(semaphore_file))
        remove(semaphore_file);
		
      if (semaphoreFile[0]) {
        /* "started" */
        semaphoreFile[0] = compose_filename(semaphoreFile[0], rootname);
	createSemaphoreFile(semaphoreFile[0]);
      }
      if (semaphoreFile[1]) {
        /* "done" */
        semaphoreFile[1] = compose_filename(semaphoreFile[1], rootname);
        if (fexists(semaphoreFile[1]))
          remove(semaphoreFile[1]);
      }
      if (semaphoreFile[2]) {
        /* "failed" */
        semaphoreFile[2] = compose_filename(semaphoreFile[2], rootname);
        if (fexists(semaphoreFile[2]))
          remove(semaphoreFile[2]);
      }
      
      /* parse the lattice file and create the beamline */
      run_conditions.lattice = compose_filename(saved_lattice, rootname);
      if (element_divisions>1)
        addDivisionSpec("*", NULL, NULL, element_divisions, 0.0);
#ifdef USE_MPE /* use the MPE library */
      if (USE_MPE) {
        int event1a, event1b;
	event1a = MPE_Log_get_event_number();
	event1b = MPE_Log_get_event_number();
	if(isMaster) 
	  MPE_Describe_state(event1a, event1b, "get_beamline", "blue");
	MPE_Log_event(event1a, 0, "start get_beamline"); /* record time spent on reading input */ 
#endif
      beamline = get_beamline(lattice, use_beamline, p_central, echo_lattice);
#ifdef  USE_MPE
	      MPE_Log_event(event1b, 0, "end get_beamline");
      }
#endif
      fprintf(stdout, "length of beamline %s per pass: %21.15e m\n", beamline->name, beamline->revolution_length);
      fflush(stdout);
      lattice = saved_lattice;
      
      /* output the magnet layout */
      if (magnets)
        output_magnets(magnets, lattice, beamline);

      delete_phase_references();    /* necessary for multi-step runs */
      reset_special_elements(beamline, 1);
      reset_driftCSR();
      last_default_order = default_order;
      run_setuped = 1;
      break;
    case GLOBAL_SETTINGS:
      processGlobalSettings(&namelist_text);
      break;
    case RUN_CONTROL:
      if (!run_setuped)
        bombElegant("run_setup must precede run_control namelist", NULL);
      vary_setup(&run_control, &namelist_text, &run_conditions, beamline);
      run_control.ready = 1;
      run_controled = 1;
      break;
    case VARY_ELEMENT:
      if (!run_controled)
        bombElegant("run_control must precede vary_element namelists", NULL);
      if (beam_type!=-1)
        bombElegant("vary_element statements must come before beam definition", NULL);
      add_varied_element(&run_control, &namelist_text, &run_conditions, beamline);
      break;
    case ERROR_CONTROL:
      if (!run_setuped || !run_controled)
        bombElegant("run_setup and run_control must precede error_control namelist", NULL);
      if (beam_type!=-1)
        bombElegant("error specifications must be completed before beam type is specified", NULL);
      error_setup(&error_control, &namelist_text, &run_conditions, beamline);
      error_controled = 1;
      break;                    
    case ERROR_ELEMENT:
      if (beam_type!=-1)
        bombElegant("error_element statements must come before beam definition", NULL);
      if (!error_controled)
        bombElegant("error_control namelist must precede error_element namelists", NULL);
      add_error_element(&error_control, &namelist_text, beamline);
      break;
    case CORRECTION_SETUP:
      if (!run_setuped)
        bombElegant("run_setup must precede correction", NULL);
      if (beam_type!=-1)
        bombElegant("beam setup (bunched_beam or sdds_beam) must follow correction setup", NULL);
      correction_setuped = 1;
      correction_setup(&correct, &namelist_text, &run_conditions, beamline); 
      delete_phase_references();
      reset_special_elements(beamline, 1);
      reset_driftCSR();
      break;
    case SET_AWE_BEAM: 
      fprintf(stdout, "This program no longer supports awe-format files.\n");
      fflush(stdout);
      fprintf(stdout, "Use awe2sdds to convert your data files, and use\n");
      fflush(stdout);
      fprintf(stdout, "the sdds_beam command instead of awe_beam.\n");
      fflush(stdout);
      break;
    case SET_BUNCHED_BEAM:
      if (!run_setuped || !run_controled)
        bombElegant("run_setup and run_control must precede bunched_beam namelist", NULL);
      setup_bunched_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, &optimize.variables,
                         &output_data, beamline, beamline->n_elems,
                         correct.mode!=-1 && 
                         (correct.track_before_and_after || correct.start_from_centroid));
      setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
      beam_type = SET_BUNCHED_BEAM;
      break;
    case SET_SDDS_BEAM: 
      if (!run_setuped || !run_controled)
        bombElegant("run_setup and run_control must precede sdds_beam namelist", NULL);
      setup_sdds_beam(&beam, &namelist_text, &run_conditions, &run_control, &error_control, 
                      &optimize.variables, &output_data, beamline, beamline->n_elems,
                      correct.mode!=-1 && 
                      (correct.track_before_and_after || correct.start_from_centroid));
      setup_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, beamline);
      beam_type = SET_SDDS_BEAM;
      break;
    case TRACK:
    case ANALYZE_MAP:
      switch (commandCode) {
      case TRACK:
        set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
        set_print_namelist_flags(0);
        if (processNamelist(&track, &namelist_text)==NAMELIST_ERROR)
          bombElegant(NULL, NULL);
        if (echoNamelists) print_namelist(stdout, &track);
        run_conditions.stopTrackingParticleLimit = stop_tracking_particle_limit;
#if USE_MPI
        if (stop_tracking_particle_limit!=-1)
          bombElegant("stop_tracking_particle_limit feature not supported in Pelegant", NULL);
#endif
        if (use_linear_chromatic_matrix && 
            !(linear_chromatic_tracking_setup_done || twiss_computed || do_twiss_output))
          bombElegant("you must compute twiss parameters or give linear_chromatic_tracking_setup to do linear chromatic tracking", NULL);
        if (longitudinal_ring_only && !(twiss_computed || do_twiss_output))
          bombElegant("you must compute twiss parameters to do longitudinal ring tracking", NULL);
        if (use_linear_chromatic_matrix && longitudinal_ring_only)
          bombElegant("can't do linear chromatic tracking and longitudinal-only tracking together", NULL);
        if (beam_type==-1)
          bombElegant("beam must be defined prior to tracking", NULL);
        break;
      case ANALYZE_MAP:
        setup_transport_analysis(&namelist_text, &run_conditions, &run_control, &error_control);
        break;
      }
      new_beam_flags = 0;
      firstPass = 1;
      while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
        /* vary_beamline asserts changes due to vary_element, error_element, and load_parameters */
        fill_double_array(starting_coord, 6, 0.0);
        correctionDone = 0;
        if (correct.mode!=-1 && (correct.track_before_and_after || correct.start_from_centroid)) {
          if (beam_type==SET_AWE_BEAM) {
            bombElegant("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
          }
          else if (beam_type==SET_SDDS_BEAM) {
            if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, 0)<0)
              break;
          }
          else
            new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, 0);
          new_beam_flags = TRACK_PREVIOUS_BUNCH;
          if (commandCode==TRACK && correct.track_before_and_after)
            track_beam(&run_conditions, &run_control, &error_control, &optimize.variables, 
                       beamline, &beam, &output_data, 
                       PRECORRECTION_BEAM, 0, &finalCharge);
        }

        /* If needed, find closed orbit, twiss parameters, moments, and response matrix, but don't do
         * any output unless requested to do so "pre-correction"
         */
        if (do_closed_orbit && !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
            !soft_failure) {
          fprintf(stdout, "Closed orbit not found---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_twiss_output && !run_twiss_output(&run_conditions, beamline, starting_coord, 0) &&
            !soft_failure) {
          fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_moments_output)
          runMomentsOutput(&run_conditions, beamline, starting_coord, 0, 1);
        if (do_response_output)
          run_response_output(&run_conditions, beamline, &correct, 0);

        if (correct.mode!=-1 || fl_do_tune_correction || do_chromatic_correction) {
          if (correct.use_actual_beam) {
            if (beam_type==SET_AWE_BEAM) {
              bombElegant("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
            }
            else if (beam_type==SET_SDDS_BEAM) {
              if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, 0)<0)
                break;
            }
            else
              new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, 0);
            new_beam_flags = TRACK_PREVIOUS_BUNCH;
          }
          for (i=failed=0; i<correction_iterations; i++) {
            if (correction_iterations>1) {
              fprintf(stdout, "\nOrbit/tune/chromaticity correction iteration %ld\n", i+1);
              fflush(stdout);
            }
            if (correct.mode!=-1 && 
                !do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, 
                               run_control.i_step, 
                               (i==0?INITIAL_CORRECTION:0)+(i==correction_iterations-1?FINAL_CORRECTION:0))) {
              fputs("warning: orbit correction failed--continuing with next step\n", stdout);
              continue;
            }
            if (fl_do_tune_correction) {
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
              if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                      run_control.i_step, i==correction_iterations-1) &&
                  !soft_failure) {
                fprintf(stdout, "Tune correction failed---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
            }
            if (do_chromatic_correction) {
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
              if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                              run_control.i_step, i==correction_iterations-1) &&
                  !soft_failure) {
                fprintf(stdout, "Chromaticity correction failed---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
            }
            correctionDone = 1;
          }
          if (failed)
            continue;
        }
        if (beam_type==SET_AWE_BEAM) {
          bombElegant("beam type of SET_AWE_BEAM in main routine--this shouldn't happen", NULL);
        }
        else if (beam_type==SET_SDDS_BEAM) {
          if (new_sdds_beam(&beam, &run_conditions, &run_control, &output_data, new_beam_flags)<0)
            break;
        }
        else if (beam_type==SET_BUNCHED_BEAM)
          new_bunched_beam(&beam, &run_conditions, &run_control, &output_data, new_beam_flags);

        /* Assert post-correction perturbations */
        perturb_beamline(&run_control, &error_control, &run_conditions, beamline); 

        /* Do post-correction output */
        if (do_closed_orbit && !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
            !soft_failure) {
          fprintf(stdout, "Closed orbit not found---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_twiss_output && !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
            !soft_failure) {
          fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_moments_output)
          runMomentsOutput(&run_conditions, beamline, starting_coord, 1, 1);
        if (do_coupled_twiss_output &&
            run_coupled_twiss_output(&run_conditions, beamline, starting_coord) &&
            !soft_failure) {
          fprintf(stdout, "Coupled twiss parameters computation failed\n");
          fflush(stdout);
        }
        if (do_response_output)
          run_response_output(&run_conditions, beamline, &correct, 1);
        if (center_on_orbit)
          center_beam_on_coords(beam.particle, beam.n_to_track, starting_coord, center_momentum_also);
	else if (offset_by_orbit)
          offset_beam_by_coords(beam.particle, beam.n_to_track, starting_coord, offset_momentum_also);
        run_matrix_output(&run_conditions, beamline);
        if (firstPass) {
          /* prevent fiducialization of RF etc. by correction etc. */
          delete_phase_references();
          reset_special_elements(beamline, 1);
          reset_driftCSR();
        }
        firstPass = 0;
        switch (commandCode) {
        case TRACK:
          /* Finally, do tracking */
          track_beam(&run_conditions, &run_control, &error_control, &optimize.variables, 
                     beamline, &beam, &output_data, 
                     (use_linear_chromatic_matrix?LINEAR_CHROMATIC_MATRIX:0)+
                     (longitudinal_ring_only?LONGITUDINAL_RING_ONLY:0)+
                     (ibs_only?IBS_ONLY_TRACKING:0), 0, &finalCharge);
          break;
        case ANALYZE_MAP:
          do_transport_analysis(&run_conditions, &run_control, &error_control, beamline, 
                                (do_closed_orbit || correct.mode!=-1?starting_coord:NULL));
          break;
        }
        if (parameters)
          dumpLatticeParameters(parameters, &run_conditions, beamline);
      }
      if (commandCode==TRACK)
        finish_output(&output_data, &run_conditions, &run_control, &error_control, &optimize.variables, 
                      beamline, beamline->n_elems, &beam, finalCharge);
      if (beam_type)
        free_beamdata(&beam);
      if (do_closed_orbit)
        finish_clorb_output();
      if (do_twiss_output)
        finish_twiss_output();
      if (do_moments_output)
        finishMomentsOutput();
      if (do_coupled_twiss_output)
        finish_coupled_twiss_output();
      if (do_response_output)
        finish_response_output();
      if (parameters)
        finishLatticeParametersFile();
#ifdef SUNOS4
      check_heap();
#endif
      switch (commandCode) {
      case TRACK:
        fprintf(stdout, "Finished tracking.\n");
        break;
      case ANALYZE_MAP:
        fprintf(stdout, "Finished transport analysis.\n");
        break;
      }
      fflush(stdout);
      /* reassert defaults for namelist run_setup */
      lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses = 
        parameters = NULL;
      combine_bunch_statistics = 0;
      random_number_seed = 987654321;
      wrap_around = 1;
      final_pass = 0;
      default_order = 2;
      concat_order = 0;
      tracking_updates = 1;
      concat_order = print_statistics = p_central = 0;
      run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
        fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_coupled_twiss_output = do_response_output = 0;
      break;
    case MATRIX_OUTPUT:
      if (!run_setuped)
        bombElegant("run_setup must precede matrix_output namelist", NULL);
      setup_matrix_output(&namelist_text, &run_conditions, beamline);
      do_matrix_output = 1;
      break;
    case TWISS_OUTPUT:
      if (!run_setuped)
        bombElegant("run_setup must precede twiss_output namelist", NULL);
      setup_twiss_output(&namelist_text, &run_conditions, beamline, &do_twiss_output,
                         run_conditions.default_order);
      if (!do_twiss_output) {
        twiss_computed = 1;
        run_twiss_output(&run_conditions, beamline, NULL, -1);
        delete_phase_references();
        reset_special_elements(beamline, 1);
        reset_driftCSR();
        finish_twiss_output();
      }
      break;
    case MOMENTS_OUTPUT:
      if (!run_setuped)
        bombElegant("run_setup must precede moments_output namelist", NULL);
      setupMomentsOutput(&namelist_text, &run_conditions, beamline, &do_moments_output,
                         run_conditions.default_order);
      if (!do_moments_output) {
        moments_computed = 1;
        runMomentsOutput(&run_conditions, beamline, NULL, -1, 1);
        delete_phase_references();
        reset_special_elements(beamline, 1);
        reset_driftCSR();
        finishMomentsOutput();
      }
      break;
    case COUPLED_TWISS_OUTPUT:
      if (!run_setuped)
        bombElegant("run_setup must precede coupled_twiss_output namelist", NULL);
      setup_coupled_twiss_output(&namelist_text, &run_conditions, beamline, &do_coupled_twiss_output,
                                 run_conditions.default_order);
      if (!do_coupled_twiss_output) {
        run_coupled_twiss_output(&run_conditions, beamline, NULL);
        delete_phase_references();
        reset_special_elements(beamline, 1);
        reset_driftCSR();
        finish_coupled_twiss_output();
      }
      break;
    case TUNE_SHIFT_WITH_AMPLITUDE:
      if (do_twiss_output)
        bombElegant("you must give tune_shift_with_amplitude before twiss_output", NULL);
      setupTuneShiftWithAmplitude(&namelist_text, &run_conditions);
      break;
    case SEMAPHORES:
      if (run_setuped)
        bombElegant("you must give the semaphores command before run_setup", NULL);
      do_semaphore_setup(semaphoreFile, &namelist_text);
      break;
    case STOP:
      lorentz_report();
      finish_load_parameters();
      if (semaphore_file)
        createSemaphoreFile(semaphore_file);
      if (semaphoreFile[1])
        createSemaphoreFile(semaphoreFile[1]);
      free_beamdata(&beam);
      printFarewell(stdout);
#if USE_MPI
      if (isSlave)
	MPI_Comm_free(&workers); 
      MPI_Group_free(&worker_group); 
      close(fd); 
      MPI_Finalize();
#endif
      exitElegant(0);
      break;
    case OPTIMIZATION_SETUP:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      do_optimization_setup(&optimize, &namelist_text, &run_conditions, beamline);
      break;
#if USE_MPI
    case PARALLEL_OPTIMIZATION_SETUP:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      do_parallel_optimization_setup(&optimize, &namelist_text, &run_conditions, beamline);
      break;
#endif
    case OPTIMIZE:
      if (beam_type==-1)
        bombElegant("beam definition must come before optimize command", NULL);
      while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
        do_optimize(&namelist_text, &run_conditions, &run_control, &error_control, beamline, &beam,
                    &output_data, &optimize, &chrom_corr_data, beam_type, do_closed_orbit,
                    do_chromatic_correction, &correct, correct.mode, &tune_corr_data, 
                    fl_do_tune_correction, do_find_aperture, do_response_output);
        if (parameters)
          dumpLatticeParameters(parameters, &run_conditions, beamline);
      }
      if (parameters)
        finishLatticeParametersFile();
      /* reassert defaults for namelist run_setup */
      lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses =
        parameters = NULL;
      combine_bunch_statistics = 0;
      random_number_seed = 987654321;
      wrap_around = 1;
      final_pass = 0;
      default_order = 2;
      concat_order = 0;
      tracking_updates = 1;
      concat_order = print_statistics = p_central = 0;
      run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
        fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_coupled_twiss_output = do_response_output = 0;
#if USE_MPI
      runInSinglePartMode = 0; /* We should set the flag to the normal parallel tracking after parallel optimization */
#endif
      break;
    case OPTIMIZATION_VARIABLE:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      add_optimization_variable(&optimize, &namelist_text, &run_conditions, beamline);
      break;
    case OPTIMIZATION_CONSTRAINT:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      add_optimization_constraint(&optimize, &namelist_text, &run_conditions, beamline);
      break;
    case OPTIMIZATION_COVARIABLE:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      add_optimization_covariable(&optimize, &namelist_text, &run_conditions, beamline);
      break;
    case OPTIMIZATION_TERM:
      if (beam_type!=-1)
        bombElegant("optimization statements must come before beam definition", NULL);
      add_optimization_term(&optimize, &namelist_text, &run_conditions, beamline);
      break;
    case SAVE_LATTICE:
      do_save_lattice(&namelist_text, &run_conditions, beamline);
      break;
    case RPN_EXPRESSION:
      run_rpn_expression(&namelist_text);
      break;
    case RPN_LOAD:
      run_rpn_load(&namelist_text, &run_conditions);
      break;
    case PROGRAM_TRACE:
      process_trace_request(&namelist_text);
      break;
    case CHROMATICITY:
      setup_chromaticity_correction(&namelist_text, &run_conditions, beamline, &chrom_corr_data);
      do_chromatic_correction = 1;
      break;
    case CORRECT_TUNES:
      setup_tune_correction(&namelist_text, &run_conditions, beamline, &tune_corr_data);
      fl_do_tune_correction = 1;
      break;
    case CLOSED_ORBIT:
      setup_closed_orbit(&namelist_text, &run_conditions, beamline);
      do_closed_orbit = 1;
      if (correction_setuped)
        fprintf(stdout, "warning: you've asked to do both closed-orbit calculation and orbit correction.\nThis may duplicate effort.\n");
      fflush(stdout);
      break;
    case FIND_APERTURE:
    case FREQUENCY_MAP:
    case MOMENTUM_APERTURE:
      switch (commandCode) {
      case FIND_APERTURE:
        setup_aperture_search(&namelist_text, &run_conditions, &run_control, &do_find_aperture);
        if (do_find_aperture) continue;
        break;
      case FREQUENCY_MAP:
        setupFrequencyMap(&namelist_text, &run_conditions, &run_control);
        break;
      case MOMENTUM_APERTURE:
        setupMomentumApertureSearch(&namelist_text, &run_conditions, &run_control);
        break;
      }
      while (vary_beamline(&run_control, &error_control, &run_conditions, beamline)) {
#if DEBUG
	fprintf(stdout, "semaphore_file = %s\n", semaphore_file?semaphore_file:NULL);
#endif  
        fill_double_array(starting_coord, 6, 0.0);
        if (correct.mode!=-1 || fl_do_tune_correction || do_chromatic_correction) {
          for (i=failed=0; i<correction_iterations; i++) {
            if (correction_iterations>1) {
              fprintf(stdout, "\nOrbit/tune/chromaticity correction iteration %ld\n", i+1);
              fflush(stdout);
            }
            if (correct.mode!=-1 && 
                !do_correction(&correct, &run_conditions, beamline, starting_coord, &beam, 
                               run_control.i_step, 
                               (i==0?INITIAL_CORRECTION:0)+(i==correction_iterations-1?FINAL_CORRECTION:0))) {
              fputs("warning: orbit correction failed--continuing with next step\n", stdout);
              continue;
            }
            if (fl_do_tune_correction) {
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
              if (!do_tune_correction(&tune_corr_data, &run_conditions, beamline, starting_coord,
                                      run_control.i_step, i==correction_iterations-1) &&
                  !soft_failure) {
                fprintf(stdout, "Tune correction failed---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
            }
            if (do_chromatic_correction) {
              if (do_closed_orbit && 
                  !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 0) &&
                  !soft_failure) {
                fprintf(stdout, "Closed orbit not found---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
              if (!do_chromaticity_correction(&chrom_corr_data, &run_conditions, beamline, starting_coord,
                                              run_control.i_step, i==correction_iterations-1) &&
                  !soft_failure) {
                fprintf(stdout, "Chromaticity correction failed---continuing to next step\n");
                fflush(stdout);
                failed = 1;
                break;
              }
            }
            correctionDone = 1;
          }
          if (failed)
            continue;
        }
        
        perturb_beamline(&run_control, &error_control, &run_conditions, beamline); 
        if (do_closed_orbit && 
            !run_closed_orbit(&run_conditions, beamline, starting_coord, &beam, 1) &&
            !soft_failure) {
          fprintf(stdout, "Closed orbit not found---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_twiss_output && 
            !run_twiss_output(&run_conditions, beamline, starting_coord, 1) &&
            !soft_failure) {
          fprintf(stdout, "Twiss parameters not defined---continuing to next step\n");
          fflush(stdout);
          continue;
        }
        if (do_coupled_twiss_output &&
            run_coupled_twiss_output(&run_conditions, beamline, starting_coord) &&
            !soft_failure) {
          fprintf(stdout, "Coupled twiss parameters calculation failed.\n");
          fflush(stdout);
        }
        run_matrix_output(&run_conditions, beamline);
        if (do_response_output)
          run_response_output(&run_conditions, beamline, &correct, 1);
        if (parameters)
          dumpLatticeParameters(parameters, &run_conditions, beamline);
        switch (commandCode) {
        case FIND_APERTURE:
          do_aperture_search(&run_conditions, &run_control, starting_coord,
			     &error_control, beamline, &apertureReturn);
          break;
        case FREQUENCY_MAP:
          doFrequencyMap(&run_conditions, &run_control, starting_coord, &error_control, beamline);
          break;
        case MOMENTUM_APERTURE:
          doMomentumApertureSearch(&run_conditions, &run_control, &error_control, beamline, 
                                   (correct.mode!=-1 || do_closed_orbit)?starting_coord:NULL);
          break;
        }
      }
      fprintf(stdout, "Finished all tracking steps.\n"); fflush(stdout);
      fflush(stdout);
      switch (commandCode) {
      case FIND_APERTURE:
        finish_aperture_search(&run_conditions, &run_control, &error_control, beamline);
        break;
      case FREQUENCY_MAP:
        finishFrequencyMap();
        break;
      case MOMENTUM_APERTURE:
        finishMomentumApertureSearch();
        break;
      }
      if (do_closed_orbit)
        finish_clorb_output();
      if (parameters)
        finishLatticeParametersFile();
      if (beam_type!=-1)
        free_beamdata(&beam);
      if (do_closed_orbit)
        finish_clorb_output();
      if (do_twiss_output)
        finish_twiss_output();
      if (do_response_output)
        finish_response_output();
#ifdef SUNOS4
      check_heap();
#endif
      switch (commandCode) {
      case FIND_APERTURE:
	fprintf(stdout, "Finished dynamic aperture search.\n");
        break;
      case FREQUENCY_MAP:
	fprintf(stdout, "Finished frequency map analysis.\n");
        break;
      case MOMENTUM_APERTURE:
	fprintf(stdout, "Finished momentum aperture search.\n");
        break;
      }
#if DEBUG
      fprintf(stdout, "semaphore_file = %s\n", semaphore_file?semaphore_file:NULL);
#endif  
      fflush(stdout);
      /* reassert defaults for namelist run_setup */
      lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses =
        parameters = NULL;
      combine_bunch_statistics = 0;
      random_number_seed = 987654321;
      wrap_around = 1;
      final_pass = 0;
      default_order = 2;
      concat_order = 0;
      tracking_updates = 1;
      concat_order = print_statistics = p_central = 0;
      run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
        fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_coupled_twiss_output = do_response_output = 0;
      break;
    case LINK_CONTROL:
      if (!run_setuped || !run_controled)
        bombElegant("run_setup and run_control must precede link_control namelist", NULL);
      element_link_control(&links, &namelist_text, &run_conditions, beamline);
      break;                    
    case LINK_ELEMENTS:
      if (!run_setuped || !run_controled)
        bombElegant("run_setup and run_control must precede link_elements namelist", NULL);
      if (!beamline)
        bombElegant("beamline not defined--can't add element links", NULL);
      add_element_links(&links, &namelist_text, beamline);
      links_present = 1;
      beamline->links = &links;
      break;
    case STEERING_ELEMENT:
      if (correction_setuped)
        bombElegant("you must define steering elements prior to giving the 'correct' namelist", NULL);
      add_steering_element(&correct, beamline, &run_conditions, &namelist_text);
      break;
    case AMPLIF_FACTORS:
      if (parameters)
        dumpLatticeParameters(parameters, &run_conditions, beamline);
      compute_amplification_factors(&namelist_text, &run_conditions, &correct, do_closed_orbit, beamline);
      break;
    case PRINT_DICTIONARY:
      set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
      set_print_namelist_flags(0);
      if (processNamelist(&print_dictionary, &namelist_text)==NAMELIST_ERROR)
        bombElegant(NULL, NULL);
      if (echoNamelists) print_namelist(stdout, &print_dictionary);
      do_print_dictionary(filename, latex_form, SDDS_form);
      break;
    case FLOOR_COORDINATES:
      if (isMaster)
        output_floor_coordinates(&namelist_text, &run_conditions, beamline);
      break;
    case CORRECTION_MATRIX_OUTPUT:
      if (!run_setuped)
        bombElegant("run setup must precede correction_matrix_output namelist", NULL);
      setup_correction_matrix_output(&namelist_text, &run_conditions, beamline, &correct,
                                     &do_response_output, 
                                     do_twiss_output+do_matrix_output+twiss_computed);
      if (!do_response_output) {
        run_response_output(&run_conditions, beamline, &correct, -1);
        delete_phase_references();
        reset_special_elements(beamline, 1);
        reset_driftCSR();
        finish_response_output();
      }
      break;
    case LOAD_PARAMETERS:
      if (!run_setuped)
        bombElegant("run_setup must precede load_parameters namelists", NULL);
      if (run_controled)
        bombElegant("load_parameters namelists must precede run_control namelist", NULL);
      if (error_controled)
        bombElegant("load_parameters namelists must precede error_control and error namelists", NULL);
      if (setup_load_parameters(&namelist_text, &run_conditions, beamline) && magnets)
        /* make sure the magnet output is right in case loading parameters changed something */
        output_magnets(magnets, lattice, beamline);
      break;
    case SUBPROCESS:
      if (isMaster)
        run_subprocess(&namelist_text, &run_conditions);
      break;
    case FIT_TRACES:
#if 0
      do_fit_trace_data(&namelist_text, &run_conditions, beamline);
      if (parameters) {
        dumpLatticeParameters(parameters, &run_conditions, beamline);
        finishLatticeParametersFile();
      }
      /* reassert defaults for namelist run_setup */
      lattice = use_beamline = acceptance = centroid = sigma = final = output = rootname = losses =
        parameters = NULL;
      combine_bunch_statistics = 0;
      random_number_seed = 987654321;
      wrap_around = 1;
      final_pass = 0;
      default_order = 2;
      concat_order = 0;
      tracking_updates = 1;
      concat_order = print_statistics = p_central = 0;
      run_setuped = run_controled = error_controled = correction_setuped = do_chromatic_correction =
        fl_do_tune_correction = do_closed_orbit = do_twiss_output = do_coupled_twiss_output = do_response_output = 0;
#endif
      break;
    case SASEFEL_AT_END:
      if (!run_setuped)
        bombElegant("run_setup must precede sasefel namelist", NULL);
      if (beam_type!=-1)
        bombElegant("sasefel namelist must precede beam definition", NULL);
      setupSASEFELAtEnd(&namelist_text, &run_conditions, &output_data);
      break;
    case ALTER_ELEMENTS:
      if (!run_setuped)
        bombElegant("run_setup must precede alter_element namelist", NULL);
      do_alter_element(&namelist_text, &run_conditions, beamline);
      break;
    case SLICE_ANALYSIS:
      if (!run_setuped)
        bombElegant("run_setup must precede slice_analysis namelist", NULL);
      if (beam_type!=-1)
        bombElegant("slice_analysis namelist must precede beam definition", NULL);
      setupSliceAnalysis(&namelist_text, &run_conditions, &output_data);
      break;
    case DIVIDE_ELEMENTS: 
      if (run_setuped)
        bombElegant("divide_elements must precede run_setup", NULL);
      setupDivideElements(&namelist_text, &run_conditions, beamline);
      break;
    case TRANSMUTE_ELEMENTS:
      if (run_setuped)
        bombElegant("transmute_elements must precede run_setup", NULL);
      setupTransmuteElements(&namelist_text, &run_conditions, beamline);
      break;
    case INSERT_SCEFFECTS:
      if (run_setuped)
        bombElegant("insert_sceffects must precede run_setup", NULL);
      setupSCEffect(&namelist_text, &run_conditions, beamline);
      break;
    case INSERT_ELEMENTS:
      if (!run_setuped)
        bombElegant("run_setup must precede insert_element namelist", NULL);
      do_insert_elements(&namelist_text, &run_conditions, beamline);
      break;
    case REPLACE_ELEMENTS:
      if (!run_setuped)
        bombElegant("run_setup must precede replace_element namelist", NULL);
      do_replace_elements(&namelist_text, &run_conditions, beamline);
      break;
    case TOUSCHEK_SCATTER:
      if (!run_setuped)
        bombElegant("run_setup must precede touschek_scatter namelist", NULL);
      TouschekEffect(&run_conditions, &run_control, &error_control, beamline, &namelist_text);
      break;
    case TWISS_ANALYSIS:
      if (do_twiss_output)
        bombElegant("twiss_analysis must come before twiss_output", NULL);
      setupTwissAnalysisRequest(&namelist_text, &run_conditions, beamline);
      break;
    case APERTURE_INPUT:
    case APERTURE_DATAX:
      readApertureInput(&namelist_text, &run_conditions);
      break;
    case LINEAR_CHROMATIC_TRACKING_SETUP:
      if (do_twiss_output)
        bombElegant("you can't give twiss_output and linear_chromatic_tracking_setup together", NULL);
      if (!run_setuped)
        bombElegant("run_setup must precede linear_chromatic_tracking_setup", NULL);
      setupLinearChromaticTracking(&namelist_text, beamline);
      linear_chromatic_tracking_setup_done = 1;
      break;
    case MODULATE_ELEMENTS:
      if (!run_setuped)
        bombElegant("run_setup must precede modulate_elements", NULL);
      addModulationElements(&(run_conditions.modulationData), &namelist_text, beamline);
      break;
    case RAMP_ELEMENTS:
      if (!run_setuped)
        bombElegant("run_setup must precede ramp_elements", NULL);
      addRampElements(&(run_conditions.rampData), &namelist_text, beamline);
      break;
    default:
      fprintf(stdout, "unknown namelist %s given.  Known namelists are:\n", namelist_text.group_name);
      fflush(stdout);
      for (i=0; i<N_COMMANDS; i++)
        fprintf(stdout, "%s\n", description[i]);
      fflush(stdout);
      exitElegant(1);
      break;
    }
#ifdef SUNOS4
    check_heap();
#endif
  }
#if DEBUG
  fprintf(stdout, "semaphore_file = %s\n", semaphore_file?semaphore_file:NULL);
#endif  
  fprintf(stdout, "End of input data encountered.\n"); fflush(stdout);
  fflush(stdout);
  lorentz_report();
  finish_load_parameters();
  if (semaphore_file)
    createSemaphoreFile(semaphore_file);
  if (semaphoreFile[1])
    createSemaphoreFile(semaphoreFile[1]);
  free_beamdata(&beam);
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
  report_stats(stdout, "statistics: ");
  fflush(stdout);
#endif
  log_exit("main");
  printFarewell(stdout);
#if USE_MPI
  close(fd);
  MPI_Finalize();
#endif
  if (load_hash)
    hdestroy(load_hash);                         /* destroy hash table */
  return(0);
}

void printFarewell(FILE *fp)
{
#if (!USE_MPI)
  fprintf(stdout, "=====================================================================================\n");
  fprintf(stdout, "Thanks for using elegant.  Please cite the following reference in your publications:\n");
  fprintf(stdout, "  M. Borland, \"elegant: A Flexible SDDS-Compliant Code for Accelerator Simulation,\"\n");
  fprintf(stdout, "  Advanced Photon Source LS-287, September 2000.\n");
  fprintf(stdout, "If you use a modified version, please indicate this in all publications.\n");
  fprintf(stdout, "=====================================================================================\n");
#else 
  fprintf(stdout, "=====================================================================================\n");
  fprintf(stdout, "Thanks for using Pelegant.  Please cite the following references in your publications:\n");
  fprintf(stdout, "  M. Borland, \"elegant: A Flexible SDDS-Compliant Code for Accelerator Simulation,\"\n");
  fprintf(stdout, "  Advanced Photon Source LS-287, September 2000.\n");
  fprintf(stdout, "  Y. Wang and M. Borland, \"Pelegant: A Parallel Accelerator Simulation Code for  \n");
  fprintf(stdout, "  Electron Generation and Tracking, Proceedings of the 12th Advanced Accelerator  \n"); 
  fprintf(stdout, "  Concepts Workshop, 2006.\n");
  fprintf(stdout, "If you use a modified version, please indicate this in all publications.\n");
  fprintf(stdout, "=====================================================================================\n");
#endif
}


double find_beam_p_central(char *input)
{
  SDDS_DATASET SDDSin;
  char s[SDDS_MAXLINE];
  double *p=NULL, psum;
  long i, rows;
#if SDDS_MPI_IO 
/* All the processes will read the wake file, but not in parallel.
   Zero the Memory when call  SDDS_InitializeInput */
  SDDSin.parallel_io = 0; 
#endif

  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, input) || !SDDS_ReadPage(&SDDSin)) {
    sprintf(s, "Problem opening beam input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if ((rows=SDDS_RowCount(&SDDSin))<=0 || !(p=SDDS_GetColumnInDoubles(&SDDSin, "p"))) {
    sprintf(s, "No data in input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  for (i=psum=0; i<rows; i++) 
    psum += p[i];
  if (SDDS_ReadPage(&SDDSin)>0)
    fprintf(stdout, "Warning: file %s has multiple pages.  Only the first is used for expand_for.\n",
            input);
    fflush(stdout);
  SDDS_Terminate(&SDDSin);
  fprintf(stdout, "Expanding about p = %21.15e\n", psum/rows);
  fflush(stdout);
  return psum/rows;
}

#ifdef SUNOS4
#include <malloc.h>

void check_heap() 
{
    struct mallinfo info;
    
    fprintf(stdout, "Performing memory heap verification..."); fflush(stdout);
    fflush(stdout);
    if (malloc_verify())
        fputs("okay.", stdout);
    else
        fputs("errors, detected.", stdout);
    fflush(stdout);
#if defined(VAX_VMS) || defined(UNIX) || defined(_WIN32)
    report_stats(stdout, "statistics: ");
    fflush(stdout);
#endif
    return;
    info.arena = info.ordblks = info.smblks = info.hblks = info.hblkhd = info.usmblks = 0;   
    info.fsmblks = info.uordblks = info.fordblks = info.keepcost = info.allocated = info.treeoverhead = 0;
    info = mallinfo();
    fprintf(stdout, "memory allocation information:\n");
    fflush(stdout);
    fprintf(stdout, "  total space in arena: %ld\n", info.arena);     
    fflush(stdout);
    fprintf(stdout, "  number of ordinary blocks: %ld\n", info.ordblks);   
    fflush(stdout);
    fprintf(stdout, "  number of small blocks: %ld\n", info.smblks);    
    fflush(stdout);
    fprintf(stdout, "  number of holding blocks: %ld\n", info.hblks);     
    fflush(stdout);
    fprintf(stdout, "  space in holding block headers: %ld\n", info.hblkhd);    
    fflush(stdout);
    fprintf(stdout, "  space in small blocks in use: %ld\n", info.usmblks);   
    fflush(stdout);
    fprintf(stdout, "  space in free small blocks: %ld\n", info.fsmblks);   
    fflush(stdout);
    fprintf(stdout, "  space in ordinary blocks in use: %ld\n", info.uordblks);  
    fflush(stdout);
    fprintf(stdout, "  space in free ordinary blocks: %ld\n", info.fordblks);  
    fflush(stdout);
    fprintf(stdout, "  cost of enabling keep option: %ld\n", info.keepcost);  
    fflush(stdout);
    fprintf(stdout, "  number of ordinary blocks allocated: %ld\n", info.allocated);
    fflush(stdout);
    fprintf(stdout, "  bytes used in maintaining the free tree: %ld\n", info.treeoverhead);
    fflush(stdout);
    }
#endif

void center_beam_on_coords(double **part, long np, double *coord, long center_dp)
{
    double sum, offset;
    long i, j, lim;
    
    if (!np)
        return;
    
    if (center_dp)
        lim = 5;
    else
        lim = 4;
    
    for (j=0; j<=lim; j++) {
        for (i=sum=0; i<np; i++)
            sum += part[i][j];
        offset = sum/np - coord[j];
        for (i=0; i<np; i++)
            part[i][j] -= offset;
        }
    }

void offset_beam_by_coords(double **part, long np, double *coord, long offset_dp)
{
    long i, j, lim;
    
    if (!np)
        return;
    
    if (offset_dp)
        lim = 5;
    else
        lim = 4;
    
    for (j=0; j<=lim; j++) {
        for (i=0; i<np; i++)
            part[i][j] += coord[j];
        }
    }

typedef struct {
  char *elementName;
  long index;
} DICTIONARY_ENTRY;

int dictionaryEntryCmp(const void *p1, const void *p2)
{
  DICTIONARY_ENTRY *de1, *de2;
  de1 = (DICTIONARY_ENTRY *)p1;
  de2 = (DICTIONARY_ENTRY *)p2;
  return strcmp(de1->elementName, de2->elementName);
}

void do_print_dictionary(char *filename, long latex_form, long SDDS_form)
{
  FILE *fp;
  long i;
  DICTIONARY_ENTRY *dictList;
  
  if (!filename)
    bombElegant("filename invalid (do_print_dictionary)", NULL);
  if (latex_form && SDDS_form)
    bombElegant("give latex_form=1 or SDDS_form=1, not both", NULL); 

  if (!(dictList = malloc(sizeof(*dictList)*N_TYPES)))
    bombElegant("memory allocation failure", NULL);
  for (i=1; i<N_TYPES; i++) {
    dictList[i-1].elementName = entity_name[i];
    dictList[i-1].index = i;
  }
  qsort((void*)dictList, N_TYPES-1, sizeof(*dictList), dictionaryEntryCmp);
  fp = fopen_e(filename, "w", 0);
#if DEBUG
  fprintf(stderr, "Opened file to print dictionary entries\n");
#endif
  if (latex_form) {
    fprintf(fp, "\\newlength{\\descwidth}\n");
    fprintf(fp, "\\setlength{\\descwidth}{2in}\n");
  } 
  if (SDDS_form) {
    fprintf(fp, "SDDS1\n");
    fprintf(fp, "&parameter name=ElementType, type=string &end\n");
    fprintf(fp, "&parameter name=ParallelCapable, type=short &end\n");
    fprintf(fp, "&column name=ParameterName, type=string &end\n");
    fprintf(fp, "&column name=Units, type=string &end\n");
    fprintf(fp, "&column name=Type, type=string &end\n");
    fprintf(fp, "&column name=Default, type=string &end\n");
    fprintf(fp, "&column name=Description, type=string &end\n");
    fprintf(fp, "&data mode=ascii &end\n");
  }
#if DEBUG
  fprintf(stderr, "Beginning loop to print dictionary entries\n");
#endif
  for (i=0; i<N_TYPES-1; i++) {
#if DEBUG
    fprintf(stderr, "Printing dictionary entry %ld (%s) of %ld\n", i, 
            dictList[i].elementName, (long)(N_TYPES-1));
#endif
    if (entity_description[dictList[i].index].flags&NO_DICT_OUTPUT)
      continue;
    print_dictionary_entry(fp, dictList[i].index, latex_form, SDDS_form);
  }
#if DEBUG
  fprintf(stderr, "Free'ing dictList\n");
#endif
  free(dictList);
#if DEBUG
  fprintf(stderr, "Closing file\n");
#endif
  fclose(fp);
#if DEBUG
  fprintf(stderr, "Exiting do_print_dictionary\n");
#endif
}

#define PRINTABLE_NULL(s) (s?s:"NULL")
char *translateUnitsToTex(char *source);
char *makeTexSafeString(char *source, long checkMath);

void print_dictionary_entry(FILE *fp, long type, long latex_form, long SDDS_form)
{
  char *type_name[3] = {"double", "long", "STRING", };
  char *description;
  long j, texLines, specialEntry, xyWaveforms;
  char buffer[16384];
  char *specialDescription = "Optionally used to assign an element to a group, with a user-defined name.  Group names will appear in the parameter output file in the column ElementGroup";
  if (latex_form) {
    fprintf(fp, "\\begin{latexonly}\n\\newpage\n\\begin{center}{\\Large\\verb|%s|}\\end{center}\n\\end{latexonly}\\subsection{%s}\n", 
            entity_name[type], entity_name[type]);
    fprintf(fp, "%s\n\\\\\n", makeTexSafeString(entity_text[type], 0));
    fprintf(fp, "Parallel capable? : %s\\\\\n", entity_description[type].flags&UNIPROCESSOR?"no":"yes");
    fprintf(fp, "\\begin{tabular}{|l|l|l|l|p{\\descwidth}|} \\hline\n");
    fprintf(fp, "Parameter Name & Units & Type & Default & Description \\\\ \\hline \n");
  } else {
    fprintf(fp, "%c***** element type %s:\n", SDDS_form?'!':'*', entity_name[type]);
    if (SDDS_form) {
      strcpy_ss(buffer, entity_name[type]);
      replace_chars(buffer, "\n\t", "  ");
      fprintf(fp, "%s\n", buffer);
      fprintf(fp, "%ld\n", (long)(entity_description[type].flags&UNIPROCESSOR?0:1));
      fprintf(fp, "%ld\n", entity_description[type].n_params+1);
    }
    else 
      fprintf(fp, "%s\n", entity_text[type]);
  }
  specialEntry = 0;
  xyWaveforms = 0;
  for (j=texLines=0; j<=entity_description[type].n_params; j++) {
    if (j==entity_description[type].n_params)
      specialEntry = 1;
    /* 35 lines fits on a latex page */
    if (latex_form && texLines>35) {
      texLines = 0;
      fprintf(fp, "\\end{tabular}\n\n");
      fprintf(fp, "\\begin{latexonly}\n\\newpage\n\\begin{center}{\\Large\\verb|%s| continued}\\end{center}\n\\end{latexonly}\n", 
              entity_name[type]);
      fprintf(fp, "%s\n\\\\\n", makeTexSafeString(entity_text[type], 0));
      fprintf(fp, "\\begin{tabular}{|l|l|l|l|p{\\descwidth}|} \\hline\n");
      fprintf(fp, "Parameter Name & Units & Type & Default & Description \\\\ \\hline \n");
    }
    if (SDDS_form) {
      fprintf(fp, "\"%s\" \"%s\" \"%s\"", 
              PRINTABLE_NULL(specialEntry ? "GROUP" : entity_description[type].parameter[j].name), 
              PRINTABLE_NULL(specialEntry ? "" : entity_description[type].parameter[j].unit),
              PRINTABLE_NULL(specialEntry ? "string" : type_name[entity_description[type].parameter[j].type-1]));
    } else if (!latex_form)
      fprintf(fp, "%20s %20s %10s", 
              PRINTABLE_NULL(specialEntry ? "GROUP" : entity_description[type].parameter[j].name), 
              PRINTABLE_NULL(specialEntry ? "" : entity_description[type].parameter[j].unit),
              PRINTABLE_NULL(specialEntry ? "string" : type_name[entity_description[type].parameter[j].type-1]));
    else {
      fprintf(fp, "%s ",
              makeTexSafeString(PRINTABLE_NULL(specialEntry ? "GROUP" : 
                                               entity_description[type].parameter[j].name), 0));
      fprintf(fp, "& %s ",
              translateUnitsToTex(PRINTABLE_NULL(specialEntry ? "" : 
                                                 entity_description[type].parameter[j].unit)));
      fprintf(fp, "& %s & ",
              makeTexSafeString(PRINTABLE_NULL(specialEntry ? "string" :
                                               type_name[entity_description[type].parameter[j].type-1]), 0));
    }
    if (!specialEntry) {
      if (entity_description[type].parameter[j].flags&PARAM_XY_WAVEFORM)
        xyWaveforms++;
      switch (entity_description[type].parameter[j].type) {
      case IS_DOUBLE:
        if (latex_form && entity_description[type].parameter[j].number==0)
          fprintf(fp, " 0.0");
        else 
          fprintf(fp, "  %.15g", entity_description[type].parameter[j].number);
        break;
      case IS_LONG:
        if (latex_form && entity_description[type].parameter[j].integer==0)
          fprintf(fp, " \\verb|0|");
        else      
          fprintf(fp, "  %-15ld", entity_description[type].parameter[j].integer);
        break;
      case IS_STRING:
        if (SDDS_form)
          fprintf(fp, " \"%s\"", 
                  PRINTABLE_NULL(entity_description[type].parameter[j].string));
        else 
	  fprintf(fp, "  %-15s", 
                  PRINTABLE_NULL(entity_description[type].parameter[j].string));
        break;
      default:
        fprintf(stdout, "Invalid parameter type for %s item of %s\n",
                PRINTABLE_NULL(entity_description[type].parameter[j].name),
                entity_name[type]);
        fflush(stdout);
        exitElegant(1);
      }
    } else {
      fprintf(fp, "NULL");
    }
    if (specialEntry)
      description = specialDescription;
    else
      description = entity_description[type].parameter[j].description;
    if (latex_form) {
      char *ptr0;
      strcpy_ss(buffer, description);
      if (strlen(ptr0 = buffer)) {
        /* add to lines counter based on estimate of wrap-around lines
           in the latex parbox. 28 is approximate number of char. in 2 in */
        texLines += strlen(ptr0) / 28;
        if (*ptr0) {
          fprintf(fp, " & %s ", 
                  makeTexSafeString(ptr0, 1));
          fprintf(fp, " \\\\ \\hline \n");
          texLines++;
        }
      } else {
        fprintf(fp, " & \\\\ \\hline \n");
        texLines++;
      }
    }
    else if (SDDS_form)
      fprintf(fp, "  \"%s\"\n", description);
    else
      fprintf(fp, "  %s\n", description);
  }
  if (latex_form) {
    char physicsFile[1024];

    fprintf(fp, "\\end{tabular}\n\n");


    sprintf(physicsFile, "%s.tex", entity_name[type]);
    str_tolower(physicsFile);
    if (fexists(physicsFile)) {
      fprintf(stderr, "Including file %s\n", physicsFile);
      fprintf(fp, "\\vspace*{0.5in}\n\\input{%s}\n", physicsFile);
    }
    if (xyWaveforms) 
      fprintf(fp, "\\vspace*{0.5in}\n\\input{xyWaveforms.tex}\n");
  }
}

#include <memory.h>

void initialize_structures(RUN *run_conditions, VARY *run_control, ERRORVAL *error_control, CORRECTION *correct, 
                           BEAM *beam, OUTPUT_FILES *output_data, OPTIMIZATION_DATA *optimize,
                           CHROM_CORRECTION *chrom_corr_data, TUNE_CORRECTION *tune_corr_data,
                           ELEMENT_LINKS *links)
{
    if (run_conditions)
        memset((void*)run_conditions, 0, sizeof(*run_conditions));
    if (run_control)
        memset((void*)run_control, 0, sizeof(*run_control));
    run_control->bunch_frequency = 2856e6;
    if (error_control)
        memset((void*)error_control, 0, sizeof(*error_control));
    if (correct)
        memset((void*)correct, 0, sizeof(*correct));
    correct->mode = correct->method = -1;
    if (beam)
        memset((void*)beam, 0, sizeof(*beam));
    if (output_data)
        memset((void*)output_data, 0, sizeof(*output_data));
    if (optimize)
        memset((void*)optimize, 0, sizeof(*optimize));
    if (chrom_corr_data)
        memset((void*)chrom_corr_data, 0, sizeof(*chrom_corr_data));
    if (tune_corr_data)
        memset((void*)tune_corr_data, 0, sizeof(*tune_corr_data));
    if (links)
        memset((void*)links, 0, sizeof(*links));
    }

char *makeTexSafeString(char *source, long checkMath)
{
  static char buffer[1024];
  long index = 0, math = 0;
  if (!source)
    return source;
  buffer[0] = 0;
  while (*source) {
    if (*source=='$' && checkMath) {
      math = !math;
      buffer[index++] = *source++;
      continue;
      }
    if (!math) {
      if (*source=='_' || *source=='^' || *source=='{' || *source=='}' || *source=='%') {
        buffer[index++] = '\\';
        buffer[index++] = *source++;
      } else if  (*source=='<' || *source=='>' || *source=='|') {
        buffer[index++] = '$';
        buffer[index++] = *source++;
        buffer[index++] = '$';
      } else
        buffer[index++] = *source++;
    }
    else {
      buffer[index++] = *source++;
    }
  }
  buffer[index] = 0;
  return buffer;
}


char *translateUnitsToTex(char *source)
{
  static char buffer[1024];
  char *ptr, *ptr1;

  if (strlen(source)==0)
    return source;
  buffer[0] = '$';
  buffer[1] = 0;
  ptr = source;
  while (*ptr) {
    if (*ptr=='$') {
      switch (*(ptr1=ptr+1)) {
      case 'a':
      case 'b':
      case 'A':
      case 'B':
        while (*ptr1 && *ptr1!='$')
          ptr1++;
        if (*ptr1 && (*(ptr1+1)=='n' || *(ptr1+1)=='N')) {
          if (*(ptr+1)=='a' || *(ptr+1)=='A')
            strcat(buffer, "^{");
          else
            strcat(buffer, "_{");
          strncat(buffer, ptr+2, (ptr1-1)-(ptr+2)+1);
          strcat(buffer, "}");
          ptr = ptr1+1;
        } else
          strncat(buffer, ptr, 1);
        break;
      default:
        fprintf(stdout, "Unrecognized $ sequence: %s\n", ptr);
        fflush(stdout);
        strncat(buffer, ptr, 1);
        break;
      }
    } else 
      strncat(buffer, ptr, 1);
    ptr++;
  }
  strcat(buffer, "$");
  return buffer;
}

void free_beamdata(BEAM *beam)
{
  if (beam->particle)
    free_czarray_2d((void**)beam->particle, beam->n_particle, 7);
  if (beam->accepted)
    free_czarray_2d((void**)beam->accepted, beam->n_particle, 7);
  if (beam->original && beam->original!=beam->particle)
    free_czarray_2d((void**)beam->original, beam->n_original, 7);
  if (beam->lostOnPass)
    free(beam->lostOnPass);
  beam->particle = beam->accepted = beam->original = NULL;
  beam->lostOnPass = NULL;
  beam->n_original = beam->n_to_track = beam->n_accepted = beam->n_saved = beam->n_particle = 0;
  beam->p0_original = beam->p0 =0.;
}  

long getTableFromSearchPath(TABLE *tab, char *file, long sampleInterval, long flags)
{
  char *filename;
  long value;
  if (!(filename=findFileInSearchPath(file)))
    return 0;
  value = get_table(tab, filename, sampleInterval, flags);
  free(filename);
  return value;
}

void do_semaphore_setup(char **semaphoreFile, 
                        NAMELIST_TEXT *nltext)
{
  
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&semaphores, nltext)==NAMELIST_ERROR)
    bombElegant("Problem processing semaphores command", NULL);
  if (echoNamelists) print_namelist(stdout, &semaphores);
 
  semaphoreFile[0] = semaphoreFile[1] = semaphoreFile[2] = NULL;
  if (writePermitted) {
    if (started)
      SDDS_CopyString(&semaphoreFile[0], started);
    if (done)
      SDDS_CopyString(&semaphoreFile[1], done);
    if (failed)
      SDDS_CopyString(&semaphoreFile[2], failed);
  }
}

void getRunControlContext (VARY *context)
{
  *context = run_control;
}

void swapParticles(double *p1, double *p2)
{
  double buffer[7];
  memcpy(buffer,     p1, sizeof(double)*7);
  memcpy(p1    ,     p2, sizeof(double)*7);
  memcpy(p2    , buffer, sizeof(double)*7);
}

void createSemaphoreFile(char *filename)
{
  FILE *fp;
  if (filename) {
    fprintf(stdout, "Creating semaphore file %s\n", filename);
  } else 
    return;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stdout, "Problem creating semaphore file %s\n", filename);
    exitElegant(1);
  }
  fclose(fp);
}

void readApertureInput(NAMELIST_TEXT *nltext, RUN *run)
{
  SDDS_DATASET SDDSin;
  char s[16384];
  long i;
  
#include "aperture_data.h"

  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&aperture_data, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &aperture_data);

  resetApertureData(&(run->apertureData));
  
  run->apertureData.initialized = 0;

  if (disable)
    return;
  
  if (!SDDS_InitializeInputFromSearchPath(&SDDSin, input)) {
    sprintf(s, "Problem opening aperture input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  
  if (!check_sdds_column(&SDDSin, "xHalfAperture", "m") ||
      !check_sdds_column(&SDDSin, "yHalfAperture", "m") ||
      !check_sdds_column(&SDDSin, "s", "m") ||
      !check_sdds_column(&SDDSin, "xCenter", "m") ||
      !check_sdds_column(&SDDSin, "yCenter", "m") ) {
    fprintf(stdout, "Necessary data quantities (s, xHalfAperture, yHalfAperture, xCenter, and yCenter) have wrong units or are not present in %s\n",
            input);
    fprintf(stdout, "Note that units must be \"m\" on all quantities\n");
    fflush(stdout);
    exitElegant(1);
  }

  if (!SDDS_ReadPage(&SDDSin)) {
    sprintf(s, "Problem reading aperture input file %s---seems to be empty", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }

  if ((run->apertureData.points=SDDS_RowCount(&SDDSin))<2) {
    fprintf(stdout, "** Warning: aperture input file %s has only %ld rows!",
            input, run->apertureData.points);
  }

  if (!(run->apertureData.s = SDDS_GetColumnInDoubles(&SDDSin, "s")) ||
      !(run->apertureData.xMax = SDDS_GetColumnInDoubles(&SDDSin, "xHalfAperture")) ||
      !(run->apertureData.yMax = SDDS_GetColumnInDoubles(&SDDSin, "yHalfAperture")) ||
      !(run->apertureData.dx = SDDS_GetColumnInDoubles(&SDDSin, "xCenter")) ||
      !(run->apertureData.dy = SDDS_GetColumnInDoubles(&SDDSin, "yCenter")) ) {
    sprintf(s, "Problem getting data from aperture input file %s", input);
    SDDS_SetError(s);
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
  }
  if (run->apertureData.s[0]!=0) {
    fprintf(stdout, "The first value of s in %s is not zero.\n", input);
    exitElegant(1);
  }
  for (i=0; i<run->apertureData.points; i++) {
    if (i && run->apertureData.s[i]<run->apertureData.s[i-1]) {
      fprintf(stdout, "s values in %s are not monotonically increasing.\n", input);
      exitElegant(1);
    }
    if (run->apertureData.xMax[i]<0 || run->apertureData.yMax[i]<0) {
      fprintf(stdout, "One or more xHalfAperture and yHalfAperture values in %s are negative.\n", input);
      exitElegant(1);
    }
  }
  run->apertureData.periodic = periodic;
  run->apertureData.persistent = persistent;
  run->apertureData.initialized = 1;

  fprintf(stdout, "\n** %ld points of aperture data read from file\n\n", run->apertureData.points);
  
  return;
}

void resetApertureData(APERTURE_DATA *apData)
{
  if (apData->initialized) {
    if (apData->s)
      free(apData->s);
    if (apData->xMax)
      free(apData->xMax);
    if (apData->yMax)
      free(apData->yMax);
    if (apData->dx)
      free(apData->dx);
    if (apData->dy)
      free(apData->dy);
    apData->s = apData->xMax = apData->yMax = apData->dx = apData->dy = NULL;    
    apData->initialized = 0;
  }
}


static long savedRandomNumberSeed ;

void seedElegantRandomNumbers(long seed, long restart)
{
  if (restart)
    seed = savedRandomNumberSeed;
  else
    savedRandomNumberSeed = seed;
  
  /* seed random number generators.  Note that random_1_elegant seeds random_2, random_3,
   * and random_4.
   * random_1_elegant is used for beamline errors.  
   * random_4 is used for beam generation 
   * random_2 is used for random scraping/sampling/scattering.  
   * random_3 is used for BPM noise.
   */
  random_1_elegant(-FABS(seed));
}

/* Routine to close output files that are associated with beamline elements
 */

void closeBeamlineOutputFiles(LINE_LIST *beamline)
{
  ELEMENT_LIST *eptr;
  CSRCSBEND *CsrCsBend;
  
  eptr = &(beamline->elem);
  while (eptr) {
    switch (eptr->type) {
    case T_WATCH:
      if (((WATCH*)(eptr->p_elem))->initialized) {
        SDDS_Terminate(&(((WATCH*)(eptr->p_elem))->SDDS_table));
        ((WATCH*)(eptr->p_elem))->initialized = 0;
      }
      break;
    case T_HISTOGRAM:
      if (((HISTOGRAM*)(eptr->p_elem))->initialized) {
        SDDS_Terminate(&(((HISTOGRAM*)(eptr->p_elem))->SDDS_table));
        ((HISTOGRAM*)(eptr->p_elem))->initialized = 0;
      }
      break;
    case T_CSRCSBEND:
      CsrCsBend = (CSRCSBEND*)(eptr->p_elem);
      if (CsrCsBend->histogramFile) 
        SDDS_Terminate(&(CsrCsBend->SDDSout));
      if (CsrCsBend->particleOutputFile)
        SDDS_Terminate(&(CsrCsBend->SDDSpart));
      break;
    case T_RFMODE:
      if (((RFMODE*)(eptr->p_elem))->record && ((RFMODE*)(eptr->p_elem))->fileInitialized)
        SDDS_Terminate(&(((RFMODE*)(eptr->p_elem))->SDDSrec));
      ((RFMODE*)(eptr->p_elem))->fileInitialized = 0;
      break;
    case T_FRFMODE:
      if (((FRFMODE*)(eptr->p_elem))->outputFile && ((FRFMODE*)(eptr->p_elem))->initialized)
        SDDS_Terminate(&(((FRFMODE*)(eptr->p_elem)))->SDDSout);
      ((FRFMODE*)(eptr->p_elem))->initialized = 0;
      break;
    case T_TRFMODE:
      if (((TRFMODE*)(eptr->p_elem))->record && ((TRFMODE*)(eptr->p_elem))->fileInitialized)
        SDDS_Terminate(&(((TRFMODE*)(eptr->p_elem))->SDDSrec));
      ((TRFMODE*)(eptr->p_elem))->fileInitialized = 0;
      break;
    case T_FTRFMODE:
      if (((FTRFMODE*)(eptr->p_elem))->outputFile && ((FTRFMODE*)(eptr->p_elem))->initialized)
        SDDS_Terminate(&(((FTRFMODE*)(eptr->p_elem)))->SDDSout);
      ((FTRFMODE*)(eptr->p_elem))->initialized = 0;
      break;
    case T_ZLONGIT:
      if (((ZLONGIT*)(eptr->p_elem))->wakes && ((ZLONGIT*)(eptr->p_elem))->SDDS_wake_initialized)
        SDDS_Terminate(&(((ZLONGIT*)(eptr->p_elem)))->SDDS_wake);
      ((ZLONGIT*)(eptr->p_elem))->SDDS_wake_initialized = 0;
      break;
    case T_ZTRANSVERSE:
      if (((ZTRANSVERSE*)(eptr->p_elem))->wakes && ((ZTRANSVERSE*)(eptr->p_elem))->SDDS_wake_initialized)
        SDDS_Terminate(&(((ZTRANSVERSE*)(eptr->p_elem)))->SDDS_wake);
      ((ZTRANSVERSE*)(eptr->p_elem))->SDDS_wake_initialized = 0;      
      break;      
    case T_TFBDRIVER:
      if (((TFBDRIVER*)(eptr->p_elem))->outputFile)
	SDDS_Terminate(&(((TFBDRIVER*)(eptr->p_elem)))->SDDSout);
      break;
    default:
      break;
    }
    eptr = eptr->succ;
  }
}

void setSigmaIndices() 
{
  long i, j, k;
  
  for (i=k=0; i<6; i++)
    for (j=i; j<6; j++, k++) {
      sigmaIndex3[i][j] = sigmaIndex3[j][i] = k;
      sigmaIndex1[k] = i;
      sigmaIndex2[k] = j;
    }
}

#define TYPE_ELECTRON 0
#define TYPE_PROTON   1
#define TYPE_MUON     2
#define TYPE_POSITRON 3
#define TYPE_OTHER    4
#define N_PARTICLE_TYPES 5
static char *particleTypeName[N_PARTICLE_TYPES] = {
  "electron", "proton", "muon", "positron", "custom"
  };
static double charge[N_PARTICLE_TYPES] = {
  -e_mks, e_mks, -e_mks, e_mks, 0
  };
static double mass[N_PARTICLE_TYPES] = {
  me_mks, 1.6726485e-27, 1.88353109e-28, me_mks, 0
  };

void process_particle_command(NAMELIST_TEXT *nltext)
{
  long code, i;
  
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&change_particle, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &change_particle);

  code = match_string(change_particle_struct.name, particleTypeName, N_PARTICLE_TYPES, EXACT_MATCH);
  
  if (code==TYPE_ELECTRON) {
    particleIsElectron = 1;
    /* ensure exact backward compatibility */
    particleMass = me_mks;
    particleCharge = e_mks;
    particleMassMV = me_mev;
    particleRadius = re_mks;
    particleRelSign = 1;
  } else if (code==TYPE_POSITRON) {
    particleIsElectron = 0;
    /* ensure exact reversal from electron */
    particleMass = me_mks;
    particleCharge = -e_mks;
    particleMassMV = me_mev;
    particleRadius = re_mks;
    particleRelSign = -1;
  } else {
    particleIsElectron = 0;
    switch (code) {
    case TYPE_PROTON:
    case TYPE_MUON:
      particleMass = mass[code];
      /* minus sign is a legacy of time when this was an electron-only code */
      particleCharge = -charge[code];  
      break;
    case TYPE_OTHER:
      if (change_particle_struct.mass_ratio<=0 || change_particle_struct.charge_ratio==0)
        bombElegant("Must have mass_ratio>0 and charge_ratio nonzero", NULL);
      particleMass = mass[TYPE_ELECTRON]*change_particle_struct.mass_ratio;
      /* minus sign is a legacy of time when this was an electron-only code */
      particleCharge = -charge[TYPE_ELECTRON]*change_particle_struct.charge_ratio;
      break;
    default:
      fprintf(stderr, "Unknown particle type.  Known types are \n");
      for (i=0; i<N_PARTICLE_TYPES; i++)
        fprintf(stderr, "%s%c", particleTypeName[i], i==(N_PARTICLE_TYPES-1)?'\n':' ');
      exitElegant(1);
      break;
    }
    particleRelSign = -SIGN(particleCharge)/SIGN(charge[TYPE_ELECTRON]);
    particleMassMV = particleMass*sqr(c_mks)/fabs(particleCharge)/1e6;
    particleRadius = sqr(particleCharge)/(4*PI*epsilon_o*particleMass*sqr(c_mks));
  }
  /* 
    printf("particleMass = %e, particleRadius = %e, particleCharge = %e, particleMassMV = %e, particleRelSign = %e\n",
    particleMass, particleRadius, particleCharge, particleMassMV, particleRelSign);
    */
  printf("************************************************\n");
  printf("* WARNING ! \n");
  printf("* Changing the particle type is not a fully tested\n");
  printf("* feature.\n");
  printf("* Please be alert for results that don't make sense!\n");
  printf("************************************************\n");
}

void processGlobalSettings(NAMELIST_TEXT *nltext)
{
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&global_settings, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &global_settings);

  inhibitFileSync = inhibit_fsync;
  echoNamelists = echo_namelists;
  if (log_file)
    freopen(log_file, "w", stdout);
  if (error_log_file)
    freopen(error_log_file, "w", stderr);
}

void bombElegant(char *error, char *usage)
{
  if (error)
    fprintf(stderr, "error: %s\n", error);
  if (usage)
    fprintf(stderr, "usage: %s\n", usage);
  if (semaphoreFile[2]) 
    createSemaphoreFile(semaphoreFile[2]);
  exit(1);
}

void exitElegant(long status)
{
  if (status && semaphoreFile[2]) 
    createSemaphoreFile(semaphoreFile[2]);
  if (!status && semaphoreFile[1])
    createSemaphoreFile(semaphoreFile[1]);
  if (!status && semaphore_file)
    createSemaphoreFile(semaphore_file);
  exit(status);
}
