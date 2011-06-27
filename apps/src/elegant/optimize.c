/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: optimize.c
 *
 * Michael Borland, 1991
 */
#include "mdb.h"
#include "track.h"
#if defined(__BORLANDC__)
#define DBL_MAX         1.7976931348623158e+308 /* max value */
#endif
#include "optimize.h"
#include "match_string.h"
#include "chromDefs.h"
#include "tuneDefs.h"
#include "correctDefs.h"

static long stopOptimization = 0;
long checkForOptimRecord(double *value, long values, long *again);
void storeOptimRecord(double *value, long values, long invalid, double result);
void rpnStoreHigherMatrixElements(VMATRIX *M, long **TijkMem, long **UijklMem, long maxOrder);

#if USE_MPI
/* Find the global minimal value and its location across all the processors */
void find_global_min_index (double *min, int *processor_ID, MPI_Comm comm);
/* Print out the best individual and value after each iteration */
void SDDS_PrintPopulations(SDDS_TABLE *popLogPtr, long iteration, double best_value, double *variable, long n_variables, double *covariable, long n_covariables);
#endif

#if !USE_MPI
long myid=0;
#endif

void do_optimization_setup(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{

    log_entry("do_optimization_setup");

    /* process the namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&optimization_setup, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &optimization_setup);

    /* check validity of input values, and copy into structure */
    if ((optimization_data->mode=match_string(mode, optimize_mode, N_OPTIM_MODES, EXACT_MATCH))<0)
        bombElegant("unknown optimization mode", NULL);
    optimization_data->equation = equation;
    if ((optimization_data->method=match_string(method, optimize_method, N_OPTIM_METHODS, EXACT_MATCH))<0)
        bombElegant("unknown optimization method", NULL);
    if ((optimization_data->tolerance=tolerance)==0)
        bombElegant("tolerance == 0", NULL);
    if ((optimization_data->n_passes=n_passes)<=0)
        bombElegant("n_passes <= 0", NULL);
    if ((optimization_data->n_evaluations=n_evaluations)<=0)
        bombElegant("n_evaluations <= 0", NULL);
    if ((optimization_data->n_restarts = n_restarts)<0)
      bombElegant("n_restarts < 0", NULL);
    if ((optimization_data->matrix_order=matrix_order)<1 ||
        matrix_order>3)
      bombElegant("matrix_order must be 1, 2, or 3", NULL);
    optimization_data->soft_failure = soft_failure;
    if (output_sparsing_factor<=0)
      output_sparsing_factor = 1;
#if USE_MPI
    if (!writePermitted)
      log_file = NULL; 
    optimization_data->random_factor = random_factor; 
#endif 
   if (log_file) {
        if (str_in(log_file, "%s"))
            log_file = compose_filename(log_file, run->rootname);
        if (strcmp(log_file, "/dev/tty")==0 || strcmp(log_file, "tt:")==0)
            optimization_data->fp_log = stdout;
        else if ((optimization_data->fp_log=fopen_e(log_file, "w", FOPEN_RETURN_ON_ERROR))==NULL)
            bombElegant("unable to open log file", NULL);
        }
    if (term_log_file) {
      if (str_in(term_log_file, "%s"))
        term_log_file = compose_filename(term_log_file, run->rootname);
    }
    optimization_data->verbose = verbose;
    if (optimization_data->mode==OPTIM_MODE_MAXIMUM && target!=-DBL_MAX)
      target = -target;
    optimization_data->target = target;
    optimization_data->simplexPassRangeFactor = simplex_pass_range_factor;
    optimization_data->simplexDivisor = simplex_divisor;
    optimization_data->includeSimplex1dScans = include_simplex_1d_scans;
    optimization_data->startFromSimplexVertex1 = start_from_simplex_vertex1;
    if ((optimization_data->restart_worst_term_factor = restart_worst_term_factor)<=0)
      bombElegant("restart_worst_term_factor <= 0", NULL);
    if ((optimization_data->restart_worst_terms=restart_worst_terms)<=0)
      bombElegant("restart_worst_terms <= 0", NULL);
    
    /* reset flags for elements that may have been varied previously */
    if (optimization_data->variables.n_variables)
        set_element_flags(beamline, optimization_data->variables.element, NULL, NULL, 
                          NULL, optimization_data->variables.n_variables,
                          PARAMETERS_ARE_STATIC, 0, 1, 0);

    /* initialize other elements of the structure */
    optimization_data->new_data_read = 0;
    optimization_data->balance_terms = balance_terms;
    optimization_data->variables.n_variables = 0;
    optimization_data->covariables.n_covariables = 0;
    optimization_data->constraints.n_constraints = 0;
    optimization_data->TijkMem = NULL;
    optimization_data->UijklMem = NULL;
    log_exit("do_optimization_setup");
    }

#if USE_MPI
void do_parallel_optimization_setup(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{ 
  runInSinglePartMode = 1;  /* All the processors will track the same particles with different parameters */

  do_optimization_setup(optimization_data, nltext, run, beamline);  
  if (optimization_data->method!=OPTIM_METHOD_HYBSIMPLEX) {
    if (population_size < n_processors) {
      fprintf (stderr, "The populaition size (%ld) can not be less than the number of processors (%d).\n",
               population_size, n_processors);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
 
  if (optimization_data->method==OPTIM_METHOD_SWARM) /* The n_restarts is used to control n_iterations for particle swarm optimization */ 
    optimization_data->n_restarts = n_iterations-1;
  optimization_data->n_iterations = n_iterations;
  optimization_data->max_no_change = max_no_change;
  optimization_data->population_size = population_size;
  optimization_data->print_all_individuals = print_all_individuals;
  if (population_log && strlen(population_log)) {
    if (str_in(population_log, "%s"))
      population_log = compose_filename(population_log, run->rootname);
  }
}
#endif

void add_optimization_variable(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long n_variables;
    OPTIM_VARIABLES *variables;
    ELEMENT_LIST *context;
    /* these are used to append a dummy name to the variables list for use with final parameters output: */
    static char *extra_name[3] = {"optimized", "optimizationFunction", "bestOptimizationFunction"};
    static char *extra_unit[3] = {"", "", ""};
    long i, extras = 3;
    char *ptr;
    
    log_entry("add_optimization_variable");

    /* process namelist text */
    /* can't use automatic defaults, because of DBL_MAX being a nonconstant object */
    set_namelist_processing_flags(0);
    set_print_namelist_flags(0);
    name = item = NULL;
    step_size = 1;
    lower_limit = -(upper_limit = DBL_MAX);
    if (processNamelist(&optimization_variable, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &optimization_variable);

    if (disable)
      return;
    
    if ((n_variables = optimization_data->variables.n_variables)==0) {
        if (optimization_data->new_data_read)
            bombElegant("improper sequencing of variation and tracking", NULL);
        optimization_data->new_data_read = 1;
        }

    variables = &(optimization_data->variables);
    variables->element = trealloc(variables->element, sizeof(*variables->element)*(n_variables+extras+1));
    variables->item = trealloc(variables->item, sizeof(*variables->item)*(n_variables+extras+1));
    variables->lower_limit = trealloc(variables->lower_limit, sizeof(*variables->lower_limit)*(n_variables+extras+1));
    variables->upper_limit = trealloc(variables->upper_limit, sizeof(*variables->upper_limit)*(n_variables+extras+1));
    variables->step = trealloc(variables->step, sizeof(*variables->step)*(n_variables+extras+1));
    variables->orig_step = trealloc(variables->orig_step, sizeof(*variables->orig_step)*(n_variables+extras+1));
    variables->varied_quan_name = trealloc(variables->varied_quan_name, sizeof(*variables->varied_quan_name)*(n_variables+extras+1));
    variables->varied_quan_unit = trealloc(variables->varied_quan_unit, sizeof(*variables->varied_quan_unit)*(n_variables+extras+1));
    variables->varied_type = trealloc(variables->varied_type, sizeof(*variables->varied_type)*(n_variables+extras+1));
    variables->varied_param = trealloc(variables->varied_param, sizeof(*variables->varied_param)*(n_variables+extras+1));
    variables->varied_quan_value = trealloc(variables->varied_quan_value, sizeof(*variables->varied_quan_value)*(n_variables+extras+1));
    variables->initial_value = trealloc(variables->initial_value, sizeof(*variables->initial_value)*(n_variables+extras+1));
    variables->memory_number = trealloc(variables->memory_number, sizeof(*variables->memory_number)*(n_variables+extras+1));

    /* check for valid input */
    if (name==NULL)
        bombElegant("element name missing in optimization_variable namelist", NULL);
    str_toupper(name);
    context = NULL;
    if (!find_element(name, &context, &(beamline->elem))) {
        fprintf(stdout, "error: cannot vary element %s--not in beamline\n", name);
        fflush(stdout);
        exitElegant(1);
        }
    cp_str(&variables->element[n_variables], name);
    variables->varied_type[n_variables] = context->type;
    if (item==NULL)
        bombElegant("item name missing in optimization_variable list", NULL);
    str_toupper(item);
    if ((variables->varied_param[n_variables] = confirm_parameter(item, context->type))<0) {
        fprintf(stdout, "error: cannot vary %s--no such parameter for %s\n",item, name);
        fflush(stdout);
        exitElegant(1);
        }
    cp_str(&variables->item[n_variables], item);
    cp_str(&variables->varied_quan_unit[n_variables], 
        entity_description[context->type].parameter[variables->varied_param[n_variables]].unit);
    if (!get_parameter_value(variables->varied_quan_value+n_variables, name, variables->varied_param[n_variables],
            context->type, beamline))
        bombElegant("unable to get initial value for parameter", NULL);
    if (lower_limit>=upper_limit)
        bombElegant("lower_limit >= upper_limit", NULL);

    variables->initial_value[n_variables] = variables->varied_quan_value[n_variables];
    if (variables->initial_value[n_variables]>upper_limit) {
      if (force_inside) {
	fprintf(stdout, "Warning: Initial value (%e) is greater than upper limit\n", variables->initial_value[n_variables]);
	variables->initial_value[n_variables] = variables->varied_quan_value[n_variables] = upper_limit;
      } else {
	fprintf(stdout, "Error: Initial value (%e) is greater than upper limit\n", variables->initial_value[n_variables]);
	exit(1);
      }
    }
    if (variables->initial_value[n_variables]<lower_limit) {
      if (force_inside) {
	fprintf(stdout, "Warning: Initial value (%e) is smaller than lower limit\n", variables->initial_value[n_variables]);
	variables->initial_value[n_variables] = variables->varied_quan_value[n_variables] = lower_limit;
      } else {
	fprintf(stdout, "Error: Initial value (%e) is smaller than lower limit\n", variables->initial_value[n_variables]);
	exit(1);
      }
    }
    variables->orig_step[n_variables]     = step_size;
    variables->varied_quan_name[n_variables]  = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+3));
    sprintf(variables->varied_quan_name[n_variables], "%s.%s", name, item);
    variables->lower_limit[n_variables] = lower_limit;
    variables->upper_limit[n_variables] = upper_limit;
    rpn_store(variables->initial_value[n_variables], NULL, 
              variables->memory_number[n_variables] = 
              rpn_create_mem(variables->varied_quan_name[n_variables], 0)); 
  
    ptr = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+4));
    sprintf(ptr, "%s.%s0", name, item);
    rpn_store(variables->initial_value[n_variables], NULL, rpn_create_mem(ptr, 0));
    free(ptr);
    
    for (i=0; i<extras; i++) {
      variables->varied_quan_name[n_variables+i+1] = extra_name[i];
      variables->varied_quan_unit[n_variables+i+1] = extra_unit[i];
    }
    
    optimization_data->variables.n_variables += 1;
    log_exit("add_optimization_variable");
    }

void add_optimization_term(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run,
                           LINE_LIST *beamline)
{
  long field_value=0, n_field_values = 0, index, nFileTerms = 0;
  char s[16834], value[100];
  char **fileTerm = NULL;
  
  /* process namelist text */
  set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
  set_print_namelist_flags(0);
  if (processNamelist(&optimization_term, nltext)==NAMELIST_ERROR)
    bombElegant(NULL, NULL);
  if (echoNamelists) print_namelist(stdout, &optimization_term);

  if (optimization_term_struct.weight==0)
    return ;
  
  if (optimization_term_struct.term==NULL && optimization_term_struct.input_file==NULL)
    bombElegant("term is invalid and no input file given", NULL);
  if (optimization_data->equation)
    bombElegant("you've already given an optimization equation, so you can't give individual terms", NULL);
  if (optimization_term_struct.field_string && optimization_term_struct.input_file)
    bombElegant("you can't use field substitution and file input together", NULL);
  if (optimization_term_struct.field_string) {
    if ((n_field_values = (optimization_term_struct.field_final_value-optimization_term_struct.field_initial_value)/optimization_term_struct.field_interval+1)<=0) 
      bombElegant("something strage about field_final_value, field_initial_value, and field_interval", NULL);
    field_value = optimization_term_struct.field_initial_value;
  } else if (optimization_term_struct.input_file) {
    SDDS_DATASET SDDSin;
    if (!optimization_term_struct.input_column)
      bombElegant("you must give input_column when giving input_file", NULL);
    if (!SDDS_InitializeInputFromSearchPath(&SDDSin, optimization_term_struct.input_file) ||
        SDDS_ReadPage(&SDDSin)!=1)
      SDDS_Bomb("problem reading optimization term input file");
    if ((nFileTerms=SDDS_RowCount(&SDDSin))==0) {
      printf("Warning: No rows loaded from file!\n");
      return ;
    }
    if (SDDS_GetNamedColumnType(&SDDSin, optimization_term_struct.input_column)!=SDDS_STRING) {
      printf("Error: column %s is nonexistent or not string type.\n", optimization_term_struct.input_column);
      exitElegant(1);
    }
    if (!(fileTerm = SDDS_GetColumn(&SDDSin, optimization_term_struct.input_column))) 
      SDDS_Bomb("error getting optimization term column from file");
    SDDS_Terminate(&SDDSin);
  } else {
    optimization_term_struct.field_initial_value = optimization_term_struct.field_final_value = 0;
    optimization_term_struct.field_interval = 1;
    n_field_values = 1;
    field_value = 0;
    nFileTerms = 0;
  }

  if (!(optimization_data->term 
        = SDDS_Realloc(optimization_data->term,
                       sizeof(*optimization_data->term)*(optimization_data->terms+n_field_values+nFileTerms))) ||
      !(optimization_data->termValue 
        = SDDS_Realloc(optimization_data->termValue,
                       sizeof(*optimization_data->termValue)*(optimization_data->terms+n_field_values+nFileTerms))) ||
      !(optimization_data->termWeight 
        = SDDS_Realloc(optimization_data->termWeight,
                       sizeof(*optimization_data->termWeight)*(optimization_data->terms+n_field_values+nFileTerms))) ||
      !(optimization_data->usersTermWeight 
        = SDDS_Realloc(optimization_data->usersTermWeight,
                       sizeof(*optimization_data->usersTermWeight)*(optimization_data->terms+n_field_values+nFileTerms))))
    bombElegant("memory allocation failure", NULL);

  if (!nFileTerms) {
    for (index=0; index<n_field_values; index++) {
      if (optimization_term_struct.field_string && strlen(optimization_term_struct.field_string)) {
        sprintf(value, "%ld", field_value);
        replaceString(s, optimization_term_struct.term, optimization_term_struct.field_string, value, -1, 0);
        field_value += optimization_term_struct.field_interval;
      } else 
        strcpy(s, optimization_term_struct.term);
      if (!SDDS_CopyString(&optimization_data->term[optimization_data->terms+index], s))
        bombElegant("memory allocation failure", NULL);
      optimization_data->termWeight[optimization_data->terms+index] = 1;
      optimization_data->usersTermWeight[optimization_data->terms+index] = optimization_term_struct.weight;
      if (optimization_term_struct.verbose)
        printf("Added term: %s\n", optimization_data->term[optimization_data->terms+index]);
    }
    optimization_data->terms += n_field_values;
  } else {
    for (index=0; index<nFileTerms; index++) {
      optimization_data->term[optimization_data->terms+index] = fileTerm[index];
      optimization_data->termWeight[optimization_data->terms+index] = 1;
      optimization_data->usersTermWeight[optimization_data->terms+index] = optimization_term_struct.weight;
      if (optimization_term_struct.verbose)
        printf("Added term: %s\n", optimization_data->term[optimization_data->terms+index]);
    }
    optimization_data->terms += nFileTerms;
    free(fileTerm);
  }
  
}

void add_optimization_covariable(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
#include "optim_covariable.h"
    long n_covariables;
    OPTIM_COVARIABLES *covariables;
    ELEMENT_LIST *context;
    char nameBuffer[100], *ptr;
    static long nameIndex = 0;
    
    log_entry("add_optimization_covariable");

    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&optimization_covariable, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &optimization_covariable);
    if (disable)
      return ;
    
    n_covariables = optimization_data->covariables.n_covariables;

    covariables = &(optimization_data->covariables);
    covariables->element = trealloc(covariables->element, sizeof(*covariables->element)*(n_covariables+1));
    covariables->item = trealloc(covariables->item, sizeof(*covariables->item)*(n_covariables+1));
    covariables->equation = trealloc(covariables->equation, sizeof(*covariables->equation)*(n_covariables+1));
    covariables->pcode = trealloc(covariables->pcode, sizeof(*covariables->pcode)*(n_covariables+1));
    covariables->varied_quan_name = trealloc(covariables->varied_quan_name, sizeof(*covariables->varied_quan_name)*(n_covariables+1));
    covariables->varied_quan_unit = trealloc(covariables->varied_quan_unit, sizeof(*covariables->varied_quan_unit)*(n_covariables+1));
    covariables->varied_type = trealloc(covariables->varied_type, sizeof(*covariables->varied_type)*(n_covariables+1));
    covariables->varied_param = trealloc(covariables->varied_param, sizeof(*covariables->varied_param)*(n_covariables+1));
    covariables->varied_quan_value = trealloc(covariables->varied_quan_value, sizeof(*covariables->varied_quan_value)*(n_covariables+1));
    covariables->memory_number = trealloc(covariables->memory_number, sizeof(*covariables->memory_number)*(n_covariables+1));

    /* check for valid input */
    if (name==NULL)
        bombElegant("element name missing in optimization_variable namelist", NULL);
    str_toupper(name);
    context = NULL;
    if (!find_element(name, &context, &(beamline->elem))) {
        fprintf(stdout, "error: cannot vary element %s--not in beamline\n", name);
        fflush(stdout);
        exitElegant(1);
        }
    cp_str(&covariables->element[n_covariables], name);
    covariables->varied_type[n_covariables] = context->type;
    if (item==NULL)
        bombElegant("item name missing in optimization_variable list", NULL);
    str_toupper(item);
    if ((covariables->varied_param[n_covariables] = confirm_parameter(item, context->type))<0) {
        fprintf(stdout, "error: cannot vary %s--no such parameter for %s\n",item, name);
        fflush(stdout);
        exitElegant(1);
        }
    cp_str(&covariables->item[n_covariables], item);
    cp_str(&covariables->varied_quan_unit[n_covariables], 
        entity_description[context->type].parameter[covariables->varied_param[n_covariables]].unit);

    if (!get_parameter_value(covariables->varied_quan_value+n_covariables, name, covariables->varied_param[n_covariables],
            context->type, beamline))
        bombElegant("unable to get initial value for parameter", NULL);

    ptr = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+4));
    sprintf(ptr, "%s.%s0", name, item);
    rpn_store(covariables->varied_quan_value[n_covariables], NULL, rpn_create_mem(ptr, 0));
    free(ptr);
    
    covariables->varied_quan_name[n_covariables]  = tmalloc(sizeof(char)*(strlen(name)+strlen(item)+3));
    sprintf(covariables->varied_quan_name[n_covariables], "%s.%s", name, item);
    cp_str(&covariables->equation[n_covariables], equation);
    rpn_store(covariables->varied_quan_value[n_covariables] = rpn(equation), NULL,
              covariables->memory_number[n_covariables] = rpn_create_mem(covariables->varied_quan_name[n_covariables], 0) );
    if (rpn_check_error()) exitElegant(1);
    fprintf(stdout, "Initial value of %s is %e %s\n", 
            covariables->varied_quan_name[n_covariables], covariables->varied_quan_value[n_covariables],
            covariables->varied_quan_unit[n_covariables]);
    fflush(stdout);

    sprintf(nameBuffer, "AOCEqn%ld", nameIndex++);
    create_udf(nameBuffer, equation);
    if (!SDDS_CopyString(&ptr, nameBuffer)) 
      SDDS_PrintErrors(stdout, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    covariables->pcode[n_covariables] = ptr;
    optimization_data->covariables.n_covariables += 1;
    log_exit("add_optimization_covariable");
    }

void add_optimization_constraint(OPTIMIZATION_DATA *optimization_data, NAMELIST_TEXT *nltext, RUN *run, LINE_LIST *beamline)
{
    long n_constraints;
    OPTIM_CONSTRAINTS *constraints;

    log_entry("add_optimization_constraint");

    constraints = &(optimization_data->constraints);
    n_constraints = constraints->n_constraints;

    constraints->quantity   = trealloc(constraints->quantity, sizeof(*constraints->quantity)*(n_constraints+1));
    constraints->index      = trealloc(constraints->index, sizeof(*constraints->index)*(n_constraints+1));
    constraints->lower      = trealloc(constraints->lower, sizeof(*constraints->lower)*(n_constraints+1));
    constraints->upper      = trealloc(constraints->upper, sizeof(*constraints->upper)*(n_constraints+1));

/*
    quantity = NULL;
    lower = upper = 0;
 */
    
    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&optimization_constraint, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &optimization_constraint);

    /* check for valid input */
    if (quantity==NULL)
        bombElegant("quantity name missing in optimization_constraint namelist", NULL);
    constraints->quantity[n_constraints] = quantity;
    if ((constraints->index[n_constraints]=get_final_property_index(quantity))<0)
        bombElegant("no match for quantity to be constrained", NULL);
    if (lower>upper)
        bombElegant("lower > upper for constraint", NULL);
    constraints->lower[n_constraints] = lower;
    constraints->upper[n_constraints] = upper;

    constraints->n_constraints += 1;
    log_exit("add_optimization_constraint");
    }


void summarize_optimization_setup(OPTIMIZATION_DATA *optimization_data)
{
    OPTIM_VARIABLES *variables;
    OPTIM_CONSTRAINTS *constraints;
    OPTIM_COVARIABLES *covariables;
    long i;

    log_entry("summarize_optimization_setup");

    constraints = &(optimization_data->constraints);
    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);

    if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "\nOptimization to be performed using method '%s' in mode '%s' with tolerance %e.\n", 
            optimize_method[optimization_data->method], optimize_mode[optimization_data->mode], optimization_data->tolerance);
#if USE_MPI
	if (runInSinglePartMode) /* For genetic optimization */
	  fprintf(optimization_data->fp_log, "    As many as %ld generations will be performed on each of %ld passes.\n",
             optimization_data->n_iterations, optimization_data->n_passes);
	else
#endif
        fprintf(optimization_data->fp_log, "    As many as %ld function evaluations will be performed on each of %ld passes.\n",
            optimization_data->n_evaluations, optimization_data->n_passes);

        if (optimization_data->terms==0)
          fprintf(optimization_data->fp_log, "    The quantity to be optimized is defined by the equation\n        %s\n    subject to %ld constraints%c\n",
                  optimization_data->equation, constraints->n_constraints,
                  constraints->n_constraints>0?':':'.');
        else
          fprintf(optimization_data->fp_log, "    The quantity to be optimized is defined by the sum of %ld terms,\nsubject to %ld constraints%c\n",
                  optimization_data->terms, constraints->n_constraints,
                  constraints->n_constraints>0?':':'.');
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(optimization_data->fp_log, "        %13.6e <= %10s <= %13.6e\n",
                constraints->lower[i], constraints->quantity[i], constraints->upper[i]);

        fprintf(optimization_data->fp_log, "    The following variables will be used in the optimization:\n");
        fprintf(optimization_data->fp_log, "    name      initial value   lower limit    upper limit     step size\n");
        fprintf(optimization_data->fp_log, "--------------------------------------------------------------------------\n");
        for (i=0; i<variables->n_variables; i++) {
            fprintf(optimization_data->fp_log, "%12s  %13.6e", variables->varied_quan_name[i], variables->initial_value[i]);
            if (variables->lower_limit[i]==variables->upper_limit[i])
                fprintf(optimization_data->fp_log, "        -              -      ");
            else
                fprintf(optimization_data->fp_log, "  %13.6e  %13.6e", variables->lower_limit[i], variables->upper_limit[i]);
            if (variables->orig_step[i]==0)
                fprintf(optimization_data->fp_log, "        -\n");
            else
                fprintf(optimization_data->fp_log, "  %13.6e\n", variables->orig_step[i]);
            }
        fputc('\n', optimization_data->fp_log);
        if (covariables->n_covariables) {
            fprintf(optimization_data->fp_log, "    The following covariables will be used in the optimization:\n");
            fprintf(optimization_data->fp_log, "    name      equation\n");
            fprintf(optimization_data->fp_log, "--------------------------------------------------------------------------\n");
            for (i=0; i<covariables->n_covariables; i++)
                fprintf(optimization_data->fp_log, "%12s  %s\n", covariables->varied_quan_name[i], covariables->equation[i]);
            fputc('\n', optimization_data->fp_log);
            }
        fputc('\n', optimization_data->fp_log);
        fflush(optimization_data->fp_log);
        }
    if (!log_file || optimization_data->fp_log!=stdout) {
        fprintf(stdout, "\nOptimization to be performed using method '%s' in mode '%s' with tolerance %e.\n", 
            optimize_method[optimization_data->method], optimize_mode[optimization_data->mode], optimization_data->tolerance);
        fflush(stdout);
        fprintf(stdout, "    As many as %ld function evaluations will be performed on each of %ld passes.\n",
            optimization_data->n_evaluations, optimization_data->n_passes);
        fflush(stdout);
         fprintf(stdout, "    The quantity to be optimized is defined by the equation\n        %s\n    subject to %ld constraints%c\n",
            optimization_data->equation, constraints->n_constraints,
            constraints->n_constraints>0?':':'.');
         fflush(stdout);
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(stdout, "        %13.6e <= %10s <= %13.6e\n",
                constraints->lower[i], constraints->quantity[i], constraints->upper[i]);
            fflush(stdout);
        fprintf(stdout, "    The following variables will be used in the optimization:\n");
        fflush(stdout);
        fprintf(stdout, "    name      initial value   lower limit    upper limit     step size\n");
        fflush(stdout);
        fprintf(stdout, "--------------------------------------------------------------------------\n");
        fflush(stdout);
        for (i=0; i<variables->n_variables; i++) {
            fprintf(stdout, "%12s  %13.6e", variables->varied_quan_name[i], variables->initial_value[i]);
            fflush(stdout);
            if (variables->lower_limit[i]==variables->upper_limit[i]) {
                fprintf(stdout, "        -              -      ");
                fflush(stdout);
            }  else {
                fprintf(stdout, "  %13.6e  %13.6e", variables->lower_limit[i], variables->upper_limit[i]);
                fflush(stdout);
            }
            if (variables->orig_step[i]==0) {
                fprintf(stdout, "        -\n");
                fflush(stdout);
            } else {
                fprintf(stdout, "  %13.6e\n", variables->orig_step[i]);
                fflush(stdout);
              }
          }
        fputc('\n', stdout);
        if (covariables->n_covariables) {
            fprintf(stdout, "    The following covariables will be used in the optimization:\n");
            fflush(stdout);
            fprintf(stdout, "    name      equation\n");
            fflush(stdout);
            fprintf(stdout, "--------------------------------------------------------------------------\n");
            fflush(stdout);
            for (i=0; i<covariables->n_covariables; i++)
                fprintf(stdout, "%12s  %s\n", covariables->varied_quan_name[i], covariables->equation[i]);
                fflush(stdout);
            fputc('\n', stdout);
            }
        fputc('\n', stdout);
        }
    log_exit("summarize_optimization_setup");
    }

/* variables needed to do tracking for optimization--this data has to be held in global
 * variables like this since the optimization function has a simple calling syntax
 */
static RUN *run;
static VARY *control;
static ERRORVAL *error;
static BEAM *beam;
static OPTIMIZATION_DATA *optimization_data;
static OUTPUT_FILES *output;
static LINE_LIST *beamline;
static CHROM_CORRECTION *chromCorrData;
static void *orbitCorrData;
static long orbitCorrMode;
static TUNE_CORRECTION *tuneCorrData;
static long beam_type_code, n_evaluations_made, n_passes_made;
static double *final_property_value;
static long final_property_values;
static double charge;
static unsigned long optim_func_flags;
static long force_output;
static long doClosedOrbit, doChromCorr, doTuneCorr, doFindAperture, doResponse;

/* structure to keep results of last N optimization function
 * evaluations, so we don't track the same thing twice.
 */
#define MAX_OPTIM_RECORDS 50
typedef struct {
  long invalid;
  double *variableValue;
  double result;
  long usedBefore;
} OPTIM_RECORD;
static long optimRecords = 0, nextOptimRecordSlot = 0, balanceTerms = 0, ignoreOptimRecords=0;
static OPTIM_RECORD optimRecord[MAX_OPTIM_RECORDS];
static double bestResult = DBL_MAX;

#if defined(UNIX)
#include <signal.h>
void traceback_handler(int signal);

void optimizationInterruptHandler(int signal)
{
  simplexMinAbort(1);
  optimAbort(1);
  fprintf(stdout, "Aborting optimization...");
}
#endif

int variableAtLimit(double value, double lower, double upper)
{
  double range = upper-lower;
  if (range==0 || (value-lower)/range<0.01 || (upper-value)/range<0.01)
    return 1;
  return 0;
}

void do_optimize(NAMELIST_TEXT *nltext, RUN *run1, VARY *control1, ERRORVAL *error1, 
                 LINE_LIST *beamline1, BEAM *beam1, OUTPUT_FILES *output1, 
                 OPTIMIZATION_DATA *optimization_data1, 
                 void *chromCorrData1, long beam_type1,
                 long doClosedOrbit1, long doChromCorr1,
                 void *orbitCorrData1, long orbitCorrMode1,
                 void *tuneCorrData1, long doTuneCorr1,
                 long doFindAperture1, long doResponse1)
{
    static long optUDFcount = 0;
    void optimization_report(double result, double *value, long pass, long n_evals, long n_dim);
    void (*optimization_report_ptr) (double result, double *value, long pass, long n_evals, long n_dim) = optimization_report ;	
    OPTIM_VARIABLES *variables;
    OPTIM_COVARIABLES *covariables;
    OPTIM_CONSTRAINTS *constraints;
    double result, lastResult;
    long i, startsLeft, i_step_saved;
#if USE_MPI
    long n_total_evaluations_made = 0;
    double scale_factor = optimization_data1->random_factor;
    int min_location = 0;
    static double *covariables_global = NULL;
#endif
    
    log_entry("do_optimize");

    optimRecords = ignoreOptimRecords = nextOptimRecordSlot = 0;
    bestResult = DBL_MAX;
    
    stopOptimization  = 0;
    run               = run1;
    control           = control1;
    error             = error1;
    beamline          = beamline1;
    beam              = beam1;
    output            = output1,
    optimization_data = optimization_data1;
    chromCorrData     = (CHROM_CORRECTION*)chromCorrData1;
    beam_type_code    = beam_type1;
    doClosedOrbit     = doClosedOrbit1;
    doChromCorr       = doChromCorr1;
    orbitCorrData     = orbitCorrData1;
    orbitCorrMode     = orbitCorrMode1;
    tuneCorrData      = (TUNE_CORRECTION*)tuneCorrData1;
    doTuneCorr        = doTuneCorr1;
    doFindAperture    = doFindAperture1;
    doResponse        = doResponse1;
    
    n_evaluations_made = 0;
    n_passes_made      = 0;
    i_step_saved = control->i_step;

    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);
    constraints = &(optimization_data->constraints);

    if (variables->n_variables==0)
        bombElegant("no variables specified for optimization", NULL);

    for (i=0; i<MAX_OPTIM_RECORDS; i++)
      optimRecord[i].variableValue = tmalloc(sizeof(*optimRecord[i].variableValue)*
                                             variables->n_variables);
    
    /* set the end-of-optimization hidden variable to 0 */
    variables->varied_quan_value[variables->n_variables] = 0;
    /* set the optimization function hidden variable to 0 */
    variables->varied_quan_value[variables->n_variables+1] = 0;
    
    if (optimization_data->equation==NULL || !strlen(optimization_data->equation)) {
      long i, length;
      if (optimization_data->terms==0)
        bombElegant("give an optimization equation or at least one optimization term", NULL);
      for (i=length=0; i<optimization_data->terms; i++)
        length += strlen(optimization_data->term[i])+4;
      if (!(optimization_data->equation = SDDS_Malloc(sizeof(*optimization_data->equation)*length)))
        bombElegant("memory allocation failue", NULL);
      optimization_data->equation[0] = '\0';
      for (i=0; i<optimization_data->terms; i++) {
        strcat(optimization_data->equation, optimization_data->term[i]);
        if (i)
          strcat(optimization_data->equation, " + ");
        else
          strcat(optimization_data->equation, " ");
      }
    }
    
    /* process namelist text */
    set_namelist_processing_flags(STICKY_NAMELIST_DEFAULTS);
    set_print_namelist_flags(0);
    if (processNamelist(&optimize, nltext)==NAMELIST_ERROR)
      bombElegant(NULL, NULL);
    if (echoNamelists) print_namelist(stdout, &optimize);

    if (summarize_setup)
        summarize_optimization_setup(optimization_data);

    for (i=0; i<variables->n_variables; i++) {
        if (!get_parameter_value(variables->varied_quan_value+i, 
				 variables->element[i], 
				 variables->varied_param[i], 
				 variables->varied_type[i], beamline))
            bombElegant("unable to get initial value for parameter", NULL);
	variables->initial_value[i] = variables->varied_quan_value[i];
	if (variables->initial_value[i]>variables->upper_limit[i]) {
	  if (force_inside) 
	    variables->varied_quan_value[i] = variables->initial_value[i] = variables->upper_limit[i];
	  else {
	    fprintf(stdout, "Initial value (%e) is greater than upper limit for %s.%s\n", 
		    variables->initial_value[i],
		    variables->element[i], variables->item[i]);
	    exitElegant(1);
	  }
	}
	if (variables->initial_value[i]<variables->lower_limit[i]) {
	  if (force_inside) 
	    variables->varied_quan_value[i] = variables->initial_value[i] = variables->lower_limit[i];
	  else {
	    fprintf(stdout, "Initial value (%e) is smaller than lower limit for %s.%s\n", 
		    variables->initial_value[i],
		    variables->element[i], variables->item[i]);
	    exitElegant(1);
	  }
	}
	fprintf(stdout, "Initial value for %s.%s is %e\n", 
		variables->element[i], variables->item[i],
		variables->initial_value[i]);
	rpn_store(variables->initial_value[i], NULL, variables->memory_number[i]);
	variables->step[i] = variables->orig_step[i];
        if (variables->step[i]==0) {
	  if (variables->lower_limit[i]==variables->upper_limit[i])
	      fprintf(stdout, "Note: step size for %s set to %e.\n", 
                        variables->varied_quan_name[i], 
                        variables->orig_step[i] = variables->step[i] = 1);
            else
	      fprintf(stdout, "Note: step size for %s set to %e.\n", variables->varied_quan_name[i],
                        variables->orig_step[i] = variables->step[i] = 
                        (variables->upper_limit[i]-variables->lower_limit[i])/10);
            fflush(stdout);
	}
    }

#if USE_MPI
    if (scale_factor != 1.0)
      for (i=0; i<variables->n_variables; i++) {
	variables->step[i] *= scale_factor;
      }
#endif 

    for (i=0; i<covariables->n_covariables; i++) {
      if (!get_parameter_value(covariables->varied_quan_value+i, 
			       covariables->element[i], 
			       covariables->varied_param[i],
			       covariables->varied_type[i],
			       beamline))
        bombElegant("unable to get initial value for parameter", NULL);
      rpn_store(covariables->varied_quan_value[i] = rpn(covariables->equation[i]), 
		NULL, covariables->memory_number[i]);
    }

    final_property_values = count_final_properties();
    final_property_value = tmalloc(sizeof(*final_property_value)*final_property_values);

    if (!output->sums_vs_z) {
        output->sums_vs_z = tmalloc(sizeof(*output->sums_vs_z));
        output->n_z_points = 0;
        }

    i = variables->n_variables;
    variables->varied_quan_value[i] = 0;   /* end-of-optimization indicator */
    control->i_step = 0;
    optim_func_flags = FINAL_SUMS_ONLY + INHIBIT_FILE_OUTPUT + SILENT_RUNNING;

    if (!optimization_data->UDFcreated) {
        char UDFname[100];
        sprintf(UDFname, "optUDF%ld", optUDFcount++);
        create_udf(UDFname, optimization_data->equation);
        cp_str(&optimization_data->UDFname, UDFname);
        optimization_data->UDFcreated = 1;
        }

    startsLeft = optimization_data->n_restarts+1;
    result = DBL_MAX;
    balanceTerms = 1;
#if defined(UNIX)
    if (optimization_data->method!=OPTIM_METHOD_POWELL)
      signal(SIGINT, optimizationInterruptHandler);
#endif
#if USE_MPI
    SDDS_PopulationSetup (population_log, &(optimization_data->popLog), &(optimization_data->variables), &(optimization_data->covariables));
#endif
    while (startsLeft-- && !stopOptimization) {
      lastResult = result;
      switch (optimization_data->method) {
      case OPTIM_METHOD_SIMPLEX:	
#if USE_MPI
      case OPTIM_METHOD_HYBSIMPLEX:
	if (optimization_data->method==OPTIM_METHOD_HYBSIMPLEX) {
	  for (i=0; i<variables->n_variables; i++)
	    variables->step[i] += (random_2(0)-0.5)*variables->orig_step[i];
	  /* Disabling the report from simplexMin routine, as it will print result from the Master only.
	     We print the best result across all the processor in a higher level routine */
	  optimization_report_ptr = NULL;
	}
#endif
	fputs("Starting simplex optimization.\n", stdout);
        if (simplexMin(&result, variables->varied_quan_value, variables->step, 
                       variables->lower_limit, variables->upper_limit, NULL, 
                       variables->n_variables, optimization_data->target, 
                       optimization_data->tolerance, optimization_function, optimization_report_ptr,
                       optimization_data->n_evaluations, optimization_data->n_passes, 12, 
                       optimization_data->simplexDivisor, 
                       optimization_data->simplexPassRangeFactor, 
                       (optimization_data->includeSimplex1dScans?0:SIMPLEX_NO_1D_SCANS)+
		       (optimization_data->verbose>1?SIMPLEX_VERBOSE_LEVEL1:0)+
		       (optimization_data->verbose>2?SIMPLEX_VERBOSE_LEVEL2:0)+
		       (optimization_data->startFromSimplexVertex1?SIMPLEX_START_FROM_VERTEX1:0))<0) {
          if (result>optimization_data->tolerance) {
            if (!optimization_data->soft_failure)
              bombElegant("optimization unsuccessful--aborting", NULL);
            else
              fputs("warning: optimization unsuccessful--continuing\n", stdout);
          }
          else
            fputs("warning: maximum number of passes reached in simplex optimization", stdout);
        }
#if USE_MPI
      /* check if the result meets the requirement for each point across all the processors here */
	if (optimization_data->method==OPTIM_METHOD_HYBSIMPLEX) {
#if MPI_DEBUG
	  fprintf (stdout, "minimal value is %g on %d\n", result, myid);
#endif
   	  find_global_min_index (&result, &min_location, MPI_COMM_WORLD);
	  MPI_Bcast(variables->varied_quan_value, variables->n_variables, MPI_DOUBLE, min_location, MPI_COMM_WORLD);

	  if (optimization_data->fp_log && optimization_data->verbose >1)
	    optimization_report (result, variables->varied_quan_value, optimization_data->n_restarts+1-startsLeft, n_evaluations_made, variables->n_variables);
	} 
#endif
        if (simplexMinAbort(0))
          stopOptimization = 1;
        break;
      case OPTIM_METHOD_POWELL:
        fputs("Starting Powell optimization.\n", stdout);
        if (powellMin(&result, variables->varied_quan_value, variables->step,
                      variables->lower_limit, variables->upper_limit,
                      variables->n_variables, optimization_data->target, 
                      optimization_data->tolerance, optimization_function, 
                      optimization_report, 
                      optimization_data->n_evaluations/optimization_data->n_passes+1,
                      optimization_data->n_evaluations*
                      (optimization_data->n_evaluations/optimization_data->n_passes+1), 3)<0) {
          if (result>optimization_data->tolerance) {
            if (!optimization_data->soft_failure)
              bombElegant("optimization unsuccessful--aborting", NULL);
            else
              fputs("warning: optimization unsuccessful--continuing\n", stdout);
          }
          else
            fputs("warning: maximum number of passes reached in powell optimization", stdout);
        }
        break;
      case OPTIM_METHOD_GRID:
        fputs("Starting grid-search optimization.", stdout);
        if (!grid_search_min(&result, variables->varied_quan_value, variables->lower_limit, 
                             variables->upper_limit, variables->step,
                             variables->n_variables, optimization_data->target,
                             optimization_function)) {
          if (!optimization_data->soft_failure)
            bombElegant("optimization unsuccessful--aborting", NULL);
          else 
            fputs("warning: optimization unsuccessful--continuing", stdout);
        }
        if (optimAbort(0))
          stopOptimization = 1;
        break;
      case OPTIM_METHOD_SAMPLE:
        fputs("Starting grid-sample optimization.", stdout);
        if (!grid_sample_min(&result, variables->varied_quan_value, variables->lower_limit, 
                             variables->upper_limit, variables->step,
                             variables->n_variables, optimization_data->target, 
                             optimization_function, optimization_data->n_evaluations*1.0,
                             random_1_elegant)) {
          if (!optimization_data->soft_failure)
            bombElegant("optimization unsuccessful--aborting", NULL);
          else
            fputs("warning: optimization unsuccessful--continuing", stdout);
        }
        if (optimAbort(0))
          stopOptimization = 1;
        break;
      case OPTIM_METHOD_RANSAMPLE:
        fputs("Starting random-sampled optimization.", stdout);
        if (!randomSampleMin(&result, variables->varied_quan_value, variables->lower_limit, 
                             variables->upper_limit, variables->n_variables, 
                             optimization_data->target, optimization_function, 
                             optimization_data->n_evaluations, random_1_elegant)) {
          if (!optimization_data->soft_failure)
            bombElegant("optimization unsuccessful--aborting", NULL);
          else
            fputs("warning: optimization unsuccessful--continuing", stdout);
        }
        if (optimAbort(0))
          stopOptimization = 1;
        break;
      case OPTIM_METHOD_RANWALK:
        fputs("Starting random-sampled optimization.", stdout);
        if (!randomWalkMin(&result, variables->varied_quan_value, 
                           variables->lower_limit, variables->upper_limit,
                           variables->step,
                           variables->n_variables, optimization_data->target,
                           optimization_function, 
                           optimization_data->n_evaluations, random_1_elegant)) {
          if (!optimization_data->soft_failure)
            bombElegant("optimization unsuccessful--aborting", NULL);
          else
            fputs("warning: optimization unsuccessful--continuing", stdout);
        }
        if (optimAbort(0))
          stopOptimization = 1;
        break;
#if USE_MPI
      case OPTIM_METHOD_GENETIC:
	fputs("Starting genetic optimization.\n", stdout);
	n_total_evaluations_made = geneticMin(&result, variables->varied_quan_value, 
					variables->lower_limit, variables->upper_limit, variables->step,
					variables->n_variables, optimization_data->target, 
					optimization_function, optimization_data->n_iterations,
					optimization_data->max_no_change,
					optimization_data->population_size, output_sparsing_factor,
					optimization_data->print_all_individuals, population_log,
					      &(optimization_data->popLog), optimization_data->verbose, variables, covariables);
        if (optimAbort(0))
          stopOptimization = 1;
	break;
      case OPTIM_METHOD_SWARM:
	if (optimization_data->n_restarts+1-startsLeft == 1)
	    fputs("Starting particle swarm optimization.\n", stdout);
	swarmMin(&result, variables->varied_quan_value, 
					    variables->lower_limit, variables->upper_limit, variables->step,
					    variables->n_variables, optimization_data->target,
					    optimization_function, optimization_data->population_size,
					    optimization_data->n_restarts+1-startsLeft, optimization_data->n_restarts);

      /* check if the result meets the requirement for each point across all the processors here */
        if (USE_MPI) {
#if MPI_DEBUG
          fprintf (stdout, "minimal value is %g for iteration %ld on %d\n", result, optimization_data->n_restarts+1-startsLeft, myid);
#endif
          find_global_min_index (&result, &min_location, MPI_COMM_WORLD);
#if MPI_DEBUG
	  fprintf (stdout, "min_location=%d\n", min_location);
#endif
          MPI_Bcast(variables->varied_quan_value, variables->n_variables, MPI_DOUBLE, min_location, MPI_COMM_WORLD);

#if MPI_DEBUG 
	  if (isSlave) {
	      fprintf (stdout, "Minimal value is %.15g after %ld iterations.\n", result, optimization_data->n_restarts+1-startsLeft);
	      fprintf(stdout, "new variable values for iteration %ld\n", optimization_data->n_restarts+1-startsLeft);
	      fflush(stdout);
	      for (i=0; i<variables->n_variables; i++)
		  fprintf(stdout, "    %10s: %23.15e\n", variables->varied_quan_name[i], variables->varied_quan_value[i]);
	      fflush(stdout);
	     
	      if (covariables->n_covariables) {
		if (optimization_data->n_restarts==startsLeft)
		  covariables_global = trealloc(covariables_global, sizeof(*covariables_global)*(covariables->n_covariables+1));	
                if (result<lastResult || (optimization_data->n_restarts==startsLeft)) {  
	          for (i=0; i<covariables->n_covariables; i++)
		    covariables_global[i] = covariables->varied_quan_value[i];	
                }
		fprintf(stdout, "new covariable values:\n");
		  for (i=0; i<covariables->n_covariables; i++) 
		    fprintf(stdout, "    %10s: %23.15e\n", covariables->varied_quan_name[i], covariables_global[i]);
	      }
	      fflush(stdout);
          }
#endif
	  if ((optimization_data->verbose>1) && optimization_data->fp_log) {
	      fprintf (optimization_data->fp_log, "Minimal value is %.15g after %ld iterations.\n", result, optimization_data->n_restarts+1-startsLeft);
	      fprintf(optimization_data->fp_log, "new variable values for iteration %ld\n", optimization_data->n_restarts+1-startsLeft);
	      fflush(optimization_data->fp_log);
	      for (i=0; i<variables->n_variables; i++)
		  fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", variables->varied_quan_name[i], variables->varied_quan_value[i]);
	      fflush(optimization_data->fp_log);
	     
	      if (covariables->n_covariables) {
		if (optimization_data->n_restarts==startsLeft)
		  covariables_global = trealloc(covariables_global, sizeof(*covariables_global)*(covariables->n_covariables+1));	
                if (result<lastResult || (optimization_data->n_restarts==startsLeft)) {  
	          for (i=0; i<covariables->n_covariables; i++)
		    covariables_global[i] = covariables->varied_quan_value[i];	
                }
		fprintf(optimization_data->fp_log, "new covariable values:\n");
		  for (i=0; i<covariables->n_covariables; i++) 
		    fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", covariables->varied_quan_name[i], covariables_global[i]);
	      }
	      fflush(optimization_data->fp_log);
          }
        }
        if (optimAbort(0))
          stopOptimization = 1;
        break;
#else
	case OPTIM_METHOD_GENETIC:
	case OPTIM_METHOD_SWARM:
	  bombElegant("This optimization method can only be used in Pelegant.", NULL);
	  break;
#endif
      default:
        bombElegant("unknown optimization method code (do_optimize())", NULL);
        break;
      }
#if USE_MPI
      if ((optimization_data->method==OPTIM_METHOD_SWARM) || (optimization_data->method==OPTIM_METHOD_HYBSIMPLEX)) { 
	/* The covariables are updated locally after each optimization_function call no matter if a better result 
	   is achieved, which might not match the optimal variables for current iteration. So we keep the global best 
	   covariable values in a separate array. */
	if (covariables->n_covariables && (optimization_data->verbose>1) && (result<lastResult || (optimization_data->n_restarts==startsLeft))) 
	  MPI_Bcast(covariables->varied_quan_value, covariables->n_covariables, MPI_DOUBLE, min_location, MPI_COMM_WORLD);
	
	if (covariables->n_covariables) {
	  if (optimization_data->n_restarts==startsLeft)
	    covariables_global = trealloc(covariables_global, sizeof(*covariables_global)*(covariables->n_covariables+1));	
	  if (result<lastResult || (optimization_data->n_restarts==startsLeft)) {  
	    for (i=0; i<covariables->n_covariables; i++)
	      covariables_global[i] = covariables->varied_quan_value[i];	
	  }
	}
	if (population_log)
	  SDDS_PrintPopulations(&(optimization_data->popLog), optimization_data->n_restarts+1-startsLeft, result, variables->varied_quan_value, variables->n_variables, covariables_global, covariables->n_covariables);
      }
#else
      /* This part looks like redundant, as this is repeated after exiting the while loop. -- Y. Wang */
      /* evaluate once more at the optimimum point to get all parameters right and to get additional output */
      force_output = 1;
      ignoreOptimRecords = 1;  /* to force re-evaluation */ 
      result = optimization_function(variables->varied_quan_value, &i);
      ignoreOptimRecords = 0;
      force_output = 0;
#endif
      if (result<=optimization_data->target) {
	if (optimization_data->verbose>1) {
	  fprintf(stdout, "Target value reached, terminating optimization\n");
	}
        break;
      }
#if USE_MPI     
      if ((optimization_data->method!=OPTIM_METHOD_SWARM) && (optimization_data->method!=OPTIM_METHOD_HYBSIMPLEX))
#endif
	if (fabs(result-lastResult)<optimization_data->tolerance) {
	  if (optimization_data->verbose>1) {
	    fprintf(stdout, "New result (%21.15e) not sufficiently different from old result (%21.15e), terminating optimization\n", result, lastResult);
	  }
	  break;
	}    
      lastResult = result;
      if (startsLeft && !stopOptimization) {
        for (i=0; i<variables->n_variables; i++) {
          variables->step[i] = variables->orig_step[i];
        }
        if (optimization_data->restart_worst_term_factor!=1 && optimization_data->terms>1) {
          long imax, imin, iworst;
          double *savedTermValue, newResult;
          if (!(savedTermValue = malloc(sizeof(*savedTermValue)*optimization_data->terms)))
            bombElegant("memory allocation failure (saving term values)", NULL);
          memcpy(savedTermValue, optimization_data->termValue, optimization_data->terms*sizeof(*savedTermValue));
          newResult = lastResult;
          for (iworst=0; iworst<optimization_data->restart_worst_terms; iworst++) {
            if (index_min_max(&imin, &imax, optimization_data->termValue, optimization_data->terms)) {
              optimization_data->termWeight[imax] *= optimization_data->restart_worst_term_factor;
              fprintf(stdout, "Adjusted weight for term from %le to %le: %s\n",
                      optimization_data->termWeight[imax]/optimization_data->restart_worst_term_factor,
                      optimization_data->termWeight[imax], optimization_data->term[imax]);
              /* just to be sure it doesn't get picked as the max again */
              optimization_data->termValue[imax] = -DBL_MAX;
            }
            else 
              break;
          }
          memcpy(optimization_data->termValue, savedTermValue, optimization_data->terms*sizeof(*savedTermValue));
          index_min_max(&imin, &imax, optimization_data->termWeight, optimization_data->terms);
          if (optimization_data->termWeight[imin])
            for (i=0; i<optimization_data->terms; i++)
              optimization_data->termWeight[i] /= optimization_data->termWeight[imin];
        }
	if (optimization_data->verbose>1) { 
          fprintf(stdout, "Redoing optimization\n");
          fflush(stdout);
	}
      }
    }

#if defined(UNIX)
    if (optimization_data->method!=OPTIM_METHOD_POWELL)
      signal(SIGINT, traceback_handler);
#endif

    fprintf(stdout, "Exited optimization loop\n");
    fflush(stdout);

#if USE_MPI
    last_optimize_function_call = 1;
    min_value_location = min_location;
#endif

    /* evaluate once more at the optimimum point to get all parameters right and to get additional output */
    optim_func_flags = 0;
    force_output = 1;
    ignoreOptimRecords = 1; /* to force re-evaluation */
    variables->varied_quan_value[variables->n_variables] = 1;   /* indicates end-of-optimization */
    result = optimization_function(variables->varied_quan_value, &i);
    ignoreOptimRecords = 0; 
    force_output = 0;
          
    for (i=0; i<MAX_OPTIM_RECORDS; i++)
      free(optimRecord[i].variableValue);

    /* change values in element definitions so that new lattice can be saved */
    change_defined_parameter_values(variables->element, variables->varied_param, variables->varied_type,
            variables->varied_quan_value, variables->n_variables);
    if (control->n_elements_to_vary)
        change_defined_parameter_values(control->element, control->varied_param, control->varied_type, 
            control->varied_quan_value, control->n_elements_to_vary);
    if (covariables->n_covariables)
        change_defined_parameter_values(covariables->element, covariables->varied_param, covariables->varied_type,
            covariables->varied_quan_value, covariables->n_covariables);

#if USE_MPI
    MPI_Reduce(&n_evaluations_made, &n_total_evaluations_made, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#if MPI_DEBUG
          if (isSlave) {
              fprintf (stdout, "Minimal value is %.15g after %ld iterations.\n", result, optimization_data->n_restarts+1-startsLeft);
              fprintf(stdout, "new variable values for iteration %ld\n", optimization_data->n_restarts+1-startsLeft);
              fflush(stdout);
              for (i=0; i<variables->n_variables; i++)
                  fprintf(stdout, "    %10s: %23.15e\n", variables->varied_quan_name[i], variables->varied_quan_value[i]);
              fflush(stdout);
          }
#endif
#endif

    /* Master only starting here ... */
    if(isMaster) {
      if (optimization_data->fp_log) {
        fprintf(optimization_data->fp_log, "Optimization results:\n  optimization function has value %.15g\n", 
                optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
        if (optimization_data->terms) {
          double sum;
          fprintf(optimization_data->fp_log, "Terms of equation: \n");
          for (i=sum=0; i<optimization_data->terms; i++) {
            rpn_clear();
            fprintf(optimization_data->fp_log, "%20s: %23.15e\n",
                    optimization_data->term[i],
                    rpn(optimization_data->term[i])*optimization_data->termWeight[i]);
            sum += rpn(optimization_data->term[i])*optimization_data->termWeight[i];
          }
        }
#if !USE_MPI
        fprintf(optimization_data->fp_log, "    A total of %ld function evaluations were made.\n", n_evaluations_made);
#else
	fprintf(optimization_data->fp_log, "    A total of %ld function evaluations were made with an average of %ld function evaluations per processor\n", n_total_evaluations_made, n_total_evaluations_made/n_processors); 
#endif
        if (constraints->n_constraints) {
	  fprintf(optimization_data->fp_log, "Constraints:\n");
	  for (i=0; i<constraints->n_constraints; i++)
	    fprintf(optimization_data->fp_log, "%10s: %23.15e\n", constraints->quantity[i], 
		    final_property_value[constraints->index[i]]);
	}
        fprintf(optimization_data->fp_log, "Optimum values of variables and changes from initial values:\n");
        for (i=0; i<variables->n_variables; i++)
	  fprintf(optimization_data->fp_log, "%10s: %23.15e  %23.15e (was %23.15e) %s\n", variables->varied_quan_name[i], 
		  variables->varied_quan_value[i], variables->varied_quan_value[i]-variables->initial_value[i],
		  variables->initial_value[i],
		  variableAtLimit(variables->varied_quan_value[i], variables->lower_limit[i],
				  variables->upper_limit[i]) ? "(near limit)" : "");
        for (i=0; i<covariables->n_covariables; i++)
	  fprintf(optimization_data->fp_log, "%10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
        fflush(optimization_data->fp_log);
      }
      if (!log_file || optimization_data->fp_log!=stdout) {
        fprintf(stdout, "Optimization results:\n    optimization function has value %.15g\n", 
                optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
        if (optimization_data->terms) {
          double sum;
          fprintf(stdout, "Terms of equation: \n");
          for (i=sum=0; i<optimization_data->terms; i++) {
            rpn_clear();
            fprintf(stdout, "%20s: %23.15e\n",
                    optimization_data->term[i],
                    optimization_data->termWeight[i]*rpn(optimization_data->term[i]));
            sum += optimization_data->termWeight[i]*rpn(optimization_data->term[i]);
          }
        }
        fflush(stdout);
#if !USE_MPI
        fprintf(stdout, "    A total of %ld function evaluations were made.\n", n_evaluations_made);
#else
	fprintf(optimization_data->fp_log, "    A total of %ld function evaluations were made with an average of %ld function evaluations per processor\n", n_total_evaluations_made, n_total_evaluations_made/n_processors);
#endif
        fflush(stdout);
        if (constraints->n_constraints) {
	  fprintf(stdout, "Constraints:\n");
	  fflush(stdout);
	  for (i=0; i<constraints->n_constraints; i++)
	    fprintf(stdout, "%10s: %23.15e\n", constraints->quantity[i], final_property_value[constraints->index[i]]);
	  fflush(stdout);
	}
        fprintf(stdout, "Optimum values of variables and changes from initial values:\n");
        fflush(stdout);
        for (i=0; i<variables->n_variables; i++)
	  fprintf(stdout, "%10s: %23.15e  %23.15e\n", variables->varied_quan_name[i], 
		  variables->varied_quan_value[i], variables->varied_quan_value[i]-variables->initial_value[i]);
        fflush(stdout);
        for (i=0; i<covariables->n_covariables; i++)
	  fprintf(stdout, "%10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
        fflush(stdout);
      }
      if (term_log_file && strlen(term_log_file) && optimization_data->terms) {
        SDDS_DATASET termLog;
        long i, iterm, ivalue;
        if (!SDDS_InitializeOutput(&termLog, SDDS_BINARY, 0, NULL, NULL, term_log_file) ||
            (iterm=SDDS_DefineColumn(&termLog, "Term", NULL, NULL, NULL, NULL, SDDS_STRING, 0))<0 ||
            (ivalue=SDDS_DefineColumn(&termLog, "Contribution", NULL, NULL, NULL, NULL, SDDS_DOUBLE, 0))<0 ||
            !SDDS_WriteLayout(&termLog) ||
            !SDDS_StartPage(&termLog, optimization_data->terms)) {
          fprintf(stdout, "Problem writing optimization term log\n");
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);
        }
        for (i=0; i<optimization_data->terms; i++) {
          if (!SDDS_SetRowValues(&termLog, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, i,
                                 iterm, optimization_data->term[i],
                                 ivalue,optimization_data->termWeight[i]*rpn(optimization_data->term[i]), 
                                 -1)) {
            fprintf(stdout, "Problem writing optimization term log\n");
            SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
            exitElegant(1);          
          }
        }
        if (!SDDS_WritePage(&termLog) || !SDDS_Terminate(&termLog)) {
          fprintf(stdout, "Problem writing optimization term log\n");
          SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors);
          exitElegant(1);          
        }
      }
    }
    /* ... Master only ending here */

    for (i=0; i<variables->n_variables; i++)
        variables->initial_value[i] = variables->varied_quan_value[i];
    control->i_step = i_step_saved;
    log_exit("do_optimize");
    }

/* Next three lines from elegant.c: */
#define SET_AWE_BEAM     5
#define SET_BUNCHED_BEAM 6
#define SET_SDDS_BEAM   33

#define N_TWISS_QUANS 88
static char *twiss_name[N_TWISS_QUANS] = {
    "betax", "alphax", "nux", "etax", "etapx", 
    "betay", "alphay", "nuy", "etay", "etapy", 
    "max.betax", "max.etax", "max.etapx", 
    "max.betay", "max.etay", "max.etapy",
    "min.betax", "min.etax", "min.etapx", 
    "min.betay", "min.etay", "min.etapy",
    "dnux/dp", "dnuy/dp", "alphac", "alphac2",
    "ave.betax", "ave.betay",
    "etaxp", "etayp",
    "waistsx", "waistsy",
    "dnux/dAx", "dnux/dAy", "dnuy/dAx", "dnuy/dAy",
    "dnux/dp2", "dnux/dp3",
    "dnuy/dp2", "dnuy/dp3",
    "etax2" , "etax3", 
    "etay2" , "etay3", 
    "nuxChromLower", "nuxChromUpper", 
    "nuyChromLower", "nuyChromUpper", 
    "dbetax/dp", "dbetay/dp", "dalphax/dp", "dalphay/dp",
    "dnux/dAx2", "dnux/dAy2", "dnuy/dAx2", "dnuy/dAy2",
    "dnux/dAxAy", "dnuy/dAxAy",
    "nuxTswaLower", "nuxTswaUpper", 
    "nuyTswaLower", "nuyTswaUpper", 
    "couplingIntegral", "emittanceRatio",
    "h21000", "h30000", "h10110", "h10020", "h10200",
    "dnux/dJx", "dnux/dJy", "dnuy/dJy",
    "h11001", "h00111", "h20001", "h00201", "h10002",
    "h22000", "h11110", "h00220", "h31000", "h40000",
    "h20110", "h11200", "h20020", "h20200", "h00310", "h00400"
    };
static long twiss_mem[N_TWISS_QUANS] = {
  -1, -1, -1, -1, -1,  
  -1, -1, -1, -1, -1,  
  -1, -1, -1,  
  -1, -1, -1, 
  -1, -1, -1,  
  -1, -1, -1, 
  -1, -1, -1, -1, 
  -1, -1, 
  -1, -1, 
  -1, -1, 
  -1, -1, -1, -1, 
  -1, -1, 
  -1, -1, 
  -1, -1,  
  -1, -1,  
  -1, -1,  
  -1, -1,  
  -1, -1, -1, -1,
  -1, -1, -1, -1,
  -1, -1,
  -1, -1, -1, -1,
  -1, -1,
  -1, -1, -1, -1, -1,
  -1, -1, -1,
  -1, -1, -1, -1, -1,
    };

static char *radint_name[13] = {
    "ex0", "Sdelta0",
    "Jx", "Jy", "Jdelta",
    "taux", "tauy", "taudelta",
    "I1", "I2", "I3", "I4", "I5",
  } ;
static long radint_mem[13] = {
  -1, -1, 
  -1, -1, -1,
  -1, -1, -1,
  -1, -1, -1, -1, -1,
} ;
static char *floorCoord_name[7] = {
  "X", "Y", "Z", "theta", "phi", "psi", "sTotal", 
};
static long floorCoord_mem[7] = {
  -1, -1, -1, -1, -1, -1, -1
};
static char *floorStat_name[6] = {
  "min.X", "min.Y", "min.Z", "max.X", "max.Y", "max.Z",
};
static long floorStat_mem[6] = {
  -1, -1, -1, -1, -1, -1,
};

int showTwissMemories(FILE *fp)
{
  long i;
  
  for (i=0; i<32; i++) {
    if (twiss_mem[i]!=-1)
      fprintf(fp, "%s = %21.15e\n", 
              twiss_name[i], rpn_recall(twiss_mem[i]));
  }
  fflush(fp);
  return 0;
}

double optimization_function(double *value, long *invalid)
{
  double rpn(char *expression);
  OPTIM_VARIABLES *variables;
  OPTIM_CONSTRAINTS *constraints;
  OPTIM_COVARIABLES *covariables;
  double conval, result=0;
  long i, iRec, recordUsedAgain;
  unsigned long unstable;
  VMATRIX *M;
  TWISS twiss_ave, twiss_min, twiss_max;
  double XYZ[3], Angle[3], XYZMin[3], XYZMax[3];
  double startingOrbitCoord[6] = {0,0,0,0,0,0};
  long rpnError = 0;
  
  log_entry("optimization_function");
  
#if DEBUG
  fprintf(stdout, "optimization_function: In optimization function\n");
  fprintf(stdout, "Beamline flags: %lx\n", beamline->flags);
  fflush(stdout);
#endif

  if (restart_random_numbers)
    seedElegantRandomNumbers(0, 1);
  
  *invalid = 0;
  unstable = 0;
  n_evaluations_made++;
  
  variables = &(optimization_data->variables);
  constraints = &(optimization_data->constraints);
  covariables = &(optimization_data->covariables);

  /* assert variable values and store in rpn memories */
  delete_phase_references();
  reset_special_elements(beamline, 1);

  if (beamline->links && beamline->links->n_links)
    reset_element_links(beamline->links, run, beamline);
  
  assert_parameter_values(variables->element, variables->varied_param, variables->varied_type,
                          value, variables->n_variables, beamline);
  for (i=0; i<variables->n_variables; i++)
    rpn_store(value[i], NULL, variables->memory_number[i]);

  /* set element flags to indicate variation of parameters */
  set_element_flags(beamline, variables->element, NULL, variables->varied_type, variables->varied_param, 
                    variables->n_variables, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);

  if (covariables->n_covariables) {
#if DEBUG
    fprintf(stdout, "optimization_function: Computing covariables\n");
    fflush(stdout);
#endif

    /* calculate values of covariables and assert these as well */
    for (i=0; i<covariables->n_covariables; i++) {
      rpn_clear();
      rpn_store(covariables->varied_quan_value[i]=rpn(covariables->pcode[i]), NULL, covariables->memory_number[i]);
      if (rpn_check_error()) exitElegant(1);
    }
    assert_parameter_values(covariables->element, covariables->varied_param, covariables->varied_type,
                            covariables->varied_quan_value, covariables->n_covariables, beamline);
    /* set element flags to indicate variation of parameters */
    set_element_flags(beamline, covariables->element, NULL, covariables->varied_type, covariables->varied_param, 
                      covariables->n_covariables, PARAMETERS_ARE_VARIED, VMATRIX_IS_VARIED, 0, 0);
  }
#if !USE_MPI || MPI_DEBUG /* The information here is from a processor locally, we print the result across all the processors after an iteration */
  if (optimization_data->verbose && optimization_data->fp_log) {
    fprintf(optimization_data->fp_log, "new variable values for evaluation %ld of pass %ld:\n", n_evaluations_made, n_passes_made+1);
    fflush(optimization_data->fp_log);
    for (i=0; i<variables->n_variables; i++)
      fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", variables->varied_quan_name[i], value[i]);
    fflush(optimization_data->fp_log);
    if (covariables->n_covariables) {
      fprintf(optimization_data->fp_log, "new covariable values:\n");
      for (i=0; i<covariables->n_covariables; i++)
        fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", covariables->varied_quan_name[i], covariables->varied_quan_value[i]);
    }
    fflush(optimization_data->fp_log);
  }
#endif
  if ((iRec=checkForOptimRecord(value, variables->n_variables, &recordUsedAgain))>=0) {
#if USE_MPI
  if (!runInSinglePartMode) /* For parallel genetic optimization, all the individuals for the first iteration are same */
#endif
    if (recordUsedAgain>20) {
      fprintf(stdout, "record used too many times---stopping optimization\n");
      stopOptimization = 1;
      simplexMinAbort(1);
    }
    if (optimization_data->verbose && optimization_data->fp_log)
      fprintf(optimization_data->fp_log, "Using previously computed value %23.15e\n\n", 
              optimRecord[iRec].result);
    *invalid = optimRecord[iRec].invalid;
    return optimRecord[iRec].result;
  }
  
  /* compute matrices for perturbed elements */
#if DEBUG
  fprintf(stdout, "optimization_function: Computing matrices\n");
  fflush(stdout);
#endif
  if (beamline->links && beamline->links->n_links)
    rebaseline_element_links(beamline->links, run, beamline);
  i = compute_changed_matrices(beamline, run) +
      assert_element_links(beamline->links, run, beamline, STATIC_LINK+DYNAMIC_LINK+LINK_ELEMENT_DEFINITION);
  if (beamline->flags&BEAMLINE_CONCAT_DONE)
    free_elements1(&(beamline->ecat));
  beamline->flags &= ~(BEAMLINE_CONCAT_CURRENT+BEAMLINE_CONCAT_DONE+
		       BEAMLINE_TWISS_CURRENT+BEAMLINE_TWISS_DONE+
		       BEAMLINE_RADINT_CURRENT+BEAMLINE_RADINT_DONE);

  if (i && beamline->matrix) {
    free_matrices(beamline->matrix);
    free(beamline->matrix);
    beamline->matrix = NULL;
  }
  
  if (optimization_data->verbose && optimization_data->fp_log) {
    fprintf(optimization_data->fp_log, "%ld matrices (re)computed\n", i);
    fflush(optimization_data->fp_log);
  }

#if DEBUG
  fprintf(stdout, "optimization_function: Generating beam\n");
  fflush(stdout);
#endif
  /* generate initial beam distribution and track it */
  switch (beam_type_code) {
  case SET_AWE_BEAM:
    bombElegant("beam type code of SET_AWE_BEAM in optimization_function--this shouldn't happen", NULL);
    break;
  case SET_BUNCHED_BEAM:
    new_bunched_beam(beam, run, control, output, 0);
    break;
  case SET_SDDS_BEAM:
    if (new_sdds_beam(beam, run, control, output, 0)<0)
      bombElegant("unable to get beam for tracking (is file empty or out of pages?)", NULL);
    break;
  default:
    bombElegant("unknown beam type code in optimization", NULL);
    break;
  }
  control->i_step++;       /* to prevent automatic regeneration of beam */
  zero_beam_sums(output->sums_vs_z, output->n_z_points+1);

  if (doClosedOrbit &&
      !run_closed_orbit(run, beamline, startingOrbitCoord, beam, 0)) {
    *invalid = 1;
    fprintf(stdout, "warning: unable to find closed orbit\n");
    fflush(stdout);
  }
  if (!*invalid && orbitCorrMode!=-1 && 
      !do_correction(orbitCorrData, run, beamline, startingOrbitCoord, beam, control->i_step, 
                     INITIAL_CORRECTION+NO_OUTPUT_CORRECTION)) {
    *invalid = 1;
    fprintf(stdout, "warning: unable to perform orbit correction\n");
    fflush(stdout);
  }
  if (!*invalid && doTuneCorr &&
      !do_tune_correction(tuneCorrData, run, beamline, startingOrbitCoord, 
                                  0, 0)) {
    *invalid = 1;
    fprintf(stdout, "warning: unable to do tune correction\n");
    fflush(stdout);
  }
  if (!*invalid && doChromCorr &&
      !do_chromaticity_correction(chromCorrData, run, beamline, startingOrbitCoord, 
                                  0, 0)) {
    *invalid = 1;
    fprintf(stdout, "warning: unable to do chromaticity correction\n");
    fflush(stdout);
  }
  if (!*invalid && doClosedOrbit &&
      !run_closed_orbit(run, beamline, startingOrbitCoord, beam, 0)) {
    *invalid = 1;
    fprintf(stdout, "warning: unable to find closed orbit\n");
    fflush(stdout);
  }

  if (!*invalid && doResponse) {
    update_response(run, beamline, orbitCorrData);
  }
  
#if DEBUG
  if (doClosedOrbit) {
    fprintf(stdout, "Closed orbit: %g, %g, %g, %g, %g, %g\n",
	    beamline->closed_orbit->centroid[0],
	    beamline->closed_orbit->centroid[1],
	    beamline->closed_orbit->centroid[2],
	    beamline->closed_orbit->centroid[3],
	    beamline->closed_orbit->centroid[4],
	    beamline->closed_orbit->centroid[5]);
  }
  fprintf(stdout, "Beamline flags: %lx\n", beamline->flags);
#endif

  if (!*invalid && beamline->flags&BEAMLINE_TWISS_WANTED) {
    if (twiss_mem[0]==-1) {
      for (i=0; i<N_TWISS_QUANS; i++)
        twiss_mem[i] = rpn_create_mem(twiss_name[i], 0);
    }
    /* get twiss mode and (beta, alpha, eta, etap) for both planes */
    if (optimization_data->verbose && optimization_data->fp_log) {
      fprintf(optimization_data->fp_log, "Computing twiss parameters for optimization\n");
      fflush(optimization_data->fp_log);
    }
    update_twiss_parameters(run, beamline, &unstable);
    run_coupled_twiss_output(run, beamline, startingOrbitCoord);
    if (unstable)
      *invalid = 1;
#if DEBUG
    fprintf(stdout, "Twiss parameters done.\n");
    fprintf(stdout, "betax=%g, alphax=%g, etax=%g\n", 
            beamline->elast->twiss->betax,
            beamline->elast->twiss->alphax,
            beamline->elast->twiss->etax);
    fprintf(stdout, "betay=%g, alphay=%g, etay=%g\n", 
            beamline->elast->twiss->betay,
            beamline->elast->twiss->alphay,
            beamline->elast->twiss->etay);
    fflush(stdout);
#endif
    /* store twiss parameters for last element */
    for (i=0; i<5; i++) {
      rpn_store(*((&beamline->elast->twiss->betax)+i)/(i==2?PIx2:1), NULL, twiss_mem[i]);
      rpn_store(*((&beamline->elast->twiss->betay)+i)/(i==2?PIx2:1), NULL, twiss_mem[i+5]);
    }
    /* store statistics */
    compute_twiss_statistics(beamline, &twiss_ave, &twiss_min, &twiss_max);
    rpn_store(twiss_max.betax, NULL, twiss_mem[10]);
    rpn_store(twiss_max.etax, NULL,  twiss_mem[11]);
    rpn_store(twiss_max.etapx, NULL, twiss_mem[12]);
    rpn_store(twiss_max.betay, NULL, twiss_mem[13]);
    rpn_store(twiss_max.etay, NULL,  twiss_mem[14]);
    rpn_store(twiss_max.etapy, NULL, twiss_mem[15]);
    rpn_store(twiss_min.betax, NULL, twiss_mem[16]);
    rpn_store(twiss_min.etax, NULL,  twiss_mem[17]);
    rpn_store(twiss_min.etapx, NULL, twiss_mem[18]);
    rpn_store(twiss_min.betay, NULL, twiss_mem[19]);
    rpn_store(twiss_min.etay, NULL,  twiss_mem[20]);
    rpn_store(twiss_min.etapy, NULL, twiss_mem[21]);
    /* chromaticity */
    rpn_store(beamline->chromaticity[0], NULL, twiss_mem[22]);
    rpn_store(beamline->chromaticity[1], NULL, twiss_mem[23]);
    /* first and second-order momentum compaction */
    rpn_store(beamline->alpha[0], NULL, twiss_mem[24]);
    rpn_store(beamline->alpha[1], NULL, twiss_mem[25]);
    /* more statistics */
    rpn_store(twiss_ave.betax, NULL, twiss_mem[26]);
    rpn_store(twiss_ave.betay, NULL, twiss_mem[27]);
    /* alternate names for etapx and etapy */
    rpn_store(beamline->elast->twiss->etapx, NULL, twiss_mem[28]);
    rpn_store(beamline->elast->twiss->etapy, NULL, twiss_mem[29]);
    /* number of waists per plane */
    rpn_store((double)beamline->waists[0], NULL, twiss_mem[30]);
    rpn_store((double)beamline->waists[1], NULL, twiss_mem[31]);
    /* amplitude-dependent tune shifts */
    rpn_store(beamline->dnux_dA[1][0], NULL, twiss_mem[32]);
    rpn_store(beamline->dnux_dA[0][1], NULL, twiss_mem[33]);
    rpn_store(beamline->dnuy_dA[1][0], NULL, twiss_mem[34]);
    rpn_store(beamline->dnuy_dA[0][1], NULL, twiss_mem[35]);
    /* higher-order chromaticities */
    rpn_store(beamline->chrom2[0], NULL, twiss_mem[36]);
    rpn_store(beamline->chrom3[0], NULL, twiss_mem[37]);
    rpn_store(beamline->chrom2[1], NULL, twiss_mem[38]);
    rpn_store(beamline->chrom3[1], NULL, twiss_mem[39]);
    /* higher-order dispersion */
    rpn_store(beamline->eta2[0], NULL, twiss_mem[40]);
    rpn_store(beamline->eta3[0], NULL, twiss_mem[41]);
    rpn_store(beamline->eta2[2], NULL, twiss_mem[42]);
    rpn_store(beamline->eta3[2], NULL, twiss_mem[43]);
    /* limits of tunes due to chromatic effects */
    rpn_store(beamline->tuneChromLower[0], NULL, twiss_mem[44]);
    rpn_store(beamline->tuneChromUpper[0], NULL, twiss_mem[45]);
    rpn_store(beamline->tuneChromLower[1], NULL, twiss_mem[46]);
    rpn_store(beamline->tuneChromUpper[1], NULL, twiss_mem[47]);
    /* derivative of beta functions with momentum offset */
    rpn_store(beamline->dbeta_dPoP[0], NULL, twiss_mem[48]);
    rpn_store(beamline->dbeta_dPoP[1], NULL, twiss_mem[49]);
    rpn_store(beamline->dalpha_dPoP[0], NULL, twiss_mem[50]);
    rpn_store(beamline->dalpha_dPoP[1], NULL, twiss_mem[51]);
    /* higher-order tune shifts with amplitude */
    rpn_store(beamline->dnux_dA[2][0], NULL, twiss_mem[52]);
    rpn_store(beamline->dnux_dA[0][2], NULL, twiss_mem[53]);
    rpn_store(beamline->dnuy_dA[2][0], NULL, twiss_mem[54]);
    rpn_store(beamline->dnuy_dA[0][2], NULL, twiss_mem[55]);
    rpn_store(beamline->dnux_dA[1][1], NULL, twiss_mem[56]);
    rpn_store(beamline->dnuy_dA[1][1], NULL, twiss_mem[57]);
    /* tune extrema due to TSWA */
    rpn_store(beamline->nuxTswaExtrema[0], NULL, twiss_mem[58]);
    rpn_store(beamline->nuxTswaExtrema[1], NULL, twiss_mem[59]);
    rpn_store(beamline->nuyTswaExtrema[0], NULL, twiss_mem[60]);
    rpn_store(beamline->nuyTswaExtrema[1], NULL, twiss_mem[61]);
    /* coupling parameters */
    rpn_store(beamline->couplingFactor[0], NULL, twiss_mem[62]);
    rpn_store(beamline->couplingFactor[2], NULL, twiss_mem[63]);
    /* geometric driving terms */
    rpn_store(beamline->drivingTerms.h21000, NULL, twiss_mem[64]);
    rpn_store(beamline->drivingTerms.h30000, NULL, twiss_mem[65]);
    rpn_store(beamline->drivingTerms.h10110, NULL, twiss_mem[66]);
    rpn_store(beamline->drivingTerms.h10020, NULL, twiss_mem[67]);
    rpn_store(beamline->drivingTerms.h10200, NULL, twiss_mem[68]);
    rpn_store(beamline->drivingTerms.dnux_dJx, NULL, twiss_mem[69]);
    rpn_store(beamline->drivingTerms.dnux_dJy, NULL, twiss_mem[70]);
    rpn_store(beamline->drivingTerms.dnuy_dJy, NULL, twiss_mem[71]);
    rpn_store(beamline->drivingTerms.h11001, NULL, twiss_mem[72]);
    rpn_store(beamline->drivingTerms.h00111, NULL, twiss_mem[73]);
    rpn_store(beamline->drivingTerms.h20001, NULL, twiss_mem[74]);
    rpn_store(beamline->drivingTerms.h00201, NULL, twiss_mem[75]);
    rpn_store(beamline->drivingTerms.h10002, NULL, twiss_mem[76]);
    rpn_store(beamline->drivingTerms.h22000, NULL, twiss_mem[77]); 
    rpn_store(beamline->drivingTerms.h11110, NULL, twiss_mem[78]);
    rpn_store(beamline->drivingTerms.h00220, NULL, twiss_mem[79]);
    rpn_store(beamline->drivingTerms.h31000, NULL, twiss_mem[80]);
    rpn_store(beamline->drivingTerms.h40000, NULL, twiss_mem[81]);
    rpn_store(beamline->drivingTerms.h20110, NULL, twiss_mem[82]);
    rpn_store(beamline->drivingTerms.h11200, NULL, twiss_mem[83]);
    rpn_store(beamline->drivingTerms.h20020, NULL, twiss_mem[84]); 
    rpn_store(beamline->drivingTerms.h20200, NULL, twiss_mem[85]);
    rpn_store(beamline->drivingTerms.h00310, NULL, twiss_mem[86]); 
    rpn_store(beamline->drivingTerms.h00400, NULL, twiss_mem[87]);

#if DEBUG
    fprintf(stdout, "Twiss parameters done.\n");
    fflush(stdout);
#endif
  }

  if (!*invalid && beamline->flags&BEAMLINE_RADINT_WANTED) {
    if (optimization_data->verbose && optimization_data->fp_log) {
      fprintf(optimization_data->fp_log, "Updating radiation integral values for optimization\n");
      fflush(optimization_data->fp_log);
    }
    if (radint_mem[0]==-1) {
      for (i=0; i<13; i++)
        radint_mem[i] = rpn_create_mem(radint_name[i], 0);
    }
    /* radiation integrals already updated by update_twiss_parameters above
       which is guaranteed to be called
       */
    rpn_store(beamline->radIntegrals.ex0>0?beamline->radIntegrals.ex0:sqrt(DBL_MAX), NULL, 
              radint_mem[0]);
    rpn_store(beamline->radIntegrals.sigmadelta, NULL, radint_mem[1]);
    rpn_store(beamline->radIntegrals.Jx, NULL, radint_mem[2]);
    rpn_store(beamline->radIntegrals.Jy, NULL, radint_mem[3]);
    rpn_store(beamline->radIntegrals.Jdelta, NULL, radint_mem[4]);
    rpn_store(beamline->radIntegrals.taux, NULL, radint_mem[5]);
    rpn_store(beamline->radIntegrals.tauy, NULL, radint_mem[6]);
    rpn_store(beamline->radIntegrals.taudelta, NULL, radint_mem[7]);
    for (i=0; i<5; i++)
      rpn_store(beamline->radIntegrals.RI[i], NULL, radint_mem[i+8]);
  }

  runMomentsOutput(run, beamline, startingOrbitCoord, 1, 0);

  if (floorCoord_mem[0]==-1) 
    for (i=0; i<7; i++)
      floorCoord_mem[i] = rpn_create_mem(floorCoord_name[i], 0);
  if (floorStat_mem[0]==-1)
    for (i=0; i<6; i++)
      floorStat_mem[i] = rpn_create_mem(floorStat_name[i], 0);
  final_floor_coordinates(beamline, XYZ, Angle, XYZMin, XYZMax);
  for (i=0; i<3; i++) {
    rpn_store(XYZ[i], NULL, floorCoord_mem[i]);
    rpn_store(Angle[i], NULL, floorCoord_mem[i+3]);
    rpn_store(XYZMin[i], NULL, floorStat_mem[i]);
    rpn_store(XYZMax[i], NULL, floorStat_mem[i+3]);
  }
  compute_end_positions(beamline);
  rpn_store(beamline->revolution_length, NULL, floorCoord_mem[6]);

  for (i=0; i<variables->n_variables; i++)
    variables->varied_quan_value[i] = value[i];

  if (!*invalid) {
    output->n_z_points = 0;
    if (output->sums_vs_z) {
      free(output->sums_vs_z);
      output->sums_vs_z = NULL;
    }
    M = accumulate_matrices(&(beamline->elem), run, NULL,
                            optimization_data->matrix_order<1?1:optimization_data->matrix_order, 0);

    if (optimization_data->verbose && optimization_data->fp_log) {
      fprintf(optimization_data->fp_log, "Tracking for optimization\n");
      fflush(optimization_data->fp_log);
    }
#if SDDS_MPI_IO
    if (isMaster && notSinglePart)
      if (beam->n_to_track_total<(n_processors-1)) {
	printf("*************************************************************************************************\n");
	printf("* Warning! The number of particles (%ld) shouldn't be less than the number of processors (%d)! *\n", beam->n_to_track, n_processors-1);
	printf("* Less number of processors are recommended!                                                    *\n");
	printf("*************************************************************************************************\n");
	MPI_Abort(MPI_COMM_WORLD, 2);    
      }
#endif
    track_beam(run, control, error, variables, beamline, beam, output, optim_func_flags, 1,
               &charge);

    if (doFindAperture) {
      double area;
      do_aperture_search(run, control, startingOrbitCoord, error, beamline, &area);
      rpn_store(area, NULL, rpn_create_mem("DaArea", 0));
      if (optimization_data->verbose)
        printf("DA area is %e\n", area);
    } 

  }
  
    if (!*invalid) {
      if (output->sasefel.active)
        storeSASEFELAtEndInRPN(&(output->sasefel));

      /* compute final parameters and store in rpn memories */
#if DEBUG
      fprintf(stdout, "optimization_function: Computing final parameters\n");
      fflush(stdout);
#endif
      if (!output->sums_vs_z)
        bombElegant("sums_vs_z element of output structure is NULL--programming error (optimization_function)", NULL);
      if ((i=compute_final_properties(final_property_value, output->sums_vs_z+output->n_z_points, 
                                      beam->n_to_track, beam->p0, M, beam->particle, 
                                      control->i_step, control->indexLimitProduct*control->n_steps,
                                      charge))
          != final_property_values) {
        fprintf(stdout, "error: compute_final_properties computed %ld quantities when %ld were expected (optimization_function)\n",
                i, final_property_values);
        fflush(stdout);
        abort();
      }
      if (isMaster || !notSinglePart) { /* Only the master will execute the block */
	rpn_store_final_properties(final_property_value, final_property_values);
	if (optimization_data->matrix_order>1)
	  rpnStoreHigherMatrixElements(M, &optimization_data->TijkMem,
				       &optimization_data->UijklMem,
				       optimization_data->matrix_order);
	free_matrices(M); free(M); M = NULL;

#if DEBUG
	fprintf(stdout, "optimization_function: Checking constraints\n");
	fflush(stdout);
#endif
	/* check constraints */
	if (optimization_data->verbose && optimization_data->fp_log && constraints->n_constraints) {
	  fprintf(optimization_data->fp_log, "    Constraints:\n");
	  fflush(optimization_data->fp_log);
	}
	for (i=0; i<constraints->n_constraints; i++) {
	  if (optimization_data->verbose && optimization_data->fp_log)
	    fprintf(optimization_data->fp_log, "    %10s: %23.15e", constraints->quantity[i], 
		    final_property_value[constraints->index[i]]);
	  if ((conval=final_property_value[constraints->index[i]])<constraints->lower[i] ||
	      conval>constraints->upper[i]) {
	    *invalid = 1;
	    if (optimization_data->verbose && optimization_data->fp_log) {
	      fprintf(optimization_data->fp_log, " ---- invalid\n\n");
	      fflush(optimization_data->fp_log);
	    }
	    log_exit("optimization_function");
	    break;
	  }
	  if (optimization_data->verbose && optimization_data->fp_log) {
	    fputc('\n', optimization_data->fp_log);
	    fflush(optimization_data->fp_log);
	  }
	}
      
#if DEBUG
	fprintf(stdout, "optimization_function: Computing rpn function\n");
	fflush(stdout);
#endif
	result = 0;
	rpn_clear();
	if (!*invalid) {
	  long i=0, terms=0;
	  double value, sum;
	  if (balanceTerms && optimization_data->balance_terms && optimization_data->terms) {
	    for (i=sum=0; i<optimization_data->terms; i++) {
	      rpn_clear();
	      if ((value=rpn(optimization_data->term[i]))!=0) {
		optimization_data->termWeight[i] = 1/fabs(value);
		terms ++;
	      }
	      else
		optimization_data->termWeight[i] = 0;
	      sum += fabs(value);
	      if (rpn_check_error()) {
		printf("Problem evaluating expression: %s\n", optimization_data->term[i]);
		rpn_clear_error();
		rpnError++;
	      }
	    }
	    if (rpnError) {
	      printf("RPN expression errors prevent balancing terms\n");
	      exitElegant(1);
	    }
	    if (terms)
	      for (i=0; i<optimization_data->terms; i++) {
		if (optimization_data->termWeight[i])
		  optimization_data->termWeight[i] *= optimization_data->usersTermWeight[i]*sum/terms;
		else
		  optimization_data->termWeight[i] = optimization_data->usersTermWeight[i]*sum/terms;
	      }
	    balanceTerms = 0;
	    fprintf(stdout, "\nOptimization terms balanced.\n");
	    fflush(stdout);
	  }
        
	  /* compute and return quantity to be optimized */
	  if (optimization_data->terms) {
	    long i;
	    for (i=result=0; i<optimization_data->terms; i++)  {
	      rpn_clear();
	      result += (optimization_data->termValue[i]=optimization_data->termWeight[i]*rpn(optimization_data->term[i]));
	      if (rpn_check_error()) {
		printf("Problem evaluating expression: %s\n", optimization_data->term[i]);
		rpn_clear_error();
		rpnError++;
	      }
	    }
	  }
	  else {
	    rpn_clear();    /* clear rpn stack */
	    result = rpn(optimization_data->UDFname);
	    if (rpn_check_error()) {
	      printf("Problem evaluating expression: %s\n", optimization_data->term[i]);
	      rpnError++;
	    }
	  }
	  if (rpnError) {
	    printf("RPN expression errors prevent optimization\n");
	    exitElegant(1);
	  }
	  if (isnan(result) || isinf(result)) {
	    *invalid = 1;
	  } else {
#if !USE_MPI /* The information here is from a processor locally, we print the result across all the processors after an iteration */
	    if (optimization_data->verbose && optimization_data->fp_log) {
	      fprintf(optimization_data->fp_log, "equation evaluates to %23.15e\n", result);
	      fflush(optimization_data->fp_log);
	      if (optimization_data->terms && !*invalid) {
		fprintf(optimization_data->fp_log, "Terms of equation: \n");
		for (i=sum=0; i<optimization_data->terms; i++) {
		  rpn_clear();
		  fprintf(optimization_data->fp_log, "%g*(%20s): %23.15e\n",
			  optimization_data->termWeight[i],
			  optimization_data->term[i],
			  optimization_data->termValue[i]);
		  sum += optimization_data->termValue[i];
		}
	      }
	      fprintf(optimization_data->fp_log, "\n\n");
	    }
#endif
	  }
	}
      }  /* End of Master only */
    
      if (*invalid) {
	result = sqrt(DBL_MAX);
	if (*invalid && optimization_data->verbose && optimization_data->fp_log) {
	  fprintf(optimization_data->fp_log, "Result is invalid\n");
	  fflush(optimization_data->fp_log);
	}
    }
      
      if (optimization_data->mode==OPTIM_MODE_MAXIMUM)
	result *= -1;
      
      /* copy the result into the "hidden" slot in the varied quantities array for output
       * to final properties file
       */
      variables->varied_quan_value[variables->n_variables+1] = 
	optimization_data->mode==OPTIM_MODE_MAXIMUM?-1*result:result;
      if (!*invalid && bestResult>result) {
	if (optimization_data->verbose && optimization_data->fp_log)
	  fprintf(optimization_data->fp_log, "** Result is new best\n");
	bestResult = result;
      }
      variables->varied_quan_value[variables->n_variables+2] = 
	optimization_data->mode==OPTIM_MODE_MAXIMUM?-1*bestResult:bestResult;

#if USE_MPI
      if (notSinglePart) /* Disable the beam output when all the processors track independently */
#endif
	if (!*invalid && (force_output || (control->i_step-2)%output_sparsing_factor==0))
	  do_track_beam_output(run, control, error, variables, beamline, beam, output, optim_func_flags,
			   charge);
  }
  
  
#if DEBUG
    fprintf(stdout, "optimization_function: Returning %le,  invalid=%ld\n", result, *invalid);
    fflush(stdout);
#endif

#if USE_MPI
    if (notSinglePart) {
      MPI_Bcast(invalid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
#endif
 
    storeOptimRecord(value, variables->n_variables, *invalid, result);
    
    log_exit("optimization_function");
    return(result);
}

long checkForOptimRecord(double *value, long values, long *again)
{
  long iRecord, iValue;
  double diff;
  if (ignoreOptimRecords)
    return -1;
  for (iRecord=0; iRecord<optimRecords; iRecord++) {
    for (iValue=0; iValue<values; iValue++) {
      diff = fabs(value[iValue]-optimRecord[iRecord].variableValue[iValue]);
      if (diff!=0)
        break;
    }
    if (iValue==values)
      break;
  }
  if (iRecord==optimRecords)  {
    *again = 0;
    return -1;
  }
  optimRecord[iRecord].usedBefore += 1;
  *again = optimRecord[iRecord].usedBefore - 1;
  return iRecord;
}

void storeOptimRecord(double *value, long values, long invalid, double result)
{
  long i;
  if (ignoreOptimRecords)
    return ;
  for (i=0; i<values; i++)
    optimRecord[nextOptimRecordSlot].variableValue[i] = value[i];
  optimRecord[nextOptimRecordSlot].usedBefore = 0;
  optimRecord[nextOptimRecordSlot].invalid = invalid;
  optimRecord[nextOptimRecordSlot].result = result;
  if (++nextOptimRecordSlot>=MAX_OPTIM_RECORDS)
    nextOptimRecordSlot = 0;
  if (++optimRecords>=MAX_OPTIM_RECORDS)
    optimRecords = MAX_OPTIM_RECORDS;
}

void optimization_report(double result, double *value, long pass, long n_evals, long n_dim)
{
    OPTIM_VARIABLES *variables;
    OPTIM_COVARIABLES *covariables;
    OPTIM_CONSTRAINTS *constraints;
    long i;

#if (!USE_MPI)
    if (!optimization_data->fp_log)
        return;
#endif 

    /* update internal values (particularly rpn) */
    ignoreOptimRecords = 1 ;  /* force reevaluation */
    result = optimization_function(value, &i);
    ignoreOptimRecords = 0;

#if (USE_MPI)    
    if (!optimization_data->fp_log)
        return;
#endif  

    variables = &(optimization_data->variables);
    covariables = &(optimization_data->covariables);
    constraints = &(optimization_data->constraints);

#if !USE_MPI
    fprintf(optimization_data->fp_log, "Optimization pass %ld completed:\n    optimization function has value %23.15e\n", 
            pass, optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
#else
    fprintf(optimization_data->fp_log, "Optimization restart %ld completed:\n    optimization function has value %23.15e\n", 
            pass, optimization_data->mode==OPTIM_MODE_MAXIMUM?-result:result);
#endif
    n_passes_made = pass;
    if (optimization_data->terms) {
      double sum;
      fprintf(optimization_data->fp_log, "Terms of equation: \n");
      for (i=sum=0; i<optimization_data->terms; i++) {
        rpn_clear();
        fprintf(optimization_data->fp_log, "%g*(%20s): %23.15e\n",
                optimization_data->termWeight[i], optimization_data->term[i],
                rpn(optimization_data->term[i])*optimization_data->termWeight[i]);
        sum += rpn(optimization_data->term[i])*optimization_data->termWeight[i];
      }
    }
    
    if (constraints->n_constraints) {
        fprintf(optimization_data->fp_log, "    Constraints:\n");
        for (i=0; i<constraints->n_constraints; i++)
            fprintf(optimization_data->fp_log, "    %10s: %23.15e\n", constraints->quantity[i], 
                    final_property_value[constraints->index[i]]);
        }

    fprintf(optimization_data->fp_log, "    Values of variables:\n");
    for (i=0; i<n_dim; i++)
      fprintf(optimization_data->fp_log, "    %10s: %23.15e %s\n", variables->varied_quan_name[i], value[i],
	      variableAtLimit(value[i], variables->lower_limit[i], variables->upper_limit[i]) ? "(near limit)" : "");
    fflush(optimization_data->fp_log);
  }


void rpnStoreHigherMatrixElements(VMATRIX *M, long **TijkMem, long **UijklMem, long maxOrder)
{
  long order, i, j, k, l, count;
  char buffer[10];

  if (!*TijkMem && maxOrder>=2) {
    for (i=count=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++)
          count++;
    if (!(*TijkMem = malloc(sizeof(**TijkMem)*count)))
      bombElegant("memory allocation failure (rpnStoreHigherMatrixElements)", NULL);
    for (i=count=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++, count++) {
          sprintf(buffer, "T%ld%ld%ld", i+1, j+1, k+1);
          (*TijkMem)[count] = rpn_create_mem(buffer, 0);
        }
  }
  if (!*UijklMem && maxOrder>=3) {
    for (i=count=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++)
          for (l=0; l<=k; l++)
            count++;
    if (!(*UijklMem = malloc(sizeof(**UijklMem)*count)))
      bombElegant("memory allocation failure (rpnStoreHigherMatrixElements)", NULL);
    for (i=count=0; i<6; i++)
      for (j=0; j<6; j++)
        for (k=0; k<=j; k++)
          for (l=0; l<=k; l++, count++) {
            sprintf(buffer, "U%ld%ld%ld%ld", i+1, j+1, k+1, l+1);
            (*UijklMem)[count] = rpn_create_mem(buffer, 0);
          }
  }
    
  for (order=2; order<=maxOrder && order<=M->order; order++) {
    switch (order) {
    case 2:
      if (!M->T)
        bombElegant("second order matrix is missing!", NULL);
      for (i=count=0; i<6; i++)
        for (j=0; j<6; j++)
          for (k=0; k<=j; k++, count++)
            rpn_store(M->T[i][j][k], NULL, (*TijkMem)[count]);
      break;
    case 3:
      if (!M->Q)
        bombElegant("third order matrix is missing!", NULL);
      for (i=count=0; i<6; i++)
        for (j=0; j<6; j++)
          for (k=0; k<=j; k++)
            for (l=0; l<=k; l++, count++)
              rpn_store(M->Q[i][j][k][l], NULL, (*UijklMem)[count]);
      break;
    default:
      break;
    }
  }
}

#if USE_MPI
void find_global_min_index (double *min, int *processor_ID, MPI_Comm comm) {
    struct {
      double val;
      int rank;
    } in, out;

    in.val = *min;
    MPI_Comm_rank(comm, &(in.rank));
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
    *min = out.val;
    *processor_ID = out.rank;
}

void SDDS_PopulationSetup(char *population_log, SDDS_TABLE *popLogPtr, OPTIM_VARIABLES *optim, OPTIM_COVARIABLES *co_optim) {
  if (isMaster) {
    if (population_log && strlen(population_log)) {
      if (!SDDS_InitializeOutput(popLogPtr, SDDS_BINARY, 1, NULL, NULL, population_log) ||
	  !SDDS_DefineSimpleParameter(popLogPtr, "Iteration", NULL, SDDS_LONG) ||
	  !SDDS_DefineSimpleParameter(popLogPtr, "OptimizationValue", NULL, SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleParameters(popLogPtr, optim->n_variables, optim->varied_quan_name, 
				    optim->varied_quan_unit, SDDS_DOUBLE) ||
	  !SDDS_DefineSimpleParameters(popLogPtr, co_optim->n_covariables,co_optim->varied_quan_name, 
				    co_optim->varied_quan_unit, SDDS_DOUBLE) ||
	  !SDDS_WriteLayout(popLogPtr)) {
	    fprintf(stdout, "Problem setting up population output file %s\n", population_log);
	    fflush(stdout);
	    SDDS_PrintErrors(stderr, SDDS_VERBOSE_PrintErrors);
	    exitElegant(1);
      }
    }
  }
}

/* Function to print populations */
void SDDS_PrintPopulations(SDDS_TABLE *popLogPtr, long iteration, double best_value, double *best_individual, long dimensions, double *covariable, long n_covariables) {
  int i;
  long offset = 2;
  
  if (isMaster && log_file) {
    if (!SDDS_StartPage(popLogPtr, 1))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    if (!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "Iteration", iteration , NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "OptimizationValue", best_value, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    for(i=0; i<dimensions; i++) {
      if (!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, offset+i, best_individual[i], -1))
	SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }

    if (n_covariables) {
      for (i=0; i<n_covariables; i++) {
	if (!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, offset+dimensions+i, covariable[i], -1))
	  SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
      }
    }

    if (!SDDS_WritePage(popLogPtr))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    SDDS_DoFSync(popLogPtr);
  }
}
#endif
