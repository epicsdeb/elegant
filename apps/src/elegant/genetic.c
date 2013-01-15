/************************************************************************* \
* Copyright (c) 2010 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: genetic.c
 * The wrapper to call pga package for genetic optimization
 * Y. Wang, 2010
 */

#include "mdb.h"
#include "track.h"

static double target_value;
static double *initial_guess = NULL;
static char *log_file = NULL;
static SDDS_TABLE *popLogPtr = NULL;
static long print_all = 0;
static int last_reduced_iter = 0;
static long sparsing_factor = 1000000;

void   N_InitString(PGAContext *ctx, int p, int pop);
int    N_StopCond(PGAContext *ctx);
void   N_PrintString(PGAContext *, FILE *, int, int);
void   N_EndOfGen   (PGAContext *);
double evaluate (PGAContext *ctx, int p, int pop);

long geneticMin(
                double *yReturn, 
                double *xGuess,
                double *xLowerLimit,
                double *xUpperLimit,
		double *xStep,
                long dimensions,
		double target,
		double (*func)(double *x, long *invalid),
		long maxIterations,
		long maxNoChange,
		long populationSize,
		long printFrequency,
		long printAllPopulations,
		char *population_log,
		SDDS_TABLE *logPtr,
		long verbose,
		long crossover_type,
		OPTIM_VARIABLES *optim,
		OPTIM_COVARIABLES *co_optim
		)
{
  /* argc and argv are not used. They are provided here to make the interface consistent with pgapack */
  int argc = 0;
  char **argv = NULL;
  
  PGAContext *ctx;
  long totalEvaluations = 0;
  int best_p;

  /* As pgapack has fixed interface for some functions, these variables are used to passing information
     between functions within this file */
  initial_guess = xGuess;
  target_value = target;
  log_file = population_log;
  print_all = printAllPopulations;
  last_reduced_iter = 0;
  sparsing_factor = printFrequency;
  
  if (!popLogPtr) {	
     popLogPtr = logPtr;
     if (population_log)
        SDDS_PopulationSetup(population_log, popLogPtr, optim, co_optim);
  }

  ctx = PGACreate(&argc, argv, PGA_DATATYPE_REAL, dimensions, PGA_MINIMIZE);

  /* Select an initial value uniformly from the given range for each variable */
  PGASetRealInitRange(ctx, xLowerLimit, xUpperLimit);

  /* Stop when the maximal iteration limit exceeded */
  PGASetMaxGAIterValue(ctx, maxIterations);

  /* Set the frequency of statistics printing with the number of iterations as the interval*/
  if (verbose > 0)
    PGASetPrintFrequencyValue(ctx, 1); 
  else  /* No intermediate result will be printed */
    PGASetPrintFrequencyValue(ctx, 1000000); 

  /* Reseed a population from the best string */
  /*PGASetRestartFlag(ctx, PGA_TRUE); */

  /* Set the frequency for restart */
  /* PGASetRestartFrequencyValue(ctx, 1000); */	
 
  /* Set the size of population */
  PGASetPopSize(ctx, populationSize);

  /* Set the crossover type */
  if (crossover_type!=PGA_CROSSOVER_TWOPT)
    PGASetCrossoverType(ctx, crossover_type); 
  if (dimensions==2) /* It is not valid for PGA_CROSSOVER_TWOPT when the dimension is 2 */
    PGASetCrossoverType(ctx, PGA_CROSSOVER_ONEPT); 

  /* The number of genes to be replaced for each iteration */
  PGASetNumReplaceValue (ctx,  PGAGetPopSize(ctx)*0.8);	
/*  PGASetNumReplaceValue (ctx,  PGAGetPopSize(ctx)-1); */

  /* String duplication will not be allowed in the population */
  PGASetNoDuplicatesFlag(ctx, PGA_TRUE);	
  
  /* Stop when no change in the best solution found */
  PGASetStoppingRuleType(ctx, PGA_STOP_NOCHANGE);
  /* Set the maximum number of iterations in which no change in the best evaluation is allowed before the GA stops. */
  PGASetMaxNoChangeValue(ctx, maxNoChange); 

  /* PGASetMutationType(ctx, PGA_MUTATION_UNIFORM); */
  /* As GA does not have a way to control the step size for each direction, we use the step of the first dimension */
  if (xGuess[0] != 0) /* GA uses relative values for mutation */
    PGASetMutationRealValue(ctx, fabs(xStep[0]/xGuess[0]));
  else /* There could be a problem if one of the initial values is 0, because the relative value is used in GA */
    PGASetMutationRealValue(ctx, xStep[0]);

  /* The allele value is suppposed to be within the initialization range by default, but that is not the case */
  PGASetMutationBoundedFlag(ctx, PGA_TRUE);

#if MPI_DEBUG
  /* Specifying a seed exlicitly to allow for reproducibility of runs */
  PGASetRandomSeed(ctx, 123456789);
#endif

  /* User defined functions */
  PGASetUserFunction(ctx, PGA_USERFUNCTION_INITSTRING, (void *)N_InitString);
  PGASetUserFunction(ctx, PGA_USERFUNCTION_STOPCOND,   (void *)N_StopCond);
  /* This function should not be defined as the user function, as it will be called for
     every individual 
     PGASetUserFunction(ctx, PGA_USERFUNCTION_PRINTSTRING,(void *)N_PrintString); */
  PGASetUserFunction(ctx, PGA_USERFUNCTION_ENDOFGEN,   (void *)N_EndOfGen);

  PGASetUp(ctx);

  /*  print the value of all fields in the context variable */
  PGAPrintContextVariable(ctx, stdout);

  PGARun(ctx, evaluate);

  /* Get the best string from genetic optimizer */
  if (isMaster) {
    best_p = PGAGetBestIndex(ctx, PGA_OLDPOP);
    PGARealGetString (ctx, xGuess, best_p, PGA_OLDPOP);
    *yReturn = PGAGetEvaluation (ctx, best_p, PGA_OLDPOP);
  }	
#if USE_MPI
  MPI_Bcast(xGuess, dimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(yReturn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  totalEvaluations = (long)(PGAGetNumReplaceValue(ctx)*PGAGetGAIterValue(ctx));

  PGADestroy(ctx);

  return(totalEvaluations);
}

double evaluate (PGAContext *ctx, int p, int pop) 
{
  static long dimensions = 0; 
  static double *GAVector = NULL, yStart;
  static long initialized = 0;
  double *string, result;
  long direction, isInvalid;
  
  if (!dimensions)
    dimensions =  (long) PGAGetStringLength (ctx);
  
  if (!GAVector)
    GAVector = tmalloc(sizeof(*GAVector)*dimensions);

  string = (double*) PGAGetIndividual(ctx, p, pop)->chrom;

  for (direction=0; direction<dimensions; direction++) {
    GAVector[direction] = string[direction]; 
#if MPI_DEBUG
    printf ("GAVector[%ld]=%.15lf\n", direction, GAVector[direction]);
#endif
  } 

  result = optimization_function(GAVector, &isInvalid);

  if (!initialized) {
    if (isInvalid) { /* We simply let it quit now. It can be improved latter if necessary */
      fprintf(stderr, "error: initial guess is invalid in geneticMin()\n"); 
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    yStart = result;
    initialized = 1;
  }

#if MPI_DEBUG
  fprintf (stdout, "Iter: %d, result:%g, isInvalid:%ld, fabs(yStart)*100:%g on %d\n", PGAGetGAIterValue(ctx), result, isInvalid, fabs(yStart)*100, myid);
  fflush (stdout);
#endif
  if (isInvalid || (result > 100*yStart))
     /* There is a problem for PGA package if the returned number is too big. So we return fabs(yStart)*100 ,
	where ystart  is the evaluation result at the initial point. The optimal result should be 
	less than or equal to the initial value for minimization */
    return 100*fabs(yStart);
  else
    return result;
}

void N_InitString(PGAContext *ctx, int p, int pop) {
  int i, dimensions;

  dimensions = PGAGetStringLength(ctx);
  for(i=0; i<dimensions; i++)
    PGASetRealAllele(ctx, p, pop, i, initial_guess[i]);
}

int N_StopCond(PGAContext *ctx) {
    int   done, best;

    done = PGACheckStoppingConditions(ctx);

    best = PGAGetBestIndex(ctx, PGA_OLDPOP);
    if ((done == PGA_FALSE) && 
	(PGAGetEvaluation(ctx, best, PGA_OLDPOP) <= target_value))
      done = PGA_TRUE;

    if ((done == PGA_FALSE) &&
	(PGAGetMutationRealValue(ctx) <= 10*DBL_EPSILON)) /* PGA will hang there if the variance is too small */
      done = PGA_TRUE;

    return(done);
}

/* Function to print strings */
void N_PrintString(PGAContext *ctx, FILE *file, int best_p, int pop) {
  int i, j, dimensions, pop_size, worst_p, median_p;
  double best_value, *best_string, *string, worst_value, median, average, sigma;
  long offset = 6;
  static double *a = NULL;
  
  if (print_all)
    pop_size = PGAGetPopSize(ctx);
  else
    pop_size = 1;
  if (!SDDS_StartPage(popLogPtr, pop_size))
    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

  if (isMaster && log_file) {
    best_string = (double*) PGAGetIndividual(ctx, best_p, PGA_NEWPOP)->chrom;
    best_value = PGAGetEvaluation (ctx, best_p, PGA_NEWPOP);
   
    worst_p = PGAGetWorstIndex(ctx, PGA_NEWPOP);
    worst_value = PGAGetEvaluation (ctx, worst_p, PGA_NEWPOP);

    PGASortPop(ctx, PGA_NEWPOP);
    median_p = PGAGetSortedPopIndex (ctx, PGAGetPopSize(ctx)/2);
    median = PGAGetEvaluation (ctx, median_p, PGA_NEWPOP);

    if (!a)
      a = (double *) tmalloc(PGAGetPopSize(ctx)*sizeof(*a));
    for (i=0; i<PGAGetPopSize(ctx); i++) /* Put all the result in an array */
      a[i] = PGAGetEvaluation (ctx, i, PGA_NEWPOP);
    /* This is not right, as they are the fitness value */
    average  = PGAMean (ctx, a, PGAGetPopSize(ctx));
    sigma = PGAStddev (ctx, a, PGAGetPopSize(ctx), average);

    if (!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "Iteration", PGAGetGAIterValue(ctx), NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "OptimizationValue", best_value, NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "WorstOptimizationValue", worst_value, NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "MedianOptimizationValue", median, NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "AverageOptimizationValue", average, NULL) ||
	!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
			    "OptimizationValueSpread", sigma, NULL))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);

    dimensions = PGAGetStringLength(ctx);
    for(i=0; i<dimensions; i++) {
      if (!SDDS_SetParameters(popLogPtr, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, offset+i, best_string[i], -1))
	SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    }

    if (print_all) {
      for(j=0; j<pop_size; j++) {
	if (!SDDS_SetRowValues(popLogPtr, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, j, 0, PGAGetEvaluation (ctx, j, PGA_NEWPOP), -1))
	  SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	string = (double*) PGAGetIndividual(ctx, j, PGA_NEWPOP)->chrom;
	for(i=0; i<dimensions; i++) {
	  if (!SDDS_SetRowValues(popLogPtr, SDDS_SET_BY_INDEX|SDDS_PASS_BY_VALUE, j, i+1, string[i], -1))
	    SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
	}
      }
    }

    if (!SDDS_WritePage(popLogPtr))
      SDDS_PrintErrors(stderr, SDDS_EXIT_PrintErrors|SDDS_VERBOSE_PrintErrors);
    
    SDDS_DoFSync(popLogPtr);
  }
}

/*  After each generation, this routine is called. */
void N_EndOfGen(PGAContext *ctx) {
    int best;
    int current_iter = PGAGetGAIterValue(ctx);

    /* If the best value is unchanged for 10 iterations, the step size will be reduced */
    if((ctx->ga.ItersOfSame % 10 == 0) && current_iter) {
      PGASetMutationRealValue(ctx, 0.1*PGAGetMutationRealValue(ctx));
      last_reduced_iter = PGAGetGAIterValue(ctx);
    } else {
      if(((current_iter-last_reduced_iter)%100 == 0) && current_iter) {
	PGASetMutationRealValue(ctx, PGAGetMutationRealValue(ctx)/0.618);
      }
    }

    if (isMaster) {
      if (log_file) {
	if ((sparsing_factor >0) && !(current_iter % sparsing_factor)) {
	  best = PGAGetBestIndex(ctx, PGA_NEWPOP);
	  N_PrintString(ctx, stdout, best, PGA_NEWPOP);
	}
      }
    }
#if MPI_DEBUG
    if (MPI_DEBUG) {
      int i;
      for (i=0; i<PGAGetPopSize(ctx); i++) {
	printf ("Result for %dth individual: %lf, update:%d\n", i, PGAGetFitness(ctx, i, PGA_NEWPOP),PGAGetEvaluationUpToDateFlag (ctx, i, PGA_NEWPOP));
	PGAPrintIndividual (ctx, stdout, i, PGA_NEWPOP);
      }
      printf ("Best index: %d\n\n", PGAGetBestIndex(ctx, PGA_NEWPOP));
      printf ("Mutation value for the genetic algorithm: %lg. The number of iterations with the same result:%d\n", PGAGetMutationRealValue(ctx),ctx->ga.ItersOfSame);
    }
#endif
}
