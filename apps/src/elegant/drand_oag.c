
/* This routine is copied from the drand.c under the SDDS directory
   for parallelizing the random number generator. */

#include "mdb.h"
#include "track.h"

double dlaran_OAG(long int *iseed);

static long initialized = 0;
static long savedRandomNumberSeed[4] = {0,0,0,0};

void seedElegantRandomNumbers(long iseed, unsigned long restart)
{
  long i, offset;
  if (restart) {
    if (!initialized)
      bombElegant("seedElegantRandomNumbers called for restart but not initialized", NULL);
  } else {
    initialized = 1;
    savedRandomNumberSeed[0] = FABS(iseed);
    savedRandomNumberSeed[2] = FABS(iseed+4);
#if (!USE_MPI)
    savedRandomNumberSeed[1] = FABS(iseed+2);
    savedRandomNumberSeed[3] = FABS(iseed+6);
#else
    switch (mpiRandomizationMode) {
    case 2:
      /* Quadratic dependence of seed on processor ID */
      savedRandomNumberSeed[1] = FABS(iseed+2*myid*myid  ); 
      savedRandomNumberSeed[3] = FABS(iseed+2*myid*myid+4);
      break;
    case 4:
      /* Use system generator to create modified seed for each processor */
      srand(FABS(iseed));
      for (i=0; i<myid; i++)
	offset = 2*(rand()/2);
      if (offset>=RAND_MAX/4)
	offset /= 4;
      iseed = FABS(iseed);
      if (iseed>=RAND_MAX/2)
	iseed /= 2;
      iseed += offset;
      savedRandomNumberSeed[1] = iseed;
      iseed += offset;
      savedRandomNumberSeed[3] = iseed;
      break;
    case 1:
      /* Original method */
      savedRandomNumberSeed[1] = FABS(iseed+2*myid  );
      savedRandomNumberSeed[3] = FABS(iseed+2*myid+4);
      break;
    case 3:
    default:
      /* Quadratic dependence of seed on processor ID, but guaranteed to have an even offset */
      savedRandomNumberSeed[1] = FABS(iseed+myid*(myid+1)  );
      savedRandomNumberSeed[3] = FABS(iseed+myid*(myid+1)+4);
      break;
    }
#endif
  }

#if USE_MPI
  if (myid==0) {
    printf("Seeding random number generators (mode=%ld)\n", mpiRandomizationMode);
    fflush(stdout);
  }
#else
  printf("Seeding random number generators\n");
  fflush(stdout);
#endif

  /* seed random number generators.  
   * random_1_elegant is used for beamline errors, same on all processors
   * random_2 is used for random scraping/sampling/scattering
   * random_3 is used for BPM noise, same on all processors
   * random_4 is used for beam generation 
   */

  if (!restart || restart&RESTART_RN_BEAMLINE) 
    random_1_elegant(-savedRandomNumberSeed[0]);
  if (!restart || restart&RESTART_RN_SCATTER)
    random_2(-savedRandomNumberSeed[1]);
  if (!restart || restart&RESTART_RN_BPMNOISE)
    random_3(-savedRandomNumberSeed[2]);
  if (!restart || restart&RESTART_RN_BEAMGEN)
    random_4(-savedRandomNumberSeed[3]);
}

double random_1_elegant(long iseed)
{

    static short initialized = 0;
    static long int seed[4] = {0,0,0,0};

    if (!initialized || iseed<0) {
        if (iseed<0)
          iseed = -iseed;
	/* random_1_elegant() is used for beamline errors, same on all processors */
        seed[3] = ((iseed & 4095)/2)*2+1;
        seed[2] = (iseed >>= 12) & 4095;
        seed[1] = (iseed >>= 12) & 4095;
        seed[0] = (iseed >>= 12) & 4095;
        initialized = 1;
        }
    if (!initialized)
        bombElegant("random_1_elegant not properly initialized", NULL);

    return dlaran_OAG(seed);

}

double dlaran_OAG(long int *iseed)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static long int it1, it2, it3, it4;


/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARAN returns a random real number from a uniform (0,1)   
    distribution.   

    Arguments   
    =========   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array   
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    Further Details   
    ===============   

    This routine uses a multiplicative congruential method with modulus   
    2**48 and multiplier 33952834046453 (see G.S.Fishman,   
    'Multiplicative congruential random number generators with modulus   
    2**b: an exhaustive analysis for b = 32 and a partial analysis for   
    b = 48', Math. Comp. 189, pp 331-344, 1990).   

    48-bit integers are stored in 4 integer array elements with 12 bits   
    per element. Hence the routine is portable across machines with   
    integers of 32 bits or more.   

    =====================================================================   


       multiply the seed by the multiplier modulo 2**48   

       Parameter adjustments */
    --iseed;

    /* Function Body */
    it4 = iseed[4] * 2549;
    it3 = it4 / 4096;
    it4 -= it3 << 12;
    it3 = it3 + iseed[3] * 2549 + iseed[4] * 2508;
    it2 = it3 / 4096;
    it3 -= it2 << 12;
    it2 = it2 + iseed[2] * 2549 + iseed[3] * 2508 + iseed[4] * 322;
    it1 = it2 / 4096;
    it2 -= it1 << 12;
    it1 = it1 + iseed[1] * 2549 + iseed[2] * 2508 + iseed[3] * 322 + iseed[4] 
	    * 494;
    it1 %= 4096;

/*     return updated seed */

    iseed[1] = it1;
    iseed[2] = it2;
    iseed[3] = it3;
    iseed[4] = it4;

/*     convert 48-bit integer to a real number in the interval (0,1) */

    ret_val = ((double) it1 + ((double) it2 + ((double) it3 + (double) it4 * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4;
    return ret_val;

}

