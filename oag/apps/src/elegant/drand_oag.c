
/* This routine is copied from the drand.c under the SDDS directory
   for parallelizing the random number generator. */

#include "mdb.h"
#include "track.h"

double dlaran_OAG(long int *iseed);

double random_1_elegant(long iseed)
{

    static short initialized = 0;
    static long int seed[4] = {0,0,0,0};

    if (!initialized || iseed<0) {
        if (iseed<0)
          iseed = -iseed;
#if (!USE_MPI)
        random_2(-(iseed+2));
        random_3(-(iseed+4));
	random_4(-(iseed+6));
#else
        /* Pelegant should give the same result as elegant if running on 2 processors. */
	random_2(-(iseed+2*myid)); 
	random_3(-(iseed+2*myid+2));
	random_4(-(iseed+2*myid+4));
#endif

        seed[3] = ((iseed & 4095)/2)*2+1;
        seed[2] = (iseed >>= 12) & 4095;
        seed[1] = (iseed >>= 12) & 4095;
        seed[0] = (iseed >>= 12) & 4095;
        initialized = 1;
        }
    if (!initialized)
        bomb("random_1_elegant not properly initialized", NULL);

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
