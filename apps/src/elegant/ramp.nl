/* file: ramp.nl
 * purpose: namelist for ramping
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist ramp_elements static
    STRING name = NULL;
    STRING item = NULL;
    STRING type = NULL;
    long start_pass = 0;
    long end_pass = LONG_MAX;
    double start_value = 0;
    double end_value = 0;
    long differential = 1;
    long multiplicative = 0;
    long start_occurence = 0;
    long end_occurence = 0;
    double exponent = 1;
    double s_start = -1;
    double s_end = -1;
    STRING before = NULL;
    STRING after = NULL;
    long verbose = 0;
#end

