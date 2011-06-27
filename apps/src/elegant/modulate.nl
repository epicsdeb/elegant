/* file: error.nl
 * purpose: namelist for random errors
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist modulate_elements static
    STRING name = NULL;
    STRING item = NULL;
    STRING type = NULL;
    STRING expression = NULL;
    STRING filename = NULL;
    STRING time_column = NULL;
    STRING amplitude_column = NULL;
    long differential = 1;
    long multiplicative = 0;
    long start_occurence = 0;
    long end_occurence = 0;
    double s_start = -1;
    double s_end = -1;
    STRING before = NULL;
    STRING after = NULL;
    long verbose = 0;
#end

