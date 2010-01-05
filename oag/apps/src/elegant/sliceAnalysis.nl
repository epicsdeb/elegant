/* file: sliceAnalysis.nl
 * contents: namelist for slice analysis vs s.
 * 
 * Michael Borland, 2002
 */
#include "namelist.h"

#namelist slice_analysis static
    STRING output = NULL;
    long n_slices = 0;
    double s_start = 0;
    double s_end = DBL_MAX;
    long final_values_only = 0;
#end


