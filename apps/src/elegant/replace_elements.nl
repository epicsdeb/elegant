/* file: replace_elements.nl
 * purpose: namelist for replace or delete elements along beamline
 * 
 * A.Xiao, 2008
 */
#include "namelist.h"

#namelist replace_elements static
        STRING name = NULL;
        STRING type = NULL;
        STRING exclude = NULL;
        long skip = 1;
        long disable = 0;
        STRING element_def = NULL;
        long total_occurrences = 0;
        long occurrence[100]={0};
        long verbose = 0;
#end
