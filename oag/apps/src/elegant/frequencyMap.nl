/* file: frequencyMap.nl
 * purpose: namelist definition for frequency map analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist frequency_map static
    STRING output = NULL;
    double xmin =  -0.1;
    double xmax =  0.1;
    double ymin =  1e-6;
    double ymax =  0.1;
    long nx = 20;
    long ny = 21;
    long verbosity = 1;
    long include_changes = 0;
#end

