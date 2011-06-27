/* file: sasefel.nl
 * contents: namelist for sase fel computations
 * 
 * Michael Borland, 1999
 */
#include "namelist.h"

#namelist sasefel static
    STRING output = NULL;
    STRING model = "Ming Xie";
    STRING beamsize_mode = "geometric mean";
    double beta = 0;
    double undulator_K = 3.1;
    double undulator_period = 0.033;
    double slice_fraction = 0.0;
    long n_slices = 0;
#end


