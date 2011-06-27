/* file: optim_covariable.nl
 * purpose: namelist for optimizing beamline
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist optimization_covariable static
    STRING name = NULL;
    STRING item = NULL;
    STRING equation = NULL;
    long disable = 0;
#end
