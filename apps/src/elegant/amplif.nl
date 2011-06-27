/* file: amplif.nl
 * purpose: namelist for requesting amplification factors
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist amplification_factors static
    STRING output = NULL;
    STRING uncorrected_orbit_function = NULL;
    STRING corrected_orbit_function = NULL;
    STRING kick_function = NULL;
    double change = 1e-3;
    STRING name = NULL;
    STRING type = NULL;
    STRING item = NULL;
    STRING plane = NULL;
    long number_to_do = -1;
    double maximum_z = 0;
#end
