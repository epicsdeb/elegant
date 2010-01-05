/* file: closed_orbit.nl
 * contents: namelist for closed_orbit
 * 
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist closed_orbit
    STRING output = NULL;
    long start_from_centroid = 1;
    long start_from_dp_centroid = 0;
    double closed_orbit_accuracy = 1e-12;
    long closed_orbit_iterations = 40;
    long fixed_length = 0;
    long start_from_recirc = 0;
    long verbosity = 0;
    double iteration_fraction = 0.9;
    long output_monitors_only = 0;
#end

