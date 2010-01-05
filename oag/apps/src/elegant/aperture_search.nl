/* file: aperture_search.nl
 * purpose: namelist definition for aperture search
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist find_aperture static
    STRING output = NULL;
    STRING search_output = NULL;
    STRING boundary = NULL;
    STRING mode = "many-particle";
    double xmin = -0.1;
    double xmax =  0.1;
    double ymin =  0.0;
    double ymax =  0.1;
    long nx  = 21;
    long n_splits = 0;
    double split_fraction = 0.5;
    double desired_resolution = 0.01;
    long ny  = 11;
    long verbosity = 0;    
    long assume_nonincreasing = 0;
    long offset_by_orbit = 0;
    long n_lines = 11;
    long optimization_mode = 0;
#end

