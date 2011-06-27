/* file: analyze.nl
 * purpose: namelist definition for transport analysis
 *
 * M.Borland, 1992 
 */
#include "namelist.h"

#namelist analyze_map static
    STRING output = NULL;
    double delta_x = 1e-6;
    double delta_xp = 1e-6;
    double delta_y = 1e-6;
    double delta_yp = 1e-6;
    double delta_s  = 1e-6;
    double delta_dp = 1e-6;
    long center_on_orbit = 0;
    long verbosity = 0;
    long n_points = 2;
#end

