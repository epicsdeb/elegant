/* file: momentumAperture.nl
 * contents: namelist for momentum aperture vs s
 * 
 * Michael Borland, 2006
 */
#include "namelist.h"

#namelist momentum_aperture static
    STRING output = NULL;
    double x_initial = 1e-5;
    double y_initial = 1e-5;
    double delta_negative_limit = -0.10;
    double delta_positive_limit = 0.10;
    double delta_negative_start = 0.0;
    double delta_positive_start = 0.0;
    double delta_step_size = 0.0025;
    long splits = 1;
    long steps_back = 1;
    long split_step_divisor = 10;
    long skip_elements = 0;
    long process_elements = 2147483647;
    double s_start = 0;
    double s_end = DBL_MAX;
    STRING include_name_pattern = NULL;
    STRING include_type_pattern = NULL;
    long fiducialize = 0;
    long verbosity = 1;        
    long soft_failure = 0;
    long output_mode = 0;
    long forbid_resonance_crossing = 0;
#end


