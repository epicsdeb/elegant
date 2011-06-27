/* file: tune.nl
 * contents: namelist for tune correction
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist correct_tunes
    STRING quadrupoles = NULL;
    double tune_x = -1;
    double tune_y = -1;
    long n_iterations = 5;
    double correction_fraction = 0.9;
    double tolerance = 0;
    long step_up_interval = 0;
    double max_correction_fraction = 0.9;
    double delta_correction_fraction = 0.1;
    STRING strength_log = NULL;
    long change_defined_values = 0;
    long use_perturbed_matrix = 0;
#end


