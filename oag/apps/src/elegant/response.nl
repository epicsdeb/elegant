/* file: response.nl 
 * purpose: namelist definition for trajectory/orbit response matrix output
 *
 * Michael Borland, 1993
 */
#include "namelist.h"

#namelist correction_matrix_output static
    STRING response[2] = {NULL, NULL};
    STRING inverse[2] = {NULL, NULL};
    long KnL_units = 0;
    long BnL_units = 0;
    long output_at_each_step = 0;
    long output_before_tune_correction = 0;
    long fixed_length = 0;
#end

