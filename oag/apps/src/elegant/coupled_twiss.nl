/* file: coupled_sigma.nl
 * contents: namelist for coupled_sigma_matrix
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist coupled_twiss_output
    STRING filename = NULL;
    long output_at_each_step = 0;
    long emittances_from_twiss_command = 1;
    double emittance_ratio = 0.01;
    double emit_x = 0;
    double sigma_dp = 0;
    long calculate_3d_coupling = 1;
    long verbosity = 0;
    long concat_order = 2;
    long output_sigma_matrix = 0;
#end

