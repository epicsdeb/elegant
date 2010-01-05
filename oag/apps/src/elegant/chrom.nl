/* file: chrom.nl
 * contents: namelist for chromaticity correction
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist chromaticity
    STRING sextupoles = NULL;
    double dnux_dp = 0;
    double dnuy_dp = 0;
    double sextupole_tweek = 1e-3;
    double correction_fraction = 0.9;
    long n_iterations = 5;
    double tolerance = 0;
    STRING strength_log = NULL;
    long change_defined_values = 0;
    double strength_limit = 0;
    long use_perturbed_matrix = 0;    
    long exit_on_failure = 0;
    long verbosity = 2;
#end


