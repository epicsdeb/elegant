/* file: vary.nl
 * purpose: namelists for varying beamline and controlling tracking
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist run_control static
    long n_steps = 1;
    double bunch_frequency = 0;
    long n_indices = 0;
    long n_passes = 1;
    long reset_rf_for_each_step = 1;
    long first_is_fiducial = 0;
    long restrict_fiducialization = 0;
#end

#namelist vary_element static
    long index_number = 0;
    long index_limit = 0;
    STRING name = NULL;
    STRING item = NULL;
    double initial = 0;
    double final = 0;
    long differential = 0;
    long geometric = 0;
    long multiplicative = 0;
    STRING enumeration_file = NULL;
    STRING enumeration_column = NULL;
#end
