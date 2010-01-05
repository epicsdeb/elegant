/* file: sdds_beam.nl
 * contents: namelist sdds beam input
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist sdds_beam
    STRING input = NULL;
    STRING input_list = NULL;
    STRING input_type = "elegant";
    STRING selection_parameter = NULL;
    STRING selection_string = NULL;
    long one_random_bunch = 0;
    long n_particles_per_ring = 0;
    long reuse_bunch = 0;
    long prebunched = 0;
    long sample_interval = 1;
    long n_tables_to_skip = 0;
    long center_transversely = 0;
    long center_arrival_time = 0;
    long reverse_t_sign = 0;
    double sample_fraction = 1;
    double p_lower = 0.0;
    double p_upper = 0.0;
    long save_initial_coordinates = 1;
#end

