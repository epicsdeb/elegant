/* file: correct.nl 
 * purpose: namelist definition for trajectory correction
 *
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist correct static
    STRING mode = "trajectory";
    STRING method = "global";
    STRING trajectory_output = NULL;
    STRING corrector_output = NULL;
    STRING statistics = NULL;
    double corrector_tweek[2] = {1e-6, 1e-6};
    double corrector_limit[2] = {0, 0};
    double correction_fraction[2] = {1, 1};
    double correction_accuracy[2] = {1e-6, 1e-6};
    long remove_smallest_SVs[2] = {0, 0};
    long keep_largest_SVs[2] = {0, 0};
    double minimum_SV_ratio[2] = {0, 0};
    long auto_limit_SVs[2] = {1, 1};
    long threading_divisor[2] = {100, 100};
    double bpm_noise[2] = {0, 0};
    double bpm_noise_cutoff[2] = {1.0, 1.0};
    STRING bpm_noise_distribution[2] = {"uniform", "uniform"};
    long verbose = 1;
    long fixed_length = 0;
    long fixed_length_matrix = 0;
    long n_xy_cycles = 1;
    long minimum_cycles = 1;
    long n_iterations = 1;
    long prezero_correctors = 1;
    long track_before_and_after = 0;
    long start_from_centroid = 1;
    long use_actual_beam = 0;
    double closed_orbit_accuracy = 1e-12;
    long closed_orbit_iterations = 40;
    double closed_orbit_iteration_fraction = 0.9;
    long use_perturbed_matrix = 0;
#end

