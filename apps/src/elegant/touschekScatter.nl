/* file: touschekScatter.nl
 * purpose: namelist for simulating Touschek scatter effects
 * 
 * Aimin Xiao, 2007
 */
#include "namelist.h"

#namelist touschek_scatter static
        double charge = 0;
        double frequency = 1;
        double emit_x = 0;
        double emit_nx = 0;
        double emit_y = 0;
        double emit_ny = 0;
        double sigma_dp = 0;
        double sigma_s = 0;
        double distribution_cutoff[3] = {3, 3, 3};
        double Momentum_Aperture_scale = 0.85;
        STRING Momentum_Aperture = NULL;
        STRING XDist = NULL;
        STRING YDist = NULL;
        STRING ZDist = NULL;
        STRING TranDist = NULL;
        STRING FullDist = NULL;
        STRING bunch = NULL;
        STRING loss = NULL;
        STRING distribution = NULL;
        STRING initial = NULL;
        STRING output = NULL;
        long nbins = 100;
        long n_simulated = 5E6;
        double ignored_portion = 0.01;
        long i_start = -1;
        long i_end = -1;
	long match_position_only = 0;
        long do_track = 0;
        long verbosity = 0;
	long overwrite_files = 1;
#end
/*
charge in C.
frequency in Hz
delta in % "energy aperture of interested"
emittance in m.rad
sigma_dp is the realative energy spread, sigma_s in m.
*/
