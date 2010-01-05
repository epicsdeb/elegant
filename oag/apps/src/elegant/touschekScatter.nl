/* file: touschekScatter.nl
 * purpose: namelist for simulating Touschek scatter effects
 * 
 * Aimin Xiao, 2007
 */
#include "namelist.h"

#namelist touschek_scatter static
        long nbins = 100;
        double charge = 0;
        double frequency = 1;
        double delta = 0;
        double p_central_mev = 0.0;
        double emittance[2]  = {0, 0};
        double sigma_dp = 0.0;
        double sigma_s = 0.0;
        double distribution_cutoff[3] = {3, 3, 3};
        STRING bunch = NULL;
        STRING loss = NULL;
        STRING distribution = NULL;
        STRING initial = NULL;
        STRING output = NULL;
        long n_simulated = 5E6;
        long i_start = 0;
        long i_end = 1;
        long do_track = 0;
        double ignored_portion = 0.05;
        long verbosity = 0;
#end
/*
charge in C.
frequency in Hz
delta in % "energy aperture of interested"
emittance in m.rad
sigma_dp is the realative energy spread, sigma_s in m.
*/
