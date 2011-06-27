/* file: randpert_multi.nl
 * contents: namelist and other stuff for randpert.c
 * 
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist perturbations static
    STRING type = "quadrupole";
    STRING name = NULL;
    STRING SDDS_output = NULL;
    STRING elegant_output = NULL;
    STRING kmult_output = NULL;
    double effective_length = 0.23;
    double bore_radius = 0.066;
    double dx_pole = 0;
    double dy_pole = 0;
    double dradius  = 0;
    double dx_split = 0;
    double dy_split = 0;
    double dphi_halves = 0;
    long n_cases = 1000;
    long n_harm = 0;
    long random_number_seed = 123456789;
#end

