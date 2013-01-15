/* file: fitTraces.nl
 * purpose: namelists for fitting multiple traces through a beamline
 * 
 * Michael Borland, 1997
 */
/*
 * $Log: not supported by cvs2svn $
 * Revision 1.5  1998/06/08 16:50:43  borland
 * Fixed problems with SVD.  Added more reporting.
 *
 * Revision 1.4  1998/04/17 22:15:37  borland
 * Use Meschach matrix library.  Supports SVD or non-SVD.  Removes BPM
 * common mode.
 *
 * Revision 1.3  1998/03/19 21:00:57  borland
 * Removed items from track.h and put them in three new files (correctDefs.h,
 * tuneDefs.h, chromDefs.h) to isolate references to matlib.h; this is to
 * allow use of the Meschach matrix library in parts of elegant.
 *
 * Revision 1.2  1997/10/20 14:57:11  borland
 * Improved trace fitting and related routines.  Added output of traces
 * after fitting.  Fixed some output-related bugs.
 *
 * Revision 1.1  1997/08/13 20:03:45  borland
 * First version.
 *
 */

#include "namelist.h"

#namelist fit_traces static
        STRING trace_data_file = NULL;
        STRING fit_parameters_file = NULL;
        long iterations = 10;
        long sub_iterations = 10;
        long use_SVD = 0;
        long SVs_to_keep = 0;
        long SVs_to_remove = 0;
        double BPM_threshold = 1.75e-6;
        double convergence_factor = 1;
        double convergence_factor_divisor = 2;
        double convergence_factor_multiplier = 1.5;
        long convergence_increase_steps = 5;
        double convergence_factor_min = 1.0;
        double convergence_factor_max = 2.0;
        double position_change_limit = 0;
        double slope_change_limit = 0;
        long trace_sub_iterations = 500;
        double trace_convergence_factor = 0.9;
        double trace_fractional_target = 1e-6;
        STRING fit_output_file = NULL;
        STRING trace_output_file = NULL;
        double target = 1e-12;
        double tolerance = 1e-12;
        long reject_BPM_common_mode = 1;
        long n_restarts = 1;
        double restart_randomization_level = 0.01;
#end

