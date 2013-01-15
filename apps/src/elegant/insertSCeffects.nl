/* file: sceff.nl
 * purpose: namelist for tracking with space charge effects
 * 
 * Aimin Xiao, 2005
 */
#include "namelist.h"

#namelist insert_sceffects static
        STRING name = NULL;
        STRING type = NULL;
        STRING exclude = NULL;
        long disable = 0;
        long clear = 0;
        STRING element_prefix = "MYSC";
        long skip = 0;
        long vertical = 0;
        long horizontal = 0;
        long longitudinal = 0;
        long nonlinear = 0;
	long uniform_distribution = 0;
#end
