/* file: run_rpnexpr.nl
 * purpose: namelist for executing rpn expressions
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist rpn_expression static
    STRING expression = NULL;
#end

#namelist rpn_load static
    STRING filename = NULL;
    STRING match_column = NULL;
    STRING match_column_value = NULL;
    long matching_row_number = -1;
    STRING match_parameter = NULL;
    STRING match_parameter_value = NULL;
    long use_row = -1;
    long use_page = -1;
    long load_parameters = 0;
    STRING tag = NULL;
#end
    
