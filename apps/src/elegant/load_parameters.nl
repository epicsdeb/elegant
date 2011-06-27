/* file: load_parameters.nl
 * purpose: namelist for loading parameters from external file
 * 
 * Michael Borland, 1993
 */
#include "namelist.h"

#namelist load_parameters
        STRING filename = NULL;
        STRING filename_list = NULL;
        STRING include_name_pattern = NULL;
        STRING include_item_pattern = NULL;
        STRING include_type_pattern = NULL;
        STRING exclude_name_pattern = NULL;
        STRING exclude_item_pattern = NULL;
        STRING exclude_type_pattern = NULL;
        long change_defined_values = 0;
        long clear_settings = 0;
        long allow_missing_files = 0;
        long allow_missing_elements = 0;
        long allow_missing_parameters = 0;
        long force_occurence_data = 0;
        long verbose = 0;
        long use_first = 0;
#end

