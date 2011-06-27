/* file: trace.nl
 * contents: namelist for program trace
 * 
 * Michael Borland, 1989
 */
#include "namelist.h"

#namelist trace static
    long traceback_on = 1;
    long trace_on = 1;
    long heap_verify_depth = 0;
    STRING filename = NULL;
    STRING memory_log = NULL;
    long record_allocation = 0;
#end


