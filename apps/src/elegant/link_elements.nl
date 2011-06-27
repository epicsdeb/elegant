/* file: link_elements.nl
 * purpose: namelist for linking elements
 * 
 * Michael Borland, 1991
 */
#include "namelist.h"

#namelist link_control
    long clear_links = 1;
    long summarize_links = 0;
    long verbosity = 0;
#endif
    
#namelist link_elements static
    STRING target = NULL;
    STRING exclude = NULL;
    STRING item = NULL;
    STRING source = NULL;
    STRING source_position = NULL;
    STRING mode = NULL;
    STRING equation = NULL;
    double minimum = -DBL_MAX;
    double maximum = DBL_MAX;
#end

#define SRC_POSITION_BEFORE 0
#define SRC_POSITION_AFTER 1
#define SRC_POSITION_ADJACENT 2
#define SRC_POSITION_NEAREST 3
#define SRC_POSITION_SAME_OCCURENCE 4
#define N_SRC_POSITIONS 5
static char *src_position_name[N_SRC_POSITIONS] = {
    "before", "after", "adjacent", "nearest", "same-occurence"
    } ;

#define N_LINK_MODES 4
static char *link_mode[N_LINK_MODES] = {
        "static", "dynamic", "post-correction", "turn-by-turn"
        } ;
static long link_mode_flag[N_LINK_MODES] = {
        STATIC_LINK, DYNAMIC_LINK, POST_CORRECTION_LINK, TURN_BY_TURN_LINK
        }  ;
#define LINK_FLAG_MASK (STATIC_LINK+DYNAMIC_LINK+POST_CORRECTION_LINK+TURN_BY_TURN_LINK)

