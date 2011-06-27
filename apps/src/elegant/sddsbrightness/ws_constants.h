/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file ws_constants.h
 * define the contants used in sddsws.c
 
 $log ws_constants.h
 *
 */
#include "constants.h"

#ifndef WS_CONSTANTS_INCLUDED
double xc, yc, kc;

#define C_EVANG (h_mks*c_mks/e_mks*1.0e10)
#define C_MM_M 1.0e-3
#define C_M_MM 1.0e3
#define C_CM_ANG 1.0e8
#define C_MRAD_RAD 1.0e-3
#define C_RAD_MRAD 1.0e3
#define C_MA_A 1.0e-3
#define C_CM_M 1.0e-2
#define PIHALF (PI/2.0)
#define TWOPI (2.0*PI)
#define SQRT3 1.732050807568877293527446
#define FINE_STRUCTURE_CONST (1.0e19*e_mks*1.0e19*e_mks/(4.0*PI*epsilon_o*1.0e38*hbar_mks*c_mks))
#define PTOT_FAC (PI/3.0*e_mks/epsilon_o/(me_mev*me_mev)*1.0e6)
#define PD_FAC (21.0/(16*PI*me_mev*me_mev)*PTOT_FAC)
#define BW 1.0e-3
#define EPSK 1.0e-4
#define CL1 (e_mks/(2*PI*me_mks*c_mks)*C_CM_M)
#define CL2 (1.5/(me_mev*me_mev)*hbar_mks/me_mks*1.0e6)
#define CL3 (3.0*FINE_STRUCTURE_CONST/(4.0*PI*PI))
#define CL4 (SQRT3*FINE_STRUCTURE_CONST/(2.0*PI))
#define ECP_FRAC 1.0e-6
#define EC1_FRAC 1.0e-3
#define EC2_FRAC 1.0e-6
#define EC3_FRAC 1.0e-4


#define WS_CONSTANTS_INCLUDED 1
#endif
