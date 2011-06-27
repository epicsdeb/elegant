/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "track.h"

void sdds_strength_output(char *outputfile, LINE_LIST *beamline, char *input);
void convert_to_patpet(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance,
                    char *header_file, char *ender_file);
void convert_to_parmela(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance, char *header_file, char *ender_file,
    double quad_aperture, double sext_aperture, double pc);
void convert_to_patricia(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance,
                    char *header_file, char *ender_file);
void convert_to_transport(char *output, LINE_LIST *beamline, long flip_k, double angle_tolerance, char *header_file, char *ender_file,
    double quad_aperture, double sext_aperture, double pc);
void convert_to_EmmaMatlab(char *outputfile, LINE_LIST *beamline, char *header_file, char *ender_file);
void convert_to_cosy(char *outputfile, LINE_LIST *beamline, long cosyOrder, double pCentral,
                     double quadBoreRadius, double sextBoreRadius);
void do_output_transport(FILE *fp, char *s);
void section(char *s, int n);
long nearly_equal(double x1, double x2, double frac);
char *quoted_label(char *s);
void do_output(FILE *fp, char *s);
