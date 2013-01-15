#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.18 2010-10-04 22:28:06 borland Exp $
#
#  Lowest Level Directroy Makefile
# $Log: not supported by cvs2svn $
# Revision 1.17  2010/03/24 14:29:27  borland
# Added modulate_elements command.
#
# Revision 1.16  2008/10/31 14:20:00  xiaoam
# Add replace_elements command to elegant
#
# Revision 1.15  2008/03/18 16:45:08  xiaoam
# Add insert_elements and touschekScatter into elegant
#
# Revision 1.14  2008/01/21 17:38:41  borland
# Added moments.c and moments.nl.
#
# Revision 1.13  2007/05/21 16:43:57  ywang25
# Added hash table functions according to Bob Jenkins's code in public domain. Customized for elegant to improve parameter loading.
#
# Revision 1.12  2007/03/30 16:55:00  soliday
# Removed a bunch of scripts and programs that have been moved to elegantTools.
#
# Revision 1.11  2007/02/08 17:09:54  ywang25
# Impoved Makefile to build Pelegant and elegant more conveniently.
#
# Revision 1.10  2006/03/23 00:05:45  borland
# Added coupled twiss parameter computation using code by V. Sajaev.
# Added momentum aperture computation to elegantRingAnalysis.
#
# Revision 1.9  2005/11/28 22:07:09  borland
# Added aperture input via an SDDS file using the new aperture_input command.
#
# Revision 1.8  2005/11/22 23:21:19  borland
# Added momentum aperture search, which necessitated adding an argument to
# do_tracking, resulting in changes in many files.
# Also improved convergence of orbit finder, adding a second iteration using
# tracking if the matrix-based method fails.
#
# Revision 1.7  2005/11/04 16:27:05  borland
# Added Xiao's code for space charge element SCMULT.
#
# Revision 1.6  2005/09/29 20:50:43  ywang25
# Modifications for parallelization.  Verified to be identical to sequential version when run in non-MPI mode.
#
# Revision 1.5  2004/04/08 16:09:36  soliday
# Build rules are now compatible with Base 3.14
#
# Revision 1.4  2002/08/14 20:23:30  soliday
# Added Open License
#
# Revision 1.3  1999/08/05 15:24:20  soliday
# Now uses Makefile.Host
#
# Revision 1.2  1997/10/20 14:57:02  borland
# Improved trace fitting and related routines.  Added output of traces
# after fitting.  Fixed some output-related bugs.
#
# Revision 1.1.1.1  1996/02/01  16:30:34  borland
# Imported files
#
#
#

TOP=../..
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/Makefile.OAG
include $(TOP)/configure/RULES

alter$(OBJ): alter.h

alter.h: ../alter.nl
	nlpp ../alter.nl alter.h

amplif$(OBJ): amplif.h

amplif.h: ../amplif.nl
	nlpp ../amplif.nl amplif.h

analyze$(OBJ): analyze.h

analyze.h: ../analyze.nl
	nlpp ../analyze.nl analyze.h

aperture_search$(OBJ): aperture_search.h

aperture_search.h: ../aperture_search.nl
	nlpp ../aperture_search.nl aperture_search.h

bunched_beam$(OBJ): bunched_beam.h

bunched_beam.h: ../bunched_beam.nl
	nlpp ../bunched_beam.nl bunched_beam.h

chrom$(OBJ): chrom.h

chrom.h: ../chrom.nl
	nlpp ../chrom.nl chrom.h

closed_orbit$(OBJ): closed_orbit.h

closed_orbit.h: ../closed_orbit.nl
	nlpp ../closed_orbit.nl closed_orbit.h

correct$(OBJ): correct.h steer_elem.h

correct.h: ../correct.nl
	nlpp ../correct.nl correct.h

coupled_twiss$(OBJ): coupled_twiss.h

coupled_twiss.h: ../coupled_twiss.nl
	nlpp ../coupled_twiss.nl coupled_twiss.h

divideElements$(OBJ): divideElements.h

divideElements.h: ../divideElements.nl
	nlpp ../divideElements.nl divideElements.h

transmuteElements$(OBJ): transmuteElements.h

transmuteElements.h: ../transmuteElements.nl
	nlpp ../transmuteElements.nl transmuteElements.h

elegant$(OBJ): elegant.h aperture_data.h

elegant.h: ../elegant.nl
	nlpp ../elegant.nl elegant.h

aperture_data.h: ../aperture_data.nl
	nlpp ../aperture_data.nl aperture_data.h

error$(OBJ): error.h

error.h: ../error.nl
	nlpp ../error.nl error.h

floor$(OBJ): floor.h

floor.h: ../floor.nl
	nlpp ../floor.nl floor.h

frequencyMap$(OBJ): frequencyMap.h

frequencyMap.h: ../frequencyMap.nl
	nlpp ../frequencyMap.nl frequencyMap.h

hashtab$(OBJ) : hashtab.c standard.h recycle.h lookupa.h hashtab.h

insertSCeffects$(OBJ): insertSCeffects.h

insertSCeffects.h: ../insertSCeffects.nl
	nlpp ../insertSCeffects.nl insertSCeffects.h

insert_elements$(OBJ): insert_elements.h

insert_elements.h: ../insert_elements.nl
	nlpp ../insert_elements.nl insert_elements.h

replace_elements$(OBJ): replace_elements.h

replace_elements.h: ../replace_elements.nl
	nlpp ../replace_elements.nl replace_elements.h

touschekScatter$(OBJ): touschekScatter.h

touschekScatter.h: ../touschekScatter.nl
	nlpp ../touschekScatter.nl touschekScatter.h

link_elements$(OBJ): link_elements.h

link_elements.h: ../link_elements.nl
	nlpp ../link_elements.nl link_elements.h

lookupa$(OBJ) : lookupa.c standard.h lookupa.h

load_parameters$(OBJ): load_parameters.h

load_parameters.h: ../load_parameters.nl
	nlpp ../load_parameters.nl load_parameters.h

matrix_output$(OBJ): matrix_output.h

matrix_output.h: ../matrix_output.nl
	nlpp ../matrix_output.nl matrix_output.h

modulate$(OBJ): modulate.h

modulate.h: ../modulate.nl
	nlpp ../modulate.nl modulate.h

momentumAperture$(OBJ): momentumAperture.h

momentumAperture.h: ../momentumAperture.nl
	nlpp ../momentumAperture.nl momentumAperture.h

moments$(OBJ): moments.h

moments.h: ../moments.nl
	nlpp ../moments.nl moments.h

optim_covariable$(OBJ): optim_covariable.h

optim_covariable.h: ../optim_covariable.nl
	nlpp ../optim_covariable.nl optim_covariable.h

optimize$(OBJ): optimize.h optim_covariable.h

optimize.h: ../optimize.nl
	nlpp ../optimize.nl optimize.h

ramp$(OBJ): ramp.h

ramp.h: ../ramp.nl
	nlpp ../ramp.nl ramp.h

recycle$(OBJ) : recycle.c standard.h recycle.h

response$(OBJ): response.h

response.h: ../response.nl
	nlpp ../response.nl response.h

run_rpnexpr$(OBJ): run_rpnexpr.h

run_rpnexpr.h: ../run_rpnexpr.nl
	nlpp ../run_rpnexpr.nl run_rpnexpr.h

sasefel$(OBJ): sasefel.h

sasefel.h: ../sasefel.nl
	nlpp ../sasefel.nl sasefel.h

get_beamline$(OBJ): save_lattice.h

save_lattice$(OBJ): save_lattice.h

save_lattice.h: ../save_lattice.nl
	nlpp ../save_lattice.nl save_lattice.h

sdds_beam$(OBJ): sdds_beam.h

sdds_beam.h: ../sdds_beam.nl
	nlpp ../sdds_beam.nl sdds_beam.h

sliceAnalysis$(OBJ): sliceAnalysis.h

sliceAnalysis.h: ../sliceAnalysis.nl
	nlpp ../sliceAnalysis.nl sliceAnalysis.h

steer_elem$(OBJ): steer_elem.h

steer_elem.h: ../steer_elem.nl
	nlpp ../steer_elem.nl steer_elem.h

subprocess$(OBJ): subprocess.h

subprocess.h: ../subprocess.nl
	nlpp ../subprocess.nl subprocess.h

trace$(OBJ): trace.h

trace.h: ../trace.nl
	nlpp ../trace.nl trace.h

tune$(OBJ): tune.h

tune.h: ../tune.nl
	nlpp ../tune.nl tune.h

twiss$(OBJ): twiss.h

twiss.h: ../twiss.nl
	nlpp ../twiss.nl twiss.h

vary$(OBJ): vary.h

vary.h: ../vary.nl
	nlpp ../vary.nl vary.h

fitTraces$(OBJ): fitTraces.h

fitTraces.h: ../fitTraces.nl
	nlpp ../fitTraces.nl fitTraces.h

elegantLocation = $(wildcard ../../bin/linux-x86/elegant elegant)
PelegantLocation = $(wildcard ../../bin/linux-x86/Pelegant Pelegant)

elegant = $(words $(notdir $(elegantLocation)))
Pelegant = $(words $(notdir $(PelegantLocation)))


PelegantNewer=0
ifeq ($(Pelegant), 1)
PelegantTime=$(shell stat --format=%Y $(PelegantLocation))
ifeq ($(elegant), 1)
elegantTime=$(shell stat --format=%Y $(elegantLocation))
PelegantNewer=$(shell rpnl "$(elegantTime) $(PelegantTime) < ? 1 : 0 $$")
else
PelegantNewer=1
endif
endif

ifeq ($(PelegantNewer), 0)
Pelegant:
	$(RM) *.o O.linux-x86/*.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) MPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

Selegant:
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

all buildInstall:
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

else
Pelegant:
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) MPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

Selegant:
	$(RM) *.o O.linux-x86/*.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

all buildInstall:
	$(RM) *.o O.linux-x86/*.o
	$(MV) $(TOP)/src/elegant/Makefile $(TOP)/src/elegant/Makefile.TMP
	$(CP) $(TOP)/src/elegant/Makefile.Pelegant $(TOP)/src/elegant/Makefile
	$(MAKE) NOMPI=1 -f Makefile
	$(RM) $(TOP)/src/elegant/Makefile
	$(MV) $(TOP)/src/elegant/Makefile.TMP $(TOP)/src/elegant/Makefile

endif



clean::
	$(RM) fitTraces.h vary.h twiss.h tune.h trace.h subprocess.h steer_elem.h sliceAnalysis.h sdds_beam.h save_lattice.h sasefel.h run_rpnexpr.h response.h optimize.h optim_covariable.h matrix_output.h load_parameters.h link_elements.h frequencyMap.h floor.h error.h elegant.h transmuteElements.h divideElements.h correct.h steer_elem.h closed_orbit.h chrom.h bunched_beam.h aperture_search.h analyze.h amplif.h alter.h insertSCeffects.h insert_elements.h touschekScatter.h momentumAperture.h aperture_data.h replace_elements.h modulate.h ramp.h



