#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.1 2007-03-30 16:50:27 soliday Exp $
#
# $Log: not supported by cvs2svn $
#
#

TOP=../../..
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/elegantTools/Makefile.OAG
include $(TOP)/configure/RULES

sddsrandmult$(OBJ): sddsrandmult.h

sddsrandmult.h: ../sddsrandmult.nl
	nlpp ../sddsrandmult.nl sddsrandmult.h

clean::
	$(RM) sddsrandmult.h 