#*************************************************************************
# Copyright (c) 2002 The University of Chicago, as Operator of Argonne
# National Laboratory.
# Copyright (c) 2002 The Regents of the University of California, as
# Operator of Los Alamos National Laboratory.
# This file is distributed subject to a Software License Agreement found
# in the file LICENSE that is included with this distribution. 
#*************************************************************************
#
# $Id: Makefile,v 1.14 2011-02-08 03:24:07 soliday Exp $
#
#
#

TOP=../../..
include $(TOP)/configure/CONFIG
include $(TOP)/src/elegant/sddsbrightness/Makefile.OAG
include $(TOP)/configure/RULES

sddsbrightness$(OBJ): ../sddsbrightness.h


ifdef WIN32
usb$(OBJ) :
	ifort /compile_only /threads /libs:dll ..\\usb.f

brighte$(OBJ) :
	ifort /compile_only /threads /libs:dll ..\\brighte.f

hunt$(OBJ) :
	ifort /compile_only /threads /libs:dll ..\\hunt.f

urgent$(OBJ):
	ifort /compile_only /threads /libs:dll ..\\urgent.f

us$(OBJ):
	ifort /compile_only /threads /libs:dll ..\\us.f
else
F77 = g77 
ifdef AIX
F77 = g77
endif
ifeq ($(OS_CLASS),Darwin)
F77 = /sw/bin/gfortran
#F77 = /usr/local/bin/gfortran -arch i386 -arch x86_64
endif
ifeq ($(HOST_ARCH),solaris)
F77 = g77
solaris_libs = -lcurses -lrt
endif
ifeq ($(HOST_ARCH),solaris-sparc)
F77 = g77
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-sparc)
F77 = g77
solaris_libs = -lcurses -lrt
ifdef RATIONAL
F77 = $(RATIONAL) f77 -g
endif
endif
ifeq ($(EPICS_HOST_ARCH),solaris-sparc64)
F77 = g77 -m64
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86)
F77 = f77
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86-gnu)
F77 = f77
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),solaris-x86_64)
F77 = f77 -m64
solaris_libs = -lcurses -lrt
endif
ifeq ($(EPICS_HOST_ARCH),linux-x86_64)
F77 = gfortran -m64 -ffixed-line-length-132
endif
ifeq ($(EPICS_HOST_ARCH),linux-x86)
F77 = gfortran -m32 -ffixed-line-length-132
endif

#RATIONAL = purify
#RATIONAL = quantify
#RATIONAL = purecov
ifdef RATIONAL
HOST_OPT = NO
STATIC_BUILD=NO
F77 := $(RATIONAL) -best-effort $(F77)
endif


ifeq ($(OS_CLASS),Linux)
ifeq ($(STATIC_BUILD),YES)
static_flags = -Wl,-Bstatic -static-libgcc
linux_libs = -lc
endif
endif



usb$(OBJ) : 
	$(F77) ../usb.f -c -o usb$(OBJ)	
brighte$(OBJ) :
	$(F77) ../brighte.f -c -o brighte$(OBJ)
hunt$(OBJ) :
	$(F77) ../hunt.f -c -o hunt$(OBJ)
urgent$(OBJ) :
	$(F77) ../urgent.f -c -o urgent$(OBJ)
us$(OBJ) :
	$(F77) ../us.f -c -o us$(OBJ)
endif

ifdef WIN32
sddsbrightness$(EXE) : $(OBJNAMES)
	ifort /exe:sddsbrightness$(EXE) /threads /libs:dll \
	$(OBJNAMES) /link /incremental:no /opt:ref /release \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbcommon.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\matlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\fftpack.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\SDDS1.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\rpnlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\namelist.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbmth.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdblib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\z.lib
sddsurgent$(EXE) : $(URGENT_OBJNAMES)
	ifort /exe:sddsurgent$(EXE) /threads /libs:dll \
	$(URGENT_OBJNAMES) /link /incremental:no /opt:ref /release \
	c:\\oag\\apps\\lib\\$(T_A)\\oagphy.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbcommon.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\matlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\fftpack.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\SDDS1.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\rpnlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\namelist.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbmth.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdblib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\z.lib
sddsfluxcurve$(EXE) : $(FLUXCURVE_OBJNAMES)
	ifort /exe:sddsfluxcurve$(EXE) /threads /libs:dll \
	$(FLUXCURVE_OBJNAMES) /link /incremental:no /opt:ref /release \
	c:\\oag\\apps\\lib\\$(T_A)\\oagphy.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbcommon.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\matlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\fftpack.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\SDDS1.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\rpnlib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\namelist.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdbmth.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\mdblib.lib \
	c:\\epics\\extensions\\lib\\$(T_A)\\z.lib
else
sddsbrightness$(EXE): $(OBJNAMES)
	$(F77) -o $@ $(OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -lz $(solaris_libs) $(linux_libs)

sddsurgent$(EXE): $(URGENT_OBJNAMES)
	$(F77) -o $@ $(URGENT_OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -loagphy -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -lz $(solaris_libs) $(linux_libs)

sddsfluxcurve$(EXE): $(FLUXCURVE_OBJNAMES)
	$(F77) -o $@ $(FLUXCURVE_OBJNAMES) $(static_flags) -L$(INSTALL_LIB) -L$(EPICS_EXTENSIONS_LIB) -loagphy -lmdbcommon -lmatlib -lfftpack -lSDDS1 -lrpnlib -lnamelist -lmdbmth -lmdblib -lz $(solaris_libs) $(linux_libs)
endif

