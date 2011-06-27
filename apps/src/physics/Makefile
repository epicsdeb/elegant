#
# $Log: Makefile,v $
# Revision 1.1  2006/01/09 22:05:21  shang
# first version
#
#
#
#  Lowest Level Directory Makefile
#
#

TOP=../..
#If epics/extensions/configure directory exists, build with it.
#Otherwise use epics/extensions/config.
ifeq (0, $(words $(notdir $(wildcard $(TOP)/configure))))
include $(TOP)/config/CONFIG_APPS
include $(TOP)/config/RULES_ARCHS
else
include $(TOP)/configure/CONFIG
include $(TOP)/src/physics/Makefile.OAG
include $(TOP)/configure/RULES
endif


