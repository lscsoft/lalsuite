## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## build everything completely static
AM_LDFLAGS = -all-static

isolated_DIR = ../FDS_isolated
isolated_SOURCES = ${isolated_DIR}/rngmed.c ${isolated_DIR}/rngmed.h \
	${isolated_DIR}/clusters.c ${isolated_DIR}/clusters.h \
	${isolated_DIR}/DopplerScan.c ${isolated_DIR}/DopplerScan.h \
	${isolated_DIR}/NeighboringGridPoints.c
bin_PROGRAMS = lalapps_ComputeFStatisticBinary lalapps_GenerateBinaryMesh

lalapps_ComputeFStatisticBinary_SOURCES = ComputeFStatisticBinary_v2.c ComputeFStatisticBinary_v2.h	\
	${isolated_SOURCES}

lalapps_GenerateBinaryMesh_SOURCES = GenerateBinaryMesh_v1.c GenerateBinaryMesh_v1.h 
