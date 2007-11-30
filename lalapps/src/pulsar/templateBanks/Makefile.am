## Process this file with automake to produce Makefile.in
# $Id$

AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_skynet getMetric testFlatPulsarMetric matchGrids getMesh FstatMetric testFlatLatticeTiling testLatticeCovering

testFlatLatticeTiling_SOURCES  = testLatticeCovering.c
testLatticeCovering_SOURCES  = testFlatLatticeTiling.c

lalapps_skynet_SOURCES = skynet.c

getMetric_SOURCES = getMetric.c

getMesh_SOURCES = getMesh.c

testFlatPulsarMetric_SOURCES = testFlatPulsarMetric.c

matchGrids_SOURCES = matchGrids.c

FstatMetric_SOURCES = FstatMetric.c

