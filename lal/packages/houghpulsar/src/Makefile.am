## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libhoughpulsar.la
libhoughpulsar_la_SOURCES = ConstructPLUT.c ParamPLUT.c NDParamPLUT.c PatchGrid.c Stereographic.c Peak2PHMD.c HoughMap.c DriveHough.c Statistics.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
