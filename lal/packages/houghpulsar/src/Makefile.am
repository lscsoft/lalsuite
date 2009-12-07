INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_LTLIBRARIES = libhoughpulsar.la
libhoughpulsar_la_SOURCES = \
	ConstructPLUT.c \
	DriveHough.c \
	HoughMap.c \
	NDParamPLUT.c \
	ParamPLUT.c \
	PatchGrid.c \
	Peak2PHMD.c \
	Statistics.c \
	Stereographic.c \
	Velocity.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
