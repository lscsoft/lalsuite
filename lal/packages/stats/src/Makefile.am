## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libstats.la
libstats_la_SOURCES = LALMoment.c LALCorrelation.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MAINTAINERCLEANFILES = LALMoment.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
EXTRA_DIST = LALMoment.m4 LALMomentBase.m4

LALMoment.c: LALMoment.m4 LALMomentBase.m4
	cd $(srcdir) && @M4@ LALMoment.m4 > LALMoment.c
@DVIMKFILE@
