## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libstd.la

libstd_la_SOURCES = LALStd.c LALVersion.c LALMalloc.c LALError.c LALGSL.c \
  XLALError.c StringToken.c StringConvert.c Grid.c

EXTRA_DIST = LALPrimer.c Grid.m4 LALCreateDestroyGrid.m4
EXTRA_TEXSOURCES = LALPrimer.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = Grid.c
dvi-local: .dvi-dep

Grid.c: Grid.m4 LALCreateDestroyGrid.m4
	cd $(srcdir) && @M4@ Grid.m4 > Grid.c
@DVIMKFILE@
