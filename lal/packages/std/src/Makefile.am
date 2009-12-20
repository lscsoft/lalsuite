INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_LTLIBRARIES = libstd.la
libstd_la_SOURCES = \
	Grid.c \
	LALComplex.c \
	LALError.c \
	LALGSL.c \
	LALMalloc.c \
	LALStd.c \
	LALString.c \
	LALVersion.c \
	StringConvert.c \
	StringToken.c \
	XLALError.c

EXTRA_DIST = \
	Grid.m4 \
	LALCreateDestroyGrid.m4 \
	LALPrimer.c

EXTRA_TEXSOURCES = LALPrimer.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = Grid.c

dvi-local: .dvi-dep

Grid.c: Grid.m4 LALCreateDestroyGrid.m4
	cd $(srcdir) && @M4@ Grid.m4 > Grid.c

@DVIMKFILE@
