## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libtdfilter.la
libtdfilter_la_SOURCES = \
        LPC.c BilinearTransform.c      CreateZPGFilter.c        IIRFilter.c \
        ButterworthTimeSeries.c  DestroyIIRFilter.c       IIRFilterVector.c \
        CreateIIRFilter.c        DestroyZPGFilter.c       IIRFilterVectorR.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
