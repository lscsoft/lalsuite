## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest
LALfctTest_SOURCES = LALfctTest.c timingfunctions.c timingfunctions.h
LALfctTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
LALfctSample_SOURCES = LALfctSample.c
LALfctSample_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
LALFCTInterfaceTest_SOURCES = LALFCTInterfaceTest.c
LALFCTInterfaceTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALfctTest LALFCTInterfaceTest
MOSTLYCLEANFILES = .dvi-dep output.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
