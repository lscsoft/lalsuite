## Process this file with automake to produce Makefile.in
#check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest
check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest LALSNRTest
LALfctTest_SOURCES = LALfctTest.c timingfunctions.c timingfunctions.h
LALfctTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALfctSample_SOURCES = LALfctSample.c
LALfctSample_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALFCTInterfaceTest_SOURCES = LALFCTInterfaceTest.c
LALFCTInterfaceTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALSNRTest_SOURCES = LALSNRTest.c
LALSNRTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALfctTest LALFCTInterfaceTest
MOSTLYCLEANFILES = .dvi-dep output.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
