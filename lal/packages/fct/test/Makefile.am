## Process this file with automake to produce Makefile.in
#check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest
if INTELFFT
check_PROGRAMS =
TESTS =
else
check_PROGRAMS = LALfctTest LALfctSample LALFCTInterfaceTest LALSNRTest
LALfctTest_SOURCES = LALfctTest.c timingfunctions.c timingfunctions.h
LALfctTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALfctSample_SOURCES = LALfctSample.c
LALfctSample_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALFCTInterfaceTest_SOURCES = LALFCTInterfaceTest.c
LALFCTInterfaceTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALSNRTest_SOURCES = LALSNRTest.c
LALSNRTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TESTS = LALfctTest LALFCTInterfaceTest
endif
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep output.dat
DISTCLEANFILES = *.tex
EXTRA_DIST = LALfctTest.c timingfunctions.c timingfunctions.h \
             LALfctSample.c LALFCTInterfaceTest.c LALSNRTest.c
dvi-local: .dvi-dep
@DVIMKFILE@
