## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALSampleTest LALTemplateTest
LALSampleTest_SOURCES = LALSampleTest.c
LALSampleTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALTemplateTest_SOURCES = LALTemplateTest.c
LALTemplateTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALSampleTest LALTemplateTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
