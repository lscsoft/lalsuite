## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALSampleTest LALTemplateTest
LALSampleTest_SOURCES = LALSampleTest.c
LALSampleTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
LALTemplateTest_SOURCES = LALTemplateTest.c
LALTemplateTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALSampleTest LALTemplateTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
