## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALSampleTest LALTemplateTest
LALSampleTest_SOURCES = LALSampleTest.c
LALSampleTest_LDADD = $(top_srcdir)/lib/liblal.la
LALTemplateTest_SOURCES = LALTemplateTest.c
LALTemplateTest_LDADD = $(top_srcdir)/lib/liblal.la
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = LALSampleTest LALTemplateTest
