## Process this file with automake to produce Makefile.in
check_PROGRAMS = CLRTest
CLRTest_SOURCES = CLRTest.c
CLRTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = CLRTest
TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
noinst_DATA = CLRindata.asc
EXTRA_DIST = $(DATA)
MOSTLYCLEANFILES = CLRoutdata.asc
