## Process this file with automake to produce Makefile.in
check_PROGRAMS = CLRTest
CLRTest_SOURCES = CLRTest.c
CLRTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = CLRTestScript.sh
EXTRA_DIST = CLRindata.asc CLRTestScript.sh NoTest
MOSTLYCLEANFILES = CLRoutdata.asc
