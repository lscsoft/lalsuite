## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = RandomTest FindRootTest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_srcdir)/lib/liblal.a
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = PrintVector.*
