## Process this file with automake to produce Makefile.in
check_PROGRAMS = OverlapTest DirichletTest CrossCorrTest
OverlapTest_SOURCES = OverlapTest.c
OverlapTest_LDADD = $(top_builddir)/lib/liblal.la
DirichletTest_SOURCES = DirichletTest.c
DirichletTest_LDADD = $(top_builddir)/lib/liblal.la
CrossCorrTest_SOURCES = CrossCorrTest.c
CrossCorrTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = OverlapTest DirichletTest CrossCorrTest
MOSTLYCLEANFILES = PrintVector.*
