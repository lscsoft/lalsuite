## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALMallocTest PrintVectorTest PrintFTSeriesTest
LALMallocTest_SOURCES = LALMallocTest.c
LALMallocTest_LDADD = $(top_builddir)/lib/liblal.la
PrintVectorTest_SOURCES = PrintVectorTest.c
PrintVectorTest_LDADD = $(top_builddir)/lib/liblal.la
PrintFTSeriesTest_SOURCES = PrintFTSeriesTest.c
PrintFTSeriesTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALMallocTest PrintVectorTest
MOSTLYCLEANFILES = *.dat *PrintVector.00*
