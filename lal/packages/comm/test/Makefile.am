## Process this file with automake to produce Makefile.in
if MPI
check_PROGRAMS  = CommTestMain
TESTS           = CommTest@MPITYPE@
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest CommTestMPICH CommTestLAM CommTestMain.c CommTestMaster.c CommTestSlave.c CommTestGlobal.h

CommTestMain_SOURCES = CommTestMain.c CommTestMaster.c CommTestSlave.c
CommTestMain_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = PI* Slave*.out Slave*.err machines killscript lamhosts CommTestSchema
CLEANFILES = CommTestMain
