## Process this file with automake to produce Makefile.in
if MPI
check_PROGRAMS  = CommTestMain
TESTS           = CommTest
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest CommTest CommTestMain.c CommTestMaster.c CommTestSlave.c CommTestGlobal.h

CC=@MPICC@

CommTestMain_SOURCES = CommTestMain.c CommTestMaster.c CommTestSlave.c
CommTestMain_LDADD = $(top_builddir)/lib/liblal.la @MPILIBS@
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = Slave*.out Slave*.err machines killscript
CLEANFILES = CommTestMain
