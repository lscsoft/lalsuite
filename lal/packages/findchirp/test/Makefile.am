## Process this file with automake to produce Makefile.in
if FRAMEANDMPI
check_PROGRAMS  = FindChirpMain
TESTS           = FindChirp@MPITYPE@
else
check_PROGRAMS =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FindChirpMPICH FindChirpLAM FindChirpMain.c FindChirpMaster.c FindChirpSlave.c FindChirpGlobal.h

FindChirpMain_SOURCES = FindChirpMain.c FindChirpMaster.c FindChirpSlave.c
FindChirpMain_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = PI* Slave*.out Slave*.err machines killscript lamhosts FindChirpSchema
CLEANFILES = FindChirpMain
