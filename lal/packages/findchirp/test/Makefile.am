## Process this file with automake to produce Makefile.in
if FRAME
if MPI
noinst_PROGRAMS = FindChirpMain
TESTS           = FindChirp
else
noinst_PROGRAMS =
TESTS           = NoTest
endif
else
noinst_PROGRAMS =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FindChirp FindChirpMain.c FindChirpMaster.c FindChirpSlave.c FindChirpGlobal.h

CC=@MPICC@

FindChirpMain_SOURCES = FindChirpMain.c FindChirpMaster.c FindChirpSlave.c
FindChirpMain_LDADD = $(top_srcdir)/lib/liblal.a @MPILIBS@
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = PI* Slave*.out Slave*.err machines killscript
