## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfindchirp.la

if MPI
libfindchirp_la_SOURCES = FindChirpExch.c \
                          FindChirpMemory.c \
                          FindChirpFilter.c \
                          FindChirpChisq.c \
                          FindChirpSPData.c \
                          FindChirpSPTemplate.c \
                          FindChirpMaster.c \
                          FindChirpSlave.c \
                          FindChirpTmplt.c
else
libfindchirp_la_SOURCES = 
endif

EXTRA_DIST = FindChirpExch.c \
             FindChirpMemory.c \
             FindChirpFilter.c \
             FindChirpChisq.c \
             FindChirpSPData.c \
             FindChirpSPTemplate.c \
             FindChirpMaster.c \
             FindChirpSlave.c \
             FindChirpTmplt.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
