## Process this file with automake to produce Makefile.in
#
# $Id$
#
noinst_LTLIBRARIES = libfindchirp.la

FINDCHIRPSOURCES = FindChirpExch.c FindChirpMemory.c FindChirpFilter.c \
  FindChirpChisq.c FindChirpSPData.c FindChirpSPTemplate.c FindChirpMaster.c \
  FindChirpSlave.c FindChirpTmplt.c FindChirpLinkedList.c
TWOINTERFFINDCHIRPSOURCES = TwoInterfFindChirpFilter.c

if MPI
libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES) $(TWOINTERFFINDCHIRPSOURCES)
EXTRA_TEXSOURCES = 
else
libfindchirp_la_SOURCES =
EXTRA_TEXSOURCES = $(FINDCHIRPSOURCES) $(TWOINTERFFINDCHIRPSOURCES)
endif

EXTRA_DIST = $(FINDCHIRPSOURCES) $(TWOINTERFFINDCHIRPSOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
