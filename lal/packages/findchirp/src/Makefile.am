## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfindchirp.la
FINDCHIRP_SOURCES = FindChirpExch.c FindChirpMemory.c FindChirpFilter.c \
  FindChirpChisq.c FindChirpSPData.c FindChirpSPTemplate.c FindChirpMaster.c \
  FindChirpSlave.c FindChirpTmplt.c FindChirpLinkedList.c
TWOINTERFFINDCHIRP_SOURCES = TwoInterfFindChirpFilter.c

if MPI
#libfindchirp_la_SOURCES = $(FINDCHIRP_SOURCES) $(TWOINTERFFINDCHIRP_SOURCES)
libfindchirp_la_SOURCES =
else
libfindchirp_la_SOURCES =
endif

EXTRA_DIST = $(FINDCHIRP_SOURCES) $(TWOINTERFFINDCHIRP_SOURCES)
EXTRA_TEXSOURCES = $(FINDCHIRP_SOURCES) $(TWOINTERFFINDCHIRP_SOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
