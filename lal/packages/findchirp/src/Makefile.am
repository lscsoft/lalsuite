## Process this file with automake to produce Makefile.in
#
# $Id$
#
noinst_LTLIBRARIES = libfindchirp.la

FINDCHIRPSOURCES = FindChirpMemory.c FindChirpFilter.c \
  FindChirpChisq.c FindChirpSPData.c FindChirpSPTemplate.c \
  FindChirpTmplt.c FindChirpLinkedList.c 
FINDCHIRPMPISOURCES = FindChirpExch.c FindChirpMaster.c FindChirpSlave.c
TWOINTERFFINDCHIRPSOURCES = TwoInterfFindChirpFilter.c

if MPI
libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES) $(FINDCHIRPMPISOURCES) \
                          $(TWOINTERFFINDCHIRPSOURCES)

EXTRA_TEXSOURCES = FindChirpView.c
else
libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES) $(TWOINTERFFINDCHIRPSOURCES)
EXTRA_TEXSOURCES = FindChirpView.c $(FINDCHIRPMPISOURCES)
endif

EXTRA_DIST = $(FINDCHIRPMPISOURCES) FindChirpView.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
