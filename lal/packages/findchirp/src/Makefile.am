## Process this file with automake to produce Makefile.in
#
# $Id$
#
noinst_LTLIBRARIES = libfindchirp.la

FINDCHIRPSOURCES = FindChirpMemory.c FindChirpFilter.c \
  FindChirpChisq.c FindChirpSPData.c FindChirpSPTemplate.c \
  FindChirpTmplt.c FindChirpLinkedList.c FindChirpSimulation.c
TWOINTERFFINDCHIRPSOURCES = TwoInterfFindChirpFilter.c \
  TwoInterfFindChirpMemory.c TwoInterfFindChirpSPData.c

libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES) $(TWOINTERFFINDCHIRPSOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
