## Process this file with automake to produce Makefile.in
#
# $Id$
#
noinst_LTLIBRARIES = libfindchirp.la

FINDCHIRPSOURCES = FindChirpMemory.c FindChirpLinkedList.c \
  FindChirpData.c FindChirpTemplate.c \
  FindChirpFilterInit.c FindChirpChisqInit.c \
  FindChirpFilter.c FindChirpChisq.c \
  FindChirpSPData.c FindChirpSPTemplate.c \
  FindChirpBCVData.c FindChirpBCVTemplate.c \
  FindChirpBCVFilter.c FindChirpBCVChisq.c \
  FindChirpBCVSpinData.c FindChirpBCVSpinTemplate.c \
  FindChirpBCVSpinFilter.c FindChirpBCVSpinChisq.c \
  FindChirpSimulation.c FindChirpFilterOutputVeto.c

TWOINTERFFINDCHIRPSOURCES = TwoInterfFindChirpFilter.c \
  TwoInterfFindChirpMemory.c TwoInterfFindChirpSPData.c \
  TwoInterfFindChirpChisq.c

libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES)

EXTRA_DIST = $(TWOINTERFFINDCHIRPSOURCES)
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
