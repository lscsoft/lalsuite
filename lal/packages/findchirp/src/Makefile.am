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
  FindChirpTDData.c FindChirpTDTemplate.c \
  FindChirpBCVData.c FindChirpBCVTemplate.c \
  FindChirpBCVFilter.c FindChirpBCVChisq.c \
  FindChirpBCVSpinData.c FindChirpBCVSpinTemplate.c \
  FindChirpBCVSpinFilter.c FindChirpBCVSpinChisq.c \
  FindChirpSimulation.c FindChirpFilterOutputVeto.c \
  FindChirpStoreEvent.c FindChirpClusterEvents.c \
  FindChirpBCVCFilter.c CoherentInspiralFilter.c \
  CoherentInspiralInput.c

libfindchirp_la_SOURCES = $(FINDCHIRPSOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
