## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libinspiral.la
libinspiral_la_SOURCES = \
LALBCVWaveform.c\
LALEOBWaveform.c\
LALEOBWaveformTemplates.c\
LALEtaTau02.c\
LALEtaTau04.c\
LALInspiralChooseModel.c\
LALInspiralFrequency3.c\
LALInspiralParameterCalc.c\
LALInspiralPhasing1.c\
LALInspiralPhasing2.c\
LALInspiralPhasing3.c\
LALInspiralPhiofVIntegrand.c\
LALInspiralSetup.c\
LALInspiralSpinningBHBinary.c\
LALInspiralStationaryPhaseApprox1.c\
LALInspiralStationaryPhaseApprox2.c\
LALInspiralTiming2.c\
LALInspiralTofV.c\
LALInspiralTofVIntegrand.c\
LALInspiralVelocity.c\
LALInspiralWaveLength.c\
LALInspiralWave.c\
LALInspiralWaveTemplates.c\
LALInspiralWave1.c\
LALInspiralWave1Templates.c\
LALInspiralWave2.c\
LALInspiralWave2Templates.c\
LALInspiralWave3.c\
LALInspiralWave3Templates.c\
LALInspiralDerivatives.c\
LALPrintTimeseries.c\
LALRungeKutta4.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
