AM_CFLAGS = @AM_CFLAGS@ $(LALSUPPORT_CFLAGS)
AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

check_PROGRAMS = \
	BCVTemplates \
	BasicInjectTest \
	CoarseTest \
	CoarseTest2 \
	FindChirpBankVetoTest \
	GenerateInspiralWaveform \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	GetOrientationEllipse \
	InjectionInterfaceTest \
	InspiralBCVSpinBankTest \
	InspiralSpinBankTest \
	LALInspiralSpinningBHBinariesTest \
	LALInspiralTest \
	LALInspiralTaylorT2Test \
	LALInspiralTaylorT3Test \
	LALInspiralTaylorT4Test \
	LALSTPNWaveformTest \
	MetricTest \
	MetricTestBCV \
	MetricTestPTF \
	PNTemplates \
	RandomInspiralSignalTest
TESTS = $(check_PROGRAMS)
#FindChirpTDTest_SOURCES = FindChirpTDTest.c FindChirpTDTest.h

# non-building test
# BCVSpinTemplates
# ChirpSpace
# FilterTest
# FindChirpACTDTest
# FindChirpBCVTest
# FindChirpSPTest
# FindChirpTmpltTest
# InspiralSpinBankwNDTemplateBankTest
# NRWaveInjectTest
# RingTest
# TestRingdown

# failing test
# LALHybridTest
# NRWaveIOTest
# SpaceCovering
# FindChirpTDTest

# requires lalstochastic
# CoherentInspiralFilterTest
# CoherentInspiralInputTest

# requires lalpulsar
# BCVTemplatesFlatMesh

EXTRA_DIST = \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml

CLEANFILES = \
	*.out \
	*.dat
