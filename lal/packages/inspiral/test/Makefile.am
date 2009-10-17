## Process this file with automake to produce Makefile.in
TESTS = NoTest.sh
EXTRA_DIST = \
	NoTest.sh \
	ChirpSpace.c \
	GenerateInspiralWaveform.c \
	LALHybridTest.c \
	LALInspiralSpinningBHBinariesTest.c \
	LALInspiralTest.c \
	LALSTPNWaveformTest.c \
	TestRingdown.c
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep
dvi-local: .dvi-dep
@DVIMKFILE@
