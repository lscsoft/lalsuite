## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a
liblal_a_SOURCES = \
ChooseModel.c \
InspiralDerivatives.c \
InspiralParameterCalc.c \
InspiralPhase.c \
InspiralSetup.c \
InspiralVelocity.c \
InspiralWave.c \
PhiofVIntegrand.c \
RungeKutta4.c \
TappRpnTdomFreq.c \
TappRpnTdomFreqPhase.c \
TappRpnTdomFreqToff.c \
TappRpnTdomTime.c \
TappRpnTdomTimeFrequency.c \
TappRpnTdomTimePhase.c \
TimeDomain2.c \
TofV.c \
TofVIntegrand.c

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

