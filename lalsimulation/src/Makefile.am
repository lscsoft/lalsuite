INCLUDES = -I$(top_builddir)/include
LDADD = liblalsimulation.la

bin_PROGRAMS = \
	lalsim-bh-qnmode \
	lalsim-bh-ringdown \
	lalsim-bh-sphwf \
	lalsim-ligo-noise \
	lalsim-sgwb

lalsimulationincludedir = $(includedir)/lal
lalsimulationinclude_HEADERS = \
	LALSimBlackHoleRingdown.h \
	LALSimBurst.h \
	LALSimIMR.h \
	LALSimInspiral.h \
	LALSimulation.h \
	LALSimInspiralSpinTaylorT4.h \
	LALSimNoise.h \
	LALSimSGWB.h \
	LALSimulationConfig.h \
	LALSimulationVCSInfo.h

noinst_HEADERS = \
	LALSimIMREOBFactorizedFlux.c \
	LALSimIMREOBFactorizedWaveform.c \
	LALSimIMREOBHybridRingdown.c \
	LALSimIMREOBNQCCorrection.c \
	LALSimIMREOBNRv2.h \
	LALSimIMREOBNewtonianMultipole.c \
	LALSimInspiralPNCoefficients.c \
	LALSimInspiraldEnergyFlux.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRPhenom.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralTaylorEt.c \
	LALSimInspiralTaylorF2ReducedSpin.c \
	LALSimInspiralTaylorT1.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
	LALSimInspiralTaylorLength.c \
	LALSimulation.c \
	LALSimInspiralSpinTaylorT4.c \
	LALSimulationVCSInfo.c \
	LALSimInspiralPSpinInspiralRD.c \
	LALSimNoisePSD.c \
	LALSimNoise.c \
	LALSimSGWB.c \
	LALSimSGWBORF.c \
	check_series_macros.h

liblalsimulation_la_LDFLAGS = -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_ligo_noise_SOURCES = ligo_noise.c
lalsim_sgwb_SOURCES = sgwb.c

CLEANFILES = \
	lalsimulation-user-env.shell \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

sysconf_DATA = \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALSimulationVCSInfo.h \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

lalsimulation-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALSIMULATION_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "$(SWIG_USER_ENV)" | $(SED) 's|~E~O~L~|\n|g' >>$@

lalsimulation-user-env.csh lalsimulation-user-env.sh: lalsimulation-user-env.shell $(srcdir)/generate_user_env.awk
	@echo "creating $@"
	@$(AWK) -f $(srcdir)/generate_user_env.awk -v package=$(PACKAGE_NAME) -v output=$@ $<

EXTRA_DIST = \
	generate_user_env.awk \
	generate_vcs_info.py \
	LALSimulationConfig.h.in \
	LALSimulationVCSInfo.h.in

LALSimulationVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALSimulationConfig.h
	@echo "creating LALSimulationVCSInfo.h"
	@$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALSimulation $(srcdir) $$PWD

@LAL_HEADER_LINK@

swig_headers = $(lalsimulationinclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
