AM_CPPFLAGS = -I$(top_builddir)/include
LDADD = liblalsimulation.la

bin_PROGRAMS = \
	lalsim-bh-qnmode \
	lalsim-bh-ringdown \
	lalsim-bh-sphwf \
	lalsim-detector-noise \
	lalsim-sgwb \
	lalsim-unicorn

lalsimulationincludedir = $(includedir)/lal
lalsimulationinclude_HEADERS = \
	LALSimBlackHoleRingdown.h \
	LALSimBurst.h \
	LALSimInspiralPrecess.h \
	LALSimInspiralWaveformFlags.h \
	LALSimInspiralTestGRParams.h \
	LALSimIMR.h \
	LALSimInspiral.h \
	LALSimulation.h \
	LALSimNoise.h \
	LALSimSGWB.h \
	LALSimulationConfig.h \
	LALSimulationLibVCSInfo.h \
	LALSimulationVCSInfo.h

noinst_HEADERS = \
	LALSimIMREOBFactorizedWaveform.c \
	LALSimIMREOBHybridRingdown.c \
	LALSimIMREOBNQCCorrection.c \
	LALSimIMREOBNRv2.h \
	LALSimIMREOBNewtonianMultipole.c \
	LALSimInspiralPNCoefficients.c \
	LALSimInspiraldEnergyFlux.c \
	LALSimIMRSpinEOBHamiltonian.c \
	LALSimIMRSpinEOBHcapNumericalDerivative.c \
	LALSimIMRSpinAlignedEOBHcapDerivative.c \
	LALSimIMRSpinEOB.h \
	LALSimIMRSpinEOBInitialConditions.c \
	LALSimIMRSpinEOBAuxFuncs.c \
	LALSimIMRSpinEOBFactorizedWaveform.c \
	LALSimIMRSpinEOBFactorizedFlux.c \
	LALSimIMREOBFactorizedFlux.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimBurstImg.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRSpinAlignedEOB.c \
	LALSimIMRPhenom.c \
	LALSimIMRPSpinInspiralRD.c \
	LALSimInspiralWaveformFlags.c \
	LALSimInspiralTestGRParams.c \
	LALSimInspiralPrecess.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralSpinTaylorT4.c \
	LALSimInspiralTaylorEt.c \
	LALSimInspiralTaylorF2.c \
	LALSimInspiralTaylorF2ReducedSpin.c \
	LALSimInspiralTaylorF2ReducedSpinMetric.c \
	LALSimInspiralTaylorF2ReducedSpinTidal.c \
	LALSimInspiralTaylorT1.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
	LALSimInspiralTaylorLength.c \
	LALSimInspiralWaveformTaper.c \
	LALSimulation.c \
	LALSimulationVCSInfo.c \
	LALSimNoisePSD.c \
	LALSimNoise.c \
	LALSimSGWB.c \
	LALSimSGWBORF.c \
	LALSimUnicorn.c \
	check_series_macros.h \
	unicorn.xpm

liblalsimulation_la_LDFLAGS = -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_sgwb_SOURCES = sgwb.c
lalsim_unicorn_SOURCES = unicorn.c
lalsim_detector_noise_SOURCES = detector_noise.c

MOSTLYCLEANFILES = \
	lalsimulation-user-env.shell \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

sysconf_DATA = \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

BUILT_SOURCES = \
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
	@printf "$(SWIG_USER_ENV)" >>$@

lalsimulation-user-env.csh lalsimulation-user-env.sh: lalsimulation-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	generate_user_env.awk \
	LALSimulationConfig.h.in

if HAVE_GIT_REPO
BUILT_SOURCES += LALSimulationVCSInfo.h
MOSTLYCLEANFILES += LALSimulationVCSInfo.h
LALSimulationVCSInfo.h .PHONY : gen-LALSimulationVCSInfo.h
gen-LALSimulationVCSInfo.h :
	$(GENERATE_VCS_INFO) LALSimulationVCSInfo.h $(srcdir)/LALSimulationVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalsimulationinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
