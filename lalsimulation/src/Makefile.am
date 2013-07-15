AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = liblalsimulation.la

NOISEPSDFILES = \
	LIGO-T0900288-v3-BHBH_20deg.txt \
	LIGO-T0900288-v3-High_Freq.txt \
	LIGO-T0900288-v3-NO_SRM.txt \
	LIGO-T0900288-v3-NSNS_Opt.txt \
	LIGO-T0900288-v3-ZERO_DET_high_P.txt \
	LIGO-T0900288-v3-ZERO_DET_low_P.txt

EOSDATAFILES = \
	LALSimNeutronStarEOS_AP4.dat \
	LALSimNeutronStarEOS_FPS.dat \
	LALSimNeutronStarEOS_SLY4.dat

bin_PROGRAMS = \
	lalsim-bh-qnmode \
	lalsim-bh-ringdown \
	lalsim-bh-sphwf \
	lalsim-detector-noise \
	lalsim-ns-eos-table \
	lalsim-ns-mass-radius \
	lalsim-ns-params \
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
	LALSimInspiralWaveformCache.h \
	LALSimulation.h \
	LALSimNeutronStar.h \
	LALSimNoise.h \
	LALSimReadData.h \
	LALSimSGWB.h \
	LALSimulationConfig.h \
	LALSimulationLibVCSInfo.h \
	LALSimulationVCSInfo.h \
	lalsimulationswig_test.i

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
	LALSimIMREOBFactorizedFlux.c \
	LALSimNeutronStarEOSPiecewisePolytrope.c \
	LALSimNeutronStarEOSTabular.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimBurstImg.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRSpinAlignedEOB.c \
	LALSimIMRPhenom.c \
	LALSimIMRPhenomC.c \
	LALSimIMRPSpinInspiralRD.c \
	LALSimPhenSpinRingDown.c \
	LALSimInspiralWaveformFlags.c \
	LALSimInspiralTestGRParams.c \
	LALSimInspiralPrecess.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralSpinTaylor.c \
	LALSimInspiralSpinTaylorF2.c \
	LALSimInspiralTaylorEt.c \
	LALSimInspiralTaylorF2.c \
	LALSimInspiralTaylorF2ReducedSpin.c \
	LALSimInspiralTaylorF2ReducedSpinMetric.c \
	LALSimInspiralTaylorF2ReducedSpinTidal.c \
	LALSimInspiralTaylorT1.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
        LALSimInspiralSpinTaylorT5.c \
	LALSimInspiralTaylorLength.c \
	LALSimInspiralWaveformCache.c \
	LALSimInspiralWaveformTaper.c \
	LALSimulation.c \
	LALSimulationVCSInfo.c \
	LALSimNeutronStarEOS.c \
	LALSimNeutronStarFamily.c \
	LALSimNeutronStarTOV.c \
	LALSimNoisePSD.c \
	LALSimNoise.c \
	LALSimReadData.c \
	LALSimSGWB.c \
	LALSimSGWBORF.c \
	LALSimUnicorn.c \
	check_series_macros.h \
	unicorn.xpm

liblalsimulation_la_LDFLAGS = -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_ns_eos_table_SOURCES = ns-eos-table.c
lalsim_ns_mass_radius_SOURCES = ns-mass-radius.c
lalsim_ns_params_SOURCES = ns-params.c
lalsim_sgwb_SOURCES = sgwb.c
lalsim_unicorn_SOURCES = unicorn.c
lalsim_detector_noise_SOURCES = detector_noise.c

MOSTLYCLEANFILES = \
	lalsimulation-user-env.shell \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

pkgdata_DATA = \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

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
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalsimulation-user-env.csh lalsimulation-user-env.sh: lalsimulation-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	config.h.in \
	generate_user_env.awk \
	LALSimulationConfig.h.in \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

if HAVE_GIT_REPO
BUILT_SOURCES += LALSimulationVCSInfo.h
MOSTLYCLEANFILES += LALSimulationVCSInfo.h
LALSimulationVCSInfo.h .PHONY : gen-LALSimulationVCSInfo.h
gen-LALSimulationVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALSimulationVCSInfo.h $(srcdir)/LALSimulationVCSInfo.h.git
endif

DISTCLEANFILES = config.h

@LAL_HEADER_LINK@

swig_headers = $(lalsimulationinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
