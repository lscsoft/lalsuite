BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lallibs_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

LDADD = liblalsimulation.la

NOISEPSDFILES = \
	LIGO-P1200087-v18-AdV_BNS_OPTIMIZED.txt \
	LIGO-P1200087-v18-AdV_DESIGN.txt \
	LIGO-P1200087-v18-AdV_EARLY_HIGH.txt \
	LIGO-P1200087-v18-AdV_EARLY_LOW.txt \
	LIGO-P1200087-v18-AdV_LATE_HIGH.txt \
	LIGO-P1200087-v18-AdV_LATE_LOW.txt \
	LIGO-P1200087-v18-AdV_MID_HIGH.txt \
	LIGO-P1200087-v18-AdV_MID_LOW.txt \
	LIGO-P1200087-v18-aLIGO_BNS_OPTIMIZED.txt \
	LIGO-P1200087-v18-aLIGO_DESIGN.txt \
	LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt \
	LIGO-P1200087-v18-aLIGO_EARLY_LOW.txt \
	LIGO-P1200087-v18-aLIGO_LATE_HIGH.txt \
	LIGO-P1200087-v18-aLIGO_LATE_LOW.txt \
	LIGO-P1200087-v18-aLIGO_MID_HIGH.txt \
	LIGO-P1200087-v18-aLIGO_MID_LOW.txt \
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

pkginclude_HEADERS = \
	LALSimBlackHoleRingdown.h \
	LALSimBurst.h \
	LALSimInspiralPrecess.h \
	LALSimInspiralWaveformFlags.h \
	LALSimInspiralTestGRParams.h \
	LALSimIMR.h \
	LALSimInspiral.h \
	LALSimInspiralSphHarmSeries.h \
	LALSimInspiralWaveformCache.h \
	LALSimulation.h \
	LALSimNeutronStar.h \
	LALSimNoise.h \
	LALSimReadData.h \
	LALSimSGWB.h \
	LALSimulationConfig.h \
	LALSimulationVCSInfo.h \
	LALSimulationVCSInfoHeader.h

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
	LALSimNeutronStarEOSTabular.c \
	LALSimIMRPhenomC_internals.c \
	LALSimIMRPhenomP.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimBurstImg.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRSpinAlignedEOB.c \
	LALSimIMRPhenom.c \
	LALSimIMRPhenomC.c \
	LALSimIMRPhenomP.c \
	LALSimIMRPSpinInspiralRD.c \
	LALSimPhenSpinRingDown.c \
	LALSimInspiralWaveformFlags.c \
	LALSimInspiralTestGRParams.c \
	LALSimInspiralPrecess.c \
	LALSimInspiral.c \
	LALSimInspiralSphHarmSeries.c \
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
nodist_liblalsimulation_la_SOURCES = LALSimulationBuildInfo.c

liblalsimulation_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_ns_eos_table_SOURCES = ns-eos-table.c
lalsim_ns_mass_radius_SOURCES = ns-mass-radius.c
lalsim_ns_params_SOURCES = ns-params.c
lalsim_sgwb_SOURCES = sgwb.c
lalsim_unicorn_SOURCES = unicorn.c
lalsim_detector_noise_SOURCES = detector_noise.c

pkgdata_DATA = \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

EXTRA_DIST = \
	config.h.in \
	LALSimulationConfig.h.in \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

DISTCLEANFILES = config.h

vcs_info_sources = LALSimulationVCSInfoHeader.h
build_info_source = LALSimulationBuildInfo.c
$(build_info_source): $(liblalsimulation_la_DEPENDENCIES)
