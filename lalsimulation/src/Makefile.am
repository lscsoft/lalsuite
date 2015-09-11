BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
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
	lalsim-burst \
	lalsim-detector-noise \
	lalsim-detector-strain \
	lalsim-inject \
	lalsim-inspiral \
	lalsim-ns-eos-table \
	lalsim-ns-mass-radius \
	lalsim-ns-params \
	lalsim-sgwb \
	lalsim-unicorn

pkginclude_HEADERS = \
	LALSimBlackHoleRingdown.h \
	LALSimBurst.h \
	LALSimIMR.h \
	LALSimInspiral.h \
	LALSimInspiralPrecess.h \
	LALSimInspiralTestGRParams.h \
	LALSimInspiralWaveformCache.h \
	LALSimInspiralWaveformFlags.h \
	LALSimNeutronStar.h \
	LALSimNoise.h \
	LALSimReadData.h \
	LALSimSGWB.h \
	LALSimSphHarmMode.h \
	LALSimSphHarmSeries.h \
	LALSimUtils.h \
	LALSimulation.h \
	LALSimulationConfig.h \
	LALSimulationVCSInfo.h \
	LALSimulationVCSInfoHeader.h

noinst_HEADERS = \
	LALSimFindAttachTime.h \
	LALSimIMREOBFactorizedFlux.c \
	LALSimIMREOBFactorizedWaveform.c \
	LALSimIMREOBHybridRingdown.c \
	LALSimIMREOBNQCCorrection.c \
	LALSimIMREOBNRv2.h \
	LALSimIMREOBNewtonianMultipole.c \
	LALSimIMRPhenomC_internals.c \
	LALSimIMRPhenomD.h \
	LALSimIMRPhenomD_internals.c \
	LALSimIMRPhenomP.c \
	LALSimIMRSEOBNRROMUtilities.c \
	LALSimIMRSpinAlignedEOBHcapDerivative.c \
	LALSimIMRSpinEOB.h \
	LALSimIMRSpinEOBAuxFuncs.c \
	LALSimIMRSpinEOBFactorizedFlux.c \
	LALSimIMRSpinEOBFactorizedWaveform.c \
	LALSimIMRSpinEOBHamiltonian.c \
	LALSimIMRSpinEOBHcapNumericalDerivative.c \
	LALSimIMRSpinEOBInitialConditions.c \
	LALSimInspiralPNCoefficients.c \
	LALSimInspiraldEnergyFlux.c \
	LALSimNeutronStarEOSPiecewisePolytrope.c \
	LALSimNeutronStarEOSTabular.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimBurstImg.c \
	LALSimFindAttachTime.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRPSpinInspiralRD.c \
	LALSimIMRPhenom.c \
	LALSimIMRPhenomC.c \
	LALSimIMRPhenomD.c \
	LALSimIMRPhenomP.c \
	LALSimIMRSEOBNRv1ROMDoubleSpin.c \
	LALSimIMRSEOBNRv1ROMEffectiveSpin.c \
	LALSimIMRSEOBNRv2ChirpTime.c \
	LALSimIMRSEOBNRv2ROMDoubleSpin.c \
	LALSimIMRSEOBNRv2ROMEffectiveSpin.c \
	LALSimIMRSpinAlignedEOB.c \
	LALSimInspiral.c \
	LALSimInspiralEccentricTD.c \
	LALSimInspiralEccentricityFD.c \
	LALSimInspiralHGimri.c \
	LALSimInspiralOptimizedCoefficientsEccentricityFD.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralPrecess.c \
	LALSimInspiralSpinDominatedWaveform.c \
	LALSimInspiralSpinTaylor.c \
	LALSimInspiralSpinTaylorF2.c \
	LALSimInspiralSpinTaylorT5.c \
	LALSimInspiralTaylorEt.c \
	LALSimInspiralTaylorF2.c \
	LALSimInspiralTaylorF2ReducedSpin.c \
	LALSimInspiralTaylorF2ReducedSpinMetric.c \
	LALSimInspiralTaylorF2ReducedSpinTidal.c \
	LALSimInspiralTaylorLength.c \
	LALSimInspiralTaylorT1.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
	LALSimInspiralTestGRParams.c \
	LALSimInspiralWaveformCache.c \
	LALSimInspiralWaveformFlags.c \
	LALSimInspiralWaveformTaper.c \
	LALSimNeutronStarEOS.c \
	LALSimNeutronStarFamily.c \
	LALSimNeutronStarTOV.c \
	LALSimNoise.c \
	LALSimNoisePSD.c \
	LALSimPhenSpinRingDown.c \
	LALSimReadData.c \
	LALSimSGWB.c \
	LALSimSGWBORF.c \
	LALSimSphHarmMode.c \
	LALSimSphHarmSeries.c \
	LALSimUnicorn.c \
	LALSimUtils.c \
	LALSimulation.c \
	LALSimulationVCSInfo.c \
	check_series_macros.h \
	check_waveform_macros.h \
	unicorn.xpm

nodist_liblalsimulation_la_SOURCES = LALSimulationBuildInfo.c

liblalsimulation_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_burst_SOURCES = burst.c
lalsim_detector_noise_SOURCES = detector_noise.c
lalsim_detector_strain_SOURCES = detector_strain.c
lalsim_inject_SOURCES = inject.c
lalsim_inspiral_SOURCES = inspiral.c
lalsim_ns_eos_table_SOURCES = ns-eos-table.c
lalsim_ns_mass_radius_SOURCES = ns-mass-radius.c
lalsim_ns_params_SOURCES = ns-params.c
lalsim_sgwb_SOURCES = sgwb.c
lalsim_unicorn_SOURCES = unicorn.c

pkgdata_DATA = \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

EXTRA_DIST = \
	LALSimulationConfig.h.in \
	config.h.in \
	lalsimburst_btlwnbexamples.svg \
	lalsimburst_sinegaussianexamples.svg \
	lalsimburst_stringcuspexamples.svg \
	lalsiminspiral_orbitelements.svg \
	lalsimulation_inject.svg \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

DISTCLEANFILES = config.h

vcs_info_sources = LALSimulationVCSInfoHeader.h
build_info_source = LALSimulationBuildInfo.c
$(build_info_source): $(liblalsimulation_la_DEPENDENCIES)
