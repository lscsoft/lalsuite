BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	mathematica_codes \
	$(END_OF_LIST)

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
	LALSimBlackHoleRingdownPrec.c \
	LALSimBlackHoleRingdownPrec.h \
	LALSimFindAttachTime.h \
	LALSimIMREOBFactorizedFlux.c \
	LALSimIMREOBFactorizedWaveform.c \
	LALSimIMREOBHybridRingdown.c \
	LALSimIMREOBHybridRingdownPrec.c \
	LALSimIMREOBNQCCorrection.c \
	LALSimIMREOBNQCTables.c \
	LALSimIMREOBNRv2.h \
	LALSimIMREOBNRv2HMROMUtilities.c \
	LALSimIMREOBNewtonianMultipole.c \
	LALSimIMRPhenomC_internals.c \
	LALSimIMRPhenomC_internals.h \
	LALSimIMRPhenomD.h \
	LALSimIMRPhenomD_internals.c \
	LALSimIMRPhenomD_internals.h \
	LALSimIMRPhenomP.c \
	LALSimIMRPhenomP.h \
	LALSimIMRSEOBNRROMUtilities.c \
	LALSimIMRSpinAlignedEOBGSLOptimizedInterpolation.c \
	LALSimIMRSpinAlignedEOBHcapDerivative.c \
	LALSimIMRSpinAlignedEOBHcapDerivativeOptimized.c \
	LALSimIMRSpinEOB.h \
	LALSimIMRSpinEOBAuxFuncs.c \
	LALSimIMRSpinEOBAuxFuncsPrec.c \
	LALSimIMRSpinEOBComputeAmpPhasefromEOMSoln.c \
	LALSimIMRSpinEOBFactorizedFlux.c \
	LALSimIMRSpinEOBFactorizedFluxOptimized.c \
	LALSimIMRSpinEOBFactorizedFluxPrec.c \
	LALSimIMRSpinEOBFactorizedWaveform.c \
	LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c \
	LALSimIMRSpinEOBFactorizedWaveformPrec.c \
	LALSimIMRSpinEOBHamiltonian.c \
	LALSimIMRSpinEOBHamiltonianOptimized.c \
	LALSimIMRSpinEOBHamiltonianPrec.c \
	LALSimIMRSpinEOBHcapExactDerivative.c \
	LALSimIMRSpinEOBHcapNumericalDerivative.c \
	LALSimIMRSpinEOBHcapNumericalDerivativePrec.c \
	LALSimIMRSpinEOBInitialConditions.c \
	LALSimIMRSpinEOBInitialConditionsPrec.c \
	LALSimInspiralPNCoefficients.c \
	LALSimInspiraldEnergyFlux.c \
	LALSimNeutronStarEOSPiecewisePolytrope.c \
	LALSimNeutronStarEOSTabular.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimBurstImg.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRSpinAlignedEOB.c \
	LALSimIMRSpinPrecEOB.c \
	LALSimIMRSEOBNRv1ROMEffectiveSpin.c \
	LALSimIMRSEOBNRv1ROMDoubleSpin.c \
	LALSimIMRSEOBNRv2ROMEffectiveSpin.c \
	LALSimIMRSEOBNRv2ROMDoubleSpin.c \
	LALSimIMRSEOBNRv2ROMDoubleSpinHI.c \
	LALSimIMRSEOBNRv2ChirpTime.c \
	LALSimIMREOBNRv2HMROM.c \
	LALSimIMRPhenom.c \
	LALSimIMRPhenomC.c \
	LALSimIMRPhenomD.c \
	LALSimIMRPhenomP.c \
	LALSimIMRPSpinInspiralRD.c \
	LALSimPhenSpinRingDown.c \
	LALSimInspiralWaveformFlags.c \
	LALSimInspiralTestGRParams.c \
	LALSimInspiralPrecess.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralSpinTaylor.c \
	LALSimInspiralSpinTaylorF2.c \
	LALSimInspiralEccentricityFD.c \
	LALSimInspiralOptimizedCoefficientsEccentricityFD.c \
	LALSimInspiralTaylorEt.c \
	LALSimInspiralTaylorF2.c \
	LALSimInspiralTaylorF2ReducedSpin.c \
	LALSimInspiralTaylorF2ReducedSpinMetric.c \
	LALSimInspiralTaylorF2ReducedSpinTidal.c \
	LALSimInspiralTaylorT1.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
	LALSimInspiralEccentricTD.c \
	LALSimInspiralSpinTaylorT5.c \
	LALSimInspiralSpinDominatedWaveform.c \
	LALSimInspiralTaylorLength.c \
	LALSimInspiralWaveformCache.c \
	LALSimInspiralWaveformTaper.c \
        LALSimIMRNRWaveforms.c \
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
	LALSimSphHarmMode.c \
	LALSimSphHarmSeries.c \
	LALSimUnicorn.c \
	LALSimUtils.c \
	check_series_macros.h \
	check_waveform_macros.h \
	unicorn.xpm \
	LALSimInspiralHGimri.c \
	LALSimFindAttachTime.c

nodist_liblalsimulation_la_SOURCES = LALSimulationBuildInfo.c

liblalsimulation_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_burst_SOURCES = burst.c
lalsim_ns_eos_table_SOURCES = ns-eos-table.c
lalsim_ns_mass_radius_SOURCES = ns-mass-radius.c
lalsim_ns_params_SOURCES = ns-params.c
lalsim_sgwb_SOURCES = sgwb.c
lalsim_unicorn_SOURCES = unicorn.c
lalsim_detector_noise_SOURCES = detector_noise.c
lalsim_detector_strain_SOURCES = detector_strain.c
lalsim_inspiral_SOURCES = inspiral.c
lalsim_inject_SOURCES = inject.c

pkgdata_DATA = \
	$(NOISEPSDFILES) \
	$(EOSDATAFILES)

EXTRA_DIST = \
	config.h.in \
	LALSimulationConfig.h.in \
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
