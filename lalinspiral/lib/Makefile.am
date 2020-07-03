BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	BBHPhenomCoeffs.h \
	CoincInspiralEllipsoid.h \
	FindChirp.h \
	FindChirpPTF.h \
	GenerateInspRing.h \
	GenerateInspiral.h \
	GeneratePPNInspiral.h \
	InspiralInjectionParams.h \
	LALEOBNRv2Waveform.h \
	LALInspiral.h \
	LALInspiralBank.h \
	LALInspiralConfig.h \
	LALInspiralVCSInfo.h \
	LALInspiralSBankOverlap.h \
	LALInspiralVCSInfoHeader.h \
	LALNoiseModelsInspiral.h \
	LALTrigScanCluster.h \
	LALSTPNWaveformErrors.h \
	LALSTPNWaveform2.h \
	LALSTPNWaveformFrameless.h \
	LALSQTPNIntegrator.h \
	LALSQTPNWaveform.h \
	LALSQTPNWaveformInterface.h \
	LIGOLwXMLInspiralRead.h \
	LIGOLwXMLRingdownRead.h \
	LIGOMetadataInspiralUtils.h \
	LIGOMetadataRingdownUtils.h \
	NRWaveIO.h \
	NRWaveInject.h \
	RingUtils.h \
	TemplateBankGeneration.h \
	TrigScanEThincaCommon.h

lib_LTLIBRARIES = liblalinspiral.la

liblalinspiral_la_SOURCES = \
	BandPassInspiralTemplate.c \
	FindChirpPTFWaveform.c \
	FindChirpPhenomWaveform.c \
	FindChirpSimulation.c \
	GenerateInspRing.c \
	GenerateInspiral.c \
	GeneratePPNAmpCorInspiral.c \
	GeneratePPNInspiral.c \
	GetErrorMatrixFromSnglInspiral.c \
	GetInspiralParams.c \
	InspiralBankGeneration.c \
	InspiralInjectionParams.c \
	InspiralSpinBank.c \
	LALBCVWaveform.c \
	LALEOBGetFactorizedWaveform.c \
	LALEOBWaveform.c \
	LALEOBNonQCCorrection.c \
	LALEOBPPWaveform.c \
	LALEtaTau02.c \
	LALEtaTau04.c \
	LALFactorizedFlux.c \
	LALHexagonVertices.c \
	LALInsidePolygon.c \
	LALInspiralAmplitude.c \
	LALInspiralAmplitudeCorrectedWave.c \
	LALInspiralBCVBank.c \
	LALInspiralBankList.c \
	LALInspiralBankUtils.c \
	LALInspiralChooseModel.c \
	LALInspiralComputeMetric.c \
	LALInspiralComputePTFMetric.c \
	LALInspiralComputeParams.c \
	LALInspiralCreateCoarseBank.c \
	LALInspiralCreateFineBank.c \
	LALInspiralDerivatives.c \
	LALInspiralEccentricity.c \
	LALInspiralFrequency3.c \
	LALInspiralGetApproximantString.c \
	LALInspiralHexagonalBank.c \
	LALInspiralHybridHexagonalBank.c \
	LALInspiralIIR.c \
	LALInspiralInit.c \
	LALInspiralLongestTemplateInBank.c \
	LALInspiralMoments.c \
	LALInspiralNextTemplate.c \
	LALInspiralNinjaInject.c \
	LALInspiralParameterCalc.c \
	LALInspiralPhasing1.c \
	LALInspiralPhasing2.c \
	LALInspiralPhasing3.c \
	LALInspiralPhiofVIntegrand.c \
	LALInspiralPolarizations.c \
	LALInspiralRingdownWave.c \
	LALInspiralSBankOverlap.c \
	LALInspiralSetParams.c \
	LALInspiralSetSearchLimits.c \
	LALInspiralSetup.c \
	LALInspiralSpinningBHBinary.c \
	LALInspiralStationaryPhaseApprox1.c \
	LALInspiralStationaryPhaseApprox2.c \
	LALInspiralTaylorEtWaveform.c \
	LALInspiralTaylorNWaveform.c \
	LALInspiralTaylorT4Waveform.c \
	LALInspiralTiming2.c \
	LALInspiralTofV.c \
	LALInspiralTofVIntegrand.c \
	LALInspiralUpdateParams.c \
	LALInspiralValidParams.c \
	LALInspiralValidTemplate.c \
	LALInspiralVelocity.c \
	LALInspiralWave.c \
	LALInspiralWave1.c \
	LALInspiralWave2.c \
	LALInspiralWave3.c \
	LALInspiralWaveCorrelate.c \
	LALInspiralWaveLength.c \
	LALInspiralWaveNormalise.c \
	LALInspiralWaveNormaliseLSO.c \
	LALInspiralWaveOverlap.c \
	LALNewtonianMultipole.c \
	LALPhenomWaveform.c \
	LALPSpinInspiralRD.c \
	LALPSpinInspiralRingdownWave.c \
	LALPrintTimeseries.c \
	LALRectangleVertices.c \
	LALRungeKutta4.c \
	LALSTPNWaveform.c \
	LALSTPNWaveform2.c \
	LALSTPNWaveformFrameless.c \
	LALSQTPNIntegrator.c \
	LALSQTPNWaveform.c \
	LALSQTPNWaveformInterface.c \
	LALTaylorF2ReducedSpin.c\
	LALTrigScanCluster.c \
	LIGOLwXMLInspiralRead.c \
	LIGOLwXMLRingdownRead.c \
	NDTemplateBank.c \
	NRWaveIO.c \
	NRWaveInject.c \
	RingUtils.c \
	SimInspiralUtils.c \
	SnglInspiralUtils.c \
	TrigScanEThincaCommon.c \
	$(END_OF_LIST)

nodist_liblalinspiral_la_SOURCES = \
	LALInspiralBuildInfoHeader.h \
	LALInspiralVCSInfo.c \
	$(END_OF_LIST)

liblalinspiral_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

DISTCLEANFILES = config.h

vcs_build_info_source = LALInspiralVCSInfo.c
vcs_info_sources = LALInspiralVCSInfoHeader.h
build_info_source = LALInspiralBuildInfoHeader.h
build_info_dependencies = $(liblalinspiral_la_DEPENDENCIES)

LDADD = liblalinspiral.la

EXTRA_DIST = \
	LALInspiralBankH4Figures.png \
	LALInspiralBankHCoarse2.png \
	LALInspiralBankHequalmass.png \
	LALInspiralBankHfine.png \
	LALInspiralConfig.h.in \
	LALInspiralHybridHexa1.png \
	LALInspiralHybridHexa2.png \
	bank.dox \
	config.h.in \
	findchirp.dox \
	inspiral.dox \
	$(END_OF_LIST)
