BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am
include $(top_srcdir)/gnuscripts/lalsuite_cuda.am

pkginclude_HEADERS = \
	BBHPhenomCoeffs.h \
	Chisq_CPU.h \
	Chisq_GPU.h \
	LALInspiralPyCBCTemplate.h \
	CoherentInspiral.h \
	CoincInspiralEllipsoid.h \
	FindChirp.h \
	FindChirpACTD.h \
	FindChirpBCV.h \
	FindChirpBCVSpin.h \
	FindChirpChisq.h \
	FindChirpDatatypes.h \
	FindChirpFilterOutputVeto.h \
	FindChirpPTF.h \
	FindChirpSP.h \
	FindChirpTD.h \
	GenerateInspRing.h \
	GenerateInspiral.h \
	GeneratePPNInspiral.h \
	GenerateRing.h \
	Inject.h \
	InspiralInjectionParams.h \
	LALEOBNRv2Waveform.h \
	LALInspiral.h \
	LALInspiralBank.h \
	LALInspiralConfig.h \
	LALInspiralVCSInfo.h \
	LALInspiralSBankOverlap.h \
	LALInspiralVCSInfoHeader.h \
	LALNoiseModelsInspiral.h \
	LALPSpinInspiralRD.h \
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
	Chisq_CPU.c \
	InjectTimeSeries.c \
	LALInspiralPyCBCTemplate.c \
	CoherentInspiralFilter.c \
	CoherentInspiralInput.c \
	CoincInspiralEllipsoid.c \
	CoincInspiralUtils.c \
	CoincRingdownUtils.c \
	FindChirpBCVCFilter.c \
	FindChirpBCVChisq.c \
	FindChirpBCVData.c \
	FindChirpBCVFilter.c \
	FindChirpBCVSpinData.c \
	FindChirpBCVSpinFilter.c \
	FindChirpBCVSpinTemplate.c \
	FindChirpBCVTemplate.c \
	FindChirpBankVeto.c \
	FindChirpChisq.c \
	FindChirpChisqInit.c \
	FindChirpData.c \
	FindChirpFilter.c \
	FindChirpFilterInit.c \
	FindChirpFilterOutputVeto.c \
	FindChirpIMRSimulation.c \
	FindChirpLinkedList.c \
	FindChirpMemory.c \
	FindChirpPTFTemplate.c \
	FindChirpPTFWaveform.c \
	FindChirpPhenomWaveform.c \
	FindChirpSPData.c \
	FindChirpSPTemplate.c \
	FindChirpSimulation.c \
	FindChirpStoreEvent.c \
	FindChirpTDData.c \
	FindChirpTDTemplate.c \
	FindChirpTemplate.c \
	GenerateInspRing.c \
	GenerateInspiral.c \
	GenerateInspiralSmooth.c \
	GeneratePPNAmpCorInspiral.c \
	GeneratePPNInspiral.c \
	GenerateRing.c \
	GetErrorMatrixFromSnglInspiral.c \
	GetInspiralParams.c \
	InspiralBankGeneration.c \
	InspiralInjectionParams.c \
	InspiralSpinBank.c \
	LALBCVWaveform.c \
	LALDeterminant3.c \
	LALDiscoverInspiralEvents.c \
	LALEOBGetFactorizedWaveform.c \
	LALEOBWaveform.c \
	LALEOBNonQCCorrection.c \
	LALEOBPPWaveform.c \
	LALEstimateEffectiveDistance.c \
	LALEtaTau02.c \
	LALEtaTau04.c \
	LALEvaluateInspiralChisqTest.c \
	LALFactorizedFlux.c \
	LALHexagonVertices.c \
	LALInsidePolygon.c \
	LALInspiralAmplitude.c \
	LALInspiralAmplitudeCorrectedWave.c \
	LALInspiralBCVBank.c \
	LALInspiralBankList.c \
	LALInspiralBankUtils.c \
	LALInspiralChooseModel.c \
	LALInspiralComputeChisq.c \
	LALInspiralComputeMetric.c \
	LALInspiralComputePTFMetric.c \
	LALInspiralComputeParams.c \
	LALInspiralComputeSNRIntegrand.c \
	LALInspiralCreateCoarseBank.c \
	LALInspiralCreateFineBank.c \
	LALInspiralDerivatives.c \
	LALInspiralEccentricity.c \
	LALInspiralFindEvents.c \
	LALInspiralFindEventsCluster.c \
	LALInspiralFindLoudestEvent.c \
	LALInspiralFrequency3.c \
	LALInspiralGetApproximantString.c \
	LALInspiralHexagonalBank.c \
	LALInspiralHybridHexagonalBank.c \
	LALInspiralIIR.c \
	LALInspiralInit.c \
	LALInspiralLongestTemplateInBank.c \
	LALInspiralMoments.c \
	LALInspiralMomentsIntegrand.c \
	LALInspiralNextTemplate.c \
	LALInspiralNinjaInject.c \
	LALInspiralParameterCalc.c \
	LALInspiralParseParameters.c \
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
	LALInverse3.c \
	LALNewtonianMultipole.c \
	LALMatrixTransform.c \
	LALPhenomWaveform.c \
	LALPSpinInspiralRD.c \
	LALPSpinInspiralRingdownWave.c \
	LALPrintTimeseries.c \
	LALRandomInspiralSignal.c \
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
	LALTruncateInvSpectrum.c \
	LIGOLwXMLInspiralRead.c \
	LIGOLwXMLRingdownRead.c \
	MultiInspiralUtils.c \
	NDTemplateBank.c \
	NRWaveIO.c \
	NRWaveInject.c \
	Playground.c \
	RingUtils.c \
	SimInspiralUtils.c \
	SimRingdownUtils.c \
	SnglInspiralUtils.c \
	SnglRingdownUtils.c \
	TrigScanEThincaCommon.c \
	$(END_OF_LIST)

nodist_liblalinspiral_la_SOURCES = \
	LALInspiralBuildInfoHeader.h \
	LALInspiralVCSInfo.c \
	$(END_OF_LIST)

if CUDA
INSPCXXSRC =
INSPCXXGENSRC = Chisq_GPU.cpp
INSPLIBCXX = liblalinspiralcxx.la
else
INSPCXXSRC =
INSPCXXGENSRC =
INSPLIBCXX =
endif

noinst_LTLIBRARIES = $(INSPLIBCXX)
liblalinspiralcxx_la_SOURCES = $(INSPCXXSRC)
nodist_liblalinspiralcxx_la_SOURCES = $(INSPCXXGENSRC)
liblalinspiral_la_LIBADD = $(INSPLIBCXX)

liblalinspiral_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

CLEANFILES = \
	Chisq_GPU.cpp

DISTCLEANFILES = config.h

vcs_build_info_source = LALInspiralVCSInfo.c
vcs_info_sources = LALInspiralVCSInfoHeader.h
build_info_source = LALInspiralBuildInfoHeader.h
build_info_dependencies = $(liblalinspiral_la_DEPENDENCIES)

LDADD = liblalinspiral.la

bin_PROGRAMS = \
	lalinspiral_version \
	$(END_OF_LIST)

lalinspiral_version_SOURCES = version.c

TESTS = \
	lalinspiral_version \
	$(END_OF_LIST)

EXTRA_DIST = \
	Chisq_GPU.cu \
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
