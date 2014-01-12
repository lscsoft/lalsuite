AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@

lalinspiralincludedir = $(includedir)/lal
lalinspiralinclude_HEADERS = \
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
	InspiralInjectionParams.h \
	LALEOBNRv2Waveform.h \
	LALInspiral.h \
	LALInspiralBank.h \
	LALInspiralBuildInfo.h \
	LALInspiralConfig.h \
	LALInspiralLibVCSInfo.h \
	LALInspiralSBankOverlap.h \
	LALInspiralVCSInfo.h \
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
	SimulateInspiral.h \
	TemplateBankGeneration.h \
	TrigScanEThincaCommon.h \
	lalinspiralswig_test.i

lib_LTLIBRARIES = liblalinspiral.la
liblalinspiral_la_SOURCES = \
	BandPassInspiralTemplate.c \
	Chisq_CPU.c \
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
	FindChirpClusterEvents.c \
	FindChirpData.c \
	FindChirpFilter.c \
	FindChirpFilterInit.c \
	FindChirpFilterOutputVeto.c \
	FindChirpIMRSimulation.c \
	FindChirpLinkedList.c \
	FindChirpMemory.c \
	FindChirpPTFFilter.c \
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
	LALInspiralVCSInfo.c \
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
	SimulateInspiral.c \
	SnglInspiralUtils.c \
	SnglRingdownUtils.c \
	TrigScanEThincaCommon.c
nodist_liblalinspiral_la_SOURCES = LALInspiralBuildInfo.c

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

MOSTLYCLEANFILES = \
	lalinspiral-user-env.shell \
	lalinspiral-user-env.sh \
	lalinspiral-user-env.csh

sysconf_DATA = \
	lalinspiral-user-env.sh \
	lalinspiral-user-env.csh

BUILT_SOURCES = \
	lalinspiral-user-env.sh \
	lalinspiral-user-env.csh

lalinspiral-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# source this file to access $(PACKAGE_NAME)" >$@
	@echo "set LALINSPIRAL_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalinspiral-user-env.csh lalinspiral-user-env.sh: lalinspiral-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED="$(SED)" -v output="$@" $<

EXTRA_DIST = \
	Chisq_GPU.cu \
	LALInspiralConfig.h.in \
	config.h.in \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALInspiralVCSInfo.h
MOSTLYCLEANFILES += LALInspiralVCSInfo.h
LALInspiralVCSInfo.h .PHONY : gen-LALInspiralVCSInfo.h
gen-LALInspiralVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALInspiralVCSInfo.h $(srcdir)/LALInspiralVCSInfo.h.git
endif

DISTCLEANFILES = config.h

MOSTLYCLEANFILES += LALInspiralBuildInfo.c

LALInspiralBuildInfo.c: Makefile $(liblalinspiral_la_DEPENDENCIES)
	$(AM_V_GEN)echo '/* $@: generated by Makefile */' >$@; \
	echo 'const char *const lalInspiralConfigureArgs = "'"$(ac_configure_args)"'";' >>$@; \
	echo 'const char *const lalInspiralConfigureDate = "'"$(configure_date)"'";' >>$@; \
	echo 'const char *const lalInspiralBuildDate = "'`date +"%Y-%m-%dT%H:%M:%S%z"`'";' >>$@

@LAL_HEADER_LINK@

# cuda build rules
include $(top_srcdir)/gnuscripts/cuda.common

swig_headers = $(lalinspiralinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
