BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lallibs_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

if LALXML
XMLHDR = LALXMLVOTableSerializersPulsar.h
XMLSRC = LALXMLVOTableSerializersPulsar.c
endif

pkginclude_HEADERS = \
	BinaryPulsarTiming.h \
	ComplexAM.h \
	ComputeFstat.h \
	ComputeSky.h \
	CWFastMath.h \
	CWMakeFakeData.h \
	DetectorStates.h \
	DopplerFullScan.h \
	DopplerScan.h \
	ExtrapolatePulsarSpins.h \
	FitToPulsar.h \
	FlatMesh.h \
	FlatPulsarMetric.h \
	GeneratePulsarSignal.h \
	GenerateSpinOrbitCW.h \
	GenerateTaylorCW.h \
	HoughMap.h \
	LALBarycenter.h \
	LALComputeAM.h \
	LALDemod.h \
	LALHough.h \
	LALInitBarycenter.h \
	LALPulsarConfig.h \
	LALPulsarVCSInfo.h \
	LALPulsarVCSInfoHeader.h \
	LALXMGRInterface.h \
	LatticeTiling.h \
	LFTandTSutils.h \
	LineRobustStats.h \
	LISAspecifics.h \
	LUT.h \
	NormalizeSFTRngMed.h \
	PHMD.h \
	ProbabilityDensity.h \
	PtoleMetric.h \
	PulsarCrossCorr.h \
	PulsarCrossCorr_v2.h \
	PulsarDataTypes.h \
	PulsarSimulateCoherentGW.h \
	PulsarTimes.h \
	SFTClean.h \
	SFTfileIO-LAL.h \
	SFTfileIO.h \
	SFTutils-LAL.h \
	SFTutils.h \
	SimulatePulsarSignal.h \
	SSBtimes.h \
	StackMetric.h \
	Statistics.h \
	SuperSkyMetrics.h \
	SynthesizeCWDraws.h \
	TransientCW_utils.h \
	TwoDMesh.h \
	UniversalDopplerMetric.h \
	Velocity.h \
	$(XMLHDR)

lib_LTLIBRARIES = liblalpulsar.la
liblalpulsar_la_SOURCES = \
	BinaryPulsarTiming.c \
	CoherentMetric.c \
	ComplexAM.c \
	ComputeFstat.c \
	ComputeSky.c \
	ConstructPLUT.c \
	CWFastMath.c \
	CWMakeFakeData.c \
	DetectorStates.c \
	DopplerFullScan.c \
	DopplerScan.c \
	DriveHough.c \
	DTEphemeris.c \
	ExtrapolatePulsarSpins.c \
	FitToPulsar.c \
	FlatMesh.c \
	FlatPulsarMetric.c \
	GenerateEllipticSpinOrbitCW.c \
	GenerateHyperbolicSpinOrbitCW.c \
	GenerateParabolicSpinOrbitCW.c \
	GeneratePulsarSignal.c \
	GenerateSpinOrbitCW.c \
	GenerateTaylorCW.c \
	GetEarthTimes.c \
	HoughMap.c \
	LALBarycenter.c \
	LALComputeAM.c \
	LALDemod.c \
	LALDemodFAST.c \
	LALInitBarycenter.c \
	LALPulsarVCSInfo.c \
	LALXMGRInterface.c \
	LatticeTiling.c \
	LFTandTSutils.c \
	LineRobustStats.c \
	LISAspecifics.c \
	NDParamPLUT.c \
	NormalizeSFTRngMed.c \
	ParamPLUT.c \
	PatchGrid.c \
	Peak2PHMD.c \
	ProbabilityDensity.c \
	ProjectMetric.c \
	PtoleMetric.c \
	PulsarCrossCorr.c \
	PulsarCrossCorr_v2.c \
	PulsarSimulateCoherentGW.c \
	SFTClean.c \
	SFTfileIO.c \
	SFTutils.c \
	SimulatePulsarSignal.c \
	SinCosLUT.c \
	SSBtimes.c \
	StackMetric.c \
	Statistics.c \
	Stereographic.c \
	SuperSkyMetrics.c \
	SynthesizeCWDraws.c \
	TBaryPtolemaic.c \
	TComp.c \
	TransientCW_utils.c \
	TSpin.c \
	TwoDMesh.c \
	TwoDMeshInternal.c \
	TwoDMeshMetrics.c \
	TwoDMeshRanges.c \
	UniversalDopplerMetric.c \
	Velocity.c \
	$(XMLSRC)
nodist_liblalpulsar_la_SOURCES = LALPulsarBuildInfo.c

EXTRA_liblalpulsar_la_SOURCES = \
	GSLHelpers.h \
	SinCosLUT.i 				\
	ComputeFstat_Demod.c 			\
	ComputeFstat_Demod_ComputeFaFb.c 	\
	ComputeFstat_DemodHL_Altivec.i		\
	ComputeFstat_DemodHL_OptC.i		\
	ComputeFstat_DemodHL_SSE.i		\
	ComputeFstat_DemodHL_Generic.i		\
	ComputeFstat_Resamp.c 			\
	SFTfileIO-LAL.c 			\
	SFTutils-LAL.c

liblalpulsar_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALPulsarConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALPulsarVCSInfoHeader.h
build_info_source = LALPulsarBuildInfo.c
$(build_info_source): $(liblalpulsar_la_DEPENDENCIES)
