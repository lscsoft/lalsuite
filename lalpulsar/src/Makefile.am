AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'

if LALXML
XMLHDR = LALXMLVOTableSerializersPulsar.h
XMLSRC = LALXMLVOTableSerializersPulsar.c
endif

lalpulsarincludedir = $(includedir)/lal
lalpulsarinclude_HEADERS = \
	BinaryPulsarTiming.h \
	ComplexAM.h \
	ComputeFstat.h \
	ComputeSky.h \
	ComputeSkyBinary.h \
	DetectorStates.h \
	DopplerFullScan.h \
	DopplerLatticeCovering.h \
	DopplerScan.h \
	ExtrapolatePulsarSpins.h \
	FitToPulsar.h \
	FlatLatticeTilingPulsar.h \
	FlatMesh.h \
	FlatPulsarMetric.h \
	FoldAmplitudes.h \
	GeneratePulsarSignal.h \
	SimulatePulsarSignal.h \
	GenerateSpinOrbitCW.h \
	GenerateTaylorCW.h \
	HoughMap.h \
	LALBarycenter.h \
	LALComputeAM.h \
	LALDemod.h \
	LALHough.h \
	LALInitBarycenter.h \
	LALPulsarConfig.h \
	LALPulsarLibVCSInfo.h \
	LALPulsarVCSInfo.h \
	LALXMGRInterface.h \
	LISAspecifics.h \
	LUT.h \
	LatticeCovering.h \
	NormalizeSFTRngMed.h \
	PHMD.h \
	ProbabilityDensity.h \
	PtoleMetric.h \
	PulsarCat.h \
	PulsarCrossCorr.h \
	PulsarCrossCorr_v2.h \
	PulsarDataTypes.h \
	PulsarSimulateCoherentGW.h \
	PulsarTimes.h \
	SFTClean.h \
	SFTfileIO.h \
	SFTutils.h \
	SFTutils-LAL.h \
	StackMetric.h \
	Statistics.h \
	SynthesizeCWDraws.h \
	TransientCW_utils.h \
	TwoDMesh.h \
	UniversalDopplerMetric.h \
	Velocity.h \
	$(XMLHDR) \
	lalpulsarswig_test.i

lib_LTLIBRARIES = liblalpulsar.la
liblalpulsar_la_SOURCES = \
	BinaryPulsarTiming.c \
	CoherentMetric.c \
	ComplexAM.c \
	ComputeFstat.c \
	ComputeSky.c \
	ComputeSkyBinary.c \
	ConstructPLUT.c \
	DTEphemeris.c \
	DetectorStates.c \
	DopplerFullScan.c \
	DopplerLatticeCovering.c \
	DopplerScan.c \
	DriveHough.c \
	ExtrapolatePulsarSpins.c \
	FitToPulsar.c \
	FlatLatticeTilingPulsar.c \
	FlatMesh.c \
	FlatPulsarMetric.c \
	FoldAmplitudes.c \
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
	LISAspecifics.c \
	LatticeCovering.c \
	NDParamPLUT.c \
	NormalizeSFTRngMed.c \
	ParamPLUT.c \
	PatchGrid.c \
	Peak2PHMD.c \
	ProbabilityDensity.c \
	ProjectMetric.c \
	PtoleMetric.c \
	PulsarCat.c \
	PulsarCatInput.c \
	PulsarCrossCorr.c \
	PulsarCrossCorr_v2.c \
	PulsarDataTypes.c \
	PulsarSimulateCoherentGW.c \
	SFTClean.c \
	SFTfileIO.c \
	SFTutils.c \
	SimulatePulsarSignal.c \
	StackMetric.c \
	Statistics.c \
	Stereographic.c \
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

EXTRA_liblalpulsar_la_SOURCES = \
	SFTutils-LAL.c

liblalpulsar_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalpulsar-user-env.shell \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

sysconf_DATA = \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

BUILT_SOURCES = \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

lalpulsar-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALPULSAR_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalpulsar-user-env.csh lalpulsar-user-env.sh: lalpulsar-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALPulsarConfig.h.in \
	config.h.in \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALPulsarVCSInfo.h
MOSTLYCLEANFILES += LALPulsarVCSInfo.h
LALPulsarVCSInfo.h .PHONY : gen-LALPulsarVCSInfo.h
gen-LALPulsarVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALPulsarVCSInfo.h $(srcdir)/LALPulsarVCSInfo.h.git
endif

DISTCLEANFILES = config.h

@LAL_HEADER_LINK@

swig_headers = $(lalpulsarinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
