INCLUDES = -I$(top_builddir)/include

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
	GenerateSpinOrbitCW.h \
	GenerateTaylorCW.h \
	HoughMap.h \
	LALComputeAM.h \
	LALDemod.h \
	LALHough.h \
	LALPulsarConfig.h \
	LALPulsarVCSInfo.h \
	LALXMGRInterface.h \
	LISAspecifics.h \
	LUT.h \
	LatticeCovering.h \
	NormalizeSFTRngMed.h \
	PHMD.h \
	PtoleMetric.h \
	PulsarCat.h \
	PulsarCrossCorr.h \
	PulsarDataTypes.h \
	PulsarTimes.h \
	Resample.h \
	SFTClean.h \
	SFTfileIO.h \
	SFTutils.h \
	StackMetric.h \
	Statistics.h \
	TwoDMesh.h \
	UniversalDopplerMetric.h \
	Velocity.h \
	$(XMLHDR)

lib_LTLIBRARIES = liblalpulsar.la
liblalpulsar_la_SOURCES = \
	ApplyResampleRules.c \
	BinaryPulsarTiming.c \
	CoherentMetric.c \
	ComplexAM.c \
	ComputeFstat.c \
	ComputeSky.c \
	ComputeSkyBinary.c \
	ConstructPLUT.c \
	CreateResampleRules.c \
	DTEphemeris.c \
	DestroyResampleRules.c \
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
	LALComputeAM.c \
	LALDemod.c \
	LALDemodFAST.c \
	LALPulsarVCSInfo.c \
	LALXMGRInterface.c \
	LISAspecifics.c \
	LatticeCovering.c \
	NDParamPLUT.c \
	NormalizeSFTRngMed.c \
	ParamPLUT.c \
	PatchGrid.c \
	Peak2PHMD.c \
	PolycoToTimingDifference.c \
	ProjectMetric.c \
	PtoleMetric.c \
	PulsarCat.c \
	PulsarCatInput.c \
	PulsarCrossCorr.c \
	PulsarDataTypes.c \
	RulesToTimingDifference.c \
	SFTClean.c \
	SFTfileIO.c \
	SFTutils.c \
	SimulatePulsarSignal.c \
	StackMetric.c \
	Statistics.c \
	Stereographic.c \
	TBaryPtolemaic.c \
	TComp.c \
	TSpin.c \
	TwoDMesh.c \
	TwoDMeshInternal.c \
	TwoDMeshMetrics.c \
	TwoDMeshRanges.c \
	UniversalDopplerMetric.c \
	Velocity.c \
	$(XMLSRC)

liblalpulsar_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

sysconf_DATA = \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALPulsarVCSInfo.h \
	lalpulsar-user-env.sh \
	lalpulsar-user-env.csh

lalpulsar-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Stochastic" >$@ ;\
	echo "LALPULSAR_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALPULSAR_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalpulsar-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Stochastic" >$@ ;\
	echo "setenv LALPULSAR_PREFIX \"$(prefix)\"" >>$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?LD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?DYLD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?PKG_CONFIG_PATH ) then" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig\"" >>$@ ;\
	echo "endif" >>$@

EXTRA_DIST = \
	LALPulsarConfig.h.in \
	LALPulsarVCSInfo.h.in \
	generate_vcs_info.py

LALPulsarVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALPulsarConfig.h
	@echo "creating LALPulsarVCSInfo.h"
	@$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALPulsar $(srcdir) $$PWD

@LAL_HEADER_LINK@
