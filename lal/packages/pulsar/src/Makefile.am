noinst_LTLIBRARIES = libpulsar.la
libpulsar_la_SOURCES = \
	ApplyResampleRules.c \
	BinaryPulsarTiming.c \
	CoherentMetric.c \
	ComplexAM.c \
	ComputeFstat.c \
	ComputeSky.c \
	ComputeSkyBinary.c \
	CreateResampleRules.c \
	DTEphemeris.c \
	DestroyResampleRules.c \
	DetectorStates.c \
	ExtrapolatePulsarSpins.c \
	FitToPulsar.c \
	FlatLatticeTilingPulsar.c \
	FlatMesh.c \
	FlatPulsarMetric.c \
	FoldAmplitudes.c \
	GetEarthTimes.c \
	LALComputeAM.c \
	LALDemod.c \
	LALDemodFAST.c \
	LISAspecifics.c \
	LatticeCovering.c \
	PolycoToTimingDifference.c \
	ProjectMetric.c \
	PtoleMetric.c \
	PulsarCat.c \
	PulsarCatInput.c \
	PulsarCrossCorr.c \
	PulsarDataTypes.c \
	RulesToTimingDifference.c \
	SFTutils.c \
	StackMetric.c \
	TBaryPtolemaic.c \
	TComp.c \
	TSpin.c \
	TwoDMesh.c \
	TwoDMeshInternal.c \
	TwoDMeshMetrics.c \
	TwoDMeshRanges.c \
	UniversalDopplerMetric.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
