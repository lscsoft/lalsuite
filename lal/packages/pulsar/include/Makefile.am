## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FlatMesh.h TwoDMesh.h PulsarTimes.h StackMetric.h 		\
	Resample.h LALComputeAM.h ComputeSky.h ComputeSkyBinary.h 		\
	LALDemod.h FoldAmplitudes.h LALBarycenter.h PtoleMetric.h 		\
	PulsarCat.h HeterodynePulsar.h FitToPulsar.h BinaryPulsarTiming.h 	\
	PulsarDataTypes.h ExtrapolatePulsarSpins.h LatticeCovering.h		\
	ComputeFstat.h SFTutils.h DetectorStates.h LISAspecifics.h		\
	FlatPulsarMetric.h ComplexAM.h UniversalDopplerMetric.h                 \
	FlatLatticeTilingPulsar.h PulsarCrossCorr.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
