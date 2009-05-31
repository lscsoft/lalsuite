## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FileIO.h PrintVector.h PrintFTSeries.h ReadFTSeries.h \
                     StreamInput.h StreamOutput.h LALCache.h \
                     LALInitBarycenter.h LALXMGRInterface.h \
                     LIGOLwXML.h LIGOLwXMLHeaders.h LIGOLwXMLInspiralHeaders.h \
		     LIGOLwXMLArray.h LIGOLwXMLRead.h ReadNoiseSpectrum.h \
		     SFTfileIO.h LogPrintf.h UserInput.h ConfigFile.h \
		     LALMathematica.h SegmentsIO.h Audio.h NRWaveIO.h	\
                     DopplerScan.h DopplerFullScan.h DopplerLatticeCovering.h \
                     VeryBasicXMLOutput.h ReadFiltersFile.h \
                     LALXML.h LALXMLVOTableCommon.h LALXMLVOTableSerializers.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
