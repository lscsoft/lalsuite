## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FileIO.h PrintVector.h PrintFTSeries.h ReadFTSeries.h \
                     StreamInput.h StreamOutput.h \
                     LALInitBarycenter.h LALXMGRInterface.h \
                     LIGOLwXML.h LIGOLwXMLHeaders.h \
                     LIGOLwXMLRead.h ReadNoiseSpectrum.h  SFTfileIO.h \
                     UserInput.h ConfigFile.h LALMathematica.h SegmentsIO.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
