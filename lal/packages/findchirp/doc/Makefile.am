## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = findchirp.tex \
                 FindChirpChisqHDoc.tex \
                 FindChirpSPHDoc.tex FindChirpSPDataCDoc.tex \
                 FindChirpEngineHDoc.tex FindChirpFilterCDoc.tex \
                 FindChirpEngineCDoc.tex FindChirpBCVHDoc.tex \
                 FindChirpBCVSpinDataCDoc.tex \
                 FindChirpDataCDoc.tex FindChirpFilterInitCDoc.tex \
                 FindChirpBCVDataCDoc.tex FindChirpBCVFilterCDoc.tex \
                 FindChirpBCVSpinHDoc.tex
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
@DOCMKFILE@
