## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = findchirp.tex \
                 FindChirpHDoc.tex FindChirpChisqHDoc.tex \
                 FindChirpSPHDoc.tex FindChirpSPDataCDoc.tex \
                 FindChirpEngineHDoc.tex FindChirpFilterCDoc.tex \
                 FindChirpEngineCDoc.tex FindChirpBCVHDoc.tex \
                 FindChirpBCVSpinDataCDoc.tex
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
