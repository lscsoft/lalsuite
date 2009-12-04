DOCSOURCEFILES = \
	support.tex \
	LALMathematicaHplot1.eps \
	LALMathematicaHplot1.pdf

MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = \
	main.dvi \
	main.pdf \
	main.ps

EXTRA_DIST = $(DOCSOURCEFILES)

distclean-local:
	-rm -fr .adoc

dvi-local: .dvi-dep

@DOCMKFILE@
