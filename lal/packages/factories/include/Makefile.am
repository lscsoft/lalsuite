## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = AVFactories.h SeqFactories.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = SeqFactories.h AVFactories.h
EXTRA_DIST = SeqFactoriesH.m4 SeqFactoriesBaseH.m4 ArraySeqFactoriesBaseH.m4 \
             AVFactoriesH.m4 VFactoriesBaseH.m4 AFactoriesBaseH.m4
all-local: include-link
dvi-local: .dvi-dep

SeqFactories.h: SeqFactoriesH.m4 SeqFactoriesBaseH.m4 ArraySeqFactoriesBaseH.m4
	cd $(srcdir) && @M4@ SeqFactoriesH.m4 > SeqFactories.h
AVFactories.h: AVFactoriesH.m4 VFactoriesBaseH.m4 AFactoriesBaseH.m4
	cd $(srcdir) && @M4@ AVFactoriesH.m4 > AVFactories.h

@INCLUDEMKFILE@
@DVIMKFILE@
