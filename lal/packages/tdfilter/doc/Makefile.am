## Process this file with automake to produce Makefile.in
TEXFILES = TimeDomainFilters.tex
DVIFILES = TimeDomainFilters.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES) \
        BandPassTestC.tex           DestroyZPGFilterC.tex \
        BandPassTimeSeriesH.tex     DestroyZPGFilterD.tex \
        BilinearTransformC.tex      Filters.tex \
        BilinearTransformD.tex      IIRFilterC.tex \
        ButterworthTimeSeriesC.tex  IIRFilterD.tex \
        ButterworthTimeSeriesD.tex  IIRFilterH.tex \
        CreateIIRFilterC.tex        IIRFilterVectorC.tex \
        CreateIIRFilterD.tex        IIRFilterVectorD.tex \
        CreateZPGFilterC.tex        IIRFilterVectorRC.tex \
        CreateZPGFilterD.tex        IIRFilterVectorRD.tex \
        DestroyIIRFilterC.tex       DestroyIIRFilterD.tex \
        ZPGFilterH.tex

dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

