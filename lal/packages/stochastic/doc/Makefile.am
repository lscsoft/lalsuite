## Process this file with automake to produce Makefile.in
TEXFILES = Overlap.tex Dirichlet.tex
DVIFILES = Overlap.dvi Dirichlet.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Overlap.dvi Overlap.pdf Dirichlet.dvi Dirichlet.pdf
EXTRA_DIST = $(TEXFILES) \
        GEO600LHO.eps GEO600LHO.pdf \
        GEO600LLO.eps GEO600LLO.pdf \
        GEO600VIRGO.eps GEO600VIRGO.pdf \
        LLOLHO.eps LLOLHO.pdf \
        TAMA300GEO600.eps TAMA300GEO600.pdf \
        TAMA300LHO.eps TAMA300LHO.pdf \
        TAMA300LLO.eps TAMA300LLO.pdf \
        TAMA300VIRGO.eps TAMA300VIRGO.pdf \
        VIRGOLHO.eps VIRGOLHO.pdf \
        VIRGOLLO.eps VIRGOLLO.pdf \
        fig1.eps fig1.pdf \
        fig2.eps fig2.pdf \
        fig3.eps fig3.pdf

dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
