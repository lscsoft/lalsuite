## Process this file with automake to produce Makefile.in
TEXFILES = \
ChooseModel.tex \
InspiralDerivatives.tex \
InspiralParameterCalc.tex \
InspiralPhase.tex \
InspiralSetup.tex \
InspiralVelocity.tex \
InspiralWave.tex \
MakeTemplate.tex \
PhiofVIntegrand.tex \
RungeKutta4.tex \
TappRpnTdomFreq.tex \
TappRpnTdomFreqPhase.tex \
TappRpnTdomFreqToff.tex \
TappRpnTdomTime.tex \
TappRpnTdomTimeFrequency.tex \
TappRpnTdomTimePhase.tex \
TimeDomain2.tex \
TofV.tex \
TofVIntegrand.tex

DVIFILES = \
ChooseModel.dvi \
InspiralDerivatives.dvi \
InspiralParameterCalc.dvi \
InspiralPhase.dvi \
InspiralSetup.dvi \
InspiralVelocity.dvi \
InspiralWave.dvi \
MakeTemplate.dvi \
PhiofVIntegrand.dvi \
RungeKutta4.dvi \
TappRpnTdomFreq.dvi \
TappRpnTdomFreqPhase.dvi \
TappRpnTdomFreqToff.dvi \
TappRpnTdomTime.dvi \
TappRpnTdomTimeFrequency.dvi \
TappRpnTdomTimePhase.dvi \
TimeDomain2.dvi \
TofV.dvi \
TofVIntegrand.dvi

SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = \
ChooseModel.dvi ChooseModel.pdf \
InspiralDerivatives.dvi InspiralDerivatives.pdf \
InspiralParameterCalc.dvi InspiralParameterCalc.pdf \
InspiralPhase.dvi InspiralPhase.pdf \
InspiralSetup.dvi InspiralSetup.pdf \
InspiralVelocity.dvi InspiralVelocity.pdf \
InspiralWave.dvi InspiralWave.pdf \
MakeTemplate.dvi MakeTemplate.pdf \
PhiofVIntegrand.dvi PhiofVIntegrand.pdf \
RungeKutta4.dvi RungeKutta4.pdf \
TappRpnTdomFreq.dvi TappRpnTdomFreq.pdf \
TappRpnTdomFreqPhase.dvi TappRpnTdomFreqPhase.pdf \
TappRpnTdomFreqToff.dvi TappRpnTdomFreqToff.pdf \
TappRpnTdomTime.dvi TappRpnTdomTime.pdf \
TappRpnTdomTimeFrequency.dvi TappRpnTdomTimeFrequency.pdf \
TappRpnTdomTimePhase.dvi TappRpnTdomTimePhase.pdf \
TimeDomain2.dvi TimeDomain2.pdf \
TofV.dvi TofV.pdf \
TofVIntegrand.dvi TofVIntegrand.pdf

EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
