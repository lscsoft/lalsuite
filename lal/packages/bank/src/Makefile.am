## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libbank.la
libbank_la_SOURCES = \
InspiralSpinBank.c\
LALDeterminant3.c\
LALInspiralComputeMetric.c\
LALInspiralComputeParams.c\
LALInspiralCreateCoarseBank.c\
LALInspiralCreateFineBank.c\
LALInspiralLongestTemplateInBank.c\
LALInspiralMoments.c\
LALInspiralMomentsIntegrand.c\
LALInspiralNextTemplate.c\
LALInspiralSetParams.c\
LALInspiralSetSearchLimits.c\
LALInspiralUpdateParams.c\
LALInspiralValidParams.c\
LALInspiralValidTemplate.c\
LALInverse3.c\
LALMatrixTransform.c\
LALRectangleVertices.c\
InspiralBankGeneration.c\
InspiralSpinBankwNDTemplateBank.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
