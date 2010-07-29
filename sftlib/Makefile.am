## Process this file with automake to produce Makefile.in
#
# $Id$
#

srcdir = @srcdir@
VPATH = @srcdir@

lib_LIBRARIES = libSFTReferenceLibrary.a

libSFTReferenceLibrary_a_SOURCES = SFTReferenceLibrary.c \
              SFTReferenceLibrary.h

bin_PROGRAMS =  SFTvalidate  \
		SFTwrite \
		SFTdumpheader \
		SFTdumpall

SFTvalidate_SOURCES =  SFTvalidate.c

SFTwrite_SOURCES =     SFTwrite.c

SFTdumpheader_SOURCES = SFTdumpheader.c

SFTdumpall_SOURCES =  SFTdumpall.c

LDADD = $(lib_LIBRARIES)	

#handling our documentation files
###docsdir=$(prefix)/share/doc/sftlib

docname         = T040164

docsources      =      README 		\
                       T040164.tex 	\
                       SCCB_2004_09_09 	\
                       $(docname).tex 	\
                       ligodcc.cls 	\
                       pdfdraftcopy.sty 

doc_DATA        =      $(docname).pdf 	\
                       $(docsources)

#rule for our pdf documentation-file
$(docname).pdf : $(docname).tex ligodcc.cls pdfdraftcopy.sty 
	cd $(srcdir); $(PDFLATEX) $(docname).tex 
	cd $(srcdir); $(PDFLATEX) $(docname).tex

# add documentation-sources to distribution
EXTRA_DIST      =      $(doc_DATA) \
                       run_validate.sh

CLEANFILES      =      SFT-test[1-7] SFT-bad[1-9] SFT-bad1[0-4] SFT-good *.aux *.log $(docname).pdf

TESTS		=      run_validate.sh
