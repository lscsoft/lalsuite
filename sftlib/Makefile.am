## Process this file with automake to produce Makefile.in
#
# $Id$
#

@SET_MAKE@

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

docsdir=$(docdir)

docs_DATA       =      README \
                       SFT_v2_specification \
                       SCCB_2004_09_09

EXTRA_DIST      =      $(docs_DATA)

CLEANFILES      =      SFT-test[1-7] SFT-bad[1-9] SFT-bad1[0-4] SFT-good

TESTS		=      run_validate.sh