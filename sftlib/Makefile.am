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

SFTvalidate_SOURCES =  SFTReferenceLibrary.c \
                       SFTReferenceLibrary.h \
                       SFTvalidate.c

SFTwrite_SOURCES =     SFTReferenceLibrary.c \
                       SFTReferenceLibrary.h \
	               SFTwrite.c

SFTdumpheader_SOURCES = SFTReferenceLibrary.c \
                        SFTReferenceLibrary.h \
		        SFTdumpheader.c

SFTdumpall_SOURCES =      SFTReferenceLibrary.c \
                       SFTReferenceLibrary.h \
		       SFTdumpall.c
docsdir=$(docdir)

docs_DATA       =      README \
                       SFT_v2_specification

EXTRA_DIST      =      $(docs_DATA)
