## Process this file with automake to produce Makefile.in

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_builddir)/src/lalapps -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

binsources = SFTReferenceLibrary.c \
             SFTReferenceLibrary.h

bin_PROGRAMS =  lalapps_SFTvalidate  \
		lalapps_SFTwrite \
		lalapps_SFTdumpheader \
		lalapps_SFTdumpall \
		lalapps_splitSFTs

lalapps_SFTvalidate_SOURCES =   SFTvalidate.c $(binsources)

lalapps_SFTwrite_SOURCES =      SFTwrite.c $(binsources)

lalapps_SFTdumpheader_SOURCES = SFTdumpheader.c $(binsources)

lalapps_SFTdumpall_SOURCES =    SFTdumpall.c $(binsources)

lalapps_splitSFTs_SOURCES =     splitSFTs.c $(binsources)

docsources      =      README 		\
                       SCCB_2004_09_09 	\
                       T040164.tex 	\
                       ligodcc.cls

EXTRA_DIST      =      $(docsources) \
                       run_validate.sh

CLEANFILES      =      SFT-test[1-7] SFT-bad[1-9] SFT-bad1[0-4] SFT-good

TESTS		=      run_validate.sh
