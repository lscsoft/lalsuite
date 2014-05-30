LDADD =
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am


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
