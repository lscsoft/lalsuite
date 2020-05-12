BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

LALFRAMEPROGRAMS = \
	lalapps_animate \
	lalapps_calfacs \
	lalapps_frinfo \
	lalapps_frread \
	lalapps_frextr \
	lalapps_frjoin \
	$(END_OF_LIST)

lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c
lalapps_frinfo_SOURCES = info.c
lalapps_frread_SOURCES = read.c
lalapps_frextr_SOURCES = extr.c
lalapps_frjoin_SOURCES = join.c

if LALMETAIO
if LALBURST
if LALINSPIRAL
LALFRAMEMETAIOBURSTPROGRAMS = lalapps_frview
lalapps_frview_SOURCES = view.c
endif
endif
endif

if CFITSIO
LALFRAMECFITSIOPROGRAMS = lalapps_xtefitstoframe
lalapps_xtefitstoframe_SOURCES = xtefitstoframe.c
endif

bin_PROGRAMS = \
	$(LALFRAMEPROGRAMS) \
	$(LALFRAMEMETAIOBURSTPROGRAMS) \
	$(LALFRAMECFITSIOPROGRAMS) \
	$(END_OF_LIST)

EXTRA_DIST = \
	animate.c \
	animate.png \
	calfacs.c \
	extr.c \
	info.c \
	join.c \
	read.c \
	$(END_OF_LIST)
