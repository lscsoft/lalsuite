include $(top_srcdir)/gnuscripts/lalapps.am

LALFRAMEPROGRAMS = \
	lalapps_animate \
	lalapps_calfacs \
	lalapps_computespec \
	lalapps_frinfo \
	lalapps_frread \
	lalapps_frextr \
	lalapps_frjoin
lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c
lalapps_computespec_SOURCES = computespec.c
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_frinfo_SOURCES = info.c
lalapps_frread_SOURCES = read.c
lalapps_frextr_SOURCES = extr.c
lalapps_frjoin_SOURCES = join.c
#lalapps_framestats_SOURCES = framestats.c

#if LALPULSAR
#LALPULSARPROGRAMS = lalapps_combinexteframes
#lalapps_combinexteframes_SOURCES = combinexteframes.c
#endif

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
	$(LALFRAMECFITSIOPROGRAMS)

EXTRA_DIST = \
	animate.c \
	animate.png \
	calfacs.c \
	combinexteframes.c \
	computespec.c \
	extr.c \
	framestats.c \
	info.c \
	join.c \
	lalDebugLevel.c \
	read.c \
	stripAdc.c \
	$(END_OF_LIST)
