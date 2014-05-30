if LALFRAME
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
endif

#if LALFRAME
#if LALPULSAR
#LALPULSARPROGRAMS = lalapps_combinexteframes
#lalapps_combinexteframes_SOURCES = combinexteframes.c
#endif
#endif

if LALFRAME
if LALMETAIO
if LALBURST
if LALINSPIRAL
LALFRAMEMETAIOBURSTPROGRAMS = lalapps_frview
lalapps_frview_SOURCES = view.c
endif
endif
endif
endif

if LALFRAME
if CFITSIO
LALFRAMECFITSIOPROGRAMS = lalapps_xtefitstoframe
lalapps_xtefitstoframe_SOURCES = xtefitstoframe.c
endif
endif

bin_PROGRAMS = \
	$(LALFRAMEPROGRAMS) \
	$(LALFRAMEMETAIOBURSTPROGRAMS) \
	$(LALFRAMECFITSIOPROGRAMS)

EXTRA_DIST = \
	stripAdc.c \
	lalDebugLevel.c \
	info.c \
	read.c \
	extr.c \
	join.c \
	animate.c \
	calfacs.c \
	combinexteframes.c \
	framestats.c \
	computespec.c

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
