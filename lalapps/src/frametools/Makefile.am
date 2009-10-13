if FRAME
FRAMEPROGRAMS = \
	lalapps_frinfo \
	lalapps_frread \
	lalapps_frextr \
	lalapps_frjoin
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_frinfo_SOURCES = info.c
lalapps_frread_SOURCES = read.c
lalapps_frextr_SOURCES = extr.c
lalapps_frjoin_SOURCES = join.c
else
FRAMEPROGRAMS =
endif

if LALFRAME
LALFRAMEPROGRAMS = \
	lalapps_animate \
	lalapps_calfacs \
	lalapps_computespec
lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c
lalapps_computespec_SOURCES = computespec.c
else
LALFRAMEPROGRAMS =
endif

if LALFRAME
if LALMETAIO
LALFRAMEMETAIOPROGRAMS = lalapps_frview
lalapps_frview_SOURCES = view.c
else
LALFRAMEMETAIOPROGRAMS =
endif
endif

bin_PROGRAMS = \
	$(FRAMEPROGRAMS) \
	$(LALFRAMEPROGRAMS) \
	$(LALFRAMEMETAIOPROGRAMS)

EXTRA_DIST = \
	stripAdc.c \
	lalDebugLevel.c \
	info.c \
	read.c \
	extr.c \
	join.c \
	animate.c \
	calfacs.c \
	computespec.c

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
