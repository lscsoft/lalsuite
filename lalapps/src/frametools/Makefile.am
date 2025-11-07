BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_animate \
	lalapps_calfacs \
	$(END_OF_LIST)

lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c

if FRAMEL

bin_PROGRAMS += \
	lalapps_frinfo \
	lalapps_frread \
	lalapps_frextr \
	lalapps_frjoin \
	$(END_OF_LIST)

lalapps_frinfo_SOURCES = info.c
lalapps_frread_SOURCES = read.c
lalapps_frextr_SOURCES = extr.c
lalapps_frjoin_SOURCES = join.c

endif

if LALMETAIO
if LALBURST
if LALINSPIRAL
bin_PROGRAMS += lalapps_frview
lalapps_frview_SOURCES = view.c
endif
endif
endif

if FRAMEL
if CFITSIO
bin_PROGRAMS += lalapps_xtefitstoframe
lalapps_xtefitstoframe_SOURCES = xtefitstoframe.c
endif
endif

# .1 man pages will be used automatically by lalsuite_help2man.am instead of the
# generated ones, so these files do not need to be explicitly added to man1_MANS

EXTRA_DIST += \
	lalapps_frextr.1 \
	lalapps_frinfo.1 \
	lalapps_frjoin.1 \
	lalapps_frread.1 \
	lalapps_mkcalfac.1 \
	lalapps_mkcalref.1 \
	lalapps_ring.1 \
	animate.png \
	$(END_OF_LIST)
