BUILT_SOURCES =
MOSTLYCLEANFILES =
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

EXTRA_DIST = \
	animate.png \
	$(END_OF_LIST)
