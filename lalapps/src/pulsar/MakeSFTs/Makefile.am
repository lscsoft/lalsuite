BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

dist_bin_SCRIPTS = lalapps_CopySFTs

bin_PROGRAMS =

if FRAME
bin_PROGRAMS += lalapps_MakeSFTs
endif

lalapps_MakeSFTs_SOURCES = MakeSFTs.c
lalapps_MakeSFTs_CPPFLAGS = $(AM_CPPFLAGS)

if FRAME
if PSS
lalapps_MakeSFTs_SOURCES += XLALPSSInterface.c
lalapps_MakeSFTs_CPPFLAGS += -DPSS_ENABLED
endif
endif

if HAVE_PYTHON
pybin_scripts = lalapps_MakeSFTDAG
pkgpython_PYTHON =
endif

EXTRA_DIST += \
	testMakeSFTs.sh \
	$(END_OF_LIST)

if FRAME
TESTS = \
	testMakeSFTs.sh \
	$(END_OF_LIST)
endif
