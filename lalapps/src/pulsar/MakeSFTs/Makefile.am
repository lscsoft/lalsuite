BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

dist_bin_SCRIPTS = lalapps_CopySFTs

bin_PROGRAMS =

if FRAMEL
bin_PROGRAMS += lalapps_MakeSFTs
endif

lalapps_MakeSFTs_SOURCES = MakeSFTs.c
lalapps_MakeSFTs_CPPFLAGS = $(AM_CPPFLAGS)

if FRAMEL
if PSS
lalapps_MakeSFTs_SOURCES += XLALPSSInterface.c
lalapps_MakeSFTs_CPPFLAGS += -DPSS_ENABLED
endif
endif

if HAVE_PYTHON
pybin_scripts = lalapps_MakeSFTDAG
pkgpython_PYTHON =
endif

# Add shell test scripts to this variable
test_scripts += testMakeSFTs.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# These tests require FRAME
if !FRAMEL
skip_tests += testMakeSFTs.sh
endif
