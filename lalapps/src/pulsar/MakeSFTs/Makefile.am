BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

dist_bin_SCRIPTS = lalapps_CopySFTs

if FRAME
frameprogs = lalapps_MakeSFTs
else
frameprogs =
endif

## PSS requires FrameLib anyway
if PSS
pss_progs =
pss_src = XLALPSSInterface.c
else
pss_progs =
pssinterface =
endif

bin_PROGRAMS = $(frameprogs) $(pss_progs)

if FRAME
lalapps_MakeSFTs_SOURCES = MakeSFTs.c $(pss_src)
endif

lalapps_MakeSFTs_CPPFLAGS = $(AM_CPPFLAGS)
if PSS
lalapps_MakeSFTs_CPPFLAGS += -DPSS_ENABLED
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
