BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalapps.am

if LALFRAME
if LALMETAIO
if LALBURST
bin_PROGRAMS = lalapps_StringSearch lalapps_StringAddFrame 
lalapps_StringSearch_SOURCES = StringSearch.c
lalapps_StringAddFrame_SOURCES = AddFrame.c
endif
endif
endif

if HAVE_PYTHON
pybin_SCRIPTS = \
	lalapps_cosmicstring_pipe \
	lalapps_string_cs_gamma
pkgpython_PYTHON = cosmicstring.py
endif

#pkgdata_DATA = S3_H1_StringDag.ini S3_StringDag.ini S4_StringDag.ini

EXTRA_DIST += cosmicstring.py $(DATA)
