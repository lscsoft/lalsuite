INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_StringSearch
lalapps_StringSearch_SOURCES = StringSearch.c
endif
endif

if HAVE_PYTHON
bin_SCRIPTS = lalapps_cosmicstring_pipe
pkgpython_DATA = cosmicstring.py
endif

pkgdata_DATA = S3_H1_StringDag.ini S3_StringDag.ini S4_StringDag.ini

EXTRA_DIST = cosmicstring.py $(DATA) lalapps_cosmicstring_pipe \
             bread_script.sh eff_script.sh \
             efficiency.m false_alarm.m string_injections.m

LDADD = $(top_builddir)/src/lalapps/liblalapps.la
