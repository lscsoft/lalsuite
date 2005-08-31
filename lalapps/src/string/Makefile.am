## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
if LALFRAME
if LALMETAIO
RINGPROG = 
STRINGPROG = lalapps_StringSearch
else
RINGPROG =
STRINGPROG =
endif
else
RINGPROG =
STRINGPROG =
endif
else
RINGPROG =
STRINGPROG =
endif
bin_SCRIPTS = lalapps_cosmicstring_pipe
bin_PROGRAMS = $(STRINGPROG)
if FRAME
if LALFRAME
if LALMETAIO
lalapps_StringSearch_SOURCES = StringSearch.c
endif
endif
endif
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = cosmicstring.py
pkgdata_DATA = S3_H1_StringDag.ini S3_StringDag.ini S4_StringDag.ini
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
CLEANFILES = python_config.sed lalapps_cosmicstring_pipe
EXTRA_DIST = StringSearch.c cosmicstring_pipe.in \
             $(SCRIPTS) $(DATA) \
             bread_script.sh eff_script.sh \
             efficiency.m false_alarm.m string_injections.m

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_cosmicstring_pipe: cosmicstring_pipe.in python_config.sed
	sed -f python_config.sed cosmicstring_pipe.in > lalapps_cosmicstring_pipe 
	chmod +x lalapps_cosmicstring_pipe
