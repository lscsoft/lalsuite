## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = ligolwbank.h ligolwbank.c inspiral.h inspiral.c \
              inspiralfrutils.c inspiralfrutils.h
TMPLTBANKSRC = ligolwbank.h ligolwbank.c tmpltbank.c \
               inspiralfrutils.c inspiralfrutils.h
if FRAME
if DATAFLOW
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
if PYTHON22
bin_SCRIPTS = inspiral_pipeline.py
else
bin_SCRIPTS = 
endif
EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) inspiral_pipeline.py

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
