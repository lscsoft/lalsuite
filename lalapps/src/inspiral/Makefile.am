## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = ligolwbank.h ligolwbank.c inspiral.h inspiral.c \
              inspiralfrutils.c inspiralfrutils.h
TMPLTBANKSRC = ligolwbank.h ligolwbank.c tmpltbank.c \
               inspiralfrutils.c inspiralfrutils.h
TRIGTOTMPLTSRC = trigtotmplt.c ligolwbank.h ligolwbank.c
PRINTINSPIRALSRC = print_inspiral.c ligolwbank.h ligolwbank.c
COHERESRC = coherent_inspiral.c ligolwbank.h ligolwbank.c inspiral.h \
            inspiralfrutils.c inspiralfrutils.h
INCASRC = inca.c
if FRAME
if DATAFLOW
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank lalapps_trigtotmplt \
               lalapps_print_inspiral lalapps_inca lalapps_coherent_inspiral
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_trigtotmplt_SOURCES = $(TRIGTOTMPLTSRC)
lalapps_print_inspiral_SOURCES = $(PRINTINSPIRALSRC)
lalapps_inca_SOURCES = $(INCASRC)
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_inspiral_pipe

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = inspiral.py

lalapps_inspiral_pipe: inspiral_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe.in > lalapps_inspiral_pipe 

EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(TRIGTOTMPLTSRC) $(PRINTINSPIRALSRC) \
             $(INCASRC) inspiral_pipe.in s2_dag.ini \
             $(COHERESRC)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
