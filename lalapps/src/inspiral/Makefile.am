## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = inspiral.h inspiral.c
TMPLTBANKSRC = tmpltbank.c tmpltbank.h
TRIGTOTMPLTSRC = trigtotmplt.c
PRINTINSPIRALSRC = print_inspiral.c
COHERESRC = coherent_inspiral.c inspiral.h
INCASRC = inca.c
if LALFRAME
if LALMETAIO
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
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_pipe_bkgr

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = inspiral.py

lalapps_inspiral_pipe: inspiral_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe.in > lalapps_inspiral_pipe 
lalapps_inspiral_pipe_bkgr: inspiral_pipe_bkgr.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe_bkgr.in > lalapps_inspiral_pipe_bkgr

EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(TRIGTOTMPLTSRC) $(PRINTINSPIRALSRC) \
             $(INCASRC) inspiral_pipe.in s2_dag.ini \
             $(COHERESRC) $(SCRIPTS) inspiral_pipe_bkgr.in

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
