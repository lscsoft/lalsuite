## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = inspiral.h inspiral.c
TMPLTBANKSRC = tmpltbank.c tmpltbank.h
COHERESRC = coherent_inspiral.c inspiral.h
INCASRC = inca.c
SIRESRC = sire.c
INSPINJSRC = inspinj.c
SPLITBANKSRC = splitbank.c
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank \
               lalapps_inca lalapps_sire \
               lalapps_inspinj lalapps_splitbank \
               lalapps_coherent_inspiral
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inca_SOURCES = $(INCASRC)
lalapps_sire_SOURCES = $(SIRESRC)
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2 \
              lalapps_inspiral_pipe_tama

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = inspiral.py

lalapps_inspiral_pipe: inspiral_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe.in > lalapps_inspiral_pipe 
lalapps_inspiral_pipe_s2: inspiral_pipe_s2.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe_s2.in > lalapps_inspiral_pipe_s2
lalapps_inspiral_pipe_tama: inspiral_pipe_tama.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g inspiral_pipe_tama.in > lalapps_inspiral_pipe_tama

pkgdata_DATA = BNSMasses.dat inspsrcs.dat inspiral_pipe.ini s2_dag.ini s2_tama.ini

EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(SIRESRC) \
             $(INCASRC) inspiral_pipe.in s2_dag.ini \
             $(COHERESRC) inspiral_pipe_s2.in inspiral_pipe.ini \
             $(INSPINJSRC) BNSMasses.dat inspsrcs.dat \
             $(SPLITBANKSRC)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
