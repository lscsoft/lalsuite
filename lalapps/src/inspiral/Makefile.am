## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = inspiral.c
TMPLTBANKSRC = tmpltbank.c
COHERESRC = coherent_inspiral.c
INCASRC = inca.c
SIRESRC = sire.c
INSPINJSRC = inspinj.c
MINJSRC = minj.c
BBHINJSRC = bbhinj.c
SPLITBANKSRC = splitbank.c
TAMA2LIGOLWSRC = tama2ligolw.c
INSPFRINJSRC = inspfrinj.c
GWF2XMLSRC = gwf2xml.c
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank \
               lalapps_inca lalapps_sire \
               lalapps_inspinj lalapps_splitbank \
               lalapps_minj lalapps_bbhinj \
	             lalapps_tama2ligolw lalapps_inspfrinj \
               lalapps_gwf2xml
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inca_SOURCES = $(INCASRC)
lalapps_sire_SOURCES = $(SIRESRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_minj_SOURCES = $(MINJSRC)
lalapps_bbhinj_SOURCES = $(BBHINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
lalapps_tama2ligolw_SOURCES = $(TAMA2LIGOLWSRC)
lalapps_inspfrinj_SOURCES = $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES = $(GWF2XMLSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2 \
              lalapps_inspiral_pipe_tama lalapps_inspiral_pipe_geo \
              lalapps_inspiral_hipe

pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = inspiral.py

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_inspiral_pipe: inspiral_pipe.in python_config.sed
	sed -f python_config.sed inspiral_pipe.in > lalapps_inspiral_pipe 
	chmod +x lalapps_inspiral_pipe
lalapps_inspiral_pipe_s2: inspiral_pipe_s2.in python_config.sed
	sed -f python_config.sed inspiral_pipe_s2.in > lalapps_inspiral_pipe_s2
	chmod +x lalapps_inspiral_pipe_s2
lalapps_inspiral_pipe_tama: inspiral_pipe_tama.in python_config.sed
	sed -f python_config.sed inspiral_pipe_tama.in > lalapps_inspiral_pipe_tama
	chmod +x lalapps_inspiral_pipe_tama
lalapps_inspiral_pipe_geo: inspiral_pipe_geo.in python_config.sed
	sed -f python_config.sed inspiral_pipe_geo.in > lalapps_inspiral_pipe_geo
	chmod +x lalapps_inspiral_pipe_geo
lalapps_inspiral_hipe: inspiral_hipe.in python_config.sed
	sed -f python_config.sed inspiral_hipe.in > lalapps_inspiral_hipe
	chmod +x lalapps_inspiral_hipe

pkgdata_DATA = BNSMasses.dat inspsrcs.dat inspiral_pipe.ini s2_dag.ini s2_tama.ini inspiral_hipe.ini 
EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(SIRESRC) $(MINJSRC) \
             $(INCASRC) inspiral_pipe.in s2_dag.ini \
             $(COHERESRC) inspiral_pipe_s2.in inspiral_pipe.ini \
             $(INSPINJSRC) BNSMasses.dat inspsrcs.dat \
             $(SPLITBANKSRC)                                                                                 

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2 \
             lalapps_inspiral_pipe_tama lalapps_inspiral_pipe_geo \
             lalapps_inspiral_hipe python_config.sed

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
