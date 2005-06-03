## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC 	= inspiral.c
TMPLTBANKSRC 	= tmpltbank.c
COHERESRC 	= coherent_inspiral.c
COHEREBANK      = coherentbank.c
INCASRC 	= inca.c
THINCASRC 	= thinca.c
LTINCASRC 	= lt_inca.c
SIRESRC 	= sire.c
INSPINJSRC 	= inspinj.c
MINJSRC 	= minj.c
BBHINJSRC 	= bbhinj.c
SPLITBANKSRC 	= splitbank.c
TRIGBANKSRC 	= trigbank.c
TAMA2LIGOLWSRC 	= tama2ligolw.c
INSPFRINJSRC 	= inspfrinj.c
GWF2XMLSRC 	= gwf2xml.c
SPININJSRC 	= spininj.c
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank \
               lalapps_inca lalapps_sire \
               lalapps_inspinj lalapps_splitbank \
               lalapps_minj lalapps_bbhinj \
	       lalapps_tama2ligolw lalapps_inspfrinj \
               lalapps_gwf2xml lalapps_lt_inca \
	       lalapps_thinca lalapps_trigbank \
	       lalapps_spininj lalapps_coherent_inspiral \
               lalapps_coherentbank

lalapps_inspiral_SOURCES 	= $(INSPIRALSRC)
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
lalapps_coherentbank_SOURCES = $(COHEREBANK)
lalapps_tmpltbank_SOURCES 	= $(TMPLTBANKSRC)
lalapps_inca_SOURCES 		= $(INCASRC)
lalapps_thinca_SOURCES 		= $(THINCASRC)
lalapps_lt_inca_SOURCES 	= $(LTINCASRC)
lalapps_sire_SOURCES 		= $(SIRESRC)
lalapps_inspinj_SOURCES 	= $(INSPINJSRC)
lalapps_minj_SOURCES 		= $(MINJSRC)
lalapps_bbhinj_SOURCES 		= $(BBHINJSRC)
lalapps_splitbank_SOURCES 	= $(SPLITBANKSRC)
lalapps_tama2ligolw_SOURCES 	= $(TAMA2LIGOLWSRC)
lalapps_trigbank_SOURCES 	= $(TRIGBANKSRC)
lalapps_spininj_SOURCES 	= $(SPININJSRC)
lalapps_inspfrinj_SOURCES 	= $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES 	= $(GWF2XMLSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2 \
              lalapps_inspiral_pipe_tama lalapps_inspiral_pipe_geo \
              lalapps_inspiral_hipe lalapps_inspiral_online_pipe

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
lalapps_inspiral_online_pipe: inspiral_online_pipe.in python_config.sed
	sed -f python_config.sed inspiral_online_pipe.in > lalapps_inspiral_online_pipe
	chmod +x lalapps_inspiral_online_pipe

pkgdata_DATA = BNSMasses.dat inspsrcs.dat inspiral_pipe.ini s2_dag.ini s2_tama.ini inspiral_hipe.ini online.ini
EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(SIRESRC) $(MINJSRC) \
             $(INCASRC) inspiral_pipe.in s2_dag.ini \
             $(COHERESRC) inspiral_pipe_s2.in inspiral_pipe.ini \
             $(INSPINJSRC) BNSMasses.dat inspsrcs.dat \
             $(SPLITBANKSRC) inspiral.py inspiral_pipe_tama.in \
	     inspiral_pipe_geo.in inspiral_online_pipe.in \
	     inspiral_hipe.in s2_tama.ini inspiral_hipe.ini online.ini

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe lalapps_inspiral_pipe_s2 \
             lalapps_inspiral_pipe_tama lalapps_inspiral_pipe_geo \
             lalapps_inspiral_hipe lalapps_inspiral_online_pipe \
             python_config.sed

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
