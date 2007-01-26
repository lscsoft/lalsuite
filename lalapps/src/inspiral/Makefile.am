## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC 	= inspiral.c
TMPLTBANKSRC 	= tmpltbank.c inspiralutils.c inspiral.h
COHERESRC 	= coherent_inspiral.c
COHEREBANK 	= coherentbank.c
INCASRC 	= inca.c
THINCASRC 	= thinca.c
SIRESRC 	= sire.c
COIRESRC 	= coire.c
SIVASRC 	= siva.c
INSPINJSRC 	= inspinj.c
MINJSRC 	= minj.c
BBHINJSRC 	= bbhinj.c
SPLITBANKSRC 	= splitbank.c
TRIGBANKSRC 	= trigbank.c
INSPFRINJSRC 	= inspfrinj.c
GWF2XMLSRC 	= gwf2xml.c
SPININJSRC 	= spininj.c
COINEXTSRC 	= coinext.c
RANDOMBANKSRC 	= randombank.c
if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank \
               lalapps_inca lalapps_sire lalapps_coire lalapps_siva \
               lalapps_inspinj lalapps_splitbank \
               lalapps_minj lalapps_bbhinj \
               lalapps_inspfrinj lalapps_gwf2xml \
               lalapps_thinca lalapps_trigbank \
               lalapps_spininj lalapps_coherent_inspiral \
               lalapps_coherentbank lalapps_coinext \
               lalapps_randombank 

lalapps_inspiral_SOURCES 	= $(INSPIRALSRC)
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
lalapps_coherentbank_SOURCES 	= $(COHEREBANK)
lalapps_tmpltbank_SOURCES 	= $(TMPLTBANKSRC)
lalapps_inca_SOURCES 		= $(INCASRC)
lalapps_thinca_SOURCES 		= $(THINCASRC)
lalapps_sire_SOURCES 		= $(SIRESRC)
lalapps_coire_SOURCES 		= $(COIRESRC)
lalapps_siva_SOURCES 		= $(SIVASRC)
lalapps_inspinj_SOURCES 	= $(INSPINJSRC)
lalapps_minj_SOURCES 		= $(MINJSRC)
lalapps_bbhinj_SOURCES 		= $(BBHINJSRC)
lalapps_splitbank_SOURCES 	= $(SPLITBANKSRC)
lalapps_trigbank_SOURCES 	= $(TRIGBANKSRC)
lalapps_spininj_SOURCES 	= $(SPININJSRC)
lalapps_inspfrinj_SOURCES 	= $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES 	= $(GWF2XMLSRC)
lalapps_coinext_SOURCES         = $(COINEXTSRC)
lalapps_randombank_SOURCES 	= $(RANDOMBANKSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif

if HAVE_PYTHON
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_hipe \
              lalapps_inspiral_online_pipe \
              lalapps_inspiral_online_pipe_inject \
              lalapps_check_hipe_times
pkgpython_DATA = inspiral.py
endif

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_inspiral_pipe: $(srcdir)/inspiral_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_pipe.in > lalapps_inspiral_pipe 
	chmod +x lalapps_inspiral_pipe
lalapps_inspiral_hipe: $(srcdir)/inspiral_hipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_hipe.in > lalapps_inspiral_hipe
	chmod +x lalapps_inspiral_hipe
lalapps_check_hipe_times: $(srcdir)/check_hipe_times.in python_config.sed
	sed -f python_config.sed $(srcdir)/check_hipe_times.in > lalapps_check_hipe_times
	chmod +x lalapps_check_hipe_times
lalapps_inspiral_online_pipe: $(srcdir)/inspiral_online_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_online_pipe.in > lalapps_inspiral_online_pipe
	chmod +x lalapps_inspiral_online_pipe
lalapps_inspiral_online_pipe_inject: $(srcdir)/inspiral_online_pipe_inject.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_online_pipe_inject.in > lalapps_inspiral_online_pipe_inject
	chmod +x lalapps_inspiral_online_pipe_inject

pkgdata_DATA = BNSMasses.dat inspsrcs.dat inspsrcs.errors inspsrcs-050Mpc.errors inspsrcs-075Mpc.errors inspsrcs-100Mpc.errors inspiral_pipe.ini inspiral_hipe.ini online.ini online_inject_bbh.ini
EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(SIRESRC) $(COIRESRC) $(SIVASRC) $(MINJSRC) \
             $(INCASRC) $(COHERESRC) $(INSPINJSRC) \
             $(SPLITBANKSRC) inspiral.py inspiral_pipe.in \
             inspiral_online_pipe.in inspiral_hipe.in \
             inspiral_online_pipe_inject.in check_hipe_times.in \
             $(DATA)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = lalapps_inspiral_pipe lalapps_inspiral_hipe \
             lalapps_inspiral_online_pipe \
             lalapps_inspiral_online_pipe_inject python_config.sed \
             lalapps_check_hipe_times

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
