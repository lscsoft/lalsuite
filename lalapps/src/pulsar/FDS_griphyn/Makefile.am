## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_narrowBandExtract
else
FRAMEPROGRAMS = 
endif
bin_PROGRAMS = $(FRAMEPROGRAMS)
bin_SCRIPTS = lalapps_GatherLFN lalapps_QueryMetadataLFN
if FRAME
lalapps_narrowBandExtract_SOURCES = narrowBandExtract.c
endif

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

EXTRA_DIST = extractSFTband.c GatherLFN.py QueryMetadataLFN.py
CLEANFILES = lalapps_GatherLFN lalapps_QueryMetadataLFN

lalapps_GatherLFN: GatherLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" GatherLFN.py > lalapps_GatherLFN
	chmod +x lalapps_GatherLFN

lalapps_QueryMetadataLFN: QueryMetadataLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" QueryMetadataLFN.py > lalapps_QueryMetadataLFN
	chmod +x lalapps_QueryMetadataLFN
