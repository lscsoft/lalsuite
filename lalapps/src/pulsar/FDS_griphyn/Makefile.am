## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_narrowBandExtract
else
FRAMEPROGRAMS = 
endif
bin_PROGRAMS = $(FRAMEPROGRAMS)
bin_SCRIPTS = lalapps_GatherLFN lalapps_QueryMetadataLFN lalapps_ClusterComputeF lalapps_QueryMetadataEphemeris
if FRAME
lalapps_narrowBandExtract_SOURCES = narrowBandExtract.c
endif

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

EXTRA_DIST = extractSFTband.c GatherLFN.py QueryMetadataLFN.py ClusterComputeF.py QueryMetadataEphemeris.py
CLEANFILES = lalapps_GatherLFN lalapps_QueryMetadataLFN lalapps_ClusterComputeF lalapps_QueryMetadataEphemeris

lalapps_GatherLFN: GatherLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" GatherLFN.py > lalapps_GatherLFN
	chmod +x lalapps_GatherLFN

lalapps_QueryMetadataLFN: QueryMetadataLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" QueryMetadataLFN.py > lalapps_QueryMetadataLFN
	chmod +x lalapps_QueryMetadataLFN

lalapps_ClusterComputeF: ClusterComputeF.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" ClusterComputeF.py > lalapps_ClusterComputeF
	chmod +x lalapps_ClusterComputeF

lalapps_QueryMetadataEphemeris: QueryMetadataEphemeris.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" QueryMetadataEphemeris.py > lalapps_QueryMetadataEphemeris
	chmod +x lalapps_QueryMetadataEphemeris
