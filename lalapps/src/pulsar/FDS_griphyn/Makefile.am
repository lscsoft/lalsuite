## Process this file with automake to produce Makefile.in
# $Id$
if LALFRAME
# this only requires LALFrame because of usage of LALFrameL.h
LALFRAMEPROGRAMS = lalapps_narrowBandExtract
else
LALFRAMEPROGRAMS =
endif
bin_PROGRAMS = $(LALFRAMEPROGRAMS)
bin_SCRIPTS = lalapps_GatherLFN lalapps_QueryMetadataLFN lalapps_ClusterComputeF lalapps_QueryMetadataEphemeris
if LALFRAME
lalapps_narrowBandExtract_SOURCES = narrowBandExtract.c
endif

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

EXTRA_DIST = extractSFTband.c GatherLFN.py QueryMetadataLFN.py ClusterComputeF.py QueryMetadataEphemeris.py
CLEANFILES = lalapps_GatherLFN lalapps_QueryMetadataLFN lalapps_ClusterComputeF lalapps_QueryMetadataEphemeris

lalapps_GatherLFN: $(srcdir)/GatherLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" $(srcdir)/GatherLFN.py > lalapps_GatherLFN
	chmod +x lalapps_GatherLFN

lalapps_QueryMetadataLFN: $(srcdir)/QueryMetadataLFN.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" $(srcdir)/QueryMetadataLFN.py > lalapps_QueryMetadataLFN
	chmod +x lalapps_QueryMetadataLFN

lalapps_ClusterComputeF: $(srcdir)/ClusterComputeF.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" $(srcdir)/ClusterComputeF.py > lalapps_ClusterComputeF
	chmod +x lalapps_ClusterComputeF

lalapps_QueryMetadataEphemeris: $(srcdir)/QueryMetadataEphemeris.py
	sed "s+@PYTHONPROG@+@PYTHON@+g" $(srcdir)/QueryMetadataEphemeris.py > lalapps_QueryMetadataEphemeris
	chmod +x lalapps_QueryMetadataEphemeris
