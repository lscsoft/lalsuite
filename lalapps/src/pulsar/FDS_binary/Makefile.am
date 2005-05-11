## Process this file with automake to produce Makefile.in
# $Id$

isolated_DIR = ../FDS_isolated

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(isolated_DIR)
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## build everything completely static
## not useable before automake-1.6
#AM_LDFLAGS = -all-static


isolated_src = $(isolated_DIR)/rngmed.c $(isolated_DIR)/DopplerScan.c $(isolated_DIR)/DopplerScan.h 

bin_PROGRAMS = lalapps_ComputeFStatisticBinary lalapps_GenerateBinaryMesh \
	lalapps_MakeSubMesh lalapps_GenerateRandomBinaryInput \
	lalapps_GenerateBinaryInput lalapps_GenerateSearchInput	\
	lalapps_CalculateSensitivity lalapps_FindCoincidence \
	lalapps_ComputeUL

EXTRA_DIST = clusters.h

lalapps_ComputeFStatisticBinary_SOURCES = ComputeFStatisticBinary_v2.c ReadWriteBinaryHeader_v1.c \
	ComputeFStatisticBinary_v2.h clusters.c ${isolated_src} ReadSourceFile_v1.c \
	GenerateBinaryMesh_v1.h

lalapps_GenerateBinaryMesh_SOURCES = GenerateBinaryMesh_v1.c OrbitalParameterConversion_v1.c \
	ReadWriteBinaryHeader_v1.c PeriapseShift_v1.c ReadSourceFile_v1.c \
	GenerateBinaryMesh_v1.h ReadSourceFile_v1.h PeriapseShift_v1.h

lalapps_MakeSubMesh_SOURCES = MakeSubMesh_v1.c ReadWriteBinaryHeader_v1.c \
	OrbitalParameterConversion_v1.c GenerateBinaryMesh_v1.h

lalapps_GenerateRandomBinaryInput_SOURCES = GenerateRandomBinaryInput_v1.c \
	ReadWriteBinaryHeader_v1.c

lalapps_GenerateBinaryInput_SOURCES = GenerateBinaryInput_v1.c ReadSourceFile_v1.c

lalapps_GenerateSearchInput_SOURCES = GenerateSearchInput_v1.c ReadWriteBinaryHeader_v1.c \
	GenerateBinaryMesh_v1.h

lalapps_CalculateSensitivity_SOURCES = CalculateSensitivity_v1.c ReadSourceFile_v1.c \
	CalculateSensitivity_v1.h

lalapps_FindCoincidence_SOURCES = FindCoincidence_v1.c ReadSourceFile_v1.c \
	ReadWriteBinaryHeader_v1.c OrbitalParameterConversion_v1.c FindCoincidence_v1.h \
	ReadSourceFile_v1.h GenerateBinaryMesh_v1.h

lalapps_ComputeUL_SOURCES = ComputeUL_v1.c ReadWriteBinaryHeader_v1.c ComputeUL_v1.h
