## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_SemiCoherentBinary lalapps_HJDtoGPSSSB \
	lalapps_ComputeFStatisticBinary lalapps_GenerateBinaryMesh \
	lalapps_MakeSubMesh lalapps_GenerateRandomBinaryInput \
	lalapps_GenerateBinaryInput lalapps_GenerateSearchInput	\
	lalapps_CalculateSensitivity lalapps_FindCoincidence \
	lalapps_ComputeUL lalapps_CombineBayesianResults

lalapps_SemiCoherentBinary_SOURCES = SemiCoherentBinary.c

lalapps_HJDtoGPSSSB_SOURCES = HJDtoGPSSSB.c

lalapps_CombineBayesianResults_SOURCES = CombineBayesianResults.c

lalapps_ComputeFStatisticBinary_SOURCES = ComputeFStatisticBinary_v2.c ReadWriteBinaryHeader_v1.c \
	ComputeFStatisticBinary_v2.h clusters.c clusters.h ReadSourceFile_v1.c \
	GenerateBinaryMesh_v1.h rngmed.c rngmed.h

lalapps_GenerateBinaryMesh_SOURCES = GenerateBinaryMesh_v1.c OrbitalParameterConversion_v1.c \
	ReadWriteBinaryHeader_v1.c PeriapseShift_v1.c ReadSourceFile_v1.c \
	GenerateBinaryMesh_v1.h ReadSourceFile_v1.h

lalapps_MakeSubMesh_SOURCES = MakeSubMesh_v1.c ReadWriteBinaryHeader_v1.c \
	OrbitalParameterConversion_v1.c GenerateBinaryMesh_v1.h

lalapps_GenerateRandomBinaryInput_SOURCES = GenerateRandomBinaryInput_v1.c \
	GenerateRandomBinaryInput_v1.h ReadWriteBinaryHeader_v1.c

lalapps_GenerateBinaryInput_SOURCES = GenerateBinaryInput_v1.c ReadSourceFile_v1.c

lalapps_GenerateSearchInput_SOURCES = GenerateSearchInput_v1.c ReadWriteBinaryHeader_v1.c \
	GenerateBinaryMesh_v1.h

lalapps_CalculateSensitivity_SOURCES = CalculateSensitivity_v1.c ReadSourceFile_v1.c \
	CalculateSensitivity_v1.h

lalapps_FindCoincidence_SOURCES = FindCoincidence_v1.c ReadSourceFile_v1.c \
	ReadWriteBinaryHeader_v1.c OrbitalParameterConversion_v1.c FindCoincidence_v1.h \
	ReadSourceFile_v1.h GenerateBinaryMesh_v1.h

lalapps_ComputeUL_SOURCES = ComputeUL_v1.c ReadWriteBinaryHeader_v1.c ComputeUL_v1.h



## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
