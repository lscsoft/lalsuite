include $(top_srcdir)/gnuscripts/topdir_common.am

ACLOCAL_AMFLAGS = -I gnuscripts

git_generated_sources += \
	python/git_version.py \
	src/lal/LALStochasticVCSInfo.h \
	$(EOL)

noswig_pkginclude_HEADERS += \
	src/lal/LALStochasticVCSInfo.h \
	$(EOL)

pkginclude_HEADERS += \
	src/lal/AstroOmega.h \
	src/lal/LALStochasticBuildInfo.h \
	src/lal/LALStochasticConfig.h \
	src/lal/LALStochasticLibVCSInfo.h \
	src/lal/LIGOLwXMLStochasticRead.h \
	src/lal/SimulatePopcorn.h \
	src/lal/SimulateSB.h \
	src/lal/StochasticCrossCorrelation.h \
	src/lal/StochasticMC.h \
	src/lal/lalstochasticswig_test.i \
	$(EOL)

pkgconfig_DATA += lalstochastic.pc
noinst_DATA += lalstochastic-uninstalled.pc

lib_LTLIBRARIES += src/liblalstochastic.la

src_liblalstochastic_la_SOURCES = \
	src/AstroOmega.c \
	src/LALStochasticVCSInfo.c \
	src/LIGOLwXMLStochasticRead.c \
	src/OverlapReductionFunction.c \
	src/SimulatePopcorn.c \
	src/SimulateSB.c \
	src/StochasticCrossCorrelation.c \
	src/StochasticInverseNoise.c \
	src/StochasticMC.c \
	src/StochasticOmegaGW.c \
	src/StochasticOptimalFilter.c \
	src/StochasticOptimalFilterNormalization.c \
	src/ZeroPadAndFFT.c \
	$(EOL)

nodist_src_liblalstochastic_la_SOURCES = \
	src/build_info.c \
	$(EOL)

src_liblalstochastic_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

LDADD += $(top_builddir)/src/liblalstochastic.la

pkgpython_PYTHON += \
	python/__init__.py \
	python/git_version.py \
	$(EOL)

TESTS += \
	test/AstroOmegaTest \
	test/CZeroPadAndFFTTest \
	test/OverlapReductionFunctionTest \
	test/SZeroPadAndFFTTest \
	test/SimulatePopcornTest \
	test/SimulateSBTest \
	test/StochasticCrossCorrelationSpectrumTest \
	test/StochasticCrossCorrelationStatisticTest \
	test/StochasticHeterodynedCrossCorrelationStatisticTest \
	test/StochasticInverseNoiseTest \
	test/StochasticOmegaGWTest \
	test/StochasticOptimalFilterNormalizationTest \
	test/StochasticOptimalFilterTest \
	$(EOL)

check_PROGRAMS += \
	$(TESTS)

MOSTLYCLEANFILES += \
	test/??_*.dat \
	test/PrintVector.0* \
	test/WHITENED-SB* \
	$(EOL)

DISTCLEANFILES += \
	config.lt \
	src/config.h \
	$(EOL)

EXTRA_DIST += \
	debian/changelog \
	debian/compat \
	debian/control \
	debian/copyright \
	debian/docs \
	debian/lalstochastic-dev.dirs \
	debian/lalstochastic-dev.install \
	debian/lalstochastic-octave.install \
	debian/lalstochastic-python.install \
	debian/lalstochastic.dirs \
	debian/lalstochastic.install \
	debian/rules \
	debian/source/format \
	lalstochastic.spec \
	test/CheckStatus.c \
	test/CheckStatus.h \
	test/H2.dat \
	test/L1.dat \
	test/pipeline \
	test/response1.in \
	test/response2.in \
	test/spectrum1.in \
	test/spectrum2.in \
	$(EOL)

include $(top_srcdir)/gnuscripts/swig_common.am
