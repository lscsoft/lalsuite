include $(top_srcdir)/gnuscripts/topdir_common.am

ACLOCAL_AMFLAGS = -I gnuscripts

git_generated_sources += \
	python/git_version.py \
	src/lal/LALXMLVCSInfo.h \
	$(EOL)

noswig_pkginclude_HEADERS += \
	src/lal/LALXMLVCSInfo.h \
	$(EOL)

pkginclude_HEADERS += \
	src/lal/LALXML.h \
	src/lal/LALXMLBuildInfo.h \
	src/lal/LALXMLConfig.h \
	src/lal/LALXMLLibVCSInfo.h \
	src/lal/LALXMLVOTableCommon.h \
	src/lal/LALXMLVOTableSerializers.h \
	src/lal/lalxmlswig_test.i \
	$(EOL)

pkgconfig_DATA += lalxml.pc
noinst_DATA += lalxml-uninstalled.pc

lib_LTLIBRARIES += src/liblalxml.la

src_liblalxml_la_SOURCES = \
	src/LALXML.c \
	src/LALXMLVCSInfo.c \
	src/LALXMLVOTableCommon.c \
	src/LALXMLVOTableSerializers.c \
	$(EOL)

nodist_src_liblalxml_la_SOURCES = \
	src/build_info.c \
	$(EOL)

src_liblalxml_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

LDADD += $(top_builddir)/src/liblalxml.la

pkgpython_PYTHON += \
	python/__init__.py \
	python/git_version.py \
	$(EOL)

TESTS += \
	test/LALXMLTest \
	$(EOL)

check_PROGRAMS += \
	$(TESTS)

pkgdata_DATA += \
	test/VOTable-1.1.xsd

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
	debian/lalxml-dev.dirs \
	debian/lalxml-dev.install \
	debian/lalxml-octave.install \
	debian/lalxml-python.install \
	debian/lalxml.dirs \
	debian/lalxml.install \
	debian/rules \
	debian/source/format \
	lalxml.spec \
	$(EOL)

include $(top_srcdir)/gnuscripts/swig_common.am
