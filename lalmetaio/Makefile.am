include $(top_srcdir)/gnuscripts/topdir_common.am

ACLOCAL_AMFLAGS = -I gnuscripts

git_generated_sources += \
	python/git_version.py \
	src/lal/LALMetaIOVCSInfo.h \
	$(EOL)

noswig_pkginclude_HEADERS += \
	src/lal/LALMetaIOVCSInfo.h \
	$(EOL)

pkginclude_HEADERS += \
	src/lal/LALMetaIOBuildInfo.h \
	src/lal/LALMetaIOConfig.h \
	src/lal/LALMetaIOLibVCSInfo.h \
	src/lal/LIGOLwXML.h \
	src/lal/LIGOLwXMLArray.h \
	src/lal/LIGOLwXMLHeaders.h \
	src/lal/LIGOLwXMLInspiralHeaders.h \
	src/lal/LIGOLwXMLRead.h \
	src/lal/LIGOMetadataTables.h \
	src/lal/LIGOMetadataUtils.h \
	src/lal/lalmetaioswig_test.i \
	$(EOL)

pkgconfig_DATA += lalmetaio.pc
noinst_DATA += lalmetaio-uninstalled.pc

lib_LTLIBRARIES += src/liblalmetaio.la

src_liblalmetaio_la_SOURCES = \
	src/CreateMetaTableDir.c \
	src/LALMetaIOVCSInfo.c \
	src/LIGOLwXML.c \
	src/LIGOLwXMLArray.c \
	src/LIGOLwXMLRead.c \
	src/LIGOMetadataUtils.c
	$(EOL)

nodist_src_liblalmetaio_la_SOURCES = \
	src/build_info.c \
	$(EOL)

src_liblalmetaio_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

LDADD += $(top_builddir)/src/liblalmetaio.la

pkgpython_PYTHON += \
	python/__init__.py \
	python/git_version.py \
	$(EOL)

test_programs += \
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
	debian/lalmetaio-dev.dirs \
	debian/lalmetaio-dev.install \
	debian/lalmetaio-octave.install \
	debian/lalmetaio-python.install \
	debian/lalmetaio.dirs \
	debian/lalmetaio.install \
	debian/rules \
	debian/source/format \
	lalmetaio.spec \
	$(EOL)

include $(top_srcdir)/gnuscripts/swig_common.am
