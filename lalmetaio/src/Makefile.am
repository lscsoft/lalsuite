INCLUDES = -I$(top_builddir)/include

lalmetaioincludedir = $(includedir)/lal
lalmetaioinclude_HEADERS = \
	LALMetaIOConfig.h \
	LALMetaIOVCSInfo.h \
	LIGOLwXML.h \
	LIGOLwXMLArray.h \
	LIGOLwXMLHeaders.h \
	LIGOLwXMLInspiralHeaders.h \
	LIGOLwXMLRead.h \
	LIGOMetadataTables.h \
	LIGOMetadataUtils.h

lib_LTLIBRARIES = liblalmetaio.la
liblalmetaio_la_SOURCES = \
	CreateMetaTableDir.c \
	LALMetaIOVCSInfo.c \
	LIGOLwXML.c \
	LIGOLwXMLArray.c \
	LIGOLwXMLRead.c \
	LIGOMetadataUtils.c

liblalmetaio_la_LDFLAGS = -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALMetaIOConfig.h.in \
	LALMetaIOVCSInfo.h.in \
	generate_vcs_info.py

CLEANFILES = \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

sysconf_DATA = \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALMetaIOVCSInfo.h \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

lalmetaio-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL MetaIO" >$@ ;\
	echo "LALMETAIO_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALMETAIO_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalmetaio-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL MetaIO" >$@ ;\
	echo "setenv LALMETAIO_PREFIX \"$(prefix)\"" >>$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?LD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?DYLD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?PKG_CONFIG_PATH ) then" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig\"" >>$@ ;\
	echo "endif" >>$@

LALMetaIOVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALMetaIOConfig.h
	@echo "creating LALMetaIOVCSInfo.h"
	$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALMetaIO $(srcdir) $$PWD

@LAL_HEADER_LINK@
