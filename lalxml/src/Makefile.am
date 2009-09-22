INCLUDES = -I$(top_builddir)/include

lalxmlincludedir = $(includedir)/lal
lalxmlinclude_HEADERS = \
	LALXML.h \
	LALXMLVOTableCommon.h \
	LALXMLVOTableSerializers.h

lib_LTLIBRARIES = liblalxml.la
liblalxml_la_SOURCES = \
	LALXML.c \
	LALXMLVOTableCommon.c \
	LALXMLVOTableSerializers.c
liblalxml_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = lalxml-user-env.sh lalxml-user-env.csh
sysconf_DATA = lalxml-user-env.sh lalxml-user-env.csh
BUILT_SOURCE = lalxml-user-env.sh lalxml-user-env.csh

lalxml-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL XML" >$@ ;\
	echo "LALXML_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALXML_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalxml-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL XML" >$@ ;\
	echo "setenv LALXML_PREFIX \"$(prefix)\"" >>$@ ;\
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

@LAL_HEADER_LINK@
