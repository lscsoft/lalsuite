INCLUDES = -I$(top_builddir)/include

lalframeincludedir = $(includedir)/lal
lalframeinclude_HEADERS = \
	Aggregation.h \
	FrameCache.h \
	FrameCalibration.h \
	FrameData.h \
	FrameStream.h \
	LALFrameConfig.h \
	LALFrameIO.h \
	LALFrameL.h \
	LALFrameVCSInfo.h

lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = \
	Aggregation.c \
	FrameCache.c \
	FrameCalibration.c \
	FrameData.c \
	FrameSeries.c \
	FrameStream.c \
	LALFrameIO.c \
	LALFrameVCSInfo.c

noinst_HEADERS = \
	FrameSeriesRead_source.c \
	FrameSeriesWrite_source.c

liblalframe_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = \
	lalframe-user-env.sh \
	lalframe-user-env.csh

sysconf_DATA = \
	lalframe-user-env.sh \
	lalframe-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALFrameVCSInfo.h \
	lalframe-user-env.sh \
	lalframe-user-env.csh

lalframe-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Frame" >$@ ;\
	echo "LALFRAME_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALFRAME_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalframe-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Frame" >$@ ;\
	echo "setenv LALFRAME_PREFIX \"$(prefix)\"" >>$@ ;\
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

EXTRA_DIST = \
	generate_vcs_info.py \
	LALFrameConfig.h.in \
	LALFrameVCSInfo.h.in

LALFrameVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALFrameConfig.h
	@echo "creating LALFrameVCSInfo.h"
	@$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALFrame $(srcdir) $$PWD

@LAL_HEADER_LINK@
