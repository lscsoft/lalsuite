INCLUDES = -I$(top_builddir)/include

lalburstincludedir = $(includedir)/lal
lalburstinclude_HEADERS = \
	CLR.h \
	EPSearch.h \
	GenerateBurst.h \
	LALBurstConfig.h \
	LALBurstVCSInfo.h \
	LALSimBurst.h \
	LIGOLwXMLBurstRead.h \
	LIGOMetadataBurstUtils.h \
	TFTransform.h \
	TSData.h \
	TSSearch.h \
	Thresholds.h \
	TimeFreq.h \
	TrackSearch.h

lib_LTLIBRARIES = liblalburst.la
liblalburst_la_SOURCES = \
	CleanAll.c \
	CreateTFPlane.c \
	CreateTimeFreqParam.c \
	CreateTimeFreqRep.c \
	DestroyTimeFreqParam.c \
	DestroyTimeFreqRep.c \
	Dwindow.c \
	EPSearch.c \
	FreqSeriesToTFPlane.c \
	GenerateBurst.c \
	HarmonicFinder.c \
	LALBurstVCSInfo.c \
	LALSimBurst.c \
	LIGOLwXMLBurstRead.c \
	RefInterference.c \
	SnglBurstUtils.c \
	TSData.c \
	TfrPswv.c \
	TfrRsp.c \
	TfrSp.c \
	TfrWv.c \
	Thresholds.c \
	TrackSearch.c \
	check_series_macros.h

CLEANFILES = \
	lalburst-user-env.sh \
	lalburst-user-env.csh

sysconf_DATA = \
	lalburst-user-env.sh \
	lalburst-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALBurstVCSInfo.h \
	lalburst-user-env.sh \
	lalburst-user-env.csh

lalburst-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Burst" >$@ ;\
	echo "LALBURST_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALBURST_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalburst-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Burst" >$@ ;\
	echo "setenv LALBURST_PREFIX \"$(prefix)\"" >>$@ ;\
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
	LALBurstConfig.h.in \
	LALBurstVCSInfo.h.in \
	git_version

LALBurstVCSInfo.h: ._vcsinfo

._vcsinfo: git_version
	@echo "creating LALBurstVCSInfo.h"
	@buildd="`pwd`" && cd $(srcdir) && sourced="`pwd`" && cd $$buildd && \
		./git_version --sed LALBurst $$sourced $$buildd

@LAL_HEADER_LINK@
