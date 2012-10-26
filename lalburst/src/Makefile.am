AM_CPPFLAGS = -I$(top_builddir)/include

lalburstincludedir = $(includedir)/lal
lalburstinclude_HEADERS = \
	CLR.h \
	EPSearch.h \
	GenerateBurst.h \
	LALBurstConfig.h \
	LALBurstLibVCSInfo.h \
	LALBurstVCSInfo.h \
	LIGOLwXMLBurstRead.h \
	LIGOMetadataBurstUtils.h \
	TFTransform.h \
	TSData.h \
	TSSearch.h \
	Thresholds.h \
	TimeFreq.h \
	TrackSearch.h \
	cs_cosmo.h \
	cs_lambda_cosmo.h

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
	cs_lambda_cosmo.c

liblalburst_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalburst-user-env.shell \
	lalburst-user-env.sh \
	lalburst-user-env.csh

sysconf_DATA = \
	lalburst-user-env.sh \
	lalburst-user-env.csh

BUILT_SOURCES = \
	lalburst-user-env.sh \
	lalburst-user-env.csh

lalburst-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALBURST_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

lalburst-user-env.csh lalburst-user-env.sh: lalburst-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALBurstConfig.h.in \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALBurstVCSInfo.h
MOSTLYCLEANFILES += LALBurstVCSInfo.h
LALBurstVCSInfo.h .PHONY : gen-LALBurstVCSInfo.h
gen-LALBurstVCSInfo.h :
	$(GENERATE_VCS_INFO) LALBurstVCSInfo.h $(srcdir)/LALBurstVCSInfo.h.git
endif

@LAL_HEADER_LINK@
