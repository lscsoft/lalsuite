INCLUDES = -I$(top_builddir)/include
LDADD = liblalsimulation.la

bin_PROGRAMS = \
	lalsim-bh-qnmode \
	lalsim-bh-ringdown \
	lalsim-bh-sphwf \
	lalsim-ligo-noise \
	lalsim-sgwb

lalsimulationincludedir = $(includedir)/lal
lalsimulationinclude_HEADERS = \
	LALSimBlackHoleRingdown.h \
	LALSimBurst.h \
	LALSimIMREOBNRv2.h \
	LALSimIMR.h \
	LALSimInspiral.h \
	LALSimInspiraldEnergyFlux.h \
	LALSimulation.h \
	LALSimInspiralSpinTaylorT4.h \
	LALSimNoise.h \
	LALSimSGWB.h \
	LALSimulationConfig.h \
	LALSimulationVCSInfo.h

noinst_HEADERS = \
	LALSimIMREOBFactorizedFlux.c \
	LALSimIMREOBFactorizedWaveform.c \
	LALSimIMREOBHybridRingdown.c \
	LALSimIMREOBNQCCorrection.c \
	LALSimIMREOBNewtonianMultipole.c \
	LALSimInspiralPNCoefficients.c

lib_LTLIBRARIES = liblalsimulation.la
liblalsimulation_la_SOURCES = \
	LALSimBlackHoleRingdown.c \
	LALSimBurst.c \
	LALSimIMREOBNRv2.c \
	LALSimIMRPhenom.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimInspiralTaylorT2.c \
	LALSimInspiralTaylorT3.c \
	LALSimInspiralTaylorT4.c \
	LALSimInspiralTofV.c \
	LALSimInspiralTofVIntegrand.c \
	LALSimInspiraldEnergyFlux.c \
	LALSimulation.c \
	LALSimInspiralSpinTaylorT4.c \
	LALSimulationVCSInfo.c \
	LALSimInspiralPSpinInspiralRD.c \
	LALSimNoisePSD.c \
	LALSimNoise.c \
	LALSimSGWB.c \
	LALSimSGWBORF.c \
	check_series_macros.h

liblalsimulation_la_LDFLAGS = -version-info $(LIBVERSION)

lalsim_bh_qnmode_SOURCES = bh_qnmode.c
lalsim_bh_sphwf_SOURCES = bh_sphwf.c
lalsim_bh_ringdown_SOURCES = bh_ringdown.c
lalsim_ligo_noise_SOURCES = ligo_noise.c
lalsim_sgwb_SOURCES = sgwb.c

CLEANFILES = \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

sysconf_DATA = \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALSimulationVCSInfo.h \
	lalsimulation-user-env.sh \
	lalsimulation-user-env.csh

lalsimulation-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Simulation" >$@ ;\
	echo "LALSIMULATION_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALSIMULATION_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalsimulation-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Simulation" >$@ ;\
	echo "setenv LALSIMULATION_PREFIX \"$(prefix)\"" >>$@ ;\
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
	LALSimulationConfig.h.in \
	LALSimulationVCSInfo.h.in

LALSimulationVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALSimulationConfig.h
	@echo "creating LALSimulationVCSInfo.h"
	@$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALSimulation $(srcdir) $$PWD

@LAL_HEADER_LINK@
