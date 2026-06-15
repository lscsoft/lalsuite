export DEB_BUILD_MAINT_OPTIONS = hardening=+all
DPKG_EXPORT_BUILDFLAGS = 1
include /usr/share/dpkg/buildflags.mk
include /usr/share/dpkg/pkg-info.mk

DH_PACKAGES := $(shell dh_listpackages)
PYTHON := /usr/bin/python3

# handle parallelism
ifneq (,$(filter parallel=%,$(DEB_BUILD_OPTIONS)))
	NUMJOBS = $(patsubst parallel=%,%,$(filter parallel=%,$(DEB_BUILD_OPTIONS)))
endif

ifneq (,$(LALSUITE_DEB_WITH_CUDA))
# build with CUDA support
# - CFLAGS: -fno-lto is needed for CUDA linking
CONFIGUREARGS_CUDA = --with-cuda=$(LALSUITE_DEB_WITH_CUDA) CFLAGS="-O3 -fno-lto"
# - ignore extra library dependencies from CUDA
DHSHLIBDEPSARGS_CUDA = --dpkg-shlibdeps-params=--ignore-missing-info
endif

%:
	dh $@ \
		--buildsystem=autoconf \
		--with=python3 \
	;

# configure the build for the 'main' python version
override_dh_auto_configure:
	dh_auto_configure -- \
		--disable-gcc-flags \
		--disable-swig-octave \
		$(CONFIGUREARGS) \
		$(CONFIGUREARGS_CUDA) \
		PYTHON=$(PYTHON) \
		;

# build for the 'main' python version
override_dh_auto_build:
	dh_auto_build --parallel

override_dh_auto_install:
	dh_auto_install --destdir=debian/tmp

ifneq (,$(LALSUITE_DEB_WITH_CUDA))
override_dh_auto_test:
	echo Tests skipped for CUDA architecture
else
override_dh_auto_test:
	dh_auto_test
endif

override_dh_missing:
	find debian -name '*.la' -delete
	dh_missing

override_dh_shlibdeps:
	dh_shlibdeps $(DHSHLIBDEPSARGS_CUDA)
	dh_numpy3

# apply debian/lintian-overrides to all binary packages
override_dh_lintian:
	perl -ne 'print "cp -f debian/lintian-overrides debian/$$1.lintian-overrides\n" if /^Package:\s*(\S+)/' debian/control | bash
	dh_lintian
