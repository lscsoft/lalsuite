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

%:
	dh $@ \
		--buildsystem=autoconf \
		--with=python3 \
	;

override_dh_auto_configure:
	# configure the build for the 'main' python version
	dh_auto_configure -- \
		--disable-gcc-flags \
		--disable-swig-octave \
		$(CONFIGUREARGS) \
		PYTHON=$(PYTHON)

override_dh_auto_build:
	# build for the 'main' python version
	dh_auto_build --parallel

override_dh_auto_install:
	dh_auto_install --destdir=debian/tmp

override_dh_auto_test:
	dh_auto_test

override_dh_shlibdeps:
	dh_shlibdeps \
	&& find debian -name '*.la' -delete \
	&& dh_numpy3 \
	;

override_dh_lintian:
	# apply debian/lintian-overrides to all binary packages
	perl -ne 'print "cp -f debian/lintian-overrides debian/$$1.lintian-overrides\n" if /^Package:\s*(\S+)/' debian/control | bash
	dh_lintian
