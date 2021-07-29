export DEB_BUILD_MAINT_OPTIONS = hardening=+all
DPKG_EXPORT_BUILDFLAGS = 1
include /usr/share/dpkg/buildflags.mk
include /usr/share/dpkg/pkg-info.mk

DH_PACKAGES := $(shell dh_listpackages)
WITH_OCTAVE := $(if $(filter $(DEB_SOURCE)-octave,$(DH_PACKAGES)),yes)
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
		$(if $(WITH_OCTAVE),,--disable-swig-octave) \
		$(CONFIGUREARGS) \
		PYTHON=$(PYTHON)

override_dh_auto_build:
	# build for the 'main' python version
	dh_auto_build --parallel

override_dh_auto_install:
	dh_auto_install --destdir=debian/tmp

override_dh_auto_test:
	dh_auto_test

override_dh_fixperms:
	dh_fixperms $(if $(WITH_OCTAVE),&& find debian -name '*.oct' | xargs chmod -x)

override_dh_strip:
	dh_strip $(if $(WITH_OCTAVE),&& find debian -name '*.oct' | xargs strip --strip-unneeded)

override_dh_shlibdeps:
	dh_shlibdeps \
	&& find debian -name '*.la' -delete \
	&& dh_numpy3 \
	$(if $(WITH_OCTAVE),&& dpkg-shlibdeps -Odebian/$(DEB_SOURCE)-octave.substvars $$(find debian -name '*.oct')) \
	;
