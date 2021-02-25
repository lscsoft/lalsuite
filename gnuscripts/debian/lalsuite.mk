export DEB_BUILD_MAINT_OPTIONS = hardening=+all
DPKG_EXPORT_BUILDFLAGS = 1
include /usr/share/dpkg/buildflags.mk
include /usr/share/dpkg/pkg-info.mk

DH_PACKAGES := $(shell dh_listpackages)
WITH_PYTHON2 := $(if $(filter python-$(DEB_SOURCE),$(DH_PACKAGES)),yes)
WITH_OCTAVE := $(if $(filter $(DEB_SOURCE)-octave,$(DH_PACKAGES)),yes)

# define 'main' python, and other python
# NOTE: to disable the 'other' python build set 'PYTHON_OTHER = '
#       (i.e. to the empty string) _before_ the 'include debian/lalsuite.mk'
#       line in debian/rules file
PYTHON := /usr/bin/python3
ifneq (,$(WITH_PYTHON2))
	PYTHON_OTHER ?= /usr/bin/python2
endif
PYTHON_OTHER_BUILD_DIR ?= _buildpyother

# handle parallelism
ifneq (,$(filter parallel=%,$(DEB_BUILD_OPTIONS)))
	NUMJOBS = $(patsubst parallel=%,%,$(filter parallel=%,$(DEB_BUILD_OPTIONS)))
endif

%:
	dh $@ \
		--buildsystem=autoconf \
		--with=python3 \
		$(if $(WITH_PYTHON2),--with=python2) \
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
ifneq ($(strip $(PYTHON_OTHER)),)
	# -- configure and build 'other' python library
	# copy over current build into new directory
	# NOTE: when this is removed, please remove
	#       'rsync' from Build-Depends
	$(RM) -r $(PYTHON_OTHER_BUILD_DIR)
	rsync -ra . $(PYTHON_OTHER_BUILD_DIR)
	# remove old makefiles
	find $(PYTHON_OTHER_BUILD_DIR) -name Makefile -delete
	# remove 'main' python compiled extensions
	$(RM) -r -v \
		$(PYTHON_OTHER_BUILD_DIR)/python/$(DEB_SOURCE)/*.la \
		$(PYTHON_OTHER_BUILD_DIR)/python/$(DEB_SOURCE)/*.lo \
		$(PYTHON_OTHER_BUILD_DIR)/python/$(DEB_SOURCE)/*.o \
		$(PYTHON_OTHER_BUILD_DIR)/python/$(DEB_SOURCE)/.libs
	# configure for the 'other' python
	dh_auto_configure --sourcedirectory $(PYTHON_OTHER_BUILD_DIR) -- \
		$(CONFIGUREARGS) \
		--disable-gcc-flags \
		--disable-doxygen \
		--disable-swig-octave \
		PYTHON=$(PYTHON_OTHER)
	# compile the new SWIG bindings, python modules, and python scripts
	$(MAKE) -j$(NUMJOBS) -C $(PYTHON_OTHER_BUILD_DIR) -C swig
	$(MAKE) -j$(NUMJOBS) -C $(PYTHON_OTHER_BUILD_DIR) -C python/$(DEB_SOURCE)
	$(MAKE) -j$(NUMJOBS) -C $(PYTHON_OTHER_BUILD_DIR) -C bin bin_PROGRAMS="" dist_bin_SCRIPTS=""
endif

override_dh_auto_install:
	dh_auto_install --destdir=debian/tmp
ifneq ($(strip $(PYTHON_OTHER)),)
	# install python other
	$(MAKE) $(MAKEARS) DESTDIR=$(CURDIR)/debian/tmp AM_UPDATE_INFO_DIR=no -C $(PYTHON_OTHER_BUILD_DIR) -C swig install-exec-am
	$(MAKE) $(MAKEARS) DESTDIR=$(CURDIR)/debian/tmp AM_UPDATE_INFO_DIR=no -C $(PYTHON_OTHER_BUILD_DIR) -C python install
endif

override_dh_auto_test:
	dh_auto_test
ifneq ($(strip $(PYTHON_OTHER)),)
	# test other python
	[ -d $(PYTHON_OTHER_BUILD_DIR)/test/python ] && $(MAKE) $(MAKEARS) -j1 VERBOSE=1 -C $(PYTHON_OTHER_BUILD_DIR) -C test/python check
endif

override_dh_fixperms:
	dh_fixperms $(if $(WITH_OCTAVE),&& find debian -name '*.oct' | xargs chmod -x)

override_dh_strip:
	dh_strip $(if $(WITH_OCTAVE),&& find debian -name '*.oct' | xargs strip --strip-unneeded)

override_dh_shlibdeps:
	dh_shlibdeps \
	&& find debian -name '*.la' -delete \
	&& dh_numpy3 \
	$(if $(WITH_OCTAVE),&& dpkg-shlibdeps -Odebian/$(DEB_SOURCE)-octave.substvars $$(find debian -name '*.oct')) \
	$(if $(WITH_PYTHON2),&& dh_numpy) \
	;
