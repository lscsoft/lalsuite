export DEB_BUILD_MAINT_OPTIONS = hardening=+all
DPKG_EXPORT_BUILDFLAGS = 1
include /usr/share/dpkg/buildflags.mk
include /usr/share/dpkg/pkg-info.mk

HAVE_PYTHON2 := $(if $(shell pyversions -r),yes)
HAVE_PYTHON3 := $(if $(shell py3versions -r),yes)
HAVE_OCTAVE := $(if $(shell grep $(DEB_SOURCE)-octave debian/control),yes)
DEB_HOST_MULTIARCH ?= $(shell dpkg-architecture -qDEB_HOST_MULTIARCH)
MAKEARGS = \
	-C {build_dir} V=1 \
	pythondir={install_dir} pyexecdir={install_dir} \
	pkgpythondir={install_dir}/$(DEB_SOURCE) \
	pkgpyexecdir={install_dir}/$(DEB_SOURCE)
export PYBUILD_SYSTEM = custom
export PYBUILD_CLEAN_ARGS = $(MAKE) $(MAKEARGS) clean || true
export PYBUILD_CONFIGURE_ARGS = cd {build_dir} && {dir}/configure \
	--prefix=/usr --sysconfdir=/etc \
	--libdir=/usr/lib/$(DEB_HOST_MULTIARCH) \
	--mandir=/usr/share/man --infodir=/usr/share/info \
	--disable-gcc-flags \
	$(CONFIGUREARGS) PYTHON=$$(which {interpreter}) pythondir={install_dir}
export PYBUILD_BUILD_ARGS = $(MAKE) $(MAKEARGS)
export PYBUILD_INSTALL_ARGS = $(MAKE) $(MAKEARGS) DESTDIR={destdir} install
export PYBUILD_TEST_ARGS = $(MAKE) $(MAKEARGS) check

%:
	dh $@ \
	$(if $(HAVE_PYTHON2),--with=python2) $(if $(HAVE_PYTHON3),--with=python3) \
	--buildsystem=pybuild

override_dh_fixperms:
	dh_fixperms $(if $(HAVE_OCTAVE),&& find debian -name '*.oct' | xargs chmod -x)

override_dh_strip:
	dh_strip $(if $(HAVE_OCTAVE),&& find debian -name '*.oct' | xargs strip --strip-unneeded)

override_dh_shlibdeps:
	dh_shlibdeps \
	&& find debian -name '*.la' -delete \
	$(if $(HAVE_OCTAVE),&& dpkg-shlibdeps -Odebian/$(DEB_SOURCE)-octave.substvars $$(find debian -name '*.oct')) \
	$(if $(HAVE_PYTHON2),&& dh_numpy) \
	$(if $(HAVE_PYTHON3),&& dh_numpy3)
