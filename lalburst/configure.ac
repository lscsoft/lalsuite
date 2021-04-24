AC_PREREQ([2.63])
AC_INIT([LALBurst],[1.5.8.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([lib/config.h lib/LALBurstConfig.h])
AC_CONFIG_SRCDIR([lib/EPSearch.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])

# set library api info
# 1. start with version information of "0:0:0" for each libtool library
# 2. update the version information only immediately before a public
#    release of your software. More frequent updates are unnecessary,
#    and only guarantee that the current interface number gets larger
#    faster.
# 3. if the library source code has changed at all since the last
#    update, then increment revision ("c:r:a" becomes "c:r+1:a").
# 4. if any interfaces have been added, removed, or changed since the
#    last update, increment current, and set revision to 0.
# 5. if any interfaces have been added since the last public release,
#    then increment age.
# 6. if any interfaces have been removed since the last public release,
#    then set age to 0.
AC_SUBST([LIBCURRENT],[6])
AC_SUBST([LIBREVISION],[3])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])
AS_CASE([${host_os}],
  [freebsd-aout|qnx|sco|sunos], [AC_SUBST([LIBMAJOR],[${LIBCURRENT}])],
  [AC_SUBST([LIBMAJOR],[$((LIBCURRENT-LIBAGE))])]
)

AC_CONFIG_FILES([ \
  Makefile \
  conda/Makefile \
  conda/meta.yaml.in \
  debian/Makefile \
  debian/control \
  debian/liblalburst${LIBMAJOR}.install:debian/liblalburst.install.in \
  doxygen/Makefile \
  doxygen/doxygen.cfg \
  include/Makefile \
  include/lal/Makefile \
  lalburst.pc \
  lalburst.spec \
  python/Makefile \
  python/lalburst/Makefile \
  python/lalburst/__init__.py \
  lib/LALBurstVCSInfo.c \
  lib/LALBurstVCSInfo.h \
  lib/Makefile \
  swig/Makefile \
  bin/Makefile \
  bin/version.c \
  test/Makefile \
  test/python/Makefile \
  test/python/test_python_imports.py \
  test/python/test_python_scripts.py \
])
AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-pax \
  dist-xz \
  no-dist-gzip \
])

# version requirements
MIN_LAL_VERSION="7.1.0"
MIN_LALMETAIO_VERSION="2.0.0"
MIN_LALSIMULATION_VERSION="2.5.0"
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_LALMETAIO_VERSION])
AC_SUBST([MIN_LALSIMULATION_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# enable Doxygen documentation
LALSUITE_ENABLE_DOXYGEN

# enable help2man documentation
LALSUITE_ENABLE_HELP2MAN

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([2.6])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for Python
LALSUITE_CHECK_PYTHON([2.6])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
LALSUITE_ADD_FLAGS([C],[${METAIO_CFLAGS}],[${METAIO_LIBS}])
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])

# check for libmetaio headers
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_LIB([LALMetaIO],[${MIN_LALMETAIO_VERSION}],[XLALLIGOLwHasTable],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_LIB([LALSimulation],[${MIN_LALSIMULATION_VERSION}],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALBurst library
AC_SUBST([LDADD], ["\$(top_builddir)/lib/liblalburst.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-DPKG_DATA_DIR='\"\$(pkgdatadir)/\"' ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

# force headers to be installed in ${includedir}/lal, ignoring library name
AC_SUBST([pkgincludedir], ["\${includedir}/lal"])

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# check whether we have the packages to pass the snglcoinc_verify.py tests
# see https://git.ligo.org/lscsoft/lalsuite/issues/208
AM_COND_IF(
  [HAVE_PYTHON],[
    AC_MSG_CHECKING(expected outcome of snglcoinc_verify.py)
    AM_CONDITIONAL(
      [SNGLCOINC_XPASS],
      [${PYTHON} -c "import scipy; import ligo.lw; assert scipy.__version__ > '0.18.999'" 2> /dev/null],
    )
    AS_IF(
      [test "x${SNGLCOINC_XPASS}" = xyes],
      [AC_MSG_RESULT(pass)],
      [AC_MSG_RESULT(fail)]
    )
  ],
  [AM_CONDITIONAL([SNGLCOINC_XPASS], [false])]
)


# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALBurst has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALBurst,
and run 'make install' to install LALBurst.
==================================================
])
