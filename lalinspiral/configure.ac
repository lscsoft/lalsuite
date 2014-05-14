AC_PREREQ([2.63])
AC_INIT([LALInspiral],[1.5.3.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALInspiralConfig.h])
AC_CONFIG_SRCDIR([src/LALInspiral.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalinspiral.spec \
  lalinspiral.pc \
  lalinspiral-uninstalled.pc \
  Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  python/Makefile \
  python/__init__.py \
  python/sbank/Makefile \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

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
AC_SUBST([LIBCURRENT],[8])
AC_SUBST([LIBREVISION],[3])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable CUDA support
LALSUITE_WITH_CUDA

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# optional lal libraries
LALSUITE_ENABLE_LALFRAME

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PATH_PYTHON([2.6],,[AC_MSG_ERROR([Python-2.6, or higher, is required])])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
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

# check for frame library
PKG_CHECK_MODULES([LIBFRAME],[libframe],[true],[false])
LALSUITE_ADD_FLAGS([C],[${LIBFRAME_CFLAGS}],[${LIBFRAME_LIBS}])
AC_CHECK_LIB([Frame],[FrLibIni],,[AC_MSG_ERROR([could not find the frame library])])

# check for frame library headers
AC_CHECK_HEADERS([FrameL.h],,[AC_MSG_ERROR([could not find the FrameL.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.12.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.12.0],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_LIB([LALMetaIO],[1.2.0],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_LIB([LALSimulation],[1.1.0],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_LIB([LALFrame],[1.1.0],[LALFrOpen],[lal/LALFrameL.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLAL LAL])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# define if CUDA is enabled
if test "${cuda}" = "true"; then
  AC_DEFINE([LALINSPIRAL_CUDA_ENABLED],[1],[Define if using cuda library])
fi

# force headers to be installed in ${includedir}/lal, ignoring library name
AC_SUBST([pkgincludedir], ["\${includedir}/lal"])

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALInspiral has now been successfully configured:

* CUDA support is $CUDA_ENABLE_VAL
* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALInspiral,
and run 'make install' to install LALInspiral.
==================================================
])
