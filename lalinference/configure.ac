AC_PREREQ([2.63])
AC_INIT([LALInference],[1.9.2.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALInferenceConfig.h])
AC_CONFIG_SRCDIR([src/LALInference.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  Makefile \
  debian/Makefile \
  debian/control \
  doxygen/Makefile \
  doxygen/doxygen.cfg \
  include/Makefile \
  include/lal/Makefile \
  lalinference.pc \
  lalinference.spec \
  python/Makefile \
  python/lalinference/Makefile \
  python/lalinference/__init__.py \
  python/lalinference/bayestar/Makefile \
  python/lalinference/io/Makefile \
  python/lalinference/rapid_pe/Makefile \
  python/lalinference/imrtgr/Makefile \
  python/lalinference/tiger/Makefile \
  python/lalinference/plot/Makefile \
  python/lalinference/popprior/Makefile \
  src/LALInferenceVCSInfo.c \
  src/LALInferenceVCSInfo.h \
  src/Makefile \
  src/version.c \
  swig/Makefile \
  test/Makefile \
  test/baseline/Makefile
])
AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-ustar \
  dist-xz \
  no-dist-gzip \
])

# version requirements
MIN_GLUE_VERSION="1.54.1"
MIN_LAL_VERSION="6.18.0"
MIN_LALMETAIO_VERSION="1.3.1"
MIN_LALFRAME_VERSION="1.4.3"
MIN_LALSIMULATION_VERSION="1.7.0"
MIN_LALINSPIRAL_VERSION="1.7.7"
MIN_LALPULSAR_VERSION="1.16.0"
MIN_LALBURST_VERSION="1.4.4"
MIN_LALXML_VERSION="1.2.4"
AC_SUBST([MIN_GLUE_VERSION])
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_LALMETAIO_VERSION])
AC_SUBST([MIN_LALFRAME_VERSION])
AC_SUBST([MIN_LALSIMULATION_VERSION])
AC_SUBST([MIN_LALINSPIRAL_VERSION])
AC_SUBST([MIN_LALPULSAR_VERSION])
AC_SUBST([MIN_LALBURST_VERSION])
AC_SUBST([MIN_LALXML_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

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
AC_SUBST([LIBCURRENT],[13])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# optional lal libraries
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALXML

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# enable Doxygen documentation
LALSUITE_ENABLE_DOXYGEN

# enable help2man documentation
LALSUITE_ENABLE_HELP2MAN

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
LALSUITE_CHECK_PYTHON([2.6])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([unistd.h])

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for OpenMP
LALSUITE_ENABLE_OPENMP

# check for HEALPix
PKG_CHECK_MODULES([CHEALPIX],[chealpix >= 3.30.0],[true],[false])
LALSUITE_ADD_FLAGS([C],[${CHEALPIX_CFLAGS}],[${CHEALPIX_LIBS}])
HEALPIX_ENABLE_VAL=ENABLED
AC_CHECK_LIB([chealpix],[nside2npix],,[HEALPIX_ENABLE_VAL=DISABLED])
AM_CONDITIONAL([HAVE_CHEALPIX],[test "x$HEALPIX_ENABLE_VAL" = "xENABLED"])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.15])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
LALSUITE_ADD_FLAGS([C],[${METAIO_CFLAGS}],[${METAIO_LIBS}])
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])

# check for libmetaio headers
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

# check for hdf5 support
HDF5_ENABLE_VAL="DISABLED"
# Note that AX_LIB_HDF5 uses h5cc to determine compiler for configure tests.
# Make sure that it uses the same compiler as we are using.
# Otherwise the configure tests may fail if we are using icc.
HDF5_CC="${CC}"; export HDF5_CC
HDF5_CLINKER="${CC}"; export HDF5_CLINKER
AX_LIB_HDF5()
unset HDF5_CC
unset HDF5_CLINKER
if test "${with_hdf5}" != "yes"; then
  AC_MSG_ERROR([could not find the HDF5 library])
fi

if test "$lalxml" = "true"; then
  # check for libxml2
  PKG_CHECK_MODULES([LIBXML2],[libxml-2.0 >= 2.6],[true],[false])
  LALSUITE_ADD_FLAGS([C],[${LIBXML2_CFLAGS}],[${LIBXML2_LIBS}])
  AC_CHECK_LIB([xml2],[xmlInitParser],,[AC_MSG_WARN([cannot find the xml2 library])])
fi

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_LIB([LALMetaIO],[${MIN_LALMETAIO_VERSION}],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_LIB([LALFrame],[${MIN_LALFRAME_VERSION}],[LALFrOpen],[lal/LALFrameL.h])
LALSUITE_CHECK_LIB([LALSimulation],[${MIN_LALSIMULATION_VERSION}],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_LIB([LALInspiral],[${MIN_LALINSPIRAL_VERSION}],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_LIB([LALPulsar],[${MIN_LALPULSAR_VERSION}],[XLALComputeFstat],[lal/ComputeFstat.h])
LALSUITE_CHECK_LIB([LALBurst],[${MIN_LALBURST_VERSION}],[XLALEPSearch],[lal/EPSearch.h])


# check for lalxml
LALSUITE_CHECK_OPT_LIB([LALXML],[${MIN_LALXML_VERSION}],[XLALXMLFilePrintElements],[lal/LALXML.h])

dnl FIXME: This doesn't work in the top-level build, because it tries
dnl to compile and link against LAL itself, which will not have been built yet.
dnl
dnl AC_MSG_CHECKING([if LAL was built with HDF5 support])
dnl AC_RUN_IFELSE(
dnl   [AC_LANG_PROGRAM(
dnl     [
dnl       #include <lal/H5FileIO.h>
dnl       #include <lal/XLALError.h>
dnl     ],
dnl     [
dnl       int ret;
dnl       (void)XLALSetSilentErrorHandler();
dnl       (void)XLALH5FileOpen(NULL, NULL);
dnl       return xlalErrno != XLAL_EFAULT;
dnl     ]
dnl   )],
dnl   [AC_MSG_RESULT([yes])],
dnl   [
dnl     AC_MSG_RESULT([no])
dnl     AC_MSG_ERROR([LAL must be built with HDF5 support in order to build LALInference])
dnl   ]
dnl )

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALInference library
AC_SUBST([LDADD], ["\$(top_builddir)/src/liblalinference.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-DPKG_DATA_DIR='\"\$(pkgdatadir)/\"' ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

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
LALInference has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* LALXML library support is $LALXML_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* HEALPix is $HEALPIX_ENABLE_VAL
* OpenMP acceleration is $OPENMP_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL
* help2man documentation is $HELP2MAN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALInference,
and run 'make install' to install LALInference.
==================================================
])
