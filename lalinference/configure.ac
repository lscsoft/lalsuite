AC_PREREQ([2.63])
AC_INIT([LALInference],[2.0.7.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([lib/config.h lib/LALInferenceConfig.h])
AC_CONFIG_SRCDIR([lib/LALInference.h])
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
AC_SUBST([LIBCURRENT],[21])
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
  debian/liblalinference${LIBMAJOR}.install:debian/liblalinference.install.in \
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
  python/lalinference/imrtgr/Makefile \
  python/lalinference/io/Makefile \
  python/lalinference/plot/Makefile \
  python/lalinference/tiger/Makefile \
  lib/LALInferenceVCSInfo.c \
  lib/LALInferenceVCSInfo.h \
  lib/Makefile \
  bin/Makefile \
  bin/mpi/Makefile \
  bin/version.c \
  swig/Makefile \
  test/Makefile \
  test/python/Makefile \
  test/python/test_python_imports.py \
  test/python/test_python_scripts.py \
])
AC_USE_SYSTEM_EXTENSIONS
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
MIN_GLUE_VERSION="1.54.1"
MIN_LAL_VERSION="7.1.0"
MIN_LALMETAIO_VERSION="2.0.0"
MIN_LALFRAME_VERSION="1.5.0"
MIN_LALSIMULATION_VERSION="2.5.0"
MIN_LALBURST_VERSION="1.5.3"
MIN_LALINSPIRAL_VERSION="2.0.0"
MIN_LALPULSAR_VERSION="3.0.0"
AC_SUBST([MIN_GLUE_VERSION])
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_LALMETAIO_VERSION])
AC_SUBST([MIN_LALFRAME_VERSION])
AC_SUBST([MIN_LALSIMULATION_VERSION])
AC_SUBST([MIN_LALBURST_VERSION])
AC_SUBST([MIN_LALINSPIRAL_VERSION])
AC_SUBST([MIN_LALPULSAR_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable MPI compilation
LALSUITE_ENABLE_MPI

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# optional lal libraries
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALBURST

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

# check for MPI compilers
bambimpi=false
if test "x$mpi" = "xtrue"; then
  AX_MPI([mpi=true],[mpi=false])
  AS_IF([test "$bambi" = "true"], [
    MPILIBS_SAVE="${MPILIBS}"
    AC_LANG([C++])
    AX_MPI([bambimpi=true],[])
    AC_LANG([Fortran])
    AX_MPI([],[bambimpi=false])
    AC_LANG([C])
    MPILIBS="${MPILIBS_SAVE}"
  ])
  AC_LANG([C])
fi
LALSUITE_ENABLE_MODULE([MPI])

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

# check for specific functions
AC_FUNC_STRNLEN
AC_CHECK_FUNC([strdup], [], [AC_MSG_ERROR([could not find the strdup function])])

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for OpenMP
LALSUITE_ENABLE_OPENMP

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.15])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_LIB([LALMetaIO],[${MIN_LALMETAIO_VERSION}],[XLALLIGOLwHasTable],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_LIB([LALFrame],[${MIN_LALFRAME_VERSION}],[XLALFrStreamOpen],[lal/LALFrStream.h])
LALSUITE_CHECK_LIB([LALSimulation],[${MIN_LALSIMULATION_VERSION}],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_LIB([LALInspiral],[${MIN_LALINSPIRAL_VERSION}],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_LIB([LALPulsar],[${MIN_LALPULSAR_VERSION}],[XLALComputeFstat],[lal/ComputeFstat.h])
LALSUITE_CHECK_LIB([LALBurst],[${MIN_LALBURST_VERSION}],[XLALEPSearch],[lal/EPSearch.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALInference library
AC_SUBST([LDADD], ["\$(top_builddir)/lib/liblalinference.la"])

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
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* OpenMP acceleration is $OPENMP_ENABLE_VAL
* MPI support is $MPI_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL
* help2man documentation is $HELP2MAN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALInference,
and run 'make install' to install LALInference.
==================================================
])
