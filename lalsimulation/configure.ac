AC_PREREQ([2.63])
AC_INIT([LALSimulation],[6.2.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([lib/config.h lib/LALSimulationConfig.h])
AC_CONFIG_SRCDIR([lib/LALSimulation.c])
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
AC_SUBST([LIBCURRENT],[39])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[2])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])
AS_CASE([${host_os}],
  [freebsd-aout|qnx|sco|sunos], [AC_SUBST([LIBMAJOR],[${LIBCURRENT}])],
  [AC_SUBST([LIBMAJOR],[$((LIBCURRENT-LIBAGE))])]
)

AC_CONFIG_FILES([ \
  Makefile \
  conda/Makefile \
  conda/meta.yaml.in \
  conda/migrations/Makefile \
  debian/Makefile \
  debian/control \
  debian/liblalsimulation${LIBMAJOR}.install:debian/liblalsimulation.install.in \
  doxygen/Makefile \
  include/Makefile \
  include/lal/Makefile \
  lalsimulation.pc \
  lalsimulation.spec \
  python/Makefile \
  python/lalsimulation/Makefile \
  python/lalsimulation/__init__.py \
  python/lalsimulation/nrfits/Makefile \
  python/lalsimulation/tilts_at_infinity/Makefile \
  python/lalsimulation/gwsignal/Makefile \
  python/lalsimulation/gwsignal/core/Makefile \
  python/lalsimulation/gwsignal/models/Makefile \
  lib/LALSimulationVCSInfo.c \
  lib/LALSimulationVCSInfo.h \
  lib/Makefile \
  lib/seobnrv3_opt_exactderivs/Makefile \
  lib/mathematica_codes/Makefile \
  bin/Makefile \
  bin/version.c \
  swig/Makefile \
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
MIN_LAL_VERSION="7.7.0"
AC_SUBST([MIN_LAL_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_PLATFORM
LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS
LALSUITE_WITH_FALLBACK_DATA_PATH

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable strict definitions
LALSUITE_ENABLE_STRICT_DEFS

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
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([])   # use default minimum Python version

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_CHECK_HEADERS([unistd.h])

# check for gethostname in unistd.h
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
          [Define if gethostname prototype is in unistd.h]),
          [AC_MSG_RESULT]([no]))

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for OpenMP
LALSUITE_ENABLE_OPENMP

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for python3-embed
PKG_CHECK_MODULES([PYTHON3_EMBED],[python-${PYTHON_VERSION}-embed],[have_python3_embed=true],[have_python3_embed=false])
if test "x$have_python3_embed" = xtrue; then

    # separate out -l and -L LIBS
    AC_SUBST([PYTHON3_EMBED_LIBADD],[])
    for flag in ${PYTHON3_EMBED_LIBS}; do
       case $flag in
          -l*) PYTHON3_EMBED_LIBADD="${PYTHON3_EMBED_LIBADD} ${flag}";;
          *) _PYTHON3_EMBED_LIBS="${_PYTHON3_EMBED_LIBS} ${flag}";;
       esac
    done
    PYTHON3_EMBED_LIBS="${_PYTHON3_EMBED_LIBS}"

    # check for NumPy
    AC_MSG_CHECKING([for NumPy])
    numpy_version=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import numpy
print(numpy.__version__)
EOD`]
    AS_IF([test $? -ne 0],[
        AC_MSG_ERROR([could not import NumPy])
    ])
    AC_MSG_RESULT([yes])

    # additional preprocessor flags for NumPy
    python_out=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import sys
import numpy
sys.stdout.write(' -I' + numpy.get_include())
EOD`]
    AS_IF([test $? -ne 0],[
        AC_MSG_ERROR([could not determine NumPy preprocessor flags])
    ])
    for flag in ${python_out}; do
        PYTHON3_EMBED_CFLAGS="${PYTHON3_EMBED_CFLAGS} ${flag}"
    done

    LALSUITE_ADD_FLAGS([C],[${PYTHON3_EMBED_CFLAGS}],[${PYTHON3_EMBED_LIBS}])

    AC_CHECK_HEADERS([Python.h],[],[
        AC_MSG_ERROR([could not find the Python.h header])
    ],[
        AC_INCLUDES_DEFAULT
    ])
    AC_CHECK_HEADER([numpy/arrayobject.h],[],[
        AC_MSG_ERROR([could not find the header "numpy/arrayobject.h"])
    ],[
        AC_INCLUDES_DEFAULT
        #include <Python.h>
    ])

    AC_DEFINE([HAVE_PYTHON3_EMBED],[1],[Define to 1 to use python3-embed])
fi

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALSimulation library
AC_SUBST([LDADD], ["\$(top_builddir)/lib/liblalsimulation.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

# force headers to be installed in ${includedir}/lal, ignoring library name
AC_SUBST([pkgincludedir], ["\${includedir}/lal"])

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# make sure local library module is in python path
LAL_PYTHON_PATH="\$(abs_top_builddir)/python:${LAL_PYTHON_PATH}"

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALSimulation has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* OpenMP acceleration is $OPENMP_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL
* help2man documentation is $HELP2MAN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALSimulation,
and run 'make install' to install LALSimulation.
==================================================
])
