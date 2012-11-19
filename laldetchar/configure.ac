AC_PREREQ([2.63])
AC_INIT([laldetchar],[0.0.0.1],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h src/LALDetCharConfig.h])
AC_CONFIG_SRCDIR([src/LALDetChar.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  laldetchar.spec \
  laldetchar.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  python/Makefile \
  python/__init__.py \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR
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
AC_SUBST([LIBCURRENT],[0])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# check for git repo
LALSUITE_CHECK_GIT_REPO

# nightly build
LALSUITE_ENABLE_NIGHTLY

# debug
LALSUITE_ENABLE_DEBUG

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# optional lal libraries
LALSUITE_ENABLE_LALFRAME

# put version/configure info in config header
laldetchar_version_major=`echo "$VERSION" | cut -d. -f1`
laldetchar_version_minor=`echo "$VERSION" | cut -d. -f2`
laldetchar_version_micro=`echo "$VERSION" | cut -d. -f3`
laldetchar_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$laldetchar_version_micro" && laldetchar_version_micro=0
test -z "$laldetchar_version_devel" && laldetchar_version_devel=0
laldetchar_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALDETCHAR_VERSION],["$VERSION"],[LALDetChar Version])
AC_DEFINE_UNQUOTED([LALDETCHAR_VERSION_MAJOR],[$laldetchar_version_major],
                   [LALDetChar Version Major Number])
AC_DEFINE_UNQUOTED([LALDETCHAR_VERSION_MINOR],[$laldetchar_version_minor],
                   [LALDetChar Version Minor Number])
AC_DEFINE_UNQUOTED([LALDETCHAR_VERSION_MICRO],[$laldetchar_version_micro],
                   [LALDetChar Version Micro Number])
AC_DEFINE_UNQUOTED([LALDETCHAR_VERSION_DEVEL],[$laldetchar_version_devel],
                   [LALDetChar Version Devel Number])
AC_DEFINE_UNQUOTED([LALDETCHAR_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALDetChar Configure Arguments])
AC_DEFINE_UNQUOTED([LALDETCHAR_CONFIGURE_DATE],["$laldetchar_configure_date"],
                   [LALDetChar Configure Date])

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for C/C++ compilers
LALSUITE_PROG_CC_CXX

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

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.7.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.7.0],[LALOpenDataFile],[lal/FileIO.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLAL LAL])
LALSUITE_SWIG_DEPENDS([lal],[true])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE([LAL_HEADER_LINK])
AC_PROG_MAKE_SET
AC_OUTPUT
echo "
==================================================
LALDetChar has now been successfully configured:

* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALDetChar,
and run 'make install' to install LALDetChar.
==================================================
"
