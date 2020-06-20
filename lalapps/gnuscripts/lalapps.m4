# lalapps.m4 - lalapps specific autoconf macros
#
# serial 18

AC_DEFUN([LALAPPS_ENABLE_CONDOR], [
  AC_ARG_ENABLE(
    [condor],
    AC_HELP_STRING([--enable-condor],[compile for use with condor @<:@default=no@:>@]),
    AS_CASE(["${enableval}"],
      [yes],[condor=true],
      [no],[condor=false],
      AC_MSG_ERROR([bad value ${enableval} for --enable-condor])
    ),
    [condor=false]
  )
  AS_IF([test "x$condor" = "xtrue"],
    AC_DEFINE([LALAPPS_CONDOR],[1],[LALApps is condor compiled])
  )
  AM_CONDITIONAL([CONDOR_ENABLED],[test "x$condor" = "xtrue"])
])

AC_DEFUN([LALAPPS_ENABLE_STATIC_BINARIES], [
  AC_REQUIRE([PKG_PROG_PKG_CONFIG])
  AC_REQUIRE([LALAPPS_ENABLE_CONDOR])
  AC_ARG_ENABLE(
    [static_binaries],
    AS_HELP_STRING([--enable-static-binaries],[build static binaries @<:@default=no, forced on for condor builds@:>@]),
    AS_CASE(["${enableval}"],
      [yes],[static_binaries=true],
      [no],[static_binaries=false],
      AC_MSG_ERROR([bad value ${enableval} for --enable-static-binaries])
    ),
    [static_binaries=false]
  )
  # force on if condor build is enabled
  AS_IF([test "x$static_binaries" != "xtrue" -a "x$condor" = "xtrue"], [
    AC_MSG_WARN([building static binaries (forced by condor)])
    static_binaries=true
  ])
  # the consequences
  AS_IF([test "x$static_binaries" = "xtrue"], [
    AC_DISABLE_SHARED
    AC_ENABLE_STATIC
    AS_IF([${PKG_CONFIG} --static --version >/dev/null 2>&1],[
      PKG_CONFIG="${PKG_CONFIG} --static"
    ],[
      AC_MSG_WARN([${PKG_CONFIG} does not support --static])
    ])
  ])
])

AC_DEFUN([LALAPPS_ENABLE_FFTW],
[AC_ARG_ENABLE(
  [fftw],
  AC_HELP_STRING([--enable-fftw],[compile code that requires FFTW3 library [default=yes]]),
  [ case "${enableval}" in
      yes) fftw=true;;
      no)  fftw=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-fftw) ;;
    esac
  ], [ fftw=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_FRAMEL],
[AC_ARG_ENABLE(
  [framel],
  AC_HELP_STRING([--enable-framel],[compile code that requires FrameL library [default=yes]]),
  [ case "${enableval}" in
      yes) framel=true;;
      no)  framel=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-framel) ;;
    esac
  ], [ framel=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_PSS],
[AC_ARG_ENABLE(
  [pss],
  AC_HELP_STRING([--enable-pss],[compile code that requires pss library [default=no]]),
  [ case "${enableval}" in
      yes) pss=true;;
      no) pss=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-pss) ;;
    esac
  ], [pss=false])
])

AC_DEFUN([LALAPPS_ENABLE_GDS],
[AC_ARG_ENABLE(
  [gds],
  AC_HELP_STRING([--enable-gds],[compile code that requires GSD library [default=no]]),
  [ case "${enableval}" in
      yes) gds=true;;
      no) gds=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-gds) ;;
    esac
  ], [gds=false])
])

AC_DEFUN([LALAPPS_CHECK_QTHREAD],
[AC_MSG_CHECKING([whether LAL has been compiled with Intel MKL and qthread])
AC_TRY_RUN([
#include <lal/LALConfig.h>
#ifdef LAL_QTHREAD
int main( void ) { return 0; }
#else
int main( void ) { return 1; }
#endif
],
AC_MSG_RESULT([yes])
[
if test x$condor != xtrue ; then
  echo "**************************************************************"
  echo "*                                                            *"
  echo "* LAL has been compiled with --enable-intelfft=condor but    *"
  echo "* but --enable-condor has not been specified when compiling  *"
  echo "* LALapps.                                                   *"
  echo "*                                                            *"
  echo "* LALapps must be configured with --condor-compile when      *"
  echo "* building LALapps against a version of LAL compiled with    *"
  echo "* the fake qthread library.                                  *"
  echo "*                                                            *"
  echo "* Reconfigure LALapps with --enable-condor or rebuild and    *"
  echo "* reinstall LAL with either --enable-intelfft=yes or         *"
  echo "* --disable-intelfft to continue.                            *"
  echo "*                                                            *"
  echo "* See the documentation in the LAL fft package information.  *"
  echo "*                                                            *"
  echo "**************************************************************"
  AC_MSG_ERROR(qthread requires condor_compile)
fi
]
,
AC_MSG_RESULT([no]),
AC_MSG_RESULT([unknown]) ) ] )
