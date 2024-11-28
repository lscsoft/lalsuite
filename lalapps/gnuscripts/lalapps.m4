# lalapps.m4 - lalapps specific autoconf macros
#
# serial 20

AC_DEFUN([LALAPPS_ENABLE_CONDOR], [
  AC_ARG_ENABLE(
    [condor],
    AS_HELP_STRING([--enable-condor],[compile for use with condor @<:@default=no@:>@]),
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
  AS_HELP_STRING([--enable-fftw],[compile code that requires FFTW3 library [default=yes]]),
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
  AS_HELP_STRING([--enable-framel],[compile code that requires FrameL library [default=yes]]),
  [ case "${enableval}" in
      yes) framel=true;;
      no)  framel=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-framel) ;;
    esac
  ], [ framel=true ] )
])
