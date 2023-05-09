# -*- mode: autoconf; -*-
# lalpulsar.m4 - LALPulsar specific macros
#
# serial 7

AC_DEFUN([LALPULSAR_ENABLE_SISTR],
[AC_ARG_ENABLE(
  [sistr],
  AS_HELP_STRING([--enable-sistr],[compile code that requires SIStr library [default=no]]),
  [ case "${enableval}" in
      yes) sistr=true;;
      no) sistr=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-sistr) ;;
    esac
  ], [sistr=false])
])

AC_DEFUN([LALPULSAR_CHECK_ALTIVEC],[
  # $0: check for Altivec support for ComputeFstat Demod hotloop
  have_altivec=no
  AS_CASE([$host_cpu],
    [powerpc*],[
      LALSUITE_PUSH_UVARS
      LALSUITE_CLEAR_UVARS
      AC_LANG_PUSH([C])
      CFLAGS="${lalsuite_uvar_CFLAGS}"
      AC_MSG_CHECKING([whether ]_AC_LANG[ compiler defines __ALTIVEC__ with CFLAGS=${CFLAGS}])
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([],[[
#if !defined(__ALTIVEC__)
#error Preprocessor macro not defined by compiler
#endif
]])
      ],[
        have_altivec=yes
      ])
      AC_MSG_RESULT([${have_altivec}])
      AC_LANG_POP([C])
      LALSUITE_POP_UVARS
    ]
  )
  AM_CONDITIONAL([HAVE_ALTIVEC],[test x"${have_altivec}" = xyes])
  AM_COND_IF([HAVE_ALTIVEC],[
    AC_DEFINE([HAVE_ALTIVEC],[1],[Define to 1 for Altivec support])
  ])
  # end $0
])

AC_DEFUN([LALPULSAR_WITH_EPHEM],[
  # $0: choose which ephemeris files to install
  AC_ARG_WITH(
    [ephem],
    AS_HELP_STRING(
      [--with-ephem],
      [choose which ephemeris files to install: no [none], minimal, yes [all; default]]
    ),
    [
      AS_CASE(["${withval}"],
        [no|minimal|yes],[:],
        [AC_MSG_ERROR([bad value ${withval} for --with-ephem])]
      )
    ],
    [with_ephem=yes]
  )
  AM_CONDITIONAL([INSTALL_MINIMAL_EPHEM],[test "X${with_ephem}" != Xno])
  AM_CONDITIONAL([INSTALL_ALL_EPHEM],[test "X${with_ephem}" = Xyes])
  # end $0
])
