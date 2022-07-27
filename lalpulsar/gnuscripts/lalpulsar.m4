# -*- mode: autoconf; -*-
# lalpulsar.m4 - LALPulsar specific macros
#
# serial 3

AC_DEFUN([LALPULSAR_ENABLE_PSS],
[AC_ARG_ENABLE(
  [pss],
  AS_HELP_STRING([--enable-pss],[compile code that requires pss library [default=no]]),
  [ case "${enableval}" in
      yes) pss=true;;
      no) pss=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-pss) ;;
    esac
  ], [pss=false])
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

AC_DEFUN([LALPULSAR_INSTALL_MINIMAL_EPHEM],[
  # $0: install only a minimal set of ephemeris files
  AC_ARG_ENABLE(
    [minimal_ephem],
    AS_HELP_STRING(
      [--enable-minimal-ephem],
      [install only a minimal set of ephemeris files [default=no; install all ephemeris files]]
    ),
    [],
    [enable_minimal_ephem=no]
  )
  AM_CONDITIONAL([INSTALL_MINIMAL_EPHEM],[test x"${enable_minimal_ephem}" = xyes])
  # end $0
])
