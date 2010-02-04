dnl lalapps.m4

AC_DEFUN([LALAPPS_WITH_EXTRA_CPPFLAGS],
[AC_ARG_WITH(
  [extra_cppflags],
  AC_HELP_STRING([--with-extra-cppflags=CPPFLAGS],[additional C preprocessor flags]),
  [ if test -n "${with_extra_cppflags}"
    then
      CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}";
    fi
  ],)
])

AC_DEFUN([LALAPPS_WITH_LAL_PREFIX],
[AC_ARG_WITH(
  [lal_prefix],
  AC_HELP_STRING([--with-lal-prefix=LAL_PREFIX],[location where to find LAL installation]),
  [ if test -n "${with_lal_prefix}"
    then
      LAL_PREFIX="${with_lal_prefix}"
    fi
  ],)
])

AC_DEFUN([LALAPPS_WITH_EXTRA_CFLAGS],
[AC_ARG_WITH(
  [extra_cflags],
  AC_HELP_STRING([--with-extra-cflags=CFLAGS],[additional C compiler flags]),
  [ if test -n "${with_extra_cflags}"
    then
      CFLAGS="$CFLAGS ${with_extra_cflags}";
    fi
  ],)
])

AC_DEFUN([LALAPPS_WITH_EXTRA_LDFLAGS],
[AC_ARG_WITH(
  [extra_ldflags],
  AC_HELP_STRING([--with-extra-ldflags=LDFLAGS],[additional linker flags]),
  [ if test -n "${with_extra_ldflags}"
    then
      LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
    fi
  ],)
])

AC_DEFUN([LALAPPS_WITH_EXTRA_LIBS],
[AC_ARG_WITH(
  [extra_libs],
  AC_HELP_STRING([--with-extra-libs=LIBS],[additional -l and -L linker flags]),
  [ if test -n "${with_extra_libs}"
    then
      LIBS="$LIBS ${with_extra_libs}";
    fi
  ],)
])

AC_DEFUN([LALAPPS_ENABLE_GCC_FLAGS],
[AC_ARG_ENABLE([gcc_flags],
  AC_HELP_STRING([--enable-gcc-flags],[turn on strict gcc warning flags (default=yes)]),
  [case "${enableval}" in
     yes) DO_ENABLE_LALAPPS_GCC_FLAGS;;
     no) ;;
     *) DO_ENABLE_LALAPPS_GCC_FLAGS;;
   esac ],
   [ DO_ENABLE_LALAPPS_GCC_FLAGS ] )
])

AC_DEFUN([DO_ENABLE_LALAPPS_GCC_FLAGS],
[
  lalapps_gcc_flags="-g3 -O4 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -Wnested-externs -Wno-format-zero-length"
])

AC_DEFUN([LALAPPS_WITH_CC],
[AC_ARG_WITH(
  [cc],
  AC_HELP_STRING([--with-cc=CC],[use the CC C compiler]),
  [ if test -n "${with_cc}"
    then
      CC="${with_cc}";
    fi
  ],)
])

AC_DEFUN([LALAPPS_ENABLE_CONDOR],
[AC_ARG_ENABLE(
  [condor],
  AC_HELP_STRING([--enable-condor],[compile for use with condor [default=no]]),
  [ case "${enableval}" in
      yes) condor=true;;
      no)  condor=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-condor) ;;
    esac
  ], [ condor=false ] )
])

AC_DEFUN([LALAPPS_ENABLE_BOINC],
[AC_ARG_ENABLE(
  [boinc],
  AC_HELP_STRING([--enable-boinc],[enable BOINC support [default=no]]),
  [ case "${enableval}" in
      yes) boinc=true;;
      no) boinc=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-boinc);;
    esac
  ], [ boinc=false ] )
AC_ARG_VAR([BOINC_PREFIX],[BOINC installation directory (optional)])
])

AC_DEFUN([LALAPPS_CHECK_BOINC],
[AC_MSG_CHECKING([whether LAL has been compiled with BOINC support])
AC_TRY_RUN([
#include <lal/LALConfig.h>
#ifdef LAL_BOINC_ENABLED
int main( void ) { return 0; }
#else
int main( void ) { return 1; }
#endif
],
AC_MSG_RESULT([yes])
[boinc=true],
AC_MSG_RESULT([no])
[boinc=false],
AC_MSG_RESULT([unknown])
[boinc=false])
])

AC_DEFUN([LALAPPS_ENABLE_FRAME],
[AC_ARG_ENABLE(
  [frame],
  AC_HELP_STRING([--enable-frame],[compile code that requires Frame library [default=yes]]),
  [ case "${enableval}" in
      yes) frame=true;;
      no)  frame=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
    esac
  ], [ frame=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_METAIO],
[AC_ARG_ENABLE(
  [metaio],
  AC_HELP_STRING([--enable-metaio],[compile code that requires metaio library [default=yes]]),
  [ case "${enableval}" in
      yes) metaio=true;;
      no)  metaio=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-metaio) ;;
    esac
  ], [ metaio=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_CFITSIO],
[AC_ARG_ENABLE(
  [cfitsio],
  AC_HELP_STRING([--enable-cfitsio],[compile code that requires cfitsio library [default=no]]),
  [ case "${enableval}" in
      yes) cfitsio=true;;
      no) cfitsio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-cfitsio) ;;
    esac
  ], [ cfitsio=false ] )
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

AC_DEFUN([LALAPPS_ENABLE_LALFRAME],
[AC_ARG_ENABLE(
  [lalframe],
  AC_HELP_STRING([--enable-lalframe],[compile code that requires lalframe library [default=yes]]),
  [ case "${enableval}" in
      yes) lalframe=true;;
      no) lalframe=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
    esac
  ], [ lalframe=true ] )
if test "$frame" = "false"; then
  lalframe=false
fi
])

AC_DEFUN([LALAPPS_ENABLE_LALMETAIO],
[AC_ARG_ENABLE(
  [lalmetaio],
  AC_HELP_STRING([--enable-lalmetaio],[compile code that requires lalmetaio library [default=yes]]),
  [ case "${enableval}" in
      yes) lalmetaio=true;;
      no) lalmetaio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-metaio) ;;
    esac
  ], [ lalmetaio=true ] )
if test "$metaio" = "false"; then
  lalmetaio=false
fi
])

AC_DEFUN([LALAPPS_ENABLE_LALBURST],
[AC_ARG_ENABLE(
  [lalburst],
  AC_HELP_STRING([--enable-lalburst],[compile code that requires lalburst library [default=yes]]),
  [ case "${enableval}" in
      yes) lalburst=true;;
      no) lalburst=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-burst) ;;
    esac
  ], [ lalburst=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_LALPULSAR],
[AC_ARG_ENABLE(
  [lalpulsar],
  AC_HELP_STRING([--enable-lalpulsar],[compile code that requires lalpulsar library [default=yes]]),
  [ case "${enableval}" in
      yes) lalpulsar=true;;
      no) lalpulsar=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalpulsar) ;;
    esac
  ], [ lalpulsar=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_LALSTOCHASTIC],
[AC_ARG_ENABLE(
  [lalstochastic],
  AC_HELP_STRING([--enable-lalstochastic],[compile code that requires lalstochastic library [default=yes]]),
  [ case "${enableval}" in
      yes) lalstochastic=true;;
      no) lalstochastic=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-stochastic) ;;
    esac
  ], [ lalstochastic=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_LALXML],
[AC_ARG_ENABLE(
  [lalxml],
  AC_HELP_STRING([--enable-lalxml],[compile code that requires lalxml library [default=no]]),
  [ case "${enableval}" in
      yes) lalxml=true;;
      no) lalxml=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalxml) ;;
    esac
  ], [ lalxml=false ] )
if test "$lalpulsar" = "false"; then
  lalxml=false
fi
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

AC_DEFUN([LALAPPS_ENABLE_NIGHTLY],
[AC_ARG_ENABLE(
  [nightly],
  AC_HELP_STRING([--enable-nightly],[nightly build [default=no]]),
  [ case "${enableval}" in
      yes) NIGHTLY_VERSION=`date +"%Y%m%d"`
           VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
      no)  NIGHTLY_VERSION="";;
      *)   NIGHTLY_VERSION="${enableval}"
           VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
    esac ],
  [ NIGHTLY_VERSION="" ] )
  AC_SUBST(NIGHTLY_VERSION)
])
