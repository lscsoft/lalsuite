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
  lalapps_gcc_flags="-g3 -O4 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -fno-common -Wnested-externs -Wno-format-zero-length"
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
  AC_HELP_STRING([--enable-metaio],[compile code that requires metaio/dataflow library [default=yes]]),
  [ case "${enableval}" in
      yes) metaio=true;;
      no)  metaio=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-metaio) ;;
    esac
  ], [ metaio=true ] )
])

AC_DEFUN([LALAPPS_ENABLE_XML],
[AC_ARG_ENABLE(
  [xml],
  AC_HELP_STRING([--enable-xml],[compile code for XML I/O [default=no]]),
  [ case "${enableval}" in
      yes) xml=true;;
      no)  xml=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-xml) ;;
    esac
  ], [ xml=false ] )
])

AC_DEFUN([LALAPPS_ENABLE_CFITSIO],
[AC_ARG_ENABLE(
  [cfitsio],
  AC_HELP_STRING([--enable-cfitsio],[compile code that requires cfitsio library [default=yes]]),
  [ case "${enableval}" in
      yes) cfitsio=true;;
      no) cfitsio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-cfitsio) ;;
    esac
  ], [ cfitsio=true ] )
])

AC_DEFUN([LALAPPS_DISABLE_FRAME],
[echo "**************************************************************"
 echo "*                                                            *"
 echo "* Frame support will be DISABLED:                            *"
 echo "* LALApps is being configured with --disable-frame settings. *"
 echo "*                                                            *"
 echo "**************************************************************"
 frame=false
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
