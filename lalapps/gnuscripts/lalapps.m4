# lalapps.m4 - lalapps specific autoconf macros
#
# serial 11

AC_DEFUN([LALAPPS_WITH_EXTRA_CPPFLAGS],
[AC_ARG_WITH(
  [extra_cppflags],
  AC_HELP_STRING([--with-extra-cppflags=CPPFLAGS],[additional C preprocessor flags]),
  AS_IF([test -n "${with_extra_cppflags}"],[CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}"])
  ,)
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

AC_DEFUN([LALAPPS_ENABLE_STATIC_BINARIES],
[AC_ARG_ENABLE(
  [static_binaries],
  AC_HELP_STRING([--enable-static-binaries],[build static binaries [default=no]]),
  [ case "${enableval}" in
      yes) static_binaries=true;;
      no)  static_binaries=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-static-binaries) ;;
    esac
  ], [ static_binaries=false ] )
if test "$condor" = "true"; then
  static_binaries=false
fi
if test "$boinc" = "true"; then
  static_binaries=false
fi
])

AC_DEFUN([LALAPPS_ENABLE_MPI],
[AC_ARG_ENABLE(
  [mpi],
  AC_HELP_STRING([--enable-mpi],[compile using MPI for supported codes [default=no]]),
  [ case "${enableval}" in
      yes) mpi=true;;
      no)  mpi=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
    esac
  ], [ mpi=false ] )
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

AC_DEFUN([LALAPPS_ENABLE_GDS],
[AC_ARG_ENABLE(
  [pss],
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

AC_DEFUN([LALAPPS_ENABLE_BAMBI],
[AC_ARG_ENABLE(
  [bambi],
  AC_HELP_STRING([--enable-bambi],[build LALInferenceBAMBI [default=no]]),
  [ case "${enableval}" in
      yes) bambi=true;;
      no)  bambi=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-bambi) ;;
    esac
  ], [ bambi=false ] )
])

AC_DEFUN([LALAPPS_CHECK_BAMBI],[
  PKG_CHECK_MODULES([BAMBI],[bambi],[
    AX_LAPACK([
      BAMBI_LIBS="$BAMBI_LIBS $LAPACK_LIBS $BLAS_LIBS $FLIBS"
      BAMBI_ENABLE_VAL="ENABLED"
      hbf=true
    ])
  ])
])
