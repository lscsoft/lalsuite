dnl lal.m4

AC_DEFUN([LAL_WITH_GCC_FLAGS],
[AC_ARG_WITH(
        [gcc_flags],   
        [  --with-gcc-flags        turn on strict gcc warning flags],
        [ if test -n "${with_gcc_flags}"
          then
            lal_gcc_flags="-g3 -O4 -ansi -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -fno-common -Wnested-externs -D__NO_STRING_INLINES"
          else
            lal_gcc_flags=""
          fi
        ], [ lal_gcc_flags="" ] )
])

AC_DEFUN([LAL_WITH_EXTRA_CPPFLAGS],
[AC_ARG_WITH(
	[extra_cppflags], 
        [  --with-extra-cppflags=CPPFLAGS  additional C preprocessor flags],
	[ if test -n "${with_extra_cppflags}"
	  then
	    CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_CFLAGS],
[AC_ARG_WITH(
	[cflags], 
        [  --with-cflags=CFLAGS        C compiler flags],
	[ if test -n "${with_cflags}"
	  then
	    CFLAGS="${with_cflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_CFLAGS],
[AC_ARG_WITH(
	[extra_cflags], 
        [  --with-extra-cflags=CFLAGS  additional C compiler flags],
	[ if test -n "${with_extra_cflags}"
	  then
	    CFLAGS="$CFLAGS ${with_extra_cflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_LDFLAGS],
[AC_ARG_WITH(
	[extra_ldflags], 
        [  --with-extra-ldflags=LDFLAGS  additional linker flags],
	[ if test -n "${with_extra_ldflags}"
	  then
	    LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_LIBS],
[AC_ARG_WITH(
	[extra_libs], 
        [  --with-extra-libs=LIBS  additional -l and -L linker flags],
	[ if test -n "${with_extra_libs}"
	  then
	    LIBS="$LIBS ${with_extra_libs}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_MPICC],
[AC_ARG_WITH(
        [mpicc], 
        [  --with-mpicc=MPICC      use the MPICC C compiler for MPI code],
        [ if test -n "${with_mpicc}"
          then
            MPICC="${with_mpicc}";
          fi
        ],)
])

AC_DEFUN([LAL_WITH_CC],
[AC_ARG_WITH(
        [cc], 
        [  --with-cc=CC            use the CC C compiler],
        [ if test -n "${with_cc}"
          then
            CC="${with_cc}";
          fi
        ],)
])

AC_DEFUN([LAL_ENABLE_FFTW3],
[AC_ARG_ENABLE(
        [fftw3],
        [  --enable-fftw3          use fftw3 [default=no] ],
        [ case "${enableval}" in
            yes) fftw3=true;;
            no)  fftw3=false;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-fftw3) ;;
          esac
        ], [ fftw3=true ] )
])

AC_DEFUN([LAL_ENABLE_FRAME],
[AC_ARG_ENABLE(
        [frame],
        [  --enable-frame          compile code that requires Frame library [default=no] ],
        [ case "${enableval}" in
            yes) frame=true;;
            no)  frame=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
          esac
        ], [ frame=false ] )
])

AC_DEFUN([LAL_ENABLE_MPI],
[AC_ARG_ENABLE(
        [mpi],
        [  --enable-mpi            compile code that requires MPI [default=no] ],
        [ case "${enableval}" in
            yes) mpi=true;;
            no)  mpi=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
          esac
        ], [ mpi=false ] )
])

AC_DEFUN([LAL_ENABLE_DATAFLOW],
[AC_ARG_ENABLE(
  [dataflow],
  [  --enable-dataflow       compile code that requires metaio/dataflow library [default=no] ],
  [ case "${enableval}" in
      yes) dataflow=true;;
      no)  dataflow=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
    esac
  ], [ dataflow=false ] )
AC_ARG_ENABLE(
  [metaio],
  [  --enable-metaio         compile code that requires metaio/dataflow library [default=no] ],
  [ case "${enableval}" in
      yes) metaio=true;;
      no)  metaio=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
    esac
  ], [ metaio=false ] )
  if test "x$metaio" = "xtrue" ; then dataflow=true ; fi
])

AC_DEFUN([LAL_ENABLE_INTELFFT],
[AC_ARG_ENABLE(
        [intelfft],
        [  --enable-intelfft       use Intel FFT libraries insted of FFTW [default=no] ],
        [ case "${enableval}" in
            yes) intelfft=true;;
            no)  intelfft=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-intelfft) ;;
          esac
        ], [ intelfft=false ] )
])

AC_DEFUN([LAL_ENABLE_DEBUG],
[AC_ARG_ENABLE(
        [debug],
        [  --enable-debug          include standard LAL debugging code [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(LAL_NDEBUG, 1, Suppress debugging code) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN([LAL_ENABLE_MACROS],
[AC_ARG_ENABLE(
        [macros],
        [  --enable-macros         use LAL macros [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(NOLALMACROS, 1, Use functions rather than macros) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN([LAL_ENABLE_PTHREAD_LOCK],
[AC_ARG_ENABLE(
        [pthread_lock],
        [  --enable-pthread-lock   use pthread mutex lock for threadsafety [default=no] ],
        [ case "${enableval}" in
            yes) lal_pthread_lock=true; AC_DEFINE(LAL_PTHREAD_LOCK, 1, Use pthread mutex lock for threadsafety) ;;
            no) lal_pthread_lock=false ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-pthread-lock) ;;
          esac
        ], [ lal_pthread_lock=false ] )
])

AC_DEFUN([LAL_INTEL_FFT_LIBS_MSG_ERROR],
[echo "**************************************************************"
 echo "* The --enable-static and --enable-shared options are        *"
 echo "* mutually exclusive with Intel FFT libraries.               *"
 echo "*                                                            *"
 echo "* Please reconfigure with:                                   *"
 echo "*                                                            *"
 echo "*   --enable-static --disable-shared                         *"
 echo "*                                                            *"
 echo "* for static libraries or                                    *"
 echo "*                                                            *"
 echo "*   --disable-static --enable-shared                         *"
 echo "*                                                            *"
 echo "* for shared libraries.                                      *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
AC_MSG_ERROR([Intel FFT must use either static or shared libraries])
])



AC_DEFUN([LAL_FFTW_MSG_ERROR],
[echo "**************************************************************"
 echo "* You must install FFTW (v >= 2.0) on your system.           *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
AC_MSG_ERROR([single precision FFTW must be properly installed.])
])

AC_DEFUN([LAL_SFFTW_WORKS],
[AC_MSG_CHECKING([whether single precison FFTW works])
AC_TRY_RUN([
#include <stdio.h>
#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif
int main() { return (sizeof(fftw_real)!=4 || fftw_sizeof_fftw_real()!=4); }
],
AC_MSG_RESULT([yes]),
AC_MSG_RESULT([no])
[
echo "**************************************************************"
echo "* FFTW does not seem to be working.                          *"
echo "* Possible problems:                                         *"
echo "*   - FFTW version < 2.0                                     *"
echo "*   - Could not find header sfftw.h, fftw.h, or fftw3.h      *"
echo "*   - FFTW was not configured with the --enable-float option *"
echo "* Consult file config.log for details                        *"
]
LAL_FFTW_MSG_ERROR,
AC_MSG_RESULT([unknown]) ) ] )

AC_DEFUN([LAL_SRFFTW_WORKS],
[AC_MSG_CHECKING([whether single precison real FFTW works])
AC_TRY_RUN([
#include <stdio.h>
#ifdef HAVE_SRFFTW_H
#include <srfftw.h>
#elif HAVE_RFFTW_H
#include <rfftw.h>
#else
#error "don't have either srfftw.h or rfftw.h"
#endif
int main() { return sizeof(fftw_real) - 4; }
],
AC_MSG_RESULT([yes]),
AC_MSG_RESULT([no])
[
echo "**************************************************************"
echo "* FFTW does not seem to be working.                          *"
echo "* Possible problems:                                         *"
echo "*   - FFTW version < 2.0                                     *"
echo "*   - Could not find header srfftw.h, rfftw.h, or fftw3.h    *"
echo "*   - FFTW was not configured with the --enable-float option *"
echo "* Consult file config.log for details                        *"
echo "**************************************************************"
]
LAL_FFTW_MSG_ERROR,
AC_MSG_RESULT([unknown]) ) ] )

AC_DEFUN([LAL_CHECK_FRAMELIB],
[ if test "${frame}" = "true"; then
        lal_check_framelib_save_LIBS="$LIBS"
        AC_CHECK_LIB(Frame, FrLibIni, ,
          [AC_MSG_ERROR([couldn't find Frame library for --enable-frame])] )
        AC_MSG_CHECKING([whether Frame library version >= 6.00])
        AC_TRY_RUN([#include "FrameL.h"
          int main() { return FRAMELIB_VERSION < 6.00 ? 1 : 0 ; }],
          AC_MSG_RESULT([yes]),
          [AC_MSG_RESULT([no])
            AC_MSG_ERROR([FrameL.h not found or FRAMELIB_VERSION < 6.00])],
          AC_MSG_RESULT([unknown]))
        LIBS="$lal_check_framelib_save_LIBS"
  fi
])

AC_DEFUN([LAL_CHECK_MPI],
[ AC_CHECK_PROGS(MPICC, mpicc hcc, $CC)
  AC_MSG_CHECKING([for mpicc flags])
  SHOWARG=""
  MPICPPFLAGS=""
  MPICFLAGS=""
  MPILDFLAGS=""
  MPITYPE=NONE
  if (($MPICC -compile_info 1>/dev/null 2>/dev/null) && ($MPICC -link_info 1>/dev/null 2>/dev/null)) ; then
    MPITYPE=mpich
    for mpiarg in `$MPICC -compile_info` ; do
      case $mpiarg in
        -D*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
        -I*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
      esac
    done
    for mpiarg in `$MPICC -link_info` ; do
      case $mpiarg in
        -L*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
        -l*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
      esac
    done
  else
    if ($MPICC -show 1>/dev/null 2>/dev/null) ; then
      MPITYPE=MPICH
      SHOWARG="-show"
    elif ($MPICC -showme 1>/dev/null 2>/dev/null) ; then
      MPITYPE=LAM
      SHOWARG="-showme"
    else
      AC_MSG_WARN([couldn't determine mpi compile flags])
    fi
    if test -n "$SHOWARG" ; then
      for mpiarg in `$MPICC $SHOWARG` ; do
        case $mpiarg in
          -D*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
          -I*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
          -L*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
          -l*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
          -pthread) MPILDFLAGS="$MPILDFLAGS -lpthread" ;;
          -Wl*) MPICFLAGS="$MPICFLAGS $mpiarg" ;;
        esac
      done
    fi
  fi
  AC_MSG_RESULT([$MPICPPFLAGS $MPICFLAGS $MPILDFLAGS])
  LIBS="$LIBS $MPILDFLAGS"
  CPPFLAGS="$CPPFLAGS $MPICPPFLAGS"
  CFLAGS="$CFLAGS $MPICFLAGS"
  AC_CHECK_HEADER(mpi.h, ,AC_MSG_ERROR([can't find mpi.h]))
  AC_MSG_CHECKING([whether mpi works])
  AC_TRY_LINK([#include <mpi.h>
    ], MPI_Finalize();,
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([mpi does not work]))
  AC_MSG_CHECKING([mpi type])
  AC_TRY_COMPILE([
    #include <mpi.h>
    #ifndef LAM_MPI
    #error "not LAM"
    #endif], , [ MPITYPE=LAM
                 AC_MSG_RESULT([lam])],
    AC_TRY_COMPILE([
      #include <mpi.h>
      #ifndef MPICH_NAME
      #erro "not MPICH"
      #endif], , [ MPITYPE=MPICH
                   AC_MSG_RESULT([mpich])],
      [
      if test $MPITYPE = MPICH ; then
      AC_MSG_RESULT([couldn't determine... guessing mpich])
      elif test $MPITYPE = LAM ; then
      AC_MSG_RESULT([couldn't determine... assuming lam])
      else
      AC_MSG_RESULT([couldn't determine.])
      AC_MSG_ERROR([mpi must be either lam or mpich])
      fi
      ]
    )
  )
  AC_SUBST(MPITYPE)dnl
])


dnl This is AC_CHECK_SIZEOF but prepends LAL.

AC_DEFUN([LAL_CHECK_SIZEOF],
[changequote(<<, >>)dnl
dnl The name to #define.
define(<<LAL_TYPE_NAME>>, translit(lal_sizeof_$1, [a-z *], [A-Z_P]))dnl
dnl The cache variable name.
define(<<LAL_CV_NAME>>, translit(lal_cv_sizeof_$1, [ *], [_p]))dnl
changequote([, ])dnl
AC_MSG_CHECKING([size of $1])
AC_CACHE_VAL(LAL_CV_NAME,
[AC_TRY_RUN([#include <stdio.h>
main()
{
  FILE *f=fopen("conftestval", "w");
  if (!f) exit(1);
  fprintf(f, "%d\n", sizeof($1));
  exit(0);
}], LAL_CV_NAME=`cat conftestval`, LAL_CV_NAME=0, ifelse([$2], , , LAL_CV_NAME=$2))
])dnl
AC_MSG_RESULT([$LAL_CV_NAME])
AC_DEFINE_UNQUOTED(LAL_TYPE_NAME, $LAL_CV_NAME)
undefine([LAL_TYPE_NAME])dnl
undefine([LAL_CV_NAME])dnl
])

