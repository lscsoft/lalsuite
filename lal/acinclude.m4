dnl acinclude.m4

AC_DEFUN(LAL_WITH_GCC_FLAGS,
[AC_ARG_WITH(
        gcc_flags,   
        [  --with-gcc-flags        turn on strict gcc warning flags],
        [ if test -n "${with_gcc_flags}"
          then
            lal_gcc_flags="-g3 -O4 -ansi -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -fno-common -Wnested-externs -D__NO_STRING_INLINES"
          else
            lal_gcc_flags=""
          fi
        ], [ lal_gcc_flags="" ] )
])

AC_DEFUN(LAL_WITH_EXTRA_CPPFLAGS,
[AC_ARG_WITH(
	extra_cppflags, 
        [  --with-extra-cppflags=CPPFLAGS  additional C preprocessor flags],
	[ if test -n "${with_extra_cppflags}"
	  then
	    CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}";
	  fi
	],)
])

AC_DEFUN(LAL_WITH_EXTRA_CFLAGS,
[AC_ARG_WITH(
	extra_cflags, 
        [  --with-extra-cflags=CFLAGS  additional C compiler flags],
	[ if test -n "${with_extra_cflags}"
	  then
	    CFLAGS="$CFLAGS ${with_extra_cflags}";
	  fi
	],)
])

AC_DEFUN(LAL_WITH_EXTRA_LDFLAGS,
[AC_ARG_WITH(
	extra_ldflags, 
        [  --with-extra-ldflags=LDFLAGS  additional linker flags],
	[ if test -n "${with_extra_ldflags}"
	  then
	    LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
	  fi
	],)
])

AC_DEFUN(LAL_WITH_EXTRA_LIBS,
[AC_ARG_WITH(
	extra_libs, 
        [  --with-extra-libs=LIBS  additional -l and -L linker flags],
	[ if test -n "${with_extra_libs}"
	  then
	    LIBS="$LIBS ${with_extra_libs}";
	  fi
	],)
])

AC_DEFUN(LAL_WITH_MPICC,
[AC_ARG_WITH(
        mpicc, 
        [  --with-mpicc=MPICC      use the MPICC C compiler for MPI code],
        [ if test -n "${with_mpicc}"
          then
            MPICC="${with_mpicc}";
          fi
        ],)
])

AC_DEFUN(LAL_WITH_CC,
[AC_ARG_WITH(
        cc, 
        [  --with-cc=CC            use the CC C compiler],
        [ if test -n "${with_cc}"
          then
            CC="${with_cc}";
          fi
        ],)
])

AC_DEFUN(LAL_ENABLE_FRAME,
[AC_ARG_ENABLE(
        frame,
        [  --enable-frame          compile code that requires Frame library [default=no] ],
        [ case "${enableval}" in
            yes) frame=true  ;;
            no)  frame=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
          esac
        ], [ frame=false ] )
])

AC_DEFUN(LAL_ENABLE_MPI,
[AC_ARG_ENABLE(
        mpi,
        [  --enable-mpi            compile code that requires MPI [default=no] ],
        [ case "${enableval}" in
            yes) mpi=true  ;;
            no)  mpi=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
          esac
        ], [ mpi=false ] )
])

AC_DEFUN(LAL_ENABLE_DEBUG,
[AC_ARG_ENABLE(
        debug,
        [  --enable-debug          include standard LAL debugging code [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(LAL_NDEBUG, 1, Suppress debugging code) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN(LAL_ENABLE_MACROS,
[AC_ARG_ENABLE(
        macros,
        [  --enable-macros         use LAL macros [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(NOLALMACROS, 1, Use functions rather than macros) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN(LAL_ENABLE_PTHREAD_LOCK,
[AC_ARG_ENABLE(
        pthread_lock,
        [  --enable-pthread-lock   use pthread mutex lock for threadsafety [default=no] ],
        [ case "${enableval}" in
            yes) lal_pthread_lock=true; AC_DEFINE(LAL_PTHREAD_LOCK, 1, Use pthread mutex lock for threadsafety) ;;
            no) lal_pthread_lock=false ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-pthread-lock) ;;
          esac
        ], [ lal_pthread_lock=false ] )
])

AC_DEFUN(LAL_FFTW_MSG_ERROR,
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

AC_DEFUN(LAL_SFFTW_WORKS,
[AC_MSG_CHECKING(whether single precison FFTW works)
AC_TRY_RUN([
#include <stdio.h>
#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif
int main() { return sizeof(fftw_real) - 4; } ],
AC_MSG_RESULT(yes),
AC_MSG_RESULT(no)
echo "**************************************************************"
echo "* FFTW does not seem to be working.                          *"
echo "* Possible problems:                                         *"
echo "*   - FFTW version < 2.0                                     *"
echo "*   - Compiler could not find header sfftw.h or fftw.h       *"
echo "*   - FFTW was not configured with the --enable-float option *"
LAL_FFTW_MSG_ERROR,
AC_MSG_RESULT(unknown) ) ] )

AC_DEFUN(LAL_SRFFTW_WORKS,
[AC_MSG_CHECKING(whether single precison real FFTW works)
AC_TRY_RUN([
#include <stdio.h>
#ifdef HAVE_SRFFTW_H
#include <srfftw.h>
#elif HAVE_RFFTW_H
#include <rfftw.h>
#else
#error "don't have either srfftw.h or rfftw.h"
#endif
int main() { return sizeof(fftw_real) - 4; } ],
AC_MSG_RESULT(yes),
AC_MSG_RESULT(no)
echo "**************************************************************"
echo "* FFTW does not seem to be working.                          *"
echo "* Possible problems:                                         *"
echo "*   - FFTW version < 2.0                                     *"
echo "*   - Compiler could not find header srfftw.h or rfftw.h     *"
echo "*   - FFTW was not configured with the --enable-float option *"
echo "**************************************************************"
LAL_FFTW_MSG_ERROR,
AC_MSG_RESULT(unknown) ) ] )

AC_DEFUN(LAL_CHECK_FRAMELIB,
[ if test "${frame}" = "true"; then
        AC_CHECK_LIB(Frame, FrLibIni, ,
          [AC_MSG_ERROR(couldn't find Frame library for --enable-frame)] , )
        AC_MSG_CHECKING([whether Frame library version >= 3.85])
        AC_TRY_RUN([#include "FrameL.h"
          int main() { return FRAMELIB_VERSION < 3.85 ? 1 : 0 ; }],
          AC_MSG_RESULT(yes),
          [AC_MSG_RESULT(no),
            AC_MSG_ERROR(FrameL.h not found or FRAMELIB_VERSION < 3.85)],
          AC_MSG_RESULT(unknown))
  fi
])

AC_DEFUN(LAL_CHECK_MPI,
[ AC_CHECK_PROGS(MPICC, mpicc hcc, $CC)
  AC_MSG_CHECKING(for mpicc flags)
  SHOWARG=""
  MPICPPFLAGS=""
  MPICFLAGS=""
  MPILDFLAGS=""
  if (($MPICC -compile_info 1>/dev/null 2>/dev/null) && ($MPICC -link_info 1>/dev/null 2>/dev/null)) ; then
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
      SHOWARG="-show"
    elif ($MPICC -showme 1>/dev/null 2>/dev/null) ; then
      SHOWARG="-showme"
    else
      AC_MSG_WARN(couldn't determine mpi compile flags)
    fi
    if test -n "$SHOWARG" ; then
      for mpiarg in `$MPICC $SHOWARG` ; do
        case $mpiarg in
          -D*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
          -I*) MPICPPFLAGS="$MPICPPFLAGS $mpiarg" ;;
          -L*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
          -l*) MPILDFLAGS="$MPILDFLAGS $mpiarg" ;;
          -Wl*) MPICFLAGS="$MPICFLAGS $mpiarg" ;;
        esac
      done
    fi
  fi
  AC_MSG_RESULT($MPICPPFLAGS $MPICFLAGS $MPILDFLAGS)
  LIBS="$LIBS $MPILDFLAGS"
  CPPFLAGS="$CPPFLAGS $MPICPPFLAGS"
  CFLAGS="$CFLAGS $MPICFLAGS"
  AC_CHECK_HEADER(mpi.h, ,AC_MSG_ERROR(can't find mpi.h))
  AC_MSG_CHECKING(whether mpi works)
  AC_TRY_LINK([#include <mpi.h>
    ], MPI_Finalize();,
    AC_MSG_RESULT(yes),
    AC_MSG_RESULT(no)
    AC_MSG_ERROR(mpi does not work))
  AC_MSG_CHECKING(mpi type)
  MPITYPE=NONE
  AC_TRY_COMPILE([
    #include <mpi.h>
    #ifndef LAM_MPI
    #error "not LAM"
    #endif], , [ MPITYPE=LAM
                 AC_MSG_RESULT(lam)],
    AC_TRY_COMPILE([
      #include <mpi.h>
      #ifndef MPICH_NAME
      #erro "not MPICH"
      #endif], , [ MPITYPE=MPICH
                   AC_MSG_RESULT(mpich)],
      AC_MSG_RESULT(couldn't determine)
      AC_MSG_ERROR(mpi must be either lam or mpich)
    )
  )
  AC_SUBST(MPITYPE)dnl
])


dnl This is AC_CHECK_SIZEOF but prepends LAL.

AC_DEFUN(LAL_CHECK_SIZEOF,
[changequote(<<, >>)dnl
dnl The name to #define.
define(<<LAL_TYPE_NAME>>, translit(lal_sizeof_$1, [a-z *], [A-Z_P]))dnl
dnl The cache variable name.
define(<<LAL_CV_NAME>>, translit(lal_cv_sizeof_$1, [ *], [_p]))dnl
changequote([, ])dnl
AC_MSG_CHECKING(size of $1)
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
AC_MSG_RESULT($LAL_CV_NAME)
AC_DEFINE_UNQUOTED(LAL_TYPE_NAME, $LAL_CV_NAME)
undefine([LAL_TYPE_NAME])dnl
undefine([LAL_CV_NAME])dnl
])



#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each `.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],[
AC_OUTPUT_COMMANDS([
test x"$AMDEP_TRUE" != x"" ||
for mf in $CONFIG_FILES; do
  case "$mf" in
  *:*) mf=`echo "$mf"|sed 's%:.*%%'` ;;
  esac
  case "$mf" in
  Makefile) dirpart=.;;
  */Makefile) dirpart=`echo "$mf" | sed -e 's|/[^/]*$||'`;;
  *) continue;;
  esac
  grep '^DEP_FILES *= *[^ #]' < "$mf" > /dev/null || continue
  # Extract the definition of DEP_FILES from the Makefile without
  # running `make'.
  DEPDIR=`sed -n -e '/^DEPDIR = / s///p' < "$mf"`
  test -z "$DEPDIR" && continue
  # When using ansi2knr, U may be empty or an underscore; expand it
  U=`sed -n -e '/^U = / s///p' < "$mf"`
  test -d "$dirpart/$DEPDIR" || mkdir "$dirpart/$DEPDIR"
  # We invoke sed twice because it is the simplest approach to
  # changing $(DEPDIR) to its actual value in the expansion.
  for file in `sed -n -e '
    /^DEP_FILES = .*\\\\$/ {
      s/^DEP_FILES = //
      :loop
        s/\\\\$//
        p
        n
        /\\\\$/ b loop
      p
    }
    /^DEP_FILES = / s/^DEP_FILES = //p' < "$mf" | \
       sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g' -e 's/\$U/'"$U"'/g'`; do
    # Make sure the directory exists.
    test -f "$dirpart/$file" && continue
    fdir=`echo "$file" | sed -e 's|/[^/]*$||'`
    $ac_aux_dir/mkinstalldirs "$dirpart/$fdir" > /dev/null 2>&1
    # echo "creating $dirpart/$file"
    echo '# dummy' > "$dirpart/$file"
  done
done
], [AMDEP_TRUE="$AMDEP_TRUE"
ac_aux_dir="$ac_aux_dir"])])




## libtool.m4 - Configure libtool for the target system. -*-Shell-script-*-
## Copyright (C) 1996-1999, 2000 Free Software Foundation, Inc.
## Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## As a special exception to the GNU General Public License, if you
## distribute this file as part of a program that contains a
## configuration script generated by Autoconf, you may include it under
## the same distribution terms that you use for the rest of that program.

# serial 40 AC_PROG_LIBTOOL
AC_DEFUN(AC_PROG_LIBTOOL,
[AC_REQUIRE([AC_LIBTOOL_SETUP])dnl

# Save cache, so that ltconfig can load it
AC_CACHE_SAVE

# Actually configure libtool.  ac_aux_dir is where install-sh is found.
CC="$CC" CFLAGS="$CFLAGS" CPPFLAGS="$CPPFLAGS" \
LD="$LD" LDFLAGS="$LDFLAGS" LIBS="$LIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB" \
DLLTOOL="$DLLTOOL" AS="$AS" OBJDUMP="$OBJDUMP" \
${CONFIG_SHELL-/bin/sh} $ac_aux_dir/ltconfig --no-reexec \
$libtool_flags --no-verify $ac_aux_dir/ltmain.sh $lt_target \
|| AC_MSG_ERROR([libtool configure failed])

# Reload cache, that may have been modified by ltconfig
AC_CACHE_LOAD

# This can be used to rebuild libtool when needed
LIBTOOL_DEPS="$ac_aux_dir/ltconfig $ac_aux_dir/ltmain.sh"

# Always use our own libtool.
LIBTOOL='$(SHELL) $(top_builddir)/libtool'
AC_SUBST(LIBTOOL)dnl

# Redirect the config.log output again, so that the ltconfig log is not
# clobbered by the next message.
exec 5>>./config.log
])

AC_DEFUN(AC_LIBTOOL_SETUP,
[AC_PREREQ(2.13)dnl
AC_REQUIRE([AC_ENABLE_SHARED])dnl
AC_REQUIRE([AC_ENABLE_STATIC])dnl
AC_REQUIRE([AC_ENABLE_FAST_INSTALL])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([AC_PROG_RANLIB])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_LD])dnl
AC_REQUIRE([AC_PROG_NM])dnl
AC_REQUIRE([AC_PROG_LN_S])dnl
dnl

case "$target" in
NONE) lt_target="$host" ;;
*) lt_target="$target" ;;
esac

# Check for any special flags to pass to ltconfig.
libtool_flags="--cache-file=$cache_file"
test "$enable_shared" = no && libtool_flags="$libtool_flags --disable-shared"
test "$enable_static" = no && libtool_flags="$libtool_flags --disable-static"
test "$enable_fast_install" = no && libtool_flags="$libtool_flags --disable-fast-install"
test "$ac_cv_prog_gcc" = yes && libtool_flags="$libtool_flags --with-gcc"
test "$ac_cv_prog_gnu_ld" = yes && libtool_flags="$libtool_flags --with-gnu-ld"
ifdef([AC_PROVIDE_AC_LIBTOOL_DLOPEN],
[libtool_flags="$libtool_flags --enable-dlopen"])
ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
[libtool_flags="$libtool_flags --enable-win32-dll"])
AC_ARG_ENABLE(libtool-lock,
  [  --disable-libtool-lock  avoid locking (might break parallel builds)])
test "x$enable_libtool_lock" = xno && libtool_flags="$libtool_flags --disable-lock"
test x"$silent" = xyes && libtool_flags="$libtool_flags --silent"

# Some flags need to be propagated to the compiler or linker for good
# libtool support.
case "$lt_target" in
*-*-irix6*)
  # Find out which ABI we are using.
  echo '[#]line __oline__ "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case "`/usr/bin/file conftest.o`" in
    *32-bit*)
      LD="${LD-ld} -32"
      ;;
    *N32*)
      LD="${LD-ld} -n32"
      ;;
    *64-bit*)
      LD="${LD-ld} -64"
      ;;
    esac
  fi
  rm -rf conftest*
  ;;

*-*-sco3.2v5*)
  # On SCO OpenServer 5, we need -belf to get full-featured binaries.
  SAVE_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -belf"
  AC_CACHE_CHECK([whether the C compiler needs -belf], lt_cv_cc_needs_belf,
    [AC_TRY_LINK([],[],[lt_cv_cc_needs_belf=yes],[lt_cv_cc_needs_belf=no])])
  if test x"$lt_cv_cc_needs_belf" != x"yes"; then
    # this is probably gcc 2.8.0, egcs 1.0 or newer; no need for -belf
    CFLAGS="$SAVE_CFLAGS"
  fi
  ;;

ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
[*-*-cygwin* | *-*-mingw*)
  AC_CHECK_TOOL(DLLTOOL, dlltool, false)
  AC_CHECK_TOOL(AS, as, false)
  AC_CHECK_TOOL(OBJDUMP, objdump, false)
  ;;
])
esac
])

# AC_LIBTOOL_DLOPEN - enable checks for dlopen support
AC_DEFUN(AC_LIBTOOL_DLOPEN, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])])

# AC_LIBTOOL_WIN32_DLL - declare package support for building win32 dll's
AC_DEFUN(AC_LIBTOOL_WIN32_DLL, [AC_BEFORE([$0], [AC_LIBTOOL_SETUP])])

# AC_ENABLE_SHARED - implement the --enable-shared flag
# Usage: AC_ENABLE_SHARED[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN(AC_ENABLE_SHARED, [dnl
define([AC_ENABLE_SHARED_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(shared,
changequote(<<, >>)dnl
<<  --enable-shared[=PKGS]  build shared libraries [default=>>AC_ENABLE_SHARED_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case "$enableval" in
yes) enable_shared=yes ;;
no) enable_shared=no ;;
*)
  enable_shared=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_shared=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_shared=AC_ENABLE_SHARED_DEFAULT)dnl
])

# AC_DISABLE_SHARED - set the default shared flag to --disable-shared
AC_DEFUN(AC_DISABLE_SHARED, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_SHARED(no)])

# AC_ENABLE_STATIC - implement the --enable-static flag
# Usage: AC_ENABLE_STATIC[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN(AC_ENABLE_STATIC, [dnl
define([AC_ENABLE_STATIC_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(static,
changequote(<<, >>)dnl
<<  --enable-static[=PKGS]  build static libraries [default=>>AC_ENABLE_STATIC_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case "$enableval" in
yes) enable_static=yes ;;
no) enable_static=no ;;
*)
  enable_static=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_static=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_static=AC_ENABLE_STATIC_DEFAULT)dnl
])

# AC_DISABLE_STATIC - set the default static flag to --disable-static
AC_DEFUN(AC_DISABLE_STATIC, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_STATIC(no)])


# AC_ENABLE_FAST_INSTALL - implement the --enable-fast-install flag
# Usage: AC_ENABLE_FAST_INSTALL[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN(AC_ENABLE_FAST_INSTALL, [dnl
define([AC_ENABLE_FAST_INSTALL_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(fast-install,
changequote(<<, >>)dnl
<<  --enable-fast-install[=PKGS]  optimize for fast installation [default=>>AC_ENABLE_FAST_INSTALL_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case "$enableval" in
yes) enable_fast_install=yes ;;
no) enable_fast_install=no ;;
*)
  enable_fast_install=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_fast_install=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_fast_install=AC_ENABLE_FAST_INSTALL_DEFAULT)dnl
])

# AC_ENABLE_FAST_INSTALL - set the default to --disable-fast-install
AC_DEFUN(AC_DISABLE_FAST_INSTALL, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_FAST_INSTALL(no)])

# AC_PROG_LD - find the path to the GNU or non-GNU linker
AC_DEFUN(AC_PROG_LD,
[AC_ARG_WITH(gnu-ld,
[  --with-gnu-ld           assume the C compiler uses GNU ld [default=no]],
test "$withval" = no || with_gnu_ld=yes, with_gnu_ld=no)
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
ac_prog=ld
if test "$ac_cv_prog_gcc" = yes; then
  # Check if gcc -print-prog-name=ld gives a path.
  AC_MSG_CHECKING([for ld used by GCC])
  ac_prog=`($CC -print-prog-name=ld) 2>&5`
  case "$ac_prog" in
    # Accept absolute paths.
changequote(,)dnl
    [\\/]* | [A-Za-z]:[\\/]*)
      re_direlt='/[^/][^/]*/\.\./'
changequote([,])dnl
      # Canonicalize the path of ld
      ac_prog=`echo $ac_prog| sed 's%\\\\%/%g'`
      while echo $ac_prog | grep "$re_direlt" > /dev/null 2>&1; do
	ac_prog=`echo $ac_prog| sed "s%$re_direlt%/%"`
      done
      test -z "$LD" && LD="$ac_prog"
      ;;
  "")
    # If it fails, then pretend we aren't using GCC.
    ac_prog=ld
    ;;
  *)
    # If it is relative, then search for the first ld in PATH.
    with_gnu_ld=unknown
    ;;
  esac
elif test "$with_gnu_ld" = yes; then
  AC_MSG_CHECKING([for GNU ld])
else
  AC_MSG_CHECKING([for non-GNU ld])
fi
AC_CACHE_VAL(ac_cv_path_LD,
[if test -z "$LD"; then
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}${PATH_SEPARATOR-:}"
  for ac_dir in $PATH; do
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$ac_prog" || test -f "$ac_dir/$ac_prog$ac_exeext"; then
      ac_cv_path_LD="$ac_dir/$ac_prog"
      # Check to see if the program is GNU ld.  I'd rather use --version,
      # but apparently some GNU ld's only accept -v.
      # Break only if it was the GNU/non-GNU ld that we prefer.
      if "$ac_cv_path_LD" -v 2>&1 < /dev/null | egrep '(GNU|with BFD)' > /dev/null; then
	test "$with_gnu_ld" != no && break
      else
	test "$with_gnu_ld" != yes && break
      fi
    fi
  done
  IFS="$ac_save_ifs"
else
  ac_cv_path_LD="$LD" # Let the user override the test with a path.
fi])
LD="$ac_cv_path_LD"
if test -n "$LD"; then
  AC_MSG_RESULT($LD)
else
  AC_MSG_RESULT(no)
fi
test -z "$LD" && AC_MSG_ERROR([no acceptable ld found in \$PATH])
AC_PROG_LD_GNU
])

AC_DEFUN(AC_PROG_LD_GNU,
[AC_CACHE_CHECK([if the linker ($LD) is GNU ld], ac_cv_prog_gnu_ld,
[# I'd rather use --version here, but apparently some GNU ld's only accept -v.
if $LD -v 2>&1 </dev/null | egrep '(GNU|with BFD)' 1>&5; then
  ac_cv_prog_gnu_ld=yes
else
  ac_cv_prog_gnu_ld=no
fi])
])

# AC_PROG_NM - find the path to a BSD-compatible name lister
AC_DEFUN(AC_PROG_NM,
[AC_MSG_CHECKING([for BSD-compatible nm])
AC_CACHE_VAL(ac_cv_path_NM,
[if test -n "$NM"; then
  # Let the user override the test.
  ac_cv_path_NM="$NM"
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}${PATH_SEPARATOR-:}"
  for ac_dir in $PATH /usr/ccs/bin /usr/ucb /bin; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/nm || test -f $ac_dir/nm$ac_exeext ; then
      # Check to see if the nm accepts a BSD-compat flag.
      # Adding the `sed 1q' prevents false positives on HP-UX, which says:
      #   nm: unknown option "B" ignored
      if ($ac_dir/nm -B /dev/null 2>&1 | sed '1q'; exit 0) | egrep /dev/null >/dev/null; then
	ac_cv_path_NM="$ac_dir/nm -B"
	break
      elif ($ac_dir/nm -p /dev/null 2>&1 | sed '1q'; exit 0) | egrep /dev/null >/dev/null; then
	ac_cv_path_NM="$ac_dir/nm -p"
	break
      else
	ac_cv_path_NM=${ac_cv_path_NM="$ac_dir/nm"} # keep the first match, but
	continue # so that we can try to find one that supports BSD flags
      fi
    fi
  done
  IFS="$ac_save_ifs"
  test -z "$ac_cv_path_NM" && ac_cv_path_NM=nm
fi])
NM="$ac_cv_path_NM"
AC_MSG_RESULT([$NM])
])

# AC_CHECK_LIBM - check for math library
AC_DEFUN(AC_CHECK_LIBM,
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
LIBM=
case "$lt_target" in
*-*-beos* | *-*-cygwin*)
  # These system don't have libm
  ;;
*-ncr-sysv4.3*)
  AC_CHECK_LIB(mw, _mwvalidcheckl, LIBM="-lmw")
  AC_CHECK_LIB(m, main, LIBM="$LIBM -lm")
  ;;
*)
  AC_CHECK_LIB(m, main, LIBM="-lm")
  ;;
esac
])

# AC_LIBLTDL_CONVENIENCE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl convenience library and INCLTDL to the include flags for
# the libltdl header and adds --enable-ltdl-convenience to the
# configure arguments.  Note that LIBLTDL and INCLTDL are not
# AC_SUBSTed, nor is AC_CONFIG_SUBDIRS called.  If DIR is not
# provided, it is assumed to be `libltdl'.  LIBLTDL will be prefixed
# with '${top_builddir}/' and INCLTDL will be prefixed with
# '${top_srcdir}/' (note the single quotes!).  If your package is not
# flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
AC_DEFUN(AC_LIBLTDL_CONVENIENCE, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  case "$enable_ltdl_convenience" in
  no) AC_MSG_ERROR([this package needs a convenience libltdl]) ;;
  "") enable_ltdl_convenience=yes
      ac_configure_args="$ac_configure_args --enable-ltdl-convenience" ;;
  esac
  LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdlc.la
  INCLTDL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
])

# AC_LIBLTDL_INSTALLABLE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl installable library and INCLTDL to the include flags for
# the libltdl header and adds --enable-ltdl-install to the configure
# arguments.  Note that LIBLTDL and INCLTDL are not AC_SUBSTed, nor is
# AC_CONFIG_SUBDIRS called.  If DIR is not provided and an installed
# libltdl is not found, it is assumed to be `libltdl'.  LIBLTDL will
# be prefixed with '${top_builddir}/' and INCLTDL will be prefixed
# with '${top_srcdir}/' (note the single quotes!).  If your package is
# not flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
# In the future, this macro may have to be called after AC_PROG_LIBTOOL.
AC_DEFUN(AC_LIBLTDL_INSTALLABLE, [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  AC_CHECK_LIB(ltdl, main,
  [test x"$enable_ltdl_install" != xyes && enable_ltdl_install=no],
  [if test x"$enable_ltdl_install" = xno; then
     AC_MSG_WARN([libltdl not installed, but installation disabled])
   else
     enable_ltdl_install=yes
   fi
  ])
  if test x"$enable_ltdl_install" = x"yes"; then
    ac_configure_args="$ac_configure_args --enable-ltdl-install"
    LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdl.la
    INCLTDL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  else
    ac_configure_args="$ac_configure_args --enable-ltdl-install=no"
    LIBLTDL="-lltdl"
    INCLTDL=
  fi
])

dnl old names
AC_DEFUN(AM_PROG_LIBTOOL, [indir([AC_PROG_LIBTOOL])])dnl
AC_DEFUN(AM_ENABLE_SHARED, [indir([AC_ENABLE_SHARED], $@)])dnl
AC_DEFUN(AM_ENABLE_STATIC, [indir([AC_ENABLE_STATIC], $@)])dnl
AC_DEFUN(AM_DISABLE_SHARED, [indir([AC_DISABLE_SHARED], $@)])dnl
AC_DEFUN(AM_DISABLE_STATIC, [indir([AC_DISABLE_STATIC], $@)])dnl
AC_DEFUN(AM_PROG_LD, [indir([AC_PROG_LD])])dnl
AC_DEFUN(AM_PROG_NM, [indir([AC_PROG_NM])])dnl

dnl This is just to silence aclocal about the macro not being used
ifelse([AC_DISABLE_FAST_INSTALL])dnl
