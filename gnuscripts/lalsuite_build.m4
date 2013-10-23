# lalsuite_build.m4 - top level build macros
#
# serial 71

AC_DEFUN([LALSUITE_ADD_CFLAGS],[
  # all flags are appended to CPPFLAGS/CFLAGS
  lalsuite_append="$1"
  # print diagnostics to config.log
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: before: CPPFLAGS=${CPPFLAGS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: before: CFLAGS=${CFLAGS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: append: ${lalsuite_append}" >&AS_MESSAGE_LOG_FD
  # CPPFLAGS gets -I and -D, CFLAGS gets everything else
  # only save unique -I flags in CPPFLAGS; first instance takes precedence
  # order non-system -I before system -I in CPPFLAGS
  lalsuite_nonsysI=""
  lalsuite_sysI=""
  lalsuite_cppflags=""
  lalsuite_cflags=""
  for arg in ${CPPFLAGS} ${CFLAGS} ${lalsuite_append}; do
    AS_CASE([${arg}],
      [-I/usr/*|-I/opt/*],[
        AS_CASE([" ${lalsuite_sysI} "],
          [*" ${arg} "*],[:],
          [lalsuite_sysI="${lalsuite_sysI} ${arg}"]
        )
      ],
      [-I*],[
        AS_CASE([" ${lalsuite_nonsysI} "],
          [*" ${arg} "*],[:],
          [lalsuite_nonsysI="${lalsuite_nonsysI} ${arg}"]
        )
      ],
      [-D*],[lalsuite_cppflags="${lalsuite_cppflags} ${arg}"],
      [lalsuite_cflags="${lalsuite_cflags} ${arg}"]
    )
  done
  CPPFLAGS="${lalsuite_nonsysI} ${lalsuite_sysI} ${lalsuite_cppflags}"
  CFLAGS="${lalsuite_cflags}"
  # print diagnostics to config.log
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: lalsuite_nonsysI=${lalsuite_nonsysI}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: lalsuite_sysI=${lalsuite_sysI}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: lalsuite_cppflags=${lalsuite_cppflags}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: lalsuite_cflags=${lalsuite_cflags}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: CPPFLAGS=${CPPFLAGS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_CFLAGS]: after: CFLAGS=${CFLAGS}" >&AS_MESSAGE_LOG_FD
])

AC_DEFUN([LALSUITE_ADD_LIBS],[
  # -l flags and non-flags are prepended to LIBS
  # all other flags are appended to LDFLAGS
  lalsuite_prepend=""
  lalsuite_append=""
  for arg in $1; do
    AS_CASE([${arg}],
      [-l*],[lalsuite_prepend="${lalsuite_prepend} ${arg}"],
      [-*],[lalsuite_append="${lalsuite_append} ${arg}"],
      [lalsuite_prepend="${lalsuite_prepend} ${arg}"]
    )
  done
  # print diagnostics to config.log
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: before: LDFLAGS=${LDFLAGS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: before: LIBS=${LIBS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: prepend: ${lalsuite_prepend}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: append: ${lalsuite_append}" >&AS_MESSAGE_LOG_FD
  # LDFLAGS gets -L and other flags, LIBS gets -l and non-flags
  # only save unique -L flags in LDFLAGS; first instance takes precedence
  # order non-system -L before system -L in LDFLAGS
  # only save unique -l flags in LIBS; last instance takes precedence
  lalsuite_nonsysL=""
  lalsuite_sysL=""
  lalsuite_ldflags=""
  lalsuite_libs_rev=""
  for arg in ${lalsuite_prepend} ${LDFLAGS} ${LIBS} ${lalsuite_append}; do
    AS_CASE([${arg}],
      [-L/usr/*|-L/opt/*],[
        AS_CASE([" ${lalsuite_sysL} "],
          [*" ${arg} "*],[:],
          [lalsuite_sysL="${lalsuite_sysL} ${arg}"]
        )
      ],
      [-L*],[
        AS_CASE([" ${lalsuite_nonsysL} "],
          [*" ${arg} "*],[:],
          [lalsuite_nonsysL="${lalsuite_nonsysL} ${arg}"]
        )
      ],
      [-l*],[lalsuite_libs_rev="${arg} ${lalsuite_libs_rev}"],
      [-*],[lalsuite_ldflags="${lalsuite_ldflags} ${arg}"],
      [lalsuite_libs_rev="${arg} ${lalsuite_libs_rev}"]
    )
  done
  lalsuite_libs=""
  for arg in ${lalsuite_libs_rev}; do
    AS_CASE([" ${lalsuite_libs} "],
      [*" ${arg} "*],[:],
      [lalsuite_libs="${arg} ${lalsuite_libs}"]
    )
  done
  LDFLAGS="${lalsuite_nonsysL} ${lalsuite_sysL} ${lalsuite_ldflags}"
  LIBS="${lalsuite_libs}"
  # print diagnostics to config.log
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: lalsuite_nonsysL=${lalsuite_nonsysL}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: lalsuite_sysL=${lalsuite_sysL}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: lalsuite_ldflags=${lalsuite_ldflags}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: lalsuite_libs=${lalsuite_libs}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: LDFLAGS=${LDFLAGS}" >&AS_MESSAGE_LOG_FD
  $as_echo "$as_me:${as_lineno-$LINENO}: [LALSUITE_ADD_LIBS]: after: LIBS=${LIBS}" >&AS_MESSAGE_LOG_FD
])

AC_DEFUN([LALSUITE_WITH_CFLAGS_LIBS],[
  AC_ARG_WITH([cflags],
    AC_HELP_STRING([--with-cflags=CFLAGS],[C preprocessor/compiler flags]),
    AS_IF([test "x${with_cflags}" != x],[
      CPPFLAGS=
      CFLAGS=
      LALSUITE_ADD_CFLAGS(${with_cflags})
    ])
  )
  AC_ARG_WITH([extra_cflags],
    AC_HELP_STRING([--with-extra-cflags=CFLAGS],[extra C preprocessor/compiler flags]),
    AS_IF([test "x${with_extra_cflags}" != x],[
      LALSUITE_ADD_CFLAGS(${with_extra_cflags})
    ])
  )
  AC_ARG_WITH([libs],
    AC_HELP_STRING([--with-libs=LIBS],[linker flags]),
    AS_IF([test "x${with_libs}" != x],[
      LDFLAGS=
      LIBS=
      LALSUITE_ADD_LIBS(${with_libs})
    ])
  )
  AC_ARG_WITH([extra_libs],
    AC_HELP_STRING([--with-extra-libs=LIBS],[extra linker flags]),
    AS_IF([test "x${with_extra_libs}" != x],[
      LALSUITE_ADD_LIBS(${with_extra_libs})
    ])
  )
])

AC_DEFUN([LALSUITE_CHECK_GIT_REPO],[
  # check for git
  AC_PATH_PROGS(GIT,[git],[false])
  # check whether building from a git repository
  have_git_repo=no
  AS_IF([test "x${GIT}" != xfalse],[
    AC_MSG_CHECKING([whether building from a git repository])
    # git log will print:
    # * the last log message, if the cwd is in a git repository
    # * nothing, if the cwd is not part of the git repo (e.g. ignored)
    # * an error msg to stderr if the cwd is not in a git repository
    git_log=`( cd "${srcdir}" && ${GIT} log --pretty=oneline -n 1 -- . ) 2>/dev/null`
    AS_IF([test "x${git_log}" != x],[have_git_repo=yes])
    AC_MSG_RESULT([${have_git_repo}])
  ])
  # conditional for git and building from a git repository
  AM_CONDITIONAL(HAVE_GIT_REPO,[test "x${have_git_repo}" = xyes])
  # command line for version information generation script
  AM_COND_IF(HAVE_GIT_REPO,[
    m4_pattern_allow([AM_DEFAULT_VERBOSITY])
    m4_pattern_allow([AM_V_GEN])
    AC_SUBST([genvcsinfo_],["\$(genvcsinfo_\$(AM_DEFAULT_VERBOSITY))"])
    AC_SUBST([genvcsinfo_0],["--am-v-gen='\$(AM_V_GEN)'"])
    GENERATE_VCS_INFO="\$(PYTHON) \$(top_srcdir)/../gnuscripts/generate_vcs_info.py --git-path='\$(GIT)' \$(genvcsinfo_\$(V))"
  ],[GENERATE_VCS_INFO=false])
  AC_SUBST(GENERATE_VCS_INFO)
])

AC_DEFUN([LALSUITE_REQUIRE_CXX],[
  # require a C++ compiler
  lalsuite_require_cxx=true
])

AC_DEFUN([LALSUITE_REQUIRE_F77],[
  # require an F77 compiler
  lalsuite_require_f77=true
])

# because we want to decide whether to run AC_PROG_CXX/AC_PROG_CXXCPP
# at ./configure run time, we must erase the following macros, which
# (in autoconf 2.64 and later) require AC_PROG_CXX/AC_PROG_CXXCPP to
# be AC_REQUIRE'd at ./configure build time, regardless of whether
# they're needed or not (which is only decided later at run time).
m4_defun([AC_LANG_COMPILER(C++)],[])
m4_defun([AC_LANG_PREPROC(C++)],[])
# Same for Fortran compilers
m4_defun([AC_LANG_COMPILER(Fortran 77)],[])
m4_defun([AC_LANG_PREPROC(Fortran 77)],[])
m4_defun([AC_LANG_COMPILER(Fortran)],[])
m4_defun([AC_LANG_PREPROC(Fortran)],[])

AC_DEFUN([LALSUITE_PROG_COMPILERS],[
  # check for C99 compiler
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_CC_C99])
  AC_REQUIRE([AC_PROG_CPP])

  # check for clang
  AS_IF([test "x$GCC" = xyes],
    [AS_IF([test "`$CC -v 2>&1 | grep -c 'clang'`" != "0"],[CLANG_CC=1])],
    [CLANG_CC=])
  AC_SUBST(CLANG_CC)

  # check for C++ compiler, if needed
  AS_IF([test "${lalsuite_require_cxx}" = true],[
    AC_PROG_CXX
    AC_PROG_CXXCPP

    # check for clang++
    AS_IF([test "x$GXX" = xyes],
      [AS_IF([test "`$CXX -v 2>&1 | grep -c 'clang'`" != "0"],[CLANG_CXX=1])],
      [CLANG_CXX=])
    AC_SUBST(CLANG_CXX)
  ],[
    CXX=
    CXXCPP=
    AM_CONDITIONAL([am__fastdepCXX],[test 1 == 0])
  ])

  # check complex numbers
  LALSUITE_CHECK_C99_COMPLEX_NUMBERS
  AS_IF([test "${lalsuite_require_cxx}" = true],[
    LALSUITE_CHECK_CXX_COMPLEX_NUMBERS
  ])

  # check for F77 compiler, if needed
  AS_IF([test "${lalsuite_require_f77}" = true],[
    AC_PROG_F77
  ],[
    F77=
  ])
])

AC_DEFUN([LALSUITE_PROG_INSTALL],[
  # check for installer
  AC_REQUIRE([AC_PROG_INSTALL])
  # add -C to preserve timestamps
  INSTALL="${INSTALL} -C"
])

AC_DEFUN([LALSUITE_USE_LIBTOOL],
[## $0: Generate a libtool script for use in configure tests
AC_REQUIRE([LT_INIT])
LT_OUTPUT
m4_append([AC_LANG(C)],
[ac_link="./libtool --mode=link --tag=CC $ac_link"
])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
[m4_append([AC_LANG(C++)],
[ac_link="./libtool --mode=link --tag=CXX $ac_link"
])])[]dnl
AC_LANG(_AC_LANG)[]dnl
]) # LALSUITE_USE_LIBTOOL

AC_DEFUN([LALSUITE_ARG_VAR],[
  AC_ARG_VAR(LALSUITE_BUILD,[Set if part of lalsuite build])
  AC_ARG_VAR(LALSUITE_SUBDIRS,[Set to subdirs configured by lalsuite])
])

AC_DEFUN([LALSUITE_MULTILIB_LIBTOOL_HACK],
[## $0: libtool incorrectly determine library path on SL6
case "${host}" in
  x86_64-*-linux-gnu*)
    case `cat /etc/redhat-release 2> /dev/null` in
      "Scientific Linux"*|"CentOS"*)
        AC_MSG_NOTICE([hacking round broken libtool multilib support on RedHat systems])
        lt_cv_sys_lib_dlsearch_path_spec="/lib64 /usr/lib64"
        ;;
    esac
    ;;
esac
]) # LALSUITE_MULTILIB_LIBTOOL_HACK

# store configure flags for 'make distcheck'
AC_DEFUN([LALSUITE_DISTCHECK_CONFIGURE_FLAGS],[
  DISTCHECK_CONFIGURE_FLAGS=
  for arg in ${ac_configure_args}; do
    case ${arg} in
      (\'--enable-*\'|\'--disable-*\')
        # save any --enable/--disable arguments
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} ${arg}";;
      (\'--with-*\'|\'--without-*\')
        # save any --with/--without arguments
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} ${arg}";;
      (\'--*\')
        # skip all other ./configure arguments
        : ;;
      (\'DISTCHECK_CONFIGURE_FLAGS=*\')
        # append value of DISTCHECK_CONFIGURE_FLAGS
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} "`expr "X${arg}" : "X'DISTCHECK_CONFIGURE_FLAGS=\(.*\)'"`;;
      (\'*=*\')
        # save any environment variables given to ./configure
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} ${arg}";;
    esac
  done
  AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)
])

AC_DEFUN([LALSUITE_ENABLE_MODULE],[
AM_CONDITIONAL([$1],[test x$$2 = xtrue])
eval $1_ENABLE_VAL="`eval test "$$2" = "true" && echo "ENABLED" || echo "DISABLED"`"
])

AC_DEFUN([LALSUITE_CHECK_LIB],[
m4_pushdef([lowercase],translit([[$1]], [A-Z], [a-z]))
m4_pushdef([uppercase],translit([[$1]], [a-z], [A-Z]))
PKG_CHECK_MODULES(uppercase,[lowercase >= $2],[lowercase="true"
  if test "x${uppercase[]_DATADIR}" = x; then
    uppercase[]_DATADIR=`${PKG_CONFIG} --variable=pkgdatadir "lowercase >= $2" 2>/dev/null`
  fi
],[lowercase="false"])
if test "$lowercase" = "true"; then
  LALSUITE_ADD_CFLAGS(${uppercase[]_CFLAGS})
  LALSUITE_ADD_LIBS(${uppercase[]_LIBS})
  LALSUITE_CHECKED_LIBS="${LALSUITE_CHECKED_LIBS} lowercase"
  if test "$LALSUITE_BUILD" = "true"; then
    AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
    lowercase="true"
  else
    AC_CHECK_LIB(lowercase,[$3],[lowercase="true"],[AC_MSG_ERROR([could not find the $1 library])])
    AC_CHECK_HEADERS([$4],,[AC_MSG_ERROR([could not find the $4 header])])
    if test "$1" != "LALSupport"; then
      LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK([$1])
    fi
    AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
  fi
else
  AC_MSG_ERROR([could not find the $1 library])
fi
m4_if(lowercase,[lalsupport],[],[
  AC_ARG_VAR(uppercase[]_DATADIR, [data directory for ]uppercase[, overriding pkg-config])
])
m4_popdef([lowercase])
m4_popdef([uppercase])
])

AC_DEFUN([LALSUITE_CHECK_OPT_LIB],[
m4_pushdef([lowercase],translit([[$1]], [A-Z], [a-z]))
m4_pushdef([uppercase],translit([[$1]], [a-z], [A-Z]))
if test "$lowercase" = "true"; then
  PKG_CHECK_MODULES(uppercase,[lowercase >= $2],[lowercase="true"
    if test "x${uppercase[]_DATADIR}" = x; then
      uppercase[]_DATADIR=`${PKG_CONFIG} --variable=pkgdatadir "lowercase >= $2" 2>/dev/null`
    fi
  ],[lowercase="false"])
  if test "$lowercase" = "true"; then
    LALSUITE_ADD_CFLAGS(${uppercase[]_CFLAGS})
    LALSUITE_ADD_LIBS(${uppercase[]_LIBS})
    LALSUITE_CHECKED_LIBS="${LALSUITE_CHECKED_LIBS} lowercase"
    if test "$LALSUITE_BUILD" = "true"; then
      AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
      lowercase="true"
    else
      AC_CHECK_LIB(lowercase,[$3],[lowercase="true"],[lowercase=false
        AC_MSG_WARN([could not find the $1 library])])
      if test "$lowercase" = true; then
        AC_CHECK_HEADERS([$4],,[lowercase=false])
        if test "$lowercase" = true; then
          if test "$1" != "LALSupport"; then
            LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK([$1])
          fi
          if test "$lowercase" = true; then
            AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
          fi
        fi
      fi
    fi
  fi
fi
LALSUITE_ENABLE_MODULE(uppercase,lowercase)
m4_if(lowercase,[lalsupport],[],[
  AC_ARG_VAR(uppercase[]_DATADIR, [data directory for ]uppercase[, overriding pkg-config])
])
m4_popdef([lowercase])
m4_popdef([uppercase])
])

AC_DEFUN([LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK],[
AC_MSG_CHECKING([whether $1 headers match the library])
lib_structure=`echo $1 | sed 's/LAL/lal/'`VCSInfo
header_structure=`echo $1 | sed 's/LAL/lal/'`HeaderVCSInfo
AC_RUN_IFELSE(
  [AC_LANG_SOURCE([[
#include <string.h>
#include <stdlib.h>
#include <lal/$1VCSInfo.h>
int main(void) { exit(XLALVCSInfoCompare(&$lib_structure, &$header_structure) ? 1 : 0); }
  ]])],
  [
    AC_MSG_RESULT(yes)
  ],
  [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Your $1 headers do not match your
library. Check config.log for details.
])
  ],
  [
    AC_MSG_WARN([cross compiling: not checking])
  ]
)
])

AC_DEFUN([LALSUITE_ENABLE_NIGHTLY],
[AC_ARG_ENABLE(
  [nightly],
  AC_HELP_STRING([--enable-nightly],[nightly build [default=no]]),
  [ case "${enableval}" in
      yes) NIGHTLY_VERSION=`date -u +"%Y%m%d"`
           VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
      no) NIGHTLY_VERSION="";;
      *) NIGHTLY_VERSION="${enableval}"
         VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
      esac ],
  [ NIGHTLY_VERSION="" ] )
  AC_SUBST(NIGHTLY_VERSION)
])

AC_DEFUN([LALSUITE_ENABLE_DEBUG],
[AC_ARG_ENABLE(
  [debug],
  AC_HELP_STRING([--enable-debug],[include standard LAL debugging code [default=yes]]),
  [AS_CASE(["${enableval}"],
    [yes],,
    [no],AC_DEFINE(LAL_NDEBUG, 1, Suppress debugging code),
    AC_MSG_ERROR(bad value for ${enableval} for --enable-debug))
  ], )
])

AC_DEFUN([LALSUITE_ENABLE_ALL_LAL],
[AC_ARG_ENABLE(
  [all_lal],
  AC_HELP_STRING([--enable-all-lal],[enable/disable compilation of all LAL libraries]),
  [ case "${enableval}" in
      yes) all_lal=true;;
      no) all_lal=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-all-lal) ;;
    esac
  ], [ all_lal= ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALFRAME],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalframe],
  AC_HELP_STRING([--enable-lalframe],[compile code that requires lalframe library [default=yes]]),
  [ case "${enableval}" in
      yes) lalframe=true;;
      no) lalframe=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalframe) ;;
    esac
  ], [ lalframe=${all_lal:-true} ] )
if test "$frame" = "false"; then
  lalframe=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALMETAIO],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalmetaio],
  AC_HELP_STRING([--enable-lalmetaio],[compile code that requires lalmetaio library [default=yes]]),
  [ case "${enableval}" in
      yes) lalmetaio=true;;
      no) lalmetaio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalmetaio) ;;
    esac
  ], [ lalmetaio=${all_lal:-true} ] )
if test "$metaio" = "false"; then
  lalmetaio=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALXML],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalxml],
  AC_HELP_STRING([--enable-lalxml],[compile code that requires lalxml library [default=no]]),
  [ case "${enableval}" in
      yes) lalxml=true;;
      no) lalxml=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalxml) ;;
    esac
  ], [ lalxml=${all_lal:-false} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALSIMULATION],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalsimulation],
  AC_HELP_STRING([--enable-lalsimulation],[compile code that requires lalsimulation library [default=yes]]),
  [ case "${enableval}" in
      yes) lalsimulation=true;;
      no) lalsimulation=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalsimulation) ;;
    esac
  ], [ lalsimulation=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALBURST],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalburst],
  AC_HELP_STRING([--enable-lalburst],[compile code that requires lalburst library [default=yes]]),
  [ case "${enableval}" in
      yes) lalburst=true;;
      no) lalburst=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalburst) ;;
    esac
  ], [ lalburst=${all_lal:-true} ] )
if test "$lalmetaio" = "false"; then
  lalburst=false
fi
if test "$lalsimulation" = "false"; then
  lalburst=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALDETCHAR],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [laldetchar],
  AC_HELP_STRING([--enable-laldetchar],[compile code that requires laldetchar library [default=no]]),
  [ case "${enableval}" in
      yes) laldetchar=true;;
      no) laldetchar=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-laldetchar) ;;
    esac
  ], [ laldetchar=${all_lal:-true} ] )
if test "$lalmetaio" = "false"; then
  laldetchar=false
fi
if test "$lalburst" = "false"; then
  laldetchar=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALINSPIRAL],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalinspiral],
  AC_HELP_STRING([--enable-lalinspiral],[compile code that requires lalinspiral library [default=yes]]),
  [ case "${enableval}" in
      yes) lalinspiral=true;;
      no) lalinspiral=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalinspiral) ;;
    esac
  ], [ lalinspiral=${all_lal:-true} ] )
if test "$lalframe" = "false"; then
  lalinspiral=false
fi
if test "$lalmetaio" = "false"; then
  lalinspiral=false
fi
if test "$lalsimulation" = "false"; then
  lalinspiral=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALPULSAR],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalpulsar],
  AC_HELP_STRING([--enable-lalpulsar],[compile code that requires lalpulsar library [default=yes]]),
  [ case "${enableval}" in
      yes) lalpulsar=true;;
      no) lalpulsar=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalpulsar) ;;
    esac
  ], [ lalpulsar=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALSTOCHASTIC],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalstochastic],
  AC_HELP_STRING([--enable-lalstochastic],[compile code that requires lalstochastic library [default=yes]]),
  [ case "${enableval}" in
      yes) lalstochastic=true;;
      no) lalstochastic=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalstochastic) ;;
    esac
  ], [ lalstochastic=${all_lal:-true} ] )
if test "$lalmetaio" = "false"; then
  lalstochastic=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALINFERENCE],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalinference],
  AC_HELP_STRING([--enable-lalinference],[compile code that requires lalinference library [default=yes]]),
  [ case "${enableval}" in
      yes) lalinference=true;;
      no) lalinference=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalinference) ;;
    esac
  ], [ lalinference=${all_lal:-true} ] )
if test "$lalmetaio" = "false"; then
  lalinference=false
fi
if test "$lalframe" = "false"; then
  lalinference=false
fi
if test "$lalinspiral" = "false"; then
  lalinference=false
fi
if test "$lalpulsar" = "false"; then
  lalinference=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALAPPS],[
  AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
  AC_ARG_ENABLE(
    [lalapps],
    AC_HELP_STRING([--enable-lalapps],[compile lalapps [default=yes]]),
    [
      case "${enableval}" in
        yes) lalapps=true ;;
        no) lalapps=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalapps) ;;
      esac
    ],[
      lalapps=${all_lal:-true}
    ]
  )
])

AC_DEFUN([LALSUITE_ENABLE_BOINC],
[AC_ARG_ENABLE(
  [boinc],
  AC_HELP_STRING([--enable-boinc],[enable BOINC support [default=no]]),
  [ case "${enableval}" in
      yes) boinc=true;;
      no) boinc=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-boinc);;
    esac
  ], [ boinc=false ] )
  AS_IF([test "${boinc}" = true],[LALSUITE_REQUIRE_CXX])
  AC_ARG_VAR([BOINC_PREFIX],[BOINC installation directory (optional)])
])

AC_DEFUN([LALSUITE_CHECK_BOINC],
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

AC_DEFUN([LALSUITE_WITH_CUDA],[
AC_ARG_WITH(
  [cuda],
  AC_HELP_STRING([--with-cuda=PATH],[specify location of CUDA [/opt/cuda]]),[
    AS_CASE([${with_cuda}],
      [no],[cuda=false],
      [yes],[cuda=true; cuda_path=/opt/cuda],
      [cuda=true; cuda_path=${with_cuda}]
    )
  ],[
    cuda=false
  ])
  AS_IF([test "${cuda}" = true],[
    LALSUITE_REQUIRE_CXX
    AC_MSG_NOTICE([Using ${with_cuda} as CUDA path])
    AS_CASE([$build_os],
      [linux*],[
        AS_IF([test "x$build_cpu" = "xx86_64"],[
          cuda_libdir=lib64
        ],[
          cuda_libdir=lib
        ])
      ],
      [cuda_libdir=lib]
    )
    CUDA_LIBS="-L${cuda_path}/${cuda_libdir} -Wl,-rpath -Wl,${cuda_path}/${cuda_libdir} -lcufft -lcudart"
    CUDA_CPPFLAGS="-I${with_cuda}/include"
    LALSUITE_ADD_LIBS(${CUDA_LIBS})
    LALSUITE_ADD_CFLAGS(${CUDA_CPPFLAGS})
    AC_SUBST(CUDA_LIBS)
    AC_SUBST(CUDA_CPPFLAGS)
    AC_PATH_PROGS(NVCC,[nvcc],[],[${cuda_path}/bin:${PATH}])
    AS_IF([test "x${NVCC}" = x],[
      AC_MSG_ERROR([could not find 'nvcc' in path])
    ])
  ])
  LALSUITE_ENABLE_MODULE([CUDA],[cuda])
])

AC_DEFUN([LALSUITE_ENABLE_FAST_GSL],
[AC_ARG_ENABLE(
  [fast_gsl],
  AC_HELP_STRING([--enable-fast-gsl],[enable fast/inline GSL code [default=no]]),
  AS_CASE(["${enableval}"],
    [yes],[AC_DEFINE([HAVE_INLINE],[1],[Define to 1 to use inline code])
           AC_DEFINE([GSL_C99_INLINE],[1],[Define to 1 to use GSL C99 inline code])
           AC_DEFINE([GSL_RANGE_CHECK_OFF],[1],[Define to 1 to turn GSL range checking off])],
    [no],,
    AC_MSG_ERROR([bad value ${enableval} for --enable-fast-gsl]))
  )
])

AC_DEFUN([LALSUITE_ENABLE_OSX_VERSION_CHECK],
[AC_ARG_ENABLE(
  [osx_version_check],
  AC_HELP_STRING([--enable-osx-version-check],[disable OS X version check [default=yes]]),
  AS_CASE(["${enableval}"],
    [yes],[osx_version_check=true],
    [no],[osx_version_check=false],
    AC_MSG_ERROR([bad value ${enableval} for --enable-osx-version-check])
  ),[osx_version_check=true])
])

AC_DEFUN([LALSUITE_OSX_VERSION_CHECK],[
LALSUITE_ENABLE_OSX_VERSION_CHECK
AS_IF([test "x${osx_version_check}" = "xtrue"],[
  AS_IF([test "x$build_vendor" = "xapple"],[
    AC_CHECK_PROGS([SW_VERS],[sw_vers])
    AS_IF([test "x$SW_VERS" != "x"],[
      AC_MSG_CHECKING([Mac OS X version])
      MACOSX_VERSION=`$SW_VERS -productVersion`
      AC_MSG_RESULT([$MACOSX_VERSION])])
    AS_CASE(["$MACOSX_VERSION"],
      [10.0*|10.1*|10.2*|10.3*],AC_MSG_ERROR([This version of Mac OS X is not supported]),
      [10.4*|10.5*|10.6*|10.7*|10.8*|10.9*],,
      AC_MSG_WARN([Unknown Mac OS X version]))
])])])

AC_DEFUN([LALSUITE_WITH_NVCC_CFLAGS],
[AC_ARG_WITH(
  [nvcc_cflags],
  AC_HELP_STRING([--with-nvcc-cflags=NVCC_CFLAGS],[NVCC compiler flags]),
  AS_IF([test -n "${with_nvcc_cflags}"],[NVCC_CFLAGS="$NVCC_CFLAGS ${with_nvcc_cflags}"]),)
])

AC_DEFUN([LALSUITE_CHECK_CUDA],
[AC_MSG_CHECKING([whether LAL has been compiled with CUDA support])
AC_TRY_RUN([
#include <lal/LALConfig.h>
#ifdef LAL_CUDA_ENABLED
int main( void ) { return 0; }
#else
int main( void ) { return 1; }
#endif
],
AC_MSG_RESULT([yes])
[cuda=true],
AC_MSG_RESULT([no])
[cuda=false],
AC_MSG_RESULT([unknown])
[cuda=false])
])

AC_DEFUN([LALSUITE_CHECK_COMPLEX_NUMBER_MEMORY_LAYOUT],[
  AC_MSG_CHECKING([memory layout of complex number type '$3'])
  AS_IF([test "$cross_compiling" = yes],[
    AC_MSG_WARN([cross compiling: not checking])
  ],[
    # compile a C file containing functions where
    # the complex number datatype is a struct
    AC_LANG_PUSH([C])
    AC_COMPILE_IFELSE([
AC_LANG_SOURCE([
AC_INCLUDES_DEFAULT
typedef struct {
  $2 re;
  $2 im;
} ComplexStruct;
const size_t zsize = sizeof(ComplexStruct);
const $2 zre = 1.414213562373095048801688724209;
const $2 zim = 3.141592653589793238462643383276;
void Function1(ComplexStruct *pz) {
  pz->re = zre;
  pz->im = zim;
}
int Function2(ComplexStruct pz) {
  return (pz.re == zre && pz.im == zim);
}
])
    ],[
      # if compilation was successful, save the compiled object
      mv -f conftest.$ac_objext conftestlink.$ac_objext
      _AS_ECHO_LOG([moved conftest.$ac_objext to conftestlink.$ac_objext])
    ],[
      AC_MSG_FAILURE([unexpected compile failure])
    ])
    AC_LANG_POP([C])
    # add the object compiled above to the objects
    # which will be linked against by the next test
    lalsuite_ccnml_LIBS=$LIBS
    LIBS="$LIBS conftestlink.$ac_objext"
    # push current language so that we can restore
    # previous linker settings at end of this test
    AC_LANG_PUSH(_AC_LANG)
    # compile a _AC_LANG file where the complex number
    # datatype is a C99/C++ complex number, and which
    # calls functions in the previously-compiled file
    # (where the complex number datatype was a struct).
    # link it against the previously-compiled object,
    # and run the resulting test program.
    AC_RUN_IFELSE([
AC_LANG_PROGRAM([
AC_INCLUDES_DEFAULT
#include <$1>
#ifdef __cplusplus
extern "C" {
#endif
extern const size_t zsize;
extern const $2 zre;
extern const $2 zim;
void Function1($3 *pz);
int Function2($3 pz);
#ifdef __cplusplus
}
#endif
],[
  $3 c;
  if (sizeof($3) != zsize) {
    return 1;
  }
  Function1(&c);
  if ($4(c) != zre || $5(c) != zim) {
    return 2;
  }
  if (!Function2(c)) {
    return 3;
  }
  return 0 /* ; */
])
    ],[
      # if test program compiled and exited
      # normally, test was successful
      AC_MSG_RESULT([compatible])
    ],[
      AC_MSG_FAILURE([memory layout of complex number type '$3' is incompatible])
    ])
    # restore previous linker settings and
    # delete remaining test object files
    LIBS=$lalsuite_ccnml_LIBS
    AC_LANG_POP(_AC_LANG)
    rm -f conftestlink.$ac_objext
  ])
])

AC_DEFUN([LALSUITE_CHECK_C99_COMPLEX_NUMBERS],[
  # check C99 complex numbers
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([complex.h],[],[AC_MSG_ERROR([could not find 'complex.h'])])
  LALSUITE_CHECK_COMPLEX_NUMBER_MEMORY_LAYOUT([complex.h],[float],[float complex],[crealf],[cimagf])
  LALSUITE_CHECK_COMPLEX_NUMBER_MEMORY_LAYOUT([complex.h],[double],[double complex],[creal],[cimag])
  AC_LANG_POP([C])
])

AC_DEFUN([LALSUITE_CHECK_CXX_COMPLEX_NUMBERS],[
  # check C++ complex numbers
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADERS([complex],[],[AC_MSG_ERROR([could not find 'complex'])])
  LALSUITE_CHECK_COMPLEX_NUMBER_MEMORY_LAYOUT([complex],[float],[std::complex<float>],[real],[imag])
  LALSUITE_CHECK_COMPLEX_NUMBER_MEMORY_LAYOUT([complex],[double],[std::complex<double>],[real],[imag])
  AC_LANG_POP([C++])
])
