AC_PREREQ([2.63])
AC_INIT([lalsuite-doxygen],[6.23.1],[lal-discuss@ligo.org])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([
  Makefile \
  gnuscripts/Makefile \
  src/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests])

AC_PATH_PROG(GIT, [git], [])
AC_PROG_LN_S

AC_SUBST([ac_configure_args])
AC_OUTPUT
