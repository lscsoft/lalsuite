AC_INIT([lalsuite],[6.5.0.2],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([
  Makefile \
  gnuscripts/Makefile \
  doxygen/Makefile \
  doxygen/doxygen-html.cfg \
  doxygen/doxygen-pdf.cfg \
  doxygen/src/Makefile
])
AM_INIT_AUTOMAKE([foreign])
AC_PREFIX_DEFAULT([/opt/lscsoft/lalsuite])

# for building doxygen
AC_SUBST([abs_lalsuite_srcdir],[`cd ${srcdir} && pwd`])
AC_PROG_LN_S
AC_CHECK_PROGS([GIT], [git], [])

LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALXML
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALSTOCHASTIC
LALSUITE_ENABLE_LALPULSAR

lalsuite_top_builddir='$(top_builddir)/..'
lalsuite_top_srcdir='$(top_srcdir)/..'
LALSUITE_BUILD="true"
LALSUITE_TOP_SRCDIR="${lalsuite_top_srcdir}"
LAL_LIBS="${lalsuite_top_builddir}/lal/lib/liblal.la"
LAL_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
export LALSUITE_BUILD LALSUITE_TOP_SRCDIR LAL_LIBS LAL_CFLAGS
LALSUPPORT_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalsupport.la"
LALSUPPORT_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
export LALSUPPORT_LIBS LALSUPPORT_CFLAGS

if test "$lalframe" = "true"; then
  LALFRAME_LIBS="${lalsuite_top_builddir}/lalframe/src/liblalframe.la"
  LALFRAME_CFLAGS="-I${lalsuite_top_builddir}/lalframe/include"
  export LALFRAME_LIBS LALFRAME_CFLAGS
fi

if test "$lalmetaio" = "true"; then
  LALMETAIO_LIBS="${lalsuite_top_builddir}/lalmetaio/src/liblalmetaio.la"
  LALMETAIO_CFLAGS="-I${lalsuite_top_builddir}/lalmetaio/include"
  export LALMETAIO_LIBS LALMETAIO_CFLAGS
fi

if test "$lalxml" = "true"; then
  LALXML_LIBS="${lalsuite_top_builddir}/lalxml/src/liblalxml.la"
  LALXML_CFLAGS="-I${lalsuite_top_builddir}/lalxml/include"
  export LALXML_LIBS LALXML_CFLAGS
fi

if test "$lalburst" = "true"; then
  LALBURST_LIBS="${lalsuite_top_builddir}/lalburst/src/liblalburst.la"
  LALBURST_CFLAGS="-I${lalsuite_top_builddir}/lalburst/include"
  export LALBURST_LIBS LALBURST_CFLAGS
fi

if test "$lalinspiral" = "true"; then
  LALINSPIRAL_LIBS="${lalsuite_top_builddir}/lalinspiral/src/liblalinspiral.la"
  LALINSPIRAL_CFLAGS="-I${lalsuite_top_builddir}/lalinspiral/include"
  export LALINSPIRAL_LIBS LALINSPIRAL_CFLAGS
fi

if test "$lalstochastic" = "true"; then
  LALSTOCHASTIC_LIBS="${lalsuite_top_builddir}/lalstochastic/src/liblalstochastic.la"
  LALSTOCHASTIC_CFLAGS="-I${lalsuite_top_builddir}/lalstochastic/include"
  export LALSTOCHASTIC_LIBS LALSTOCHASTIC_CFLAGS
fi

if test "$lalpulsar" = "true"; then
  LALPULSAR_LIBS="${lalsuite_top_builddir}/lalpulsar/src/liblalpulsar.la"
  LALPULSAR_CFLAGS="-I${lalsuite_top_builddir}/lalpulsar/include"
  export LALPULSAR_LIBS LALPULSAR_CFLAGS
fi

AC_CONFIG_SUBDIRS([lal])

if test "$lalframe" = "true"; then
  AC_CONFIG_SUBDIRS([lalframe])
fi

if test "$lalmetaio" = "true"; then
  AC_CONFIG_SUBDIRS([lalmetaio])
fi

if test "$lalxml" = "true"; then
  AC_CONFIG_SUBDIRS([lalxml])
fi

if test "$lalburst" = "true"; then
  AC_CONFIG_SUBDIRS([lalburst])
fi

if test "$lalinspiral" = "true"; then
  AC_CONFIG_SUBDIRS([lalinspiral])
fi

if test "$lalstochastic" = "true"; then
  AC_CONFIG_SUBDIRS([lalstochastic])
fi

if test "$lalpulsar" = "true"; then
  AC_CONFIG_SUBDIRS([lalpulsar])
fi

AC_CONFIG_SUBDIRS([lalapps])

AC_OUTPUT
