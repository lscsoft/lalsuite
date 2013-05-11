AC_PREREQ([2.63])
AC_INIT([lalsuite],[6.16.1],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([
  Makefile \
  gnuscripts/Makefile \
  doxygen/Makefile \
  doxygen/src/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests])

LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# for building doxygen
AC_PROG_LN_S
AC_CHECK_PROGS([GIT], [git], [])

LALSUITE_ENABLE_ALL_LAL
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALXML
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALDETCHAR
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALSTOCHASTIC
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALINFERENCE
LALSUITE_ENABLE_LALAPPS

lalsuite_top_builddir='$(top_builddir)/..'
lalsuite_top_srcdir='$(top_srcdir)/..'

LALSUITE_BUILD="true"
LALSUITE_SUBDIRS="${ac_subdirs_all}"
export LALSUITE_BUILD LALSUITE_SUBDIRS

LAL_LIBS="-L${libdir} ${lalsuite_top_builddir}/lal/lib/lal/liblal.la"
LAL_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
LAL_DATADIR="${lalsuite_top_srcdir}/lal/test"
export LAL_LIBS LAL_CFLAGS LAL_DATADIR

LALSUPPORT_LIBS="-L${libdir} ${lalsuite_top_builddir}/lal/lib/lalsupport/src/liblalsupport.la"
LALSUPPORT_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
export LALSUPPORT_LIBS LALSUPPORT_CFLAGS

if test "$lalframe" = "true"; then
  LALFRAME_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalframe/src/liblalframe.la"
  LALFRAME_CFLAGS="-I${lalsuite_top_builddir}/lalframe/include"
  LALFRAME_DATADIR="${lalsuite_top_srcdir}/lalframe/test"
  export LALFRAME_LIBS LALFRAME_CFLAGS LALFRAME_DATADIR
fi

if test "$lalmetaio" = "true"; then
  LALMETAIO_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalmetaio/src/liblalmetaio.la"
  LALMETAIO_CFLAGS="-I${lalsuite_top_builddir}/lalmetaio/include"
  LALMETAIO_DATADIR="${lalsuite_top_srcdir}/lalmetaio/test"
  export LALMETAIO_LIBS LALMETAIO_CFLAGS LALMETAIO_DATADIR
fi

if test "$lalxml" = "true"; then
  LALXML_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalxml/src/liblalxml.la"
  LALXML_CFLAGS="-I${lalsuite_top_builddir}/lalxml/include"
  LALXML_DATADIR="${lalsuite_top_srcdir}/lalxml/test"
  export LALXML_LIBS LALXML_CFLAGS LALXML_DATADIR
fi

if test "$lalsimulation" = "true"; then
  LALSIMULATION_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalsimulation/src/liblalsimulation.la"
  LALSIMULATION_CFLAGS="-I${lalsuite_top_builddir}/lalsimulation/include"
  LALSIMULATION_DATADIR="${lalsuite_top_srcdir}/lalsimulation/test"
  export LALSIMULATION_LIBS LALSIMULATION_CFLAGS LALSIMULATION_DATADIR
fi

if test "$lalburst" = "true"; then
  LALBURST_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalburst/src/liblalburst.la"
  LALBURST_CFLAGS="-I${lalsuite_top_builddir}/lalburst/include"
  LALBURST_DATADIR="${lalsuite_top_srcdir}/lalburst/test"
  export LALBURST_LIBS LALBURST_CFLAGS LALBURST_DATADIR
fi

if test "$laldetchar" = "true"; then
  LALDETCHAR_LIBS="-L${libdir} ${lalsuite_top_builddir}/laldetchar/src/liblaldetchar.la"
  LALDETCHAR_CFLAGS="-I${lalsuite_top_builddir}/laldetchar/include"
  LALDETCHAR_DATADIR="${lalsuite_top_srcdir}/laldetchar/test"
  export LALDETCHAR_LIBS LALDETCHAR_CFLAGS LALDETCHAR_DATADIR
fi

if test "$lalinspiral" = "true"; then
  LALINSPIRAL_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalinspiral/src/liblalinspiral.la"
  LALINSPIRAL_CFLAGS="-I${lalsuite_top_builddir}/lalinspiral/include"
  LALINSPIRAL_DATADIR="${lalsuite_top_srcdir}/lalinspiral/test"
  export LALINSPIRAL_LIBS LALINSPIRAL_CFLAGS LALINSPIRAL_DATADIR
fi

if test "$lalstochastic" = "true"; then
  LALSTOCHASTIC_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalstochastic/src/liblalstochastic.la"
  LALSTOCHASTIC_CFLAGS="-I${lalsuite_top_builddir}/lalstochastic/include"
  LALSTOCHASTIC_DATADIR="${lalsuite_top_srcdir}/lalstochastic/test"
  export LALSTOCHASTIC_LIBS LALSTOCHASTIC_CFLAGS LALSTOCHASTIC_DATADIR
fi

if test "$lalpulsar" = "true"; then
  LALPULSAR_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalpulsar/src/liblalpulsar.la"
  LALPULSAR_CFLAGS="-I${lalsuite_top_builddir}/lalpulsar/include"
  LALPULSAR_DATADIR="${lalsuite_top_srcdir}/lalpulsar/test"
  export LALPULSAR_LIBS LALPULSAR_CFLAGS LALPULSAR_DATADIR
fi

if test "$lalinference" = "true"; then
  LALINFERENCE_LIBS="-L${libdir} ${lalsuite_top_builddir}/lalinference/src/liblalinference.la"
  LALINFERENCE_CFLAGS="-I${lalsuite_top_builddir}/lalinference/include"
  LALINFERENCE_DATADIR="${lalsuite_top_srcdir}/lalinference/test"
  export LALINFERENCE_LIBS LALINFERENCE_CFLAGS LALINFERENCE_DATADIR
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

if test "$lalsimulation" = "true"; then
  AC_CONFIG_SUBDIRS([lalsimulation])
fi

if test "$lalburst" = "true"; then
  AC_CONFIG_SUBDIRS([lalburst])
fi

if test "$laldetchar" = "true"; then
  AC_CONFIG_SUBDIRS([laldetchar])
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

if test "$lalinference" = "true"; then
  AC_CONFIG_SUBDIRS([lalinference])
fi

if test "$lalapps" = "true"; then
  AC_CONFIG_SUBDIRS([lalapps])
fi

AC_SUBST([ac_configure_args])
AC_OUTPUT
