AC_INIT([lalsuite],[6.4],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_SRCDIR([lal/configure.ac])
AC_CONFIG_AUX_DIR([lal/misc])
AC_CONFIG_MACRO_DIR([lal/misc])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([lalsuite.spec])
AM_INIT_AUTOMAKE([foreign])
AC_PREFIX_DEFAULT([/opt/lscsoft/lalsuite])

LAL_ENABLE_METAIO

AM_PROG_LIBTOOL

lalsuite_top_builddir='$(top_builddir)/..'
lalsuite_top_srcdir='$(top_srcdir)/..'
LALSUITE_BUILD="true"
LAL_LIBS="${lalsuite_top_builddir}/lal/lib/liblal.la"
LAL_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
LAL_TOP_SRCDIR="${lalsuite_top_srcdir}/lal"
export LALSUITE_BUILD LAL_LIBS LAL_CFLAGS LAL_TOP_SRCDIR
LALSUPPORT_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalsupport.la"
export LALSUPPORT_LIBS
LALMETAIO_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalmetaio.la"
export LALMETAIO_LIBS
LALFRAME_LIBS="${lalsuite_top_builddir}/lalframe/src/liblalframe.la"
LALFRAME_CFLAGS="-I${lalsuite_top_builddir}/lalframe/include"
export LALFRAME_LIBS LALFRAME_CFLAGS
LALSTOCHASTIC_LIBS="${lalsuite_top_builddir}/lalstochastic/src/liblalstochastic.la"
LALSTOCHASTIC_CFLAGS="-I${lalsuite_top_builddir}/lalstochastic/include"
export LALSTOCHASTIC_LIBS LALSTOCHASTIC_CFLAGS
LALXML_LIBS="${lalsuite_top_builddir}/lalxml/src/liblalxml.la"
LALXML_CFLAGS="-I${lalsuite_top_builddir}/lalxml/include"
export LALXML_LIBS LALXML_CFLAGS

if test "$metaio" = "true" ; then
	LIBS_save="$LIBS"
	CPPFLAGS_save="$CPPFLAGS"
	PKG_CHECK_MODULES(METAIO, libmetaio, [ true ], [ false ])
	LIBS="$METAIO_LIBS $LIBS"
	CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
	AC_SEARCH_LIBS(MetaioOpen, metaio, , [AC_MSG_WARN(disabling metaio support)]
		[metaio="false"])
	CPPFLAGS="$CPPFLAGS_save"
	LIBS="$LIBS_save"
	LIBS_save="$LIBS"
	CPPFLAGS_save="$CPPFLAGS"
fi

if test "$metaio" = "true" ; then
	LALMETAIO_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalmetaio.la"
	export LALMETAIO_LIBS
fi

AC_CONFIG_SUBDIRS([lal])
AC_CONFIG_SUBDIRS([lalframe])
AC_CONFIG_SUBDIRS([lalxml])
AC_CONFIG_SUBDIRS([lalstochastic])
AC_CONFIG_SUBDIRS([lalapps])

AC_OUTPUT
