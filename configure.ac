AC_INIT([lalsuite],[6.2],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_SRCDIR([lal/configure.ac])
AC_CONFIG_AUX_DIR([lal/misc])
AC_CONFIG_MACRO_DIR([lal/misc])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([lalsuite.spec])
AM_INIT_AUTOMAKE([foreign])
AC_PREFIX_DEFAULT([/opt/lscsoft/lalsuite])

LAL_ENABLE_FRAME
LAL_ENABLE_METAIO
LAL_ENABLE_XML

AM_PROG_LIBTOOL

lalsuite_top_builddir='$(top_builddir)/..'
LALSUITE_BUILD="true"
LAL_LIBS="${lalsuite_top_builddir}/lal/lib/liblal.la"
LAL_CFLAGS="-I${lalsuite_top_builddir}/lal/include"
export LALSUITE_BUILD LAL_LIBS LAL_CFLAGS
LALSUPPORT_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalsupport.la"
export LALSUPPORT_LIBS

if test "$frame" = "true" ; then
	LIBS_save="$LIBS"
	CPPFLAGS_save="$CPPFLAGS"
	PKG_CHECK_MODULES(FRAME, libframe, [ true ], [ false ])
	LIBS="$FRAME_LIBS $LIBS"
	CPPFLAGS="$CPPFLAGS $FRAME_CFLAGS"
	AC_SEARCH_LIBS(FrLibIni, Frame, , [AC_MSG_WARN(disabling framelib support)]
		[frame="false"])
	CPPFLAGS="$CPPFLAGS_save"
	LIBS="$LIBS_save"
fi
if test "$metaio" = "true" ; then
	LIBS_save="$LIBS"
	CPPFLAGS_save="$CPPFLAGS"
	PKG_CHECK_MODULES(METAIO, libmetaio, [ true ], [ false ])
	LIBS="$METAIO_LIBS $LIBS"
	CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
	AC_SEARCH_LIBS(MetaioOpen, metaio dataflow, , [AC_MSG_WARN(disabling metaio support)]
		[metaio="false"])
	CPPFLAGS="$CPPFLAGS_save"
	LIBS="$LIBS_save"
	LIBS_save="$LIBS"
	CPPFLAGS_save="$CPPFLAGS"
fi
if test "$xml" = "true" ; then
	PKG_CHECK_MODULES(LIBXML2, libxml-2.0 >= 2.6, [ true ], [ false ])
	LIBS="$LIBXML2_LIBS $LIBS"
	CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
	AC_SEARCH_LIBS(xmlInitParser, xml2, , [AC_MSG_WARN(disabling xml support)]
		[xml="false"])
	CPPFLAGS="$CPPFLAGS_save"
	LIBS="$LIBS_save"
fi

if test "$frame" = "true" ; then
	LALFRAME_LIBS="${lalsuite_top_builddir}/lal/packages/framedata/src/liblalframe.la"
	export LALFRAME_LIBS

fi
if test "$metaio" = "true" ; then
	LALMETAIO_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalmetaio.la"
	export LALMETAIO_LIBS
fi
if test "$xml" = "true" ; then
	LALXML_LIBS="${lalsuite_top_builddir}/lal/packages/support/src/liblalxml.la"
	export LALXML_LIBS
fi

AC_CONFIG_SUBDIRS([lal])
AC_CONFIG_SUBDIRS([lalapps])

AC_OUTPUT
