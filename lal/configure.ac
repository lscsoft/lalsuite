dnl Process this file with autoconf to produce a configure script.
AC_INIT([lal],[5.3],[lal-discuss@gravity.phys.uwm.edu])
AM_INIT_AUTOMAKE([])
AM_CONFIG_HEADER([include/config.h] [include/lal/LALConfig.h])
AH_TOP([
#ifndef CONFIG_H__LAL
#define CONFIG_H__LAL
])
AH_TEMPLATE([LAL_SIZEOF_DOUBLE], [The number of bytes in a double.])
AH_TEMPLATE([LAL_SIZEOF_FLOAT], [The number of bytes in a float.])
AH_TEMPLATE([LAL_SIZEOF_INT], [The number of bytes in a int.])
AH_TEMPLATE([LAL_SIZEOF_LONG], [The number of bytes in a long.])
AH_TEMPLATE([LAL_SIZEOF_LONG_LONG], [The number of bytes in a long long.])
AH_TEMPLATE([LAL_SIZEOF_SHORT], [The number of bytes in a short.])
AH_BOTTOM([
#endif /* CONFIG_H__LAL */
])

AC_CANONICAL_HOST

LAL_WITH_CC
LAL_WITH_MPICC
LAL_WITH_CFLAGS
LAL_WITH_EXTRA_CFLAGS
LAL_WITH_EXTRA_CPPFLAGS
LAL_WITH_EXTRA_LDFLAGS
LAL_WITH_EXTRA_LIBS
LAL_WITH_GCC_FLAGS

LAL_ENABLE_INTELFFT
LAL_ENABLE_FRAME
LAL_ENABLE_METAIO
LAL_ENABLE_XML
LAL_ENABLE_MPI

LAL_ENABLE_DEBUG
LAL_ENABLE_MACROS
LAL_ENABLE_PTHREAD_LOCK
LAL_ENABLE_NIGHTLY

AC_PREFIX_DEFAULT([/opt/lscsoft/lal])

lal_version_major=`echo "$VERSION" | cut -d. -f1`
lal_version_minor=`echo "$VERSION" | cut -d. -f2`
lal_version_micro=`echo "$VERSION" | cut -d. -f3-`
test -z "$lal_version_micro" && lal_version_micro=0
lal_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
# This sed script is more complex than it needs to be in order to prevent CVS
# from corrupting it.
lal_cvs_tag=`echo '$Name$' | sed -e 's? \\$\$??' -e 's?^\\$Name: ??'`
AC_DEFINE_UNQUOTED(LAL_VERSION, "$VERSION", [LAL Version])
AC_DEFINE_UNQUOTED(LAL_VERSION_MAJOR, $lal_version_major, [LAL Version Major Number] )
AC_DEFINE_UNQUOTED(LAL_VERSION_MINOR, $lal_version_minor, [LAL Version Minor Number] )
AC_DEFINE_UNQUOTED(LAL_VERSION_MICRO, $lal_version_micro, [LAL Version Micro Number] )
AC_DEFINE_UNQUOTED(LAL_CONFIGURE_ARGS, "$ac_configure_args", [LAL Configure Arguments])
AC_DEFINE_UNQUOTED(LAL_CONFIGURE_DATE, "$lal_configure_date", [LAL Configure Date])
AC_DEFINE_UNQUOTED(LAL_CVS_TAG, "$lal_cvs_tag", [LAL CVS Tag])
test "$prefix"="NONE" && lal_prefix="$ac_default_prefix" || lal_prefix="$prefix"
AC_DEFINE_UNQUOTED(LAL_PREFIX, "$lal_prefix", [LAL Install Prefix])

lal_pthread_available="yes"
ACX_PTHREAD( , [lal_pthread_available="no"])
if test "$lal_pthread_lock" = "true" ; then
  if test "$lal_pthread_available" = "yes" ; then
    LIBS="$PTHREAD_LIBS $LIBS"
    CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC"
  else
    AC_MSG_ERROR([don't know how to compile posix threads])
  fi
fi

dnl Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AM_PROG_LIBTOOL
AC_CHECK_PROGS(M4, gm4 m4, m4)
AC_CHECK_PROGS(LATEX, pdflatex latex, echo)
AC_CHECK_PROGS(MKIND, makeindex, echo)
AC_CHECK_PROGS(DVIPS, dvips, echo)
if test -z "$TEXFLAGS"
then
AC_MSG_CHECKING(if $LATEX accepts --interaction=batchmode)
cat > conflatextest.tex <<\CEOF
\documentclass{article}
\begin{document}
\end{document}
CEOF
if $LATEX conflatextest.tex --interaction=batchmode >/dev/null 2>&1
then
TEXFLAGS="--interaction=batchmode"
AC_MSG_RESULT(yes)
else
TEXFLAGS=""
AC_MSG_RESULT(no)
fi
rm -f conflatextest.*
fi
AC_SUBST(TEXFLAGS)

dnl Set needed compiler flags
AC_MSG_CHECKING([for additional compiler flags])
case "${host_os}" in
  solaris*) CFLAGS="$CFLAGS -D_POSIX_C_SOURCE=200101L"
  AC_MSG_RESULT([-D_POSIX_SOURCE=200101L]);;
  *) AC_MSG_RESULT(none);;
esac

dnl Checks for libraries.
LALSUPPORT_LIBADD=""
lal_save_LIBS="$LIBS"
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(z, compress, [ LALSUPPORT_LIBADD="$LALSUPPORT_LIBADD -lz" ])
dnl if solaris, check for sunmath library
case "${host_os}" in
  solaris*) AC_CHECK_LIB(sunmath, sincosp);;
esac
PKG_CHECK_MODULES(GSL, gsl, [ true ], [ false ])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$LIBS $GSL_LIBS"
AC_CHECK_LIB(gslcblas, main)
AC_CHECK_LIB(gsl, gsl_strerror, , [AC_MSG_ERROR(could not find GSL library)])
LAL_CHECK_GSL_VERSION(1.9)

dnl Checks for fft libraries.
if test "${intelfft}" = "false" ; then
  fftw3="true"
  PKG_CHECK_MODULES(FFTW3, fftw3 fftw3f, [ true ], [ false ])
  CPPFLAGS="$CPPFLAGS $FFTW3_CFLAGS"
  LIBS="$LIBS $FFTW3_LIBS"
  AC_CHECK_LIB(fftw3f, fftwf_execute_dft, , [AC_MSG_ERROR(could not find fftw3f library)], -lm)
  AC_CHECK_LIB(fftw3, fftw_execute_dft, , [AC_MSG_ERROR(could not find fftw3 library)], -lm)
else
  AC_MSG_WARN(Using Intel FFT routines)
  if test "x${qthread}" = "xtrue" ; then
    LAL_INTEL_MKL_QTHREAD_WARNING
  else
    AC_CHECK_LIB(pthread, pthread_create, , [AC_MSG_ERROR(could not find pthread library)])
  fi
  AC_CHECK_LIB(guide, __kmp_abort, , [AC_MSG_ERROR(could not find guide library included with Intel MKL)])
  AC_CHECK_HEADERS(pthread.h, , [AC_MSG_ERROR(could not find pthread header)])
  if test "${enable_static}" = "yes"; then
    if test "${enable_shared}" = "yes" ; then
      LAL_INTEL_FFT_LIBS_MSG_ERROR
    else
      AC_CHECK_LIB(mkl_ia32, DftiCreateDescriptor, , [AC_MSG_ERROR(could not find the Intel FFT library)])
    fi  
  else
    AC_CHECK_LIB(mkl, DftiCreateDescriptor, , [AC_MSG_ERROR(could not find the Intel FFT library)])
  fi
  AC_CHECK_HEADERS(mkl_fft.h, , [AC_MSG_ERROR(could not find Intel FFT header)])
fi

dnl Checks for optional FrameL, metaio, and xml2 libraries.
if test "${metaio}" = "true"; then
  PKG_CHECK_MODULES(METAIO, libmetaio, [ true ], [ false ])
  lal_pre_metaio_LIBS="$LIBS"
  LIBS="$LIBS $METAIO_LIBS"
  AC_SEARCH_LIBS(MetaioOpen, metaio dataflow, ,
                 [AC_MSG_WARN(disabling metaio support)]
                 [metaio="false"]
                 [LIBS="$lal_pre_metaio_LIBS"])
fi
if test "${frame}" = "true"; then
  PKG_CHECK_MODULES(FRAME, libframe, [ true ], [ false ])
  lal_pre_frame_LIBS="$LIBS"
  LIBS="$LIBS $FRAME_LIBS"
  AC_SEARCH_LIBS(FrLibIni, Frame, ,
                 [AC_MSG_WARN(disabling frame library)]
                 [frame="false"]
                 [LIBS="$lal_pre_frame_LIBS"])
fi
if test "${xml}" = "true"; then
  PKG_CHECK_MODULES(LIBXML2, libxml-2.0 >= 2.6, [ true ], [ false ])
  lal_pre_xml2_LIBS="$LIBS"
  LIBS="$LIBS $LIBXML2_LIBS"
  AC_SEARCH_LIBS(xmlInitParser, xml2, ,
                 [AC_MSG_WARN(disabling xml support)]
                 [xml="false"]
                 [LIBS="$lal_pre_xml2_LIBS"])
fi


dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(sys/time.h unistd.h getopt.h malloc.h zlib.h regex.h glob.h)
case "${host_os}" in
  solaris*) AC_CHECK_HEADERS(sunmath.h);;
esac
AC_CHECK_HEADERS(gsl/gsl_errno.h, , [AC_MSG_ERROR(could not find gsl/gsl_errno.h)])
if test "${intelfft}" = "false" ; then
  AC_CHECK_HEADERS(fftw3.h, , [AC_MSG_ERROR(could not find fftw3.h)])
fi

AC_MSG_CHECKING(for gethostname prototype in unistd.h)
AC_EGREP_HEADER(gethostname, unistd.h,AC_MSG_RESULT(yes)
AC_DEFINE(HAVE_GETHOSTNAME_PROTOTYPE, 1, [Define if gethostname prototype is in unistd.h]), AC_MSG_RESULT(no))

if test "${metaio}" = "true"; then
  lal_pre_metaio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
  AC_CHECK_HEADERS(metaio.h, ,
                   [AC_MSG_WARN(disabling metaio support)]
                   [metaio="false"]
                   [CPPFLAGS="$lal_pre_metaio_CPPFLAGS"])
fi
if test "${frame}" = "true"; then
  lal_pre_frame_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $FRAME_CFLAGS"
  AC_CHECK_HEADERS(FrameL.h, ,
                   [AC_MSG_WARN(disabling frame library support)]
                   [frame="false"]
                   [CPPFLAGS="$lal_pre_frame_CPPFLAGS"])
fi
if test "${xml}" = "true"; then
  lal_pre_xml2_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
  AC_CHECK_HEADERS(libxml/xmlversion.h, ,
                   [AC_MSG_WARN(disabling xml support)]
                   [xml="false"]
                   [CPPFLAGS="$lal_pre_xml2_CPPFLAGS"])
fi
if test "${lal_pthread_lock}" = "true"; then
  AC_CHECK_HEADERS(pthread.h, break)
fi


dnl Checks for optional MPI libraries and header files.
if test "${mpi}" = "true"; then
  AC_CHECK_PROGS(MPICC, mpicc hcc, $CC)
  LAL_CHECK_MPI_FLAGS
  lal_pre_mpi_LIBS="$LIBS"
  lal_pre_mpi_CPPFLAGS="$CPPFLAGS"
  lal_pre_mpi_CFLAGS="$CFLAGS"
  LIBS="$LIBS $MPI_LDFLAGS"
  CPPFLAGS="$CPPFLAGS $MPI_CPPFLAGS"
  CFLAGS="$CFLAGS $MPI_CFLAGS"
  AC_SEARCH_LIBS(MPI_Init, mpi, , [AC_MSG_ERROR([can't find MPI library])])
  AC_CHECK_HEADERS(mpi.h, , [AC_MSG_ERROR([can't find mpi.h header])])
  LIBS="$lal_pre_mpi_LIBS"
  CPPFLAGS="$lal_pre_mpi_CPPFLAGS"
  CFLAGS="$lal_pre_mpi_CFLAGS"
else
  MPI_CFLAGS=
  MPI_CPPFLAGS=
  MPI_LDFLAGS=
fi


dnl Define variables describing what was enabled.
if test "${fftw3}" = "true"; then
  AC_DEFINE(LAL_FFTW3_ENABLED, 1, [Define if using fftw3 library])
fi
if test "${metaio}" = "true"; then
  AC_DEFINE(LAL_METAIO_ENABLED, 1, [Define if using MetaIo library])
fi
if test "${frame}" = "true"; then
  AC_DEFINE(LAL_FRAME_ENABLED, 1, [Define if using frame library])
fi
if test "${xml}" = "true"; then
  AC_DEFINE(LAL_XML_ENABLED, 1, [Define if using xml2 library])
fi
if test "${mpi}" = "true"; then
  AC_DEFINE(LAL_MPI_ENABLED, 1, [Define if using MPI])
fi


dnl Record the libraries needed for linking LAL.
LALLIBADD="$LIBS"
LIBS="$lal_save_LIBS"
AC_SUBST(LALLIBADD)
LAL_LIBADD="$LIBS $LALLIBADD"
AC_SUBST(LAL_LIBADD)
AC_SUBST(LALSUPPORT_LIBADD)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(LAL_MPI_ENABLED)


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_STRINGIZE
AC_TYPE_SIZE_T
AC_HEADER_TIME
LAL_CHECK_SIZEOF(short, 2)
LAL_CHECK_SIZEOF(int, 4)
LAL_CHECK_SIZEOF(long, 4)
LAL_CHECK_SIZEOF(long long, 8)
LAL_CHECK_SIZEOF(float, 4)
LAL_CHECK_SIZEOF(double, 8)

dnl Checks for library functions.
AC_CHECK_FUNCS(getopt)
AC_CHECK_FUNCS(getopt_long)
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(vsnprintf)
AC_CHECK_FUNCS(putenv)

dnl Set package documentation install path
pkgdocdir="\${prefix}/share/doc/${PACKAGE}-${VERSION}"
AC_SUBST(pkgdocdir)

dnl Setup automake conditionals
AM_CONDITIONAL(FRAME, test x$frame = xtrue)
AM_CONDITIONAL(MPI, test x$mpi = xtrue)
AM_CONDITIONAL(FRAMEANDMPI, test x$frame = xtrue -a x$mpi = xtrue)
AM_CONDITIONAL(METAIO, test x$metaio = xtrue)
AM_CONDITIONAL(XML, test x$xml = xtrue)
AM_CONDITIONAL(INTELFFT, test x$intelfft = xtrue)
AM_CONDITIONAL(QTHREAD, test x$qthread = xtrue)


dnl Configure subdirectory lalapps if present
if test -d $srcdir/lalapps ; then
  lal_build_dir='$(top_builddir)/..'
  LAL_BUILD="true"
  LAL_LIBS=$lal_build_dir/lib/liblal.la
  LAL_CFLAGS=-I$lal_build_dir/include
  LALSUPPORT_LIBS=$lal_build_dir/packages/support/src/liblalsupport.la
  LALSUPPORT_CFLAGS=-I$lal_build_dir/include
  if test x$frame = xtrue ; then
    LALFRAME_LIBS=$lal_build_dir/packages/framedata/src/liblalframe.la
    LALFRAME_CFLAGS=-I$lal_build_dir/include
  else
    LALFRAME_LIBS=""
    LALFRAME_CFLAGS=""
  fi
  if test x$metaio = xtrue ; then
    LALMETAIO_LIBS=$lal_build_dir/packages/support/src/liblalmetaio.la
    LALMETAIO_CFLAGS=-I$lal_build_dir/include
  else
    LALMETAIO_LIBS=""
    LALMETAIO_CFLAGS=""
  fi
  if test x$xml = xtrue ; then
    LALXML_LIBS=$lal_build_dir/packages/support/src/liblalxml.la
    LALXML_CFLAGS=-I$lal_build_dir/include
  else
    LALXML_LIBS=""
    LALXML_CFLAGS=""
  fi
  export LAL_BUILD LAL_LIBS LAL_CFLAGS LALSUPPORT_LIBS LALSUPPORT_CFLAGS LALFRAME_LIBS LALFRAME_CFLAGS LALMETAIO_LIBS LALMETAIO_CFLAGS LALXML_LIBS LALXML_CFLAGS
  AC_CONFIG_SUBDIRS(lalapps)
fi


dnl Add gcc specific flags --- must be done at end or tests will break!
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS $lal_gcc_flags"
fi
LIBLISTBASE='$(top_builddir)/packages/std/src/libstd.la $(top_builddir)/packages/tools/src/libtools.la $(top_builddir)/packages/sample/src/libsample.la $(top_builddir)/packages/hello/src/libhello.la $(top_builddir)/packages/factories/src/libfactories.la $(top_builddir)/packages/vectorops/src/libvectorops.la $(top_builddir)/packages/utilities/src/libutilities.la $(top_builddir)/packages/stats/src/libstats.la $(top_builddir)/packages/inject/src/libinject.la $(top_builddir)/packages/date/src/libdate.la $(top_builddir)/packages/tdfilter/src/libtdfilter.la $(top_builddir)/packages/window/src/libwindow.la $(top_builddir)/packages/fft/src/libfft.la $(top_builddir)/packages/clremoval/src/libclremoval.la $(top_builddir)/packages/stochastic/src/libstochastic.la $(top_builddir)/packages/inspiral/src/libinspiral.la $(top_builddir)/packages/noisemodels/src/libnoisemodels.la $(top_builddir)/packages/bank/src/libbank.la $(top_builddir)/packages/ring/src/libring.la $(top_builddir)/packages/findchirp/src/libfindchirp.la $(top_builddir)/packages/pulsar/src/libpulsar.la $(top_builddir)/packages/houghpulsar/src/libhoughpulsar.la $(top_builddir)/packages/burstsearch/src/libburstsearch.la $(top_builddir)/packages/tracksearch/src/libtracksearch.la $(top_builddir)/packages/tfclusters/src/libtfclusters.la $(top_builddir)/packages/slopefilters/src/libslopefilters.la $(top_builddir)/packages/block/src/libblock.la $(top_builddir)/packages/timefreq/src/libtimefreq.la '
LIBLISTMPI='$(top_builddir)/packages/comm/src/libcomm.la '
LIBLISTFRAME='$(top_builddir)/packages/framedata/src/libframedata.la '
LIBLISTMPIFRAME=''
if test x$mpi = xfalse; then
  LIBLISTMPI=
  LIBLISTMPIFRAME=
fi
if test x$frame = xfalse; then
  LIBLISTFRAME=
  LIBLISTMPIFRAME=
fi
PKGLIST="std tools sample hello factories vectorops utilities stats inject date tdfilter window fft clremoval stochastic inspiral noisemodels bank ring findchirp pulsar houghpulsar burstsearch tracksearch tfclusters slopefilters block timefreq support comm framedata  "
AC_SUBST(PKGLIST)
AC_SUBST(LIBLISTBASE)
AC_SUBST(LIBLISTMPI)
AC_SUBST(LIBLISTFRAME)
AC_SUBST(LIBLISTMPIFRAME)
INCLUDEMKFILE=$srcdir/misc/include.mk
DOCMKFILE=$srcdir/misc/doc.mk
DVIMKFILE=$srcdir/misc/dvi.mk
AC_SUBST_FILE(INCLUDEMKFILE)
AC_SUBST_FILE(DOCMKFILE)
AC_SUBST_FILE(DVIMKFILE)
AC_OUTPUT_COMMANDS([chmod +x lal-config])
AC_CONFIG_FILES([lal.spec])
AC_CONFIG_FILES([lal.pc])
AC_CONFIG_FILES([lalframe.pc])
AC_CONFIG_FILES([lalmetaio.pc])
AC_CONFIG_FILES([lalxml.pc])
AC_CONFIG_FILES([lalmpi.pc])
AC_CONFIG_FILES([lalsupport.pc])
AC_CONFIG_FILES([lal-config])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([debian/Makefile])
AC_CONFIG_FILES([doc/Makefile])
AC_CONFIG_FILES([doc/laldoc/Makefile])
AC_CONFIG_FILES([doc/autodoc/Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([include/lal/Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([misc/Makefile])
AC_CONFIG_FILES([packages/Makefile])
AC_CONFIG_FILES([packages/std/Makefile])
AC_CONFIG_FILES([packages/std/doc/Makefile])
AC_CONFIG_FILES([packages/std/include/Makefile])
AC_CONFIG_FILES([packages/std/src/Makefile])
AC_CONFIG_FILES([packages/std/test/Makefile])
AC_CONFIG_FILES([packages/tools/Makefile])
AC_CONFIG_FILES([packages/tools/doc/Makefile])
AC_CONFIG_FILES([packages/tools/include/Makefile])
AC_CONFIG_FILES([packages/tools/src/Makefile])
AC_CONFIG_FILES([packages/tools/test/Makefile])
AC_CONFIG_FILES([packages/sample/Makefile])
AC_CONFIG_FILES([packages/sample/doc/Makefile])
AC_CONFIG_FILES([packages/sample/include/Makefile])
AC_CONFIG_FILES([packages/sample/src/Makefile])
AC_CONFIG_FILES([packages/sample/test/Makefile])
AC_CONFIG_FILES([packages/hello/Makefile])
AC_CONFIG_FILES([packages/hello/doc/Makefile])
AC_CONFIG_FILES([packages/hello/include/Makefile])
AC_CONFIG_FILES([packages/hello/src/Makefile])
AC_CONFIG_FILES([packages/hello/test/Makefile])
AC_CONFIG_FILES([packages/factories/Makefile])
AC_CONFIG_FILES([packages/factories/doc/Makefile])
AC_CONFIG_FILES([packages/factories/include/Makefile])
AC_CONFIG_FILES([packages/factories/src/Makefile])
AC_CONFIG_FILES([packages/factories/test/Makefile])
AC_CONFIG_FILES([packages/vectorops/Makefile])
AC_CONFIG_FILES([packages/vectorops/doc/Makefile])
AC_CONFIG_FILES([packages/vectorops/include/Makefile])
AC_CONFIG_FILES([packages/vectorops/src/Makefile])
AC_CONFIG_FILES([packages/vectorops/test/Makefile])
AC_CONFIG_FILES([packages/utilities/Makefile])
AC_CONFIG_FILES([packages/utilities/doc/Makefile])
AC_CONFIG_FILES([packages/utilities/include/Makefile])
AC_CONFIG_FILES([packages/utilities/src/Makefile])
AC_CONFIG_FILES([packages/utilities/test/Makefile])
AC_CONFIG_FILES([packages/stats/Makefile])
AC_CONFIG_FILES([packages/stats/doc/Makefile])
AC_CONFIG_FILES([packages/stats/include/Makefile])
AC_CONFIG_FILES([packages/stats/src/Makefile])
AC_CONFIG_FILES([packages/stats/test/Makefile])
AC_CONFIG_FILES([packages/inject/Makefile])
AC_CONFIG_FILES([packages/inject/doc/Makefile])
AC_CONFIG_FILES([packages/inject/include/Makefile])
AC_CONFIG_FILES([packages/inject/src/Makefile])
AC_CONFIG_FILES([packages/inject/test/Makefile])
AC_CONFIG_FILES([packages/date/Makefile])
AC_CONFIG_FILES([packages/date/doc/Makefile])
AC_CONFIG_FILES([packages/date/include/Makefile])
AC_CONFIG_FILES([packages/date/src/Makefile])
AC_CONFIG_FILES([packages/date/test/Makefile])
AC_CONFIG_FILES([packages/tdfilter/Makefile])
AC_CONFIG_FILES([packages/tdfilter/doc/Makefile])
AC_CONFIG_FILES([packages/tdfilter/include/Makefile])
AC_CONFIG_FILES([packages/tdfilter/src/Makefile])
AC_CONFIG_FILES([packages/tdfilter/test/Makefile])
AC_CONFIG_FILES([packages/window/Makefile])
AC_CONFIG_FILES([packages/window/doc/Makefile])
AC_CONFIG_FILES([packages/window/include/Makefile])
AC_CONFIG_FILES([packages/window/src/Makefile])
AC_CONFIG_FILES([packages/window/test/Makefile])
AC_CONFIG_FILES([packages/fft/Makefile])
AC_CONFIG_FILES([packages/fft/doc/Makefile])
AC_CONFIG_FILES([packages/fft/include/Makefile])
AC_CONFIG_FILES([packages/fft/src/Makefile])
AC_CONFIG_FILES([packages/fft/test/Makefile])
AC_CONFIG_FILES([packages/clremoval/Makefile])
AC_CONFIG_FILES([packages/clremoval/doc/Makefile])
AC_CONFIG_FILES([packages/clremoval/include/Makefile])
AC_CONFIG_FILES([packages/clremoval/src/Makefile])
AC_CONFIG_FILES([packages/clremoval/test/Makefile])
AC_CONFIG_FILES([packages/stochastic/Makefile])
AC_CONFIG_FILES([packages/stochastic/doc/Makefile])
AC_CONFIG_FILES([packages/stochastic/include/Makefile])
AC_CONFIG_FILES([packages/stochastic/src/Makefile])
AC_CONFIG_FILES([packages/stochastic/test/Makefile])
AC_CONFIG_FILES([packages/inspiral/Makefile])
AC_CONFIG_FILES([packages/inspiral/doc/Makefile])
AC_CONFIG_FILES([packages/inspiral/include/Makefile])
AC_CONFIG_FILES([packages/inspiral/src/Makefile])
AC_CONFIG_FILES([packages/inspiral/test/Makefile])
AC_CONFIG_FILES([packages/noisemodels/Makefile])
AC_CONFIG_FILES([packages/noisemodels/doc/Makefile])
AC_CONFIG_FILES([packages/noisemodels/include/Makefile])
AC_CONFIG_FILES([packages/noisemodels/src/Makefile])
AC_CONFIG_FILES([packages/noisemodels/test/Makefile])
AC_CONFIG_FILES([packages/bank/Makefile])
AC_CONFIG_FILES([packages/bank/doc/Makefile])
AC_CONFIG_FILES([packages/bank/include/Makefile])
AC_CONFIG_FILES([packages/bank/src/Makefile])
AC_CONFIG_FILES([packages/bank/test/Makefile])
AC_CONFIG_FILES([packages/ring/Makefile])
AC_CONFIG_FILES([packages/ring/doc/Makefile])
AC_CONFIG_FILES([packages/ring/include/Makefile])
AC_CONFIG_FILES([packages/ring/src/Makefile])
AC_CONFIG_FILES([packages/ring/test/Makefile])
AC_CONFIG_FILES([packages/findchirp/Makefile])
AC_CONFIG_FILES([packages/findchirp/doc/Makefile])
AC_CONFIG_FILES([packages/findchirp/include/Makefile])
AC_CONFIG_FILES([packages/findchirp/src/Makefile])
AC_CONFIG_FILES([packages/findchirp/test/Makefile])
AC_CONFIG_FILES([packages/pulsar/Makefile])
AC_CONFIG_FILES([packages/pulsar/doc/Makefile])
AC_CONFIG_FILES([packages/pulsar/include/Makefile])
AC_CONFIG_FILES([packages/pulsar/src/Makefile])
AC_CONFIG_FILES([packages/pulsar/test/Makefile])
AC_CONFIG_FILES([packages/houghpulsar/Makefile])
AC_CONFIG_FILES([packages/houghpulsar/doc/Makefile])
AC_CONFIG_FILES([packages/houghpulsar/include/Makefile])
AC_CONFIG_FILES([packages/houghpulsar/src/Makefile])
AC_CONFIG_FILES([packages/houghpulsar/test/Makefile])
AC_CONFIG_FILES([packages/burstsearch/Makefile])
AC_CONFIG_FILES([packages/burstsearch/doc/Makefile])
AC_CONFIG_FILES([packages/burstsearch/include/Makefile])
AC_CONFIG_FILES([packages/burstsearch/src/Makefile])
AC_CONFIG_FILES([packages/burstsearch/test/Makefile])
AC_CONFIG_FILES([packages/tracksearch/Makefile])
AC_CONFIG_FILES([packages/tracksearch/doc/Makefile])
AC_CONFIG_FILES([packages/tracksearch/include/Makefile])
AC_CONFIG_FILES([packages/tracksearch/src/Makefile])
AC_CONFIG_FILES([packages/tracksearch/test/Makefile])
AC_CONFIG_FILES([packages/tfclusters/Makefile])
AC_CONFIG_FILES([packages/tfclusters/doc/Makefile])
AC_CONFIG_FILES([packages/tfclusters/include/Makefile])
AC_CONFIG_FILES([packages/tfclusters/src/Makefile])
AC_CONFIG_FILES([packages/tfclusters/test/Makefile])
AC_CONFIG_FILES([packages/slopefilters/Makefile])
AC_CONFIG_FILES([packages/slopefilters/doc/Makefile])
AC_CONFIG_FILES([packages/slopefilters/include/Makefile])
AC_CONFIG_FILES([packages/slopefilters/src/Makefile])
AC_CONFIG_FILES([packages/slopefilters/test/Makefile])
AC_CONFIG_FILES([packages/block/Makefile])
AC_CONFIG_FILES([packages/block/doc/Makefile])
AC_CONFIG_FILES([packages/block/include/Makefile])
AC_CONFIG_FILES([packages/block/src/Makefile])
AC_CONFIG_FILES([packages/block/test/Makefile])
AC_CONFIG_FILES([packages/timefreq/Makefile])
AC_CONFIG_FILES([packages/timefreq/doc/Makefile])
AC_CONFIG_FILES([packages/timefreq/include/Makefile])
AC_CONFIG_FILES([packages/timefreq/src/Makefile])
AC_CONFIG_FILES([packages/timefreq/test/Makefile])
AC_CONFIG_FILES([packages/support/Makefile])
AC_CONFIG_FILES([packages/support/doc/Makefile])
AC_CONFIG_FILES([packages/support/include/Makefile])
AC_CONFIG_FILES([packages/support/src/Makefile])
AC_CONFIG_FILES([packages/support/test/Makefile])
AC_CONFIG_FILES([packages/comm/Makefile])
AC_CONFIG_FILES([packages/comm/doc/Makefile])
AC_CONFIG_FILES([packages/comm/include/Makefile])
AC_CONFIG_FILES([packages/comm/src/Makefile])
AC_CONFIG_FILES([packages/comm/test/Makefile])
AC_CONFIG_FILES([packages/framedata/Makefile])
AC_CONFIG_FILES([packages/framedata/doc/Makefile])
AC_CONFIG_FILES([packages/framedata/include/Makefile])
AC_CONFIG_FILES([packages/framedata/src/Makefile])
AC_CONFIG_FILES([packages/framedata/test/Makefile])
AC_OUTPUT
frameenabled="`test x${frame} = xtrue && echo "ENABLED" || echo "DISABLED"`"
metaioenabled="`test x${metaio} = xtrue && echo "ENABLED" || echo "DISABLED"`"
xmlenabled="`test x${xml} = xtrue && echo "ENABLED" || echo "DISABLED"`"
mpienabled="`test x${mpi} = xtrue && echo "ENABLED" || echo "DISABLED"`"
echo "
================================================================

        LAL has now been successfully configured:

                LAL Frame library support is $frameenabled
                LAL MetaIo library support is $metaioenabled
                LAL XML library support is $xmlenabled
                LAL MPI support is $mpienabled

        Now run 'make' to build LAL
 	and run 'make install' to install LAL
        LAL will be installed under the directory:

                $prefix

================================================================="
