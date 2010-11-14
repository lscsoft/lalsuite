#!/bin/bash

# $Id: eah_build2.sh,v 1.76 2010/07/19 17:12:34 bema Exp $

# simple failure
fail() {
    echo
    echo "*******************************************"
    echo "Automated installation has failed" 1>&2
    echo "*******************************************"
    echo

    if test -n "$LOGFILE" -a -r "$LOGFILE" ; then
	echo "Transcript of failure is in ${LOGFILE}"
	echo "Final fifteen lines are:"
	tail -15 "$LOGFILE"
    fi

    exit 1
} ## fail()

log_and_show() {
    echo "$@" >> "$LOGFILE"
    echo "$@" >&2
}

log_and_do() {
    echo "$@" >> "$LOGFILE"
    "$@" >> "$LOGFILE" 2>&1 || fail
}

log_and_dont_fail() {
    echo "$@" >> "$LOGFILE"
    "$@" >> "$LOGFILE" 2>&1
}

download() {
    echo curl "$1/$2 > $2" >> "$LOGFILE"
    curl "$1/$2" > "$2" 2>> "$LOGFILE"
}

eah_build2_loc="`echo $PWD/$0 | sed 's%/[^/]*$%%'`"

boinc_rev=-r22561
#previous:-r22503 -r22363 -r21777 -r'{2008-12-01}'

for i; do
    case "$i" in
	--win32)
	    build_win32=true ;;
	--static)
	    shared_copt="--disable-shared" ;;
	--rebuild)
	    rebuild_boinc=true
	    rebuild_lal=true
	    rebuild=true ;;
	--rebuild-lal)
	    rebuild_lal=true ;;
	--rebuild-boinc)
	    rebuild_boinc=true ;;
	--rebuild-binutils)
	    rebuild_binutils=true ;;
	--release)
	    rebuild_binutils=true
	    rebuild_boinc=true
	    rebuild_lal=true
	    rebuild=true
	    release=true
	    CFLAGS="-O3 $CFLAGS"
	    LDFLAGS="-static-libgcc $LDFLAGS"
	    shared_copt="--disable-shared"  ;;
	--norebuild) # dangerous, for testing only!
	    rebuild_binutils=""
	    rebuild_boinc=""
	    rebuild_lal=""
	    rebuild="" ;;
	--64)
	    CPPFLAGS="-m64 $CPPFLAGS"
	    CXXFLAGS="-m64 $CXXFLAGS"
	    CFLAGS="-m64 $CFLAGS"
	    LDFLAGS="-m64 $LDFLAGS" ;;
	--sse)
	    CPPFLAGS="-DENABLE_SSE_EXCEPTIONS $CPPFLAGS"
	    CFLAGS="-msse -march=pentium3 $CFLAGS"
	    fftw_copts=--enable-sse
	    acc="_sse";;
	--sse2)
	    CPPFLAGS="-DENABLE_SSE_EXCEPTIONS $CPPFLAGS"
	    CFLAGS="-msse -msse2 -mfpmath=sse -march=pentium-m $CFLAGS"
	    fftw_copts="--enable-sse --enable-sse2"
	    acc="_sse2";;
	--altivec)
	    CPPFLAGS="-maltivec -faltivec $CPPFLAGS"
	    CFLAGS="-fast -mcpu=G4 -maltivec -faltivec $CFLAGS"
	    CXXFLAGS="-mcpu=G4 $CXXFLAGS"
	    fftw_copts=--enable-altivec
	    acc="_altivec";;
	--cuda)
	    cuda=true
	    acc="_cuda" ;;
	--panther)
	    export MACOSX_DEPLOYMENT_TARGET=10.3
	    export SDKROOT="/Developer/SDKs/MacOSX10.3.9.sdk"
	    pflags="-arch ppc -D_NONSTD_SOURCE -isystem $SDKROOT"
	    CPPFLAGS="$pflags $CPPFLAGS"
	    CFLAGS="$pflags -Wno-long-double $CFLAGS"
	    CXXFLAGS="$pflags -Wno-long-double $CXXFLAGS"
	    LDFLAGS="$pflags -Wl,-syslibroot,$SDKROOT $LDFLAGS"
	    export RELEASE_LDADD="./libstdc++.a"
	    export CC=gcc-3.3 CXX=g++-3.3 ;;
	--with-ssl=*)
	    WITH_SSL="$i"
	    ssldir=`echo "$i" | sed 's/--with-ssl=//'`
	    CPPFLAGS="$CPPFLAGS -I$ssldir/include"
	    LIBS="$LIBS -L$ssldir/lib" ;;
	--check)
	    check=true ;;
	--check-only)
	    check=true
	    check_only=true ;;
	--check-app=*)
	    check=true
	    check_only=true
	    check_app=`echo $PWD/$i | sed 's/--check-app=//;s%.*//%/%'`;;
	--boinc-rev=*)
	    boinc_rev="`echo $i | sed s/^--boinc-rev=//`";;
	--help)
	    echo "$0 builds Einstein@home Applications of LALApps HierarchicalSearch codes"
	    echo "  --win32           cros-compile a Win32 App (requires MinGW, target i586-mingw32msvc-gcc)"
	    echo "  --rebuild         build FFTW, gsl, BOINC and LAL from source even if they are found by pkg-config"
	    echo "  --rebuild-lal     rebuild lalsuite"
	    echo "  --rebuild-boinc   rebuild BOINC"
	    echo "  --static          try to link statically"
	    echo "  --64              build 64Bit (add -m64 to  CPPFLAGS, CXXFLAGS, CFLAGS and LDFLAGS)"
	    echo "  --panther         build to run on Mac OS 10.3.9"
	    echo "  --cuda            build an App that uses CUDA"
	    echo "  --sse             build an App that uses SSE"
	    echo "  --sse2            build an App that uses SSE2"
	    echo "  --altivec         build an App that uses AltiVec"
	    echo "  --boinc-rev=<rev> specify a BOINC SVN trunk revision to use"
	    echo "  --with-ssl=<path> gets paased to BOINC configure"
	    echo "  --check           test the newly built HierarchSearchGC App"
	    echo "  --check-only      only test the already built HierarchSearchGC App"
	    echo "  --check-app=<app> only test the app specified, not necessarily the one just built"
	    echo "  --release         use some dark magic to make the App most compatible and add remote debugging."
	    echo "                    Implies --static and --rebuild and even more dirty hacks on Linux to work on Woody"
	    echo "  --norebuild       disables --rebuild on --release. DANGEROUS! Use only for testing the build script"
	    echo "  --help            show this message and exit"
	    exit ;;
	*) echo "unknown option '$i', try $0 --help"; exit ;;
    esac
done

EAH="$PWD/EinsteinAtHome"
LOGFILE="$EAH/build.log"
SOURCE="$EAH/source"
BUILD="$EAH/build$acc"
INSTALL="$EAH/install$acc"

if [ ."$build_win32" = ."true" ] ; then
    export CC=i586-mingw32msvc-gcc
    export CXX=i586-mingw32msvc-g++
    export AR=i586-mingw32msvc-ar
    export RANLIB=i586-mingw32msvc-ranlib
    CPPFLAGS="-DMINGW_WIN32 -DWIN32 -D_WIN32 -D_WIN32_WINDOWS=0x0410 $CPPFLAGS"
    # -include $INSTALL/include/win32_hacks.h
    cross_copt=--host=i586-pc-mingw32
    ext=".exe"
    wine=`which wine`
    if [ ".$wine" = "." -a ".$check" = ".true" ]; then
        log_and_show "WARNING: 'wine' not found, disabling check as it won't work"
        check=false
    fi
    if [ ".$cuda" = ".true" ] ; then
	test ".$WINEPREFIX" = "." &&
	WINEPREFIX="$HOME/.wine"
	CPPFLAGS="-I$WINEPREFIX/drive_c/CUDA/include $CPPFLAGS"
	export CUDART="$WINEPREFIX/drive_c/CUDA/lib/cudart.lib"
	export NVCC="$PWD/nvcc-wine-wrapper.sh"
    elif [ ".$release" = ".true" ] ; then
	CPPFLAGS="-DHAVE_EXCHNDL -I$INSTALL/include/bfd $CPPFLAGS"
	CFLAGS="-gstabs3 $CFLAGS"
	CXXFLAGS="-gstabs3 $CXXFLAGS"
	export RELEASE_DEPS="exchndl.o"
	export RELEASE_LDADD="exchndl.o -lbfd -liberty -lintl"
	build_binutils=true
    fi
else
    case `uname -s` in
	Darwin)
	    LDFLAGS="-framework Carbon -framework AppKit -framework IOKit -framework CoreFoundation $LDFLAGS" ;;
	Linux)
	    LDFLAGS="-lpthread $LDFLAGS"
	    if [ x"$release" = x"true" ]; then
		CPPFLAGS="-DEXT_STACKTRACE -I$INSTALL/include/bfd $CPPFLAGS"
		export RELEASE_DEPS="erp_execinfo_plus.o libstdc++.a libz.a"
		export RELEASE_LDADD="erp_execinfo_plus.o -lbfd -liberty"
		build_binutils=true
		enable_linux_compatibility_workarounds=true
	    fi ;;
    esac
fi

# if --pather and not --altivec, make sure the binary runs on G3
test ."$MACOSX_DEPLOYMENT_TARGET" = ."10.3" -a ."$acc" = ."" &&
CFLAGS="-mcpu=G3 $CFLAGS" &&
CXXFLAGS="-mcpu=G3 $CXXFLAGS"

if [ ".$cuda" = ".true" -a ."$build_win32" = ."true" ]; then
    export CFLAGS="-g0 $CFLAGS"
else
    export CFLAGS="-g $CFLAGS"
fi
export CPPFLAGS="-DBOINC_APIV6 -D__NO_CTYPE -DUSE_BOINC -DEAH_BOINC -I$INSTALL/include $CPPFLAGS"
export LDFLAGS="-L$INSTALL/lib $LDFLAGS"
export LD_LIBRARY_PATH="$INSTALL/lib:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$INSTALL/lib:$DYLD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$INSTALL/lib/pkgconfig:$PKG_CONFIG_PATH"
export BOINC_PREFIX="$INSTALL"


# make sure the E@H directory exists (for logging)
mkdir -p "$EAH" || fail

echo " " >> "$LOGFILE"
log_and_show "==========================================="
log_and_show "Build start `date`"

echo "$0" "$@" >> "$LOGFILE"
echo CC="\"$CC\"" >> "$LOGFILE"
echo CXX="\"$CXX\"" >> "$LOGFILE"
echo AR="\"$AR\"" >> "$LOGFILE"
echo RANLIB="\"$RANLIB\"" >> "$LOGFILE"
echo CFLAGS="\"$CFLAGS\"" >> "$LOGFILE"
echo CPPFLAGS="\"$CPPFLAGS\"" >> "$LOGFILE"
echo LDFLAGS="\"$LDFLAGS\"" >> "$LOGFILE"
echo LD_LIBRARY_PATH="\"$LD_LIBRARY_PATH\"" >> "$LOGFILE"
echo DYLD_LIBRARY_PATH="\"$DYLD_LIBRARY_PATH\"" >> "$LOGFILE"
echo PKG_CONFIG_PATH="\"$PKG_CONFIG_PATH\"" >> "$LOGFILE"
echo BOINC_PREFIX="\"$BOINC_PREFIX\"" >> "$LOGFILE"
echo RELEASE_DEPS="\"$RELEASE_DEPS\"" >> "$LOGFILE"
echo RELEASE_LDADD="\"$RELEASE_LDADD\"" >> "$LOGFILE"

gsl=gsl-1.9
fftw=fftw-3.2.2
binutils=binutils-2.19

if ! [ .$check_only = .true ]; then

log_and_do rm -rf "$BUILD"
test -n "$rebuild" &&
  log_and_do rm -rf "$INSTALL"

log_and_do mkdir -p "$SOURCE" "$BUILD/$fftw" "$BUILD/$gsl" "$BUILD/lal" "$BUILD/lalpulsar" "$BUILD/lalapps" "$BUILD/boinc" "$INSTALL/include"

if [ ."$build_win32" = ."true" ] ; then
    echo '#define bzero(a,b) memset(a,0,b)' > "$INSTALL/include/win32_hacks.h"
    echo '#define index(s,c) strchr(s,c)'  >> "$INSTALL/include/win32_hacks.h"
fi



log_and_do cd "$SOURCE"

if test -z "$rebuild" && pkg-config --exists gsl; then
    log_and_show "using existing gsl"
else
    log_and_show "retrieving $gsl"
    download http://www.aei.mpg.de/~repr/EaH_packages $gsl.tar.gz
    log_and_do tar xzf "$gsl.tar.gz"
fi

if test -z "$rebuild" && pkg-config --exists fftw3 fftw3f; then
    log_and_show "using existing fftw"
else
    log_and_show "retrieving $fftw"
    download ftp://ftp.fftw.org/pub/fftw $fftw.tar.gz
    log_and_do tar xzf "$fftw.tar.gz"
fi

if test -n "$build_binutils" -a -n "$rebuild_binutils"; then
    log_and_show "retrieving $binutils"
#    download http://www.aei.mpg.de/~repr/EaH_packages $binutils.tar.gz
    download ftp://ftp.fu-berlin.de/unix/gnu/binutils $binutils.tar.gz
    log_and_do rm -rf "$binutils"
    log_and_do tar xzf "$binutils.tar.gz"
#    log_and_do sh -c "grep -v '^ *SUBDIRS *=' $binutils/bfd/Makefile.am > $binutils/bfd/Makefile.tmp"
#    log_and_do mv $binutils/bfd/Makefile.tmp $binutils/bfd/Makefile.am
fi

if test -z "$rebuild_boinc" -a -d "$SOURCE/boinc" ; then
    log_and_show "using existing boinc source"
else
    log_and_show "retrieving boinc"
    log_and_do svn co "$boinc_rev" http://boinc.berkeley.edu/svn/trunk/boinc
fi

if test \! -d lalsuite/.git ; then
    log_and_do rm -rf lalsuite
    log_and_do ln -s "$eah_build2_loc/../../../../../.." lalsuite
fi

if test -z "$rebuild" && pkg-config --exists fftw3 fftw3f; then
    log_and_show "using existing fftw"
else
    log_and_show "compiling fftw"
    log_and_do cd "$BUILD/$fftw"
    log_and_do "$SOURCE/$fftw/configure" $fftw_copts "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    log_and_do make
    log_and_do make install
    log_and_do "$SOURCE/$fftw/configure" $fftw_copts --enable-single "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    log_and_do make
    log_and_do make install
fi

if test -z "$rebuild" && pkg-config --exists gsl; then
    log_and_show "using existing gsl"
else
    log_and_show "compiling gsl"
    log_and_do cd "$BUILD/$gsl"
    log_and_do "$SOURCE/$gsl/configure" "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    log_and_do make
    log_and_do make install
fi

if test -n "$build_binutils"; then
  if test -z "$rebuild_binutils"; then
    log_and_show "using existing gsl"
  else
    log_and_show "compiling binutils"
    log_and_do mkdir -p "$BUILD/$binutils"
    log_and_do cd "$BUILD/$binutils"
    log_and_do "$SOURCE/$binutils/configure" "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    if [ ".$enable_linux_compatibility_workarounds" = ".true" ]; then
        log_and_dont_fail make -k
        log_and_dont_fail make -k install
    else
        log_and_do make
        log_and_do make install
    fi
    # some post-build installation due to targets missing in the library
    log_and_do cd "$SOURCE/$binutils"
    log_and_do mkdir -p "$INSTALL/include/bfd"
    log_and_do cp -r include/* bfd/*.h "$BUILD/$binutils/binutils/config.h" "$INSTALL/include/bfd"
    if [ ."$build_win32" = ."true" ] ; then
	log_and_do cd "$BUILD/$binutils"
	log_and_do cp "intl/libintl.a" "$INSTALL/lib"
        # patch a few headers
	( cd "$INSTALL/include/bfd" &&
	    patch -N -p0 <<EOF
diff -ur include.org/coff/internal.h include/coff/internal.h
--- coff/internal.h	2008-10-06 17:29:08.000000000 +0200
+++ coff/internal.h	2008-10-06 17:31:26.000000000 +0200
@@ -98,11 +98,6 @@
 #define F_DLL           (0x2000)
 
 /* Extra structure which is used in the optional header.  */
-typedef struct _IMAGE_DATA_DIRECTORY 
-{
-  bfd_vma VirtualAddress;
-  long    Size;
-}  IMAGE_DATA_DIRECTORY;
 #define PE_EXPORT_TABLE			0
 #define PE_IMPORT_TABLE			1
 #define PE_RESOURCE_TABLE		2
Only in include/coff: internal.h~
diff -ur include.org/libcoff.h include/libcoff.h
--- libcoff.h	2008-10-06 17:29:08.000000000 +0200
+++ libcoff.h	2008-10-06 17:33:28.000000000 +0200
@@ -256,8 +256,10 @@
   /* Symbol type.  */
   unsigned short type;
 
+#ifndef __cplusplus
   /* Symbol class.  */
   unsigned char class;
+#endif
 
   /* Number of auxiliary entries.  */
   char numaux;
@@ -396,8 +398,10 @@
   /* Next type with the same name.  */
   struct coff_debug_merge_type *next;
 
+#ifndef __cplusplus
   /* Class of type.  */
   int class;
+#endif
 
   /* Symbol index where this type is defined.  */
   long indx;
Only in include: libcoff.h~
EOF
)
    fi
  fi
fi

if test -z "$rebuild_boinc" && test -r "$INSTALL/lib/libboinc_api.a" ; then
    log_and_show "using installed boinc"
else
    log_and_show "compiling boinc"
    if [ ."$build_win32" = ."true" ] ; then
	log_and_do cd "$BUILD/boinc"
	export BOINC_SRC="$SOURCE/boinc" BOINC_PREFIX="$INSTALL"
	log_and_dont_fail make -f "$SOURCE/boinc/lib/Makefile.mingw" uninstall
	log_and_do make -f "$SOURCE/boinc/lib/Makefile.mingw" clean
	log_and_do make -f "$SOURCE/boinc/lib/Makefile.mingw"
	log_and_do make -f "$SOURCE/boinc/lib/Makefile.mingw" install
    else
	log_and_do cd "$SOURCE/boinc"
	log_and_do ./_autosetup
	log_and_do cd "$BUILD/boinc"
	log_and_do "$SOURCE/boinc/configure" --disable-server --disable-manager --disable-client "$WITH_SSL" "$shared_copt" "$cross_copt" --prefix="$INSTALL" # --target=powerpc-apple-darwin7.9.0
	if [ .$MACOSX_DEPLOYMENT_TARGET = .10.3 -a -r "$SDKROOT/usr/lib/gcc/darwin/3.3/libstdc++.a" ]; then
	    log_and_do sed -i~ "s%-lstdc++%$SDKROOT/usr/lib/gcc/darwin/3.3/libstdc++.a%g" `find . -name Makefile`
	    log_and_do sed -i~ 's/#define HAVE_SYS_STATVFS_H 1/#undef HAVE_SYS_STATVFS_H/' config.h
	fi
	log_and_dont_fail make uninstall
	log_and_do make
	log_and_do make install
    fi
fi

if test -z "$rebuild_lal" && pkg-config --exists lal; then
    log_and_show "using existing lal"
else
    log_and_show "compiling LAL"
    log_and_do cd "$SOURCE/lalsuite/lal"
    if [ ".$MACOSX_DEPLOYMENT_TARGET" = ".10.3" ]; then
        log_and_do sed -i~ s/-mmacosx-version-min=10.4// configure.ac
    fi
    log_and_do ./00boot
    log_and_do cd "$BUILD/lal"
    log_and_do "$SOURCE/lalsuite/lal/configure" --disable-gcc-flags --disable-debug --enable-boinc --disable-silent-rules "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    log_and_do make
    log_and_do make install
fi

if test -z "$rebuild_lal" && pkg-config --exists lalpulsar; then
    log_and_show "using existing lalpulsar"
else
    log_and_show "compiling LALPulsar"
    log_and_do cd "$SOURCE/lalsuite/lalpulsar"
    if [ ".$MACOSX_DEPLOYMENT_TARGET" = ".10.3" ]; then
        log_and_do sed -i~ s/-mmacosx-version-min=10.4// configure.ac
    fi
    log_and_do ./00boot
    log_and_do cd "$BUILD/lalpulsar"
    log_and_do "$SOURCE/lalsuite/lalpulsar/configure" --disable-gcc-flags --disable-debug --enable-boinc --disable-silent-rules "$shared_copt" "$cross_copt" --prefix="$INSTALL"
    log_and_dont_fail make uninstall
    log_and_do make
    log_and_do make install
fi

# work around a bug in current LALApps build
#log_and_do cp -R "$SOURCE/lalsuite/lalmetaio/src/LIGOMetadataTables.h" "$INSTALL/include/lal"

log_and_show "configuring LALApps"
log_and_do cd "$SOURCE/lalsuite/lalapps"
if [ ."$MACOSX_DEPLOYMENT_TARGET" = ."10.3" ]; then
    log_and_do sed -i~ s/-mmacosx-version-min=10.4// configure.ac
    if grep ^eah_HierarchSearchGCT_manual: src/pulsar/GCT/Makefile.am >/dev/null; then
	:
    else
	echo '
## pretty dirty hack, particulary for compiling an E@H Mac OS 10.3 App on 10.4
eah_HierarchSearchGCT_manual: eah_HierarchSearchGCT-HierarchSearchGCT.o \
    eah_HierarchSearchGCT-GCTtoplist.o \
    eah_HierarchSearchGCT-HeapToplist.o \
    eah_HierarchSearchGCT-Fstat_v3.o \
    eah_HierarchSearchGCT-ComputeFstat_RS.o \
    eah_HierarchSearchGCT-hs_boinc_extras.o \
    eah_HierarchSearchGCT-hs_boinc_options.o \
    eah_HierarchSearchGCT-win_lib.o \
    eah_HierarchSearchGCT-LocalComputeFstat.o \
    libstdc++.a
	$(CXX) $(CPPFLAGS) $(AC_CPPFLAGS) -o $@ $^ $(LDFLAGS) $(AC_LDFLAGS) $(LDADD) $(AC_LDADD) $(LIBS) $(AC_LIBS) libstdc++.a
' >> src/pulsar/GCT/Makefile.am
    fi
fi

log_and_do ./00boot
if [ ."$build_win32" = ."true" ] ; then
    sed -i~ 's/test  *"${boinc}"  *=  *"true"/test "true" = "true"/' configure
fi
log_and_do cd "$BUILD/lalapps"
log_and_do "$SOURCE/lalsuite/lalapps/configure" --disable-gcc-flags --disable-debug --enable-boinc --disable-silent-rules "$shared_copt" "$cross_copt" --prefix="$INSTALL"

log_and_show "building Apps"

log_and_do cd "$BUILD/lalapps/src/lalapps"
# Windows may need a bit of a hack here, just collecting some infos in comments for now
if [ ."$build_win32" = ."true" ] ; then
    echo '/**/' > processtable.c
fi
if [ ! .$MACOSX_DEPLOYMENT_TARGET = .10.3 ] ; then
    log_and_do make LALAppsVCSInfo.h liblalapps.la
else
    log_and_do make LALAppsVCSInfo.h LALAppsVCSInfo.o lalapps.o
    log_and_do ar cru liblalapps.la lalapps.o LALAppsVCSInfo.o
fi

if [ ! .$MACOSX_DEPLOYMENT_TARGET = .10.3 ] ; then
    log_and_do cd "$BUILD/lalapps/src/pulsar/hough/src2"
    log_and_dont_fail make gitID
    if [ ".$cuda" = ".true" ] ; then
	log_and_do make "eah_HierarchicalSearch$acc$ext"
	log_and_do cp "eah_HierarchicalSearch$acc$ext" "$EAH/eah_HierarchicalSearch$acc$ext"
    else
	log_and_do make "eah_HierarchicalSearch$ext"
	log_and_do cp "eah_HierarchicalSearch$ext" "$EAH/eah_HierarchicalSearch$acc$ext"
    fi
fi

log_and_do cd "$BUILD/lalapps/src/pulsar/GCT"
log_and_dont_fail make gitID
if [  .$MACOSX_DEPLOYMENT_TARGET = .10.3 -a -r "$SDKROOT/usr/lib/gcc/darwin/3.3/libstdc++.a" ] ; then
  log_and_do rm -f libstdc++.a
  log_and_do ln -s "$SDKROOT/usr/lib/gcc/darwin/3.3/libstdc++.a"
  log_and_do make "eah_HierarchSearchGCT_manual"
  log_and_do cp "eah_HierarchSearchGCT_manual" "$EAH/eah_HierarchSearchGCT$acc"
else
  log_and_do make "eah_HierarchSearchGCT$ext"
  log_and_do cp "eah_HierarchSearchGCT$ext" "$EAH/eah_HierarchSearchGCT$acc$ext"
fi

if [ ! .$MACOSX_DEPLOYMENT_TARGET = .10.3 ] ; then
    log_and_do cd "$BUILD/lalapps/src/pulsar/Injections"
    log_and_do make eah_Makefakedata_v4$ext
    log_and_do cp eah_Makefakedata_v4$ext "$EAH"
    log_and_do cd "$BUILD/lalapps/src/pulsar/FDS_isolated"
    log_and_do make eah_PredictFStat$ext
    log_and_do cp eah_PredictFStat$ext "$EAH"
fi

log_and_show "==========================================="
log_and_show "Einstein@home Apps were built, should be in"
log_and_show "$EAH"
log_and_show "==========================================="

fi # check-only

if [ .$check = .true ]; then
    if [ -z "$check_app" ]; then
        check_app="../eah_HierarchSearchGCT$acc$ext"
    fi
    log_and_show "Running test"
    log_and_do cd "$EAH"
    log_and_do rm -rf test
    log_and_do mkdir test
    log_and_do cd test
    log_and_do cp ../eah_Makefakedata_v4$ext lalapps_Makefakedata_v4$ext
    log_and_do cp ../eah_PredictFStat$ext lalapps_PredictFStat$ext
    log_and_do cp "$INSTALL"/share/lalpulsar/*05-09.dat .
    PATH=".:$PATH" LAL_DATA_PATH="$PWD" log_and_do ../source/lalsuite/lalapps/src/pulsar/GCT/testHS.sh $wine "$check_app" --Dterms=8
    log_and_show "==========================================="
    log_and_show "Test passed"
    log_and_show "==========================================="
fi
