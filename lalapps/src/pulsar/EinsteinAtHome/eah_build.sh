#!/bin/bash


## function definitions follow here, the 'main-code' starts at MAIN

## ---------- some portability checks/adjustments [stolen from configure] ----------
## 'echo -n' is not portable..
case `echo "testing\c"; echo 1,2,3`,`echo -n testing; echo 1,2,3` in
  *c*,-n*) ECHO_N= ECHO_C='
' ECHO_T='	' ;;
  *c*,*  ) ECHO_N=-n ECHO_C= ECHO_T= ;;
  *)       ECHO_N= ECHO_C='\c' ECHO_T= ;;
esac
##----------


# simple failure
fail() 
{
    echo
    echo "********************************************************************************"
    echo "Automated installation has failed" 1>&2
    echo "********************************************************************************"

    if test -n "$LOGFILE" ; then
	echo "Transcript of failure is in ${LOGFILE}"
	echo "Final fifteen lines are:"
	tail -15 "$LOGFILE"
    fi

    log_close

    exit 1
} ## fail()

log_popup()
{
    if [ "${eah_show_log}" = X ]; then 
	xterm ${XTERMGEOM} -T `basename "$LOGFILE"`\@`hostname` -e tail -f "$LOGFILE" >& /dev/null& 
    elif [ "${eah_show_log}" = rxvt ]; then 
	rxvt -T `basename "$LOGFILE"`\@`hostname` -e tail -f "$LOGFILE" >& /dev/null& 
    elif [ "${eah_show_log}" = mac ]; then 
        rm -rf /tmp/follow /tmp/follow.command
	echo '#!/bin/sh
	      tail -f /tmp/follow & pid=$! ; while [ -r /tmp/follow ] ; do sleep 1 ; done ; kill $pid' > /tmp/follow.command
	chmod +x /tmp/follow.command
	ln -s "$LOGFILE" /tmp/follow
	open ${eah_mac_term}
    fi
} # log_popup()

log_start()
{
    LOGFILE="$1"
    rm -f "$LOGFILE"
    touch "$LOGFILE"

    log_popup
} # log_start()

log_close()
{
    if [ "${eah_show_log}" = X -o "${eah_show_log}" = rxvt ]; then 
	sleep 1
	kill $! >& /dev/null
    elif [ "${eah_show_log}" = mac ]; then 
        rm -f /tmp/follow
	sleep 2
    fi
} # log_close()

log_and_do()
{
    echo "$@" >> "$LOGFILE"
    "$@" >> "$LOGFILE" 2>&1 || fail
}

log_and_dont_fail()
{
    echo "$@" >> "$LOGFILE"
    "$@" >> "$LOGFILE" 2>&1
}

## ----------------------------------------------------------------------
## Check that given command $1 has version >= $2.$3
## return 0 if ok, 1 too old or not found  (-> shell conventions).
## ----------------------------------------------------------------------
check_version()
{
    local foundit
    local version_ok

    ## get current version of $1
    echo $ECHO_N "Checking version of '$1' >= $2.$3... $ECHO_C"
    fullpath=`type -p $1`;
    if [ -x "$fullpath" ]; then
	foundit=yes;
    fi

    if [ "$foundit" != yes ]; then 
	echo "Didn't find application";
	ac_major=0; ac_minor=0; 
    else
	cmdline="$fullpath --version 2> /dev/null";
	version=`$cmdline`;
	if [ -n "${version}" ]; then
	    ## HACK: the pre-pended space "x " is important for the following regexp 
	    ## to always work even if the version-number comes back without text prepended
	    version="x ${version}"
	else
	    version=0
	fi
	ac_major=`echo $version | sed 's/.* v*\([0-9]*\)[.]\([0-9]*\).*/\1/g'`
	ac_minor=`echo $version | sed 's/.* v*\([0-9]*\)[.]\([0-9]*\).*/\2/g'`
    fi

    if [ "${ac_major}" -gt $2 ]; then
	version_ok=yes
    elif [ "${ac_major}" -eq $2 ] && [ "${ac_minor}" -ge $3 ]; then
	version_ok=yes
    else
	version_ok=no
    fi

    if [ "${version_ok}" = "yes" ]; then
	echo "ok. (${ac_major}.${ac_minor})"
	return 0;
    else
	echo "failed. (${ac_major}.${ac_minor})"
	return 1;
    fi

} ## check_version()

##----------------------------------------------------------------------
## wrapper to wget: can use either wget, curl or lynx with same cmd-syntax
## $1: url 
## $2: pkg-name
##----------------------------------------------------------------------
wget_url()
{
    ## if target-package is found locally, remove it first
    if [ -n "${eah_WGET}" -o -n "${eah_CURL}" ]; then
	if [ -r "$2" ]; then
	    rm $2
	fi
    fi

    ## use wget or curl/lynx according to their download-syntax
    if [ -n "${eah_WGET}" ]; then
	echo $ECHO_N "Download package '$2' using '${eah_WGET}' ... $ECHO_C"  2>&1 | tee -a "$LOGFILE" 
	${eah_WGET} $1/$2 >> "$LOGFILE" 2>&1 || fail
    elif [ -n "${eah_CURL}" ]; then
	echo $ECHO_N "Download package '$2' using '${eah_CURL}' ... $ECHO_C "  2>&1 | tee -a "$LOGFILE" 
	${eah_CURL} $1/$2 > $2 2> "$LOGFILE" || fail
    else
	echo "Error... no usable download-program found (wget,curl,lynx)"	 2>&1 | tee -a "$LOGFILE"
	fail
    fi
    echo "done."	 2>&1 | tee -a "$LOGFILE"


} ## wget_url

## ----------------------------------------------------------------------
## Helper function: update/install an 'auxiliary' package (if necessary)
## 
## The package $1.tar.gz is downloaded from $EAH_PKGS_URL and installed into
## $AUX_INSTALL. The build takes place in $AUX_SOURCES, and a file
## $AUX_BUILD/$1.built is created to signal a successful installation.
##
## the build-process can be controlled with configure options "$2 $3 $4..."
##
## ----------------------------------------------------------------------
build_lsc_aux()
{
    if [ -z "$1" ]; then
	echo "Failure: build_aux() called without package-name"	 2>&1 | tee -a "$LOGFILE"
	fail
    fi
    aux_name=$1
    aux_tgz=${aux_name}.tar.gz
    aux_built=${AUX_INSTALL}/${aux_name}.built

    cd ${AUX_SOURCES} >> "$LOGFILE" 2>&1 || fail 

    ## if --rebuild-aux was given: download & rebuild this package from scratch
    if [ "${eah_rebuild_aux}" = "yes" ]; then
	if [ -r "${aux_built}" ]; then rm -f ${aux_built}; fi
    fi

    if [ ! -r "${aux_built}" ]; then
	# some special packages
	if [ "${aux_name}" = "libtool-2.2" ]; then
	    wget_url ftp://ftp.gnu.org/gnu/libtool "${aux_name}.tar.gz"
	elif [ "${aux_name}" = "fftw-3.1.2" ]; then
	    wget_url ftp://ftp.fftw.org/pub/fftw "${aux_name}.tar.gz"
	elif [ "${aux_name}" = "fftw-3.1.2_f" ]; then
	    cp fftw-3.1.2.tar.gz "${aux_name}.tar.gz"
	else
	    wget_url ${EAH_PKGS_URL} ${aux_tgz}
	fi

	echo $ECHO_N "Unpacking into '`basename ${AUX_SOURCES}`'... $ECHO_C"	 2>&1 | tee -a "$LOGFILE"
	cd ${AUX_SOURCES}
	cat ${aux_tgz} | gunzip | tar xf - >> "$LOGFILE" 2>&1  || fail
	echo "done."  	 2>&1 | tee -a "$LOGFILE"

	echo $ECHO_N "Building ${aux_name}... $ECHO_C"	 2>&1 | tee -a "$LOGFILE"
	cd ${AUX_SOURCES}/${aux_name} || fail
	shift 
	if [ -x ./configure ]; then
	    if [ "$eah_host_compile" = "yes" ]; then
		eah_configure="./configure"
	    elif [ "$eah_64bit" = "yes" -o "$eah_32bit" = "yes" ]; then
		eah_configure="./configure CPPFLAGS='$CPPGLAGS' CFLAGS='$CFLAGS' CXXFLAGS='$CXXFLAGS' LDFLAGS='$LDFLAGS' LIBS='$LIBS'"
	    else
		eah_configure="./configure"
	    fi
	else
	    echo "Failed: no configure-script found! "	 2>&1 | tee -a "$LOGFILE"
	    fail
	fi
	eah_next="${eah_configure} --prefix=${AUX_INSTALL} $@"
	echo ${eah_next} >> "$LOGFILE" 2>&1
	eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
	make >> "$LOGFILE" 2>&1 || fail
	make install >> "$LOGFILE" 2>&1 || fail
	touch ${aux_built}
	echo "done."	 2>&1 | tee -a "$LOGFILE"
    else
	echo
	echo "Using previously built version of ${aux_name} in `basename ${AUX_SOURCES}`." 2>&1 | tee -a "$LOGFILE"
	echo "Use '--rebuild-aux' to force rebuild of auxiliary packages" 2>&1 | tee -a "$LOGFILE"
	echo
    fi

} ## build_lsc_aux()


##----------------------------------------------------------------------
## checkout/update the given package-CVS from given cvs-server and 
## put it in the right places. Don't do update if eah_no_update=yes
## 
## you can optionally specify a list of modules (directories), in which case
## *ONLY* these directories will be checked out/updated (no subdirs!)
##
## $1: cvs-server URL
## $2: package-name
## $3: cvs-modules to get (default: all)
## $4: additional cvs update-options
## ----------------------------------------------------------------------
get_pkg_cvs()
{
    local fresh_cvs=no
    local pkg_name=$2
    local pkg_cvs=${pkg_name}-CVS
    local url=$1
    local modules="$3"
    local cvs_co_opts
    local cvs_update_opts

    if [ -z "${modules}" ]; then 
	co_modules=${pkg_name}
	cvs_co_opts="$4"
	cvs_update_opts="-d -P $4"
    else
	cvs_co_opts="-l $4"
	cvs_update_opts="-d -l $4"
	co_modules=
	for i in ${modules}; do
	    co_modules="${co_modules} ${pkg_name}/$i"
	done
    fi

    ## have we got the cvs-tree already in source- and build-directories?
    if [ ! -d ${BUILD_LOCATION}/extra_sources/${pkg_cvs} ]; then
	if [ ! -d ${SOURCE_LOCATION}/${pkg_cvs} ]; then
	    echo $ECHO_N "Initial checkout of ${pkg_cvs}... $ECHO_C"
	    cd ${SOURCE_LOCATION} >> "$LOGFILE" 2>&1 || fail
	    next_cmd="cvs -z3 -d ${url} co ${cvs_co_opts} ${co_modules}"
	    echo ${next_cmd} >> "$LOGFILE" 2>&1 
	    eval ${next_cmd} >> "$LOGFILE" 2>&1 || fail
	    mv -f ${pkg_name} ${pkg_cvs} >> "$LOGFILE" 2>&1 || fail
	    fresh_cvs=yes
	    echo "done."
	fi
	echo $ECHO_N "Copying ${pkg_cvs} into '$BUILD_LOCATION/extra_sources'... $ECHO_C"
	${COPYRECURSIVE} ${SOURCE_LOCATION}/${pkg_cvs} ${BUILD_LOCATION}/extra_sources >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi

    ## update cvs-tree in build-directory, unless it's 'fresh' or eah_no_update=yes
    if [ "$fresh_cvs" != yes ] && [ "${eah_no_update}" != yes ]; then
	echo $ECHO_N "Updating ${pkg_cvs} tree in '${BUILD_LOCATION}/extra_sources'... $ECHO_C"
	cd ${BUILD_LOCATION}/extra_sources/${pkg_cvs} >> "$LOGFILE" 2>&1 || fail
	next_cmd="cvs -z3 update ${cvs_update_opts} ${modules}"
	echo ${next_cmd} >> "$LOGFILE" 2>&1
	eval ${next_cmd} >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi

} ## get_pkg_cvs()


##----------------------------------------------------------------------
## download and unpack the given .tar.gz package from given download-server
## (this gets downloaded directly into BUILD_LOCATION/extra_sources)
## 
## $1: download-server URL
## $2: package-name
## ----------------------------------------------------------------------
get_pkg_tgz()
{
    local url=$1
    local pkg_name=$2
    local pkg_tgz=${pkg_name}.tar.gz

    ## NOTE: in order not to miss upgrades of the tar-packages, 
    ## we *ALWAYS* download and unpack the tarball (directly into BUILD_LOCATION),
    ## because we don't have any version-information
    echo $ECHO_N "Downloading ${pkg_tgz}... $ECHO_C"
    cd ${BUILD_LOCATION}/extra_sources >> "$LOGFILE" 2>&1 || fail
    wget_url $url $pkg_tgz
    echo "done."
	
    echo $ECHO_N "Unpacking ${pkg_tgz} into $BUILD_LOCATION/extra_sources... $ECHO_C"
    cat ${aux_tgz} | gunzip | tar xf - >> "$LOGFILE" 2>&1  || fail
    echo "done."
    
} ## get_pkg_tgz()


##----------------------------------------------------------------------
## Bootstrap the build-environment:
## try to figure out if the build environment is sufficiently new and 
## complete, or if we need to download and build something.
##
## STEP 1: build 0th order prerequisites: automake,autoconf,libtool, pkgconfig
##----------------------------------------------------------------------
step1() 
{
    echo "----------------------------------------------------------------------"
    echo "STEP 1: autoconf, automake, libtool, pkgconfig"
    echo "----------------------------------------------------------------------"
    LOGFILE=${BUILD_LOCATION}/step1.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"
    
    log_popup

    ## this checks for compilation tools on the host, not for libs for the target
    eah_host_compile=yes

    ## some sorry systems don't have proper GNU-make... build one
    if [ "$eah_force_aux" == "yes" ] || ! check_version make 3 79; then
	build_lsc_aux "make-3.80"
    fi
    ## FreeBSD's m4 seems to be broken? Download a fresh one
    if [ "$eah_force_aux" == "yes" ] || ! check_version m4 1 4; then
	build_lsc_aux "m4-1.4.1"
    fi
    if [ "$eah_force_aux" == "yes" ] || ! check_version autoconf 2 58; then
	build_lsc_aux "autoconf-2.59"
    fi
    if [ "$eah_force_aux" == "yes" ] || ! check_version automake 1 9; then
	build_lsc_aux "automake-1.9.6"
    fi
    if [ "$eah_force_aux" == "yes" ] || ! check_version libtool 1 4; then
	build_lsc_aux "libtool-1.5.6"
    fi
    if [ "$eah_force_aux" == "yes" ] || ! check_version pkg-config 0 15; then
	build_lsc_aux "pkgconfig-0.15.0"
	( cd "${AUX_INSTALL}/bin" &&
	    test -x pkg-config.dSYM &&
	    ! test -x pkg-config &&
	    ln -s pkg-config.dSYM pkg-config )
    fi

    eah_host_compile=no

    log_close

} # step1()


step1a()
{
    echo "----------------------------------------------------------------------"
    echo "STEP 1a: MinGW"
    echo "----------------------------------------------------------------------"

    aux_name=MinGW
    aux_built=${AUX_INSTALL}/${aux_name}.built
    LOGFILE=${BUILD_LOCATION}/step1a.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"
    
    log_popup

    ## this checks for compilation tools on the host, not for libs for the target
    eah_host_compile=yes

    if [ "${eah_rebuild_aux}" = "yes" ]; then
	if [ -r "${aux_built}" ]; then rm -f ${aux_built}; fi
    fi

    if [ ! -r "${aux_built}" ]; then

        # check for lex and yacc - bintools compilation will fail weirdly if absent
	passed=true
	if ! lex --version >/dev/null 2>/dev/null ; then
	    echo ERROR: lex missing - you probably want to install flex >> "$LOGFILE"
	    echo ERROR: lex missing - you probably want to install flex
	    passed=false
	fi
	if ! yacc --version >/dev/null 2>/dev/null ; then
	    echo ERROR: yacc missing - you probably want to install bison >> "$LOGFILE"
	    echo ERROR: yacc missing - you probably want to install bison
	    passed=false
	fi
	$passed || fail

	mkdir -p ${SOURCE_LOCATION}/MinGW
	cd ${SOURCE_LOCATION}/MinGW
	cvs -z3 -d ':pserver:anonymous@mingw.cvs.sourceforge.net:/cvsroot/mingw' checkout -P xscripts >> "$LOGFILE" 2>&1
	cd xscripts || fail
	chmod +x x86-mingw32-build.sh >> "$LOGFILE" || fail
	sed -i.dist '/assume  *DOWNLOAD_HOST/d;/assume  *WORKING_DIR/d;/assume  *PACKAGE_DIR/d;/assume  *INSTALL_DIR/d;/option  *GCC_LANGUAGE_OPTIONS/d' x86-mingw32-build.sh.conf
	echo "
          assume DOWNLOAD_HOST http://mesh.dl.sourceforge.net/mingw
          assume WORKING_DIR   ${SOURCE_LOCATION}/MinGW/build
          assume PACKAGE_DIR   ${SOURCE_LOCATION}/MinGW/packages
          assume INSTALL_DIR   ${AUX_INSTALL}
          option GCC_LANGUAGE_OPTIONS c++
          option GCC_LANGUAGE_OPTIONS objc
        " >> x86-mingw32-build.sh.conf

	echo "Building MinGW (this will take quite a while)..." | tee -a "$LOGFILE"
	./x86-mingw32-build.sh --unattended --no-post-clean $TARGET_HOST >> "$LOGFILE" 2>&1 || fail
	( cd "${AUX_INSTALL}/bin" && ln -s "${TARGET_HOST}-ar" ar )
	touch ${aux_built}
	echo "done."	 2>&1 | tee -a "$LOGFILE"

    else
	echo
	echo "Using previously built version of ${aux_name} in `basename ${AUX_SOURCES}`." 2>&1 | tee -a "$LOGFILE"
	echo "Use '--rebuild-aux' to force rebuild of auxiliary packages" 2>&1 | tee -a "$LOGFILE"
	echo
    fi

    eah_host_compile=no

    log_close
} ## step1a


step1b() 
{
    echo "----------------------------------------------------------------------"
    echo "STEP 1b: zlib, binutils"
    echo "----------------------------------------------------------------------"
    LOGFILE=${BUILD_LOCATION}/step1b.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"

    log_popup

    if [ "$eah_win32_cross_build" = "yes" ]; then
	aux_name=zlib-1.2.3
	aux_tgz=${aux_name}.tar.gz
	aux_built=${AUX_INSTALL}/${aux_name}.built
	cd ${AUX_SOURCES} >> "$LOGFILE" 2>&1 || fail 

        ## if --rebuild-aux was given: download & rebuild this package from scratch
	if [ "${eah_rebuild_aux}" = "yes" ]; then
	    if [ -r "${aux_built}" ]; then rm -f ${aux_built}; fi
	fi

	if [ ! -r "${aux_built}" ]; then
	    #wget_url "http://www.zlib.net" "${aux_name}.tar.gz" >> "$LOGFILE" 2>&1 || fail
	    log_and_do wget_url "http://www.gzip.org/zlib" "${aux_name}.tar.gz"
	    echo tar -xzvf zlib-1.2.3.tar.gz  >> "$LOGFILE"
	    tar -xzvf zlib-1.2.3.tar.gz  >> "$LOGFILE" 2>&1 || fail
	    echo cd zlib-1.2.3 >> "$LOGFILE"
	    cd zlib-1.2.3 >> "$LOGFILE" 2>&1 || fail
	    echo rm -f gcc ar >> "$LOGFILE"
	    rm -f gcc ar >> "$LOGFILE" 2>&1 || fail
	    echo ln -s "${AUX_INSTALL}/bin/${TARGET_HOST}-gcc" gcc >> "$LOGFILE"
	    ln -s "${AUX_INSTALL}/bin/${TARGET_HOST}-gcc" gcc >> "$LOGFILE" 2>&1 || fail
	    echo ln -s "${AUX_INSTALL}/bin/${TARGET_HOST}-ar"  ar  >> "$LOGFILE"
	    ln -s "${AUX_INSTALL}/bin/${TARGET_HOST}-ar"  ar  >> "$LOGFILE" 2>&1 || fail
	    echo PATH="$PWD:$PATH" make -f win32/Makefile.gcc libz.a >> "$LOGFILE"
	    PATH="$PWD:$PATH" make -f win32/Makefile.gcc libz.a >> "$LOGFILE" 2>&1 || fail
	    echo cp zlib.h zconf.h "${AUX_INSTALL}/include" >> "$LOGFILE"
	    cp zlib.h zconf.h "${AUX_INSTALL}/include" >> "$LOGFILE" 2>&1 || fail
	    echo cp libz.a "${AUX_INSTALL}/lib" >> "$LOGFILE"
	    cp libz.a "${AUX_INSTALL}/lib" >> "$LOGFILE" 2>&1 || fail
	    touch "${aux_built}"
	else
	    echo
	    echo "Using previously built version of ${aux_name} in `basename ${AUX_SOURCES}`." 2>&1 | tee -a "$LOGFILE"
	    echo "Use '--rebuild-aux' to force rebuild of auxiliary packages" 2>&1 | tee -a "$LOGFILE"
	    echo
	fi
    fi

    aux_name="binutils-2.18"
    aux_built=${AUX_INSTALL}/${aux_name}.built
    if [ "${eah_rebuild_aux}" = "yes" ]; then
	rm -rf "${aux_built}" "${AUX_SOURCES}/${aux_name}"
    fi
	
    if [ ! -r "${aux_built}" ]; then
	build_lsc_aux "$aux_name" "${CROSS_CONFIG_OPTS}"

        # some post-build installation due to targets missing in the library
	( cd "${AUX_SOURCES}/${aux_name}" &&
	    mkdir -p "${AUX_INSTALL}/include/bfd" &&
	    cp -r include/* bfd/*.h "${AUX_INSTALL}/include/bfd"
	) >> "$LOGFILE" 2>&1 || fail
        # patch a few headers
	if [ "$eah_win32_cross_build" = "yes" ]; then
	    mkdir -p "${AUX_INSTALL}/lib" &&
	      cp "${AUX_SOURCES}/${aux_name}/intl/libintl.a" "${AUX_INSTALL}/lib"
	    ( cd  "${AUX_INSTALL}/include/bfd" &&
		patch -N -p0 < ${eah_here}/binutils-MinGW-cross-build.patch
	    ) >> "$LOGFILE" 2>&1 || fail
	fi
    else
	echo
	echo "Using previously built version of ${aux_name} in `basename ${AUX_SOURCES}`." 2>&1 | tee -a "$LOGFILE"
	echo "Use '--rebuild-aux' to force rebuild of auxiliary packages" 2>&1 | tee -a "$LOGFILE"
	echo
    fi

    log_close

} ## step1b()


##----------------------------------------------------------------------
## STEP 2: FFTW
##----------------------------------------------------------------------
step2() 
{
    echo "----------------------------------------------------------------------"
    echo "STEP 2: FFTW"
    echo "----------------------------------------------------------------------"
    LOGFILE=${BUILD_LOCATION}/step2.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"

    log_popup

    ## build double-precision version
    build_lsc_aux "fftw-3.0.1" "--disable-shared ${CROSS_CONFIG_OPTS}"

    ## build single-precision version
    build_lsc_aux "fftw-3.0.1_f" "--disable-shared --enable-single ${CROSS_CONFIG_OPTS}"

    log_close

} ## step2()

##----------------------------------------------------------------------
## STEP 3
##----------------------------------------------------------------------
step3()
{
    echo "----------------------------------------------------------------------"
    echo "STEP 3: GSL"
    echo "----------------------------------------------------------------------"
    LOGFILE=${BUILD_LOCATION}/step3.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"

    log_popup

    config_opts="--disable-shared ${CROSS_CONFIG_OPTS}"

    if [ "${eah_os}" = darwin -a "${build_cpu}" = i686 ]; then
      config_opts="$config_opts --build=i386"
    fi
    if [ "$eah_win32_cygwin_build" = "yes" ]; then
      config_opts="$config_opts $eah_CPPFLAGS $eah_CFLAGS $eah_CXXFLAGS $eah_LDFLAGS"
    fi
    build_lsc_aux "gsl-1.9" "$config_opts"

    log_close

} ## step3()

##----------------------------------------------------------------------
## STEP 4
##----------------------------------------------------------------------
step4() 
{
    echo "----------------------------------------------------------------------"
    echo "STEP4: BOINC"
    echo "----------------------------------------------------------------------"
    LOGFILE="${BUILD_LOCATION}/step4.log"
    rm -f "$LOGFILE"
    touch "$LOGFILE"

    log_popup

    local boinc_src

    ## ---------- get cvs/svn-archive

    if [ -n "$BOINC_SVN_TAG" ]; then
        if [ "$BOINC_SVN_TAG" = "server_stable" ]; then
	    boinc_src=boinc-svn-server_stable
	    boinc_url="http://boinc.berkeley.edu/svn/branches/server_stable"
	    boinc_tag="HEAD"
	elif [ "$BOINC_SVN_TAG" = "boinc_core_release_6_2" ]; then
	    boinc_src=boinc-svn-62
	    boinc_url="http://boinc.berkeley.edu/svn/branches/boinc_core_release_6_2"
	    boinc_tag="HEAD"
	else
	    boinc_src=boinc-svn
	    boinc_url="http://boinc.berkeley.edu/svn/trunk/boinc"
	    boinc_tag="$BOINC_SVN_TAG"
	fi


	if ! [ -d "${BUILD_LOCATION}/extra_sources/$boinc_src" ]; then
	    if ! [ -d "${SOURCE_LOCATION}/$boinc_src" ]; then
		echo svn co -r "$boinc_tag" "$boinc_url" "${SOURCE_LOCATION}/$boinc_src" >> "$LOGFILE" 2>&1
		echo -n Checking out BOINC svn ...
		echo svn co -r "$boinc_tag" "$boinc_url" "${SOURCE_LOCATION}/$boinc_src" >> "$LOGFILE"
		svn co -r "$boinc_tag" "$boinc_url" "${SOURCE_LOCATION}/$boinc_src" >> "$LOGFILE" 2>&1 || fail
		if ! [ -d "${SOURCE_LOCATION}/$boinc_src/zip" ]; then
		    echo svn co 'http://boinc.berkeley.edu/svn/trunk/boinc/zip@18044' "${SOURCE_LOCATION}/$boinc_src/zip" >> "$LOGFILE"
		    svn co 'http://boinc.berkeley.edu/svn/trunk/boinc/zip@18044' "${SOURCE_LOCATION}/$boinc_src/zip" >> "$LOGFILE" 2>&1 || fail
		fi
		echo \ done.
	    fi
	    echo -n Copying BOINC svn ...
	    echo cp -R -p -f "${SOURCE_LOCATION}/$boinc_src" "${BUILD_LOCATION}/extra_sources/" >> "$LOGFILE"
	    cp -R -p -f "${SOURCE_LOCATION}/$boinc_src" "${BUILD_LOCATION}/extra_sources/" >> "$LOGFILE" 2>&1 || fail
	    echo \ done.
	fi
	if [ "${eah_no_update}" != yes ]; then
	    echo INFO: using subversion with tag "$BOINC_SVN_TAG" >> "$LOGFILE" 2>&1
	    
	    echo -n Updating BOINC svn ...
	    cd "${BUILD_LOCATION}/extra_sources/${boinc_src}" >> "$LOGFILE" 2>&1 || fail
	    echo svn up -r "$boinc_tag" >> "$LOGFILE"
	    svn up -r "$boinc_tag"  >> "$LOGFILE" 2>&1
	    echo \ done.
	fi
	if ! [ -r "${BUILD_LOCATION}/extra_sources/$boinc_src/zip/Makefile.am" ]; then
	    echo cp -R -p -f "${SOURCE_LOCATION}/$boinc_src/zip" "${BUILD_LOCATION}/extra_sources/$boinc_src" >> "$LOGFILE"
	    cp -R -p -f "${SOURCE_LOCATION}/$boinc_src/zip" "${BUILD_LOCATION}/extra_sources/$boinc_src" >> "$LOGFILE" 2>&1 || fail
	fi
    elif [ "${eah_boinc_rpdist}" = yes ]; then
	get_pkg_tgz "${EAH_PKGS_URL}" boinc-rpdist
	boinc_src=boinc-rpdist
	echo done.
    else
	get_pkg_cvs "${BOINC_CVS}" "boinc" "" "${boinc_update_opts}"
	boinc_src=boinc-CVS
    fi

    cd "${BUILD_LOCATION}/extra_sources/${boinc_src}" >> "$LOGFILE" 2>&1 || fail

    ## ---------- run patch script if told to
    if [ -n "$BOINC_PATCH" ]; then
        echo $ECHO_N "Patching BOINC... $ECHO_C"
        $BOINC_PATCH >> "$LOGFILE" 2>&1 || fail
    fi

  if [ "$eah_win32_cross_build" = "yes" ] || [ "$eah_win32_cygwin_build" = "yes" ]; then
    ## -------------------------
    ## WIN32 MinGW (Cross-)Build
    ## -------------------------

    patch -N -p0 < "${eah_here}/boinc_zip-MinGW.patch"   >> "$LOGFILE" 2>&1 #|| fail
    makefile="$PWD/lib/Makefile.mingw"
    if ! [ -r "$makefile" ]; then
	makefile="${eah_here}/Makefile.mingw"      
	patch -N -p0 < "${eah_here}/boinc-MinGW-build.patch" >> "$LOGFILE" 2>&1 #|| fail
    fi

    if [ "$eah_win32_cross_build" = "yes" ]; then 
      boinc_make_opts="-f $makefile BOINC_SRC=${PWD} BOINC_PREFIX=${BUILD_INSTALL} CC=${TARGET_HOST}-gcc CXX=${TARGET_HOST}-g++ RANLIB=${TARGET_HOST}-ranlib ${eah_CPPFLAGS}"
    elif [ "$eah_win32_cygwin_build" = "yes" ]; then
      boinc_make_opts="-f $makefile BOINC_SRC=${PWD} BOINC_PREFIX=${BUILD_INSTALL} NOCYGWIN=-mno-cygwin ${eah_CPPFLAGS}"
      # NOTE: adding -mno-cygwin undefines __CYGWIN__ macro on Cygwin's MinGW/cpp
    fi

    ## IMPORTANT: make sure we uninstall any previous libs/headers
    echo $ECHO_N "Uninstalling previous BOINC... $ECHO_C"
    eah_next="make ${boinc_make_opts} uninstall"
    echo ${eah_next} >> "$LOGFILE" 2>&1
    eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    echo "done."

    ## clean if asked for
    if [ "${eah_recompile}" = yes ]; then
	echo $ECHO_N "Cleaning source-tree... $ECHO_C"
	eah_next="make ${boinc_make_opts} clean"
	echo ${eah_next} >> "$LOGFILE" 2>&1
	eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi

    ## use the (patched) Makefile.mingw
    echo $ECHO_N "Building BOINC... $ECHO_C"
    eah_next="make ${boinc_make_opts} all libboinc_zip.a"
    echo ${eah_next} >> "$LOGFILE" 2>&1
    eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    echo "done."

    ## install
    echo $ECHO_N "Installing BOINC... $ECHO_C"
    eah_next="make ${boinc_make_opts} install install-zip"
    echo ${eah_next} >> "$LOGFILE" 2>&1
    eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    if [ -r version.h ]; then
      test -d "${BUILD_INSTALL}/include/BOINC" &&
	cp version.h "${BUILD_INSTALL}/include/BOINC"
      test -d "${BUILD_INSTALL}/include/boinc" &&
	cp version.h "${BUILD_INSTALL}/include/boinc"
    fi
    echo "done."

  else # $eah_win32_cross_build

    ## -------------------------------
    ## (ordinary) autoconf BOINC build
    ## -------------------------------

    # compile an own jpeglib
    aux_name="jpeg-6b"
    aux_built="${AUX_INSTALL}/${aux_name}.built"
    if [ "${eah_rebuild_aux}" = "yes" ]; then
	rm -rf "${aux_built}" "${AUX_SOURCES}/${aux_name}"
    fi

    if [ ! -r "${aux_built}" ]; then
	log_and_do cd "${AUX_SOURCES}"
	wget_url "$EAH_PKGS_URL" "jpegsrc.v6b.tar.gz"
	log_and_do tar -xzf jpegsrc.v6b.tar.gz
	echo -n Building libjpeg ...
	log_and_do cd jpeg-6b
	log_and_do ./configure --prefix="${AUX_INSTALL}" CPPFLAGS="$CPPGLAGS" CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS" LIBS="$LIBS"
	log_and_do make
	log_and_do make install-lib
	touch "${aux_built}"
	echo \ done.
    fi

    log_and_do cd "${BUILD_LOCATION}/extra_sources/${boinc_src}"

    # fix the Makefile.am to build boinc_zip if too new
    if ! grep 'API_SUBDIRS *=.* zip' Makefile.am >/dev/null; then
	sed -i+ 's/\(API_SUBDIRS *=.*\)/\1 zip/' Makefile.am 
	sed -i+ 's%AC_CONFIG_FILES(\[%AC_CONFIG_FILES(\[\
                 zip/Makefile\
                 zip/zip/Makefile\
                 zip/unzip/Makefile%' configure.ac
    fi

    ## ---------- configure boinc
    echo $ECHO_N "Configuring BOINC... $ECHO_C"
    eah_next="./_autosetup" 
    echo ${eah_next} >> "$LOGFILE" 2>&1 
    ${eah_next} >> "$LOGFILE" 2>&1     || fail

    test -n "$BOINCC" && boincc="CC=${BOINCC}"
    boinc_configure_args=`echo "${eah_configure_args}" | sed 's/--disable-debug//g'`

    sdkopt=""
    if [ "$eah_macos_sdk_for_boinc_only" = "yes" ]; then
	OLD_SDKROOT="$SDKROOT"
	OLD_MACOSX_DEPLOYMENT_TARGET="$MACOSX_DEPLOYMENT_TARGET"
	export SDKROOT="/Developer/SDKs/MacOSX10.4u.sdk"
	export MACOSX_DEPLOYMENT_TARGET="10.4"
	sdkopt="$sdkopt LDFLAGS='-isysroot /Developer/SDKs/MacOSX10.4u.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -arch i386 $LDFLAGS'"
	sdkopt="$sdkopt CPPFLAGS='-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CPPFLAGS'"
	sdkopt="$sdkopt CFLAGS='-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CFLAGS'"
	sdkopt="$sdkopt CXXFLAGS='-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CXXFLAGS'"
    fi
    eah_next="./configure --prefix=${BUILD_INSTALL} $boincc $sdkopt --disable-server --disable-client --disable-manager --disable-shared --enable-static ${with_ssl} ${boinc_configure_args}"
    echo ${eah_next} >> "$LOGFILE" 2>&1
    eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    echo "done."

    ## IMPORTANT: make sure we uninstall any previous libs/headers to avoid
    ## mixing old and new compile+installs
    echo $ECHO_N "Uninstalling previous BOINC-build ... $ECHO_C"
    echo "make uninstall" >> "$LOGFILE" 2>&1
    make uninstall >> "$LOGFILE" 2>&1
    echo "done."

    ## handle recompile option
    if [ "${eah_recompile}" = yes ]; then
	echo $ECHO_N "Cleaning source-tree... $ECHO_C"
	make clean >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi

    ## make a brave attempt at using boinc's build-infrastructure (good luck!)
    echo $ECHO_N "Building BOINC-libs... $ECHO_C"
    make >> "$LOGFILE" 2>&1 || fail
    make install >> "$LOGFILE" 2>&1 || fail
    if [ -r version.h ]; then
      test -d "${BUILD_INSTALL}/include/BOINC" &&
	cp version.h "${BUILD_INSTALL}/include/BOINC"
      test -d "${BUILD_INSTALL}/include/boinc" &&
	cp version.h "${BUILD_INSTALL}/include/boinc"
    fi
    if [ "$eah_macos_sdk_for_boinc_only" = "yes" ]; then
	SDKROOT="$OLD_SDKROOT"
	MACOSX_DEPLOYMENT_TARGET="$OLD_MACOSX_DEPLOYMENT_TARGET"
    fi
    echo "done."

  fi # $eah_win32_cross_build

  log_close

} ## step4()


##----------------------------------------------------------------------
## STEP 5
##----------------------------------------------------------------------
step5()
{
    echo "----------------------------------------------------------------------"
    echo "STEP5: LAL (LSC Algorithm library)"
    echo "----------------------------------------------------------------------"
    log_start "${BUILD_LOCATION}/step5.log"

    ## options to be used to configure LAL
    local lal_config_opts="--enable-frame=no --enable-metaio=no --enable-mpi=no --disable-shared ${eah_configure_args} --enable-boinc ${CROSS_CONFIG_OPTS}"

    log_and_do cd "${SOURCE_LOCATION}"

    if [ -d lalsuite/.git ]; then
        echo -n Using lalsuite repo ...
    else
        echo -n Creating lalsuite link ...
        rm -rf lalsuite
        ln -s ../../../../../../.. lalsuite
    fi
    echo \ done.

    ## configure
    log_and_do cd lalsuite/lal
    echo $ECHO_N "Configuring LAL... $ECHO_C"
    if [ "${eah_recompile}" = yes -o ! -r configure ] ; then
	log_and_do ./00boot
	log_and_do cd ../lalpulsar
	log_and_do ./00boot
    fi

    log_and_do mkdir -p "${BUILD_LOCATION}/lalsuite/lal"
    log_and_do cd "${BUILD_LOCATION}/lalsuite/lal"

    if [ "${eah_recompile}" = yes -o ! -r config.status ]; then
	eah_next="'${SOURCE_LOCATION}/lalsuite/lal/configure' --prefix=${BUILD_INSTALL} ${lal_config_opts}"
	echo ${eah_next} >> "$LOGFILE" 2>&1
	eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    fi
    echo "done."
	
    ## make
    if [ "${eah_recompile}" = yes ]; then
	echo $ECHO_N "Cleaning source-tree... $ECHO_C"
	make clean >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi
    echo $ECHO_N "Building LAL ... $ECHO_C"
    make >> "$LOGFILE" 2>&1 || fail
    make install >> "$LOGFILE" 2>&1 || fail
    echo "done."

    log_and_do mkdir -p ../lalpulsar
    log_and_do cd ../lalpulsar

    ## configure
    echo $ECHO_N "Configuring LALPulsar... $ECHO_C"
    if [ "${eah_recompile}" = yes -o ! -r config.status ]; then
	eah_next="'${SOURCE_LOCATION}/lalsuite/lalpulsar/configure' --prefix=${BUILD_INSTALL} ${lal_config_opts}"
	echo ${eah_next} >> "$LOGFILE" 2>&1
	eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    fi
    echo "done."
	
    ## make
    if [ "${eah_recompile}" = yes ]; then
	echo $ECHO_N "Cleaning source-tree... $ECHO_C"
	make clean >> "$LOGFILE" 2>&1 || fail
	echo "done."
    fi
    echo $ECHO_N "Building LALPulsar ... $ECHO_C"
    make >> "$LOGFILE" 2>&1 || fail
    make install >> "$LOGFILE" 2>&1 || fail
    echo "done."

    log_close

} ## step5()


##----------------------------------------------------------------------
## STEP 6
##----------------------------------------------------------------------
step6()
{
    echo "----------------------------------------------------------------------"
    echo "STEP6: LALApps and Einstein@Home"
    echo "----------------------------------------------------------------------"
    LOGFILE=${BUILD_LOCATION}/step6.log
    rm -f  "$LOGFILE"
    touch "$LOGFILE"

    log_popup

    ## use the lalsuite checked out in step 5 (LAL)
    cd "${SOURCE_LOCATION}" >> "$LOGFILE" 2>&1 || fail
    if [ ! -d lalsuite/lalapps ]; then
	echo "could not find $PWD/lalsuite/lalapps - run step 5 once first!" >> "$LOGFILE" 2>&1
	fail
    fi
    log_and_do cd lalsuite/lalapps

    echo $ECHO_N "Configuring LALApps... $ECHO_C"
    log_and_do ./00boot

    log_and_do mkdir -p "${BUILD_LOCATION}/lalsuite/lalapps"
    log_and_do cd "${BUILD_LOCATION}/lalsuite/lalapps"

    eah_next="PATH='.:$PATH' '${SOURCE_LOCATION}/lalsuite/lalapps/configure'"
    eah_next="$eah_next --prefix=${BUILD_INSTALL} ${eah_configure_args} ${CROSS_CONFIG_OPTS}"
    eah_next="$eah_next --disable-gcc-flags --disable-frame --disable-metaio --enable-boinc --disable-silent-rules"
    echo ${eah_next} >> "$LOGFILE" 2>&1 || fail
    eval ${eah_next} >> "$LOGFILE" 2>&1 || fail
    echo "done."

    echo $ECHO_N "Building Einstein@Home... $ECHO_C"
    log_and_do cd src/lalapps
    log_and_do make LALAppsVCSInfo.h liblalapps.la
    log_and_do cd ../pulsar/hough/src2
    log_and_do make eah_HierarchicalSearch${eah_target_ext}
    echo "done."

    log_close
} ## step6()


## ======================================================================
##
## MAIN
##
## ======================================================================

## ---------- before we have a build-specific logfile, use a local one
LOGFILE=eah_build.log
touch "$LOGFILE"
echo "----------------------------------------------------------------------" >> "$LOGFILE"
echo `date` >> "$LOGFILE"
## ----------

SHELL=/bin/sh
TODAY=`date +%Y_%m_%d`
eah_here=`pwd`

# arguments for starting up xterm
XTERMGEOM="-geometry 120x30+60+10"
eah_startdir=`pwd`

## what kind of 'wget' is installed:
if [ -x "`type -p wget`" ]; then
    eah_WGET=wget
elif [ -x "`type -p curl`" ]; then
    eah_CURL=curl
elif [ -x "`type -p lynx`" ]; then
    eah_CURL="lynx -source"
else
    echo "Sorry, no download-program found (wget, curl or lynx), please install one of these"
    fail
fi

## url's to get sources from
LSCSOFT_CVS=":pserver:anonymous@gravity.phys.uwm.edu:2402/usr/local/cvs/lscsoft"
BOINC_CVS=":pserver:anonymous@alien.ssl.berkeley.edu:/home/cvs/cvsroot"
EAH_PKGS_URL="http://www.aei.mpg.de/~repr/EaH_packages/"  

## set CVS tags - keep from environment, from config file or builtin default
if [ -r "eah_config.sh" ]; then
  if [ -z "${LALAPPS_TAG}" ]; then
    LALAPPS_TAG="`sed -n 's/^LALAPPS_TAG=//p' eah_config.sh`"
  fi
  if [ -z "${LAL_TAG}" ]; then
    LAL_TAG="`sed -n 's/^LAL_TAG=//p' eah_config.sh`"
  fi
  if [ -z "${BOINC_TAG}" ]; then
    BOINC_TAG="`sed -n 's/^BOINC_TAG=//p' eah_config.sh`"
  fi
  if [ -z "${BOINC_SVN_TAG}" ]; then
    BOINC_SVN_TAG="`sed -n 's/"//g;s/^BOINC_SVN_TAG=//p' eah_config.sh`"
  fi
  if [ -z "${EAH_TAG}" ]; then
    EAH_TAG="`sed -n 's/^EAH_TAG=//p' eah_config.sh`"
  fi
  if [ -z "${APP_NAME}" ]; then
    APP_NAME="`sed -n 's/^APP_NAME=//p' eah_config.sh`"
  fi
fi
if [ -z "${LALAPPS_TAG}" ]; then
    LALAPPS_TAG=-A
fi
if [ -z "${LAL_TAG}" ]; then
    LAL_TAG=-A
fi
if [ -z "${EAH_TAG}" ]; then
    EAH_TAG=-A
fi
if [ -z "${BOINC_TAG}" ]; then
    BOINC_TAG="-D'Mon Feb  7 00:42:56 CET 2005'"
fi
if [ -z "${APP_NAME}" ]; then
    APP_NAME=einstein
fi



## ========================================
## get some info about the build-machine
## ========================================

## we're using automake's 'config.guess' to determine the build-machine,
## so get that script first:
if [ ! -r ./config.guess ]; then
    wget_url ${EAH_PKGS_URL} config.guess
fi

eah_config_guess=`$SHELL ./config.guess`
if [ -z "${eah_config_guess}" ]; then
    echo "ERROR: failed to determine canonical system-name" 
    fail
fi
## extract components (taken from configure)
build_cpu=`echo $eah_config_guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\1/'`
build_vendor=`echo $eah_config_guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\2/'`
build_os_full=`echo $eah_config_guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\3/'`

## short form for internal use
eah_os=`echo ${build_os_full} | sed 's/^\([a-zA-Z]*\).*$/\1/'`

echo eah_os: $eah_os

## BOINC's canonical names doesn't contain OS - version-numbers
build_os_short=`echo ${build_os_full} | sed 's/^\([a-zA-Z-]*\)[0-9.]*$/\1/'`
eah_system_name=${build_cpu}-${build_vendor}-${build_os_short}

## ==============================
## portable way of recursive copying
COPYRECURSIVE="cp -Rp"

## print the usage-string
usage() 
{
    echo
    echo "Usage: $0"
    echo
    echo "  -h, --help                print this help-string"
    echo "      --nox                 don't try to open log-xterms"
    echo "      --rebuild-aux         force rebuild of auxiliary packages (pkgconfig,gsl,fftw)"
    echo "      --force-aux           force building own automake, autoconf etc. in case the local ones don't do what they should (e.g. MacPorts)"
    echo "      --no-update           don't cvs-update any sources [default = update-all]"
    echo "      --recompile           force recompile of main code-components (boinc,lal,Einstein@Home)"
    echo "      --optimize            build a fully optimized binary"
    echo "      --debug               add debugging switches to build"
    echo "      --profile             add profiling switches to build [implies --recompile]"
    echo "      --release <version>   build a portable release-binary <version> [implies --optimize,--recompile]"
    echo "      --new-boinc           build with latest BOINC SVN-version [default: \$BOINC_TAG ]"
    echo "      --no-graphics         force building of non-graphical Einstein@Home application"
    echo "      --no-dlopen           do not try to use dlopen for dynamical screensaver-lib under Linux"
    echo "      --with-ssl=<path>     pass path to ssl installation to BOINC configure"
    echo "      --sdk                 use a specific system SDK (currently only 10.4-i386 and 10.3-ppc work for MacOS)"
    echo "      --xcode3x86           turn on hack for Xcode3 under MacOS 10.5 Intel (equal to --sdk 10.4-i386)"
    echo "      --win32-cross         cross-compile a win32 App using MinGW"
    echo "      --win32-cygwin        compile a win32 App on Cygwin"
    echo "      --n0-git-repo         use the Einstein@home git repository @AEI-H instead of the canonical UWM lalsuite"
    echo "      --cuda                build a CUDA app (highly experimental)"
    echo "      --opencl              build a OpenCL app (highly experimental)"
    echo "      --64                  build a 64 bit app"
    echo "      --32                  force building a 32 bit app (on a 64 bit machine)"
    echo
    exit 0

} ## usage()

## user-options
if [ `uname -s` = Darwin ]; then
  eah_show_log=mac
elif uname -s | grep CYGWIN >/dev/null; then
  eah_show_log=rxvt
else
  eah_show_log=X
fi
eah_rebuild_aux=no
eah_no_update=no
eah_optimize=no
eah_debug=no
## don't put anything as default into release-string!
eah_release=
eah_recompile=no
eah_cvs_lal=yes
eah_new_lal=no
eah_new_boinc=no
eah_boinc_rpdist=no
eah_graphics=yes
eah_win32_cross_build=no
eah_cuda=no
eah_opencl=no
eah_64bit=no
eah_32bit=no
n0_git_repo=no
eah_host_compile=no
eah_force_aux=no
## only try to use dynamic screensaver-lib under Linux (can be turned off by user even then)
if [ "${eah_os}" = linux ]; then
    eah_dlopen=yes
else
    eah_dlopen=no
fi

cmdline_log="$0 $@"  ## keep for log-file
## -----------------------------------
## read through command-line options
while [ -n "$1" ]; do
    case "$1" in
	--nox) 
	    eah_show_log=no
	    ;;
	--rebuild-aux)
	    eah_rebuild_aux=yes
	    ;;
	--force-aux)
	    eah_force_aux=yes
	    ;;
	--xcode3x86)
	    eah_sdk=10.4-i386
	    ;;
	--n0-git-repo)
	    n0_git_repo=yes
	    ;;
	--cuda)
	    eah_cuda=yes
	    ;;
	--opencl)
	    eah_opencl=yes
	    ;;
	--64)
	    eah_64bit=yes
	    ;;
	--32)
	    eah_32bit=yes
	    ;;
	--sdk)
	    shift
	    eah_sdk="$1"
	    ;;
	-h | --help) 
	    usage
	    ;;
	--no-update)
	    eah_no_update=yes
	    ;;
	--optimize)
	    eah_optimize=yes
	    ;;
	--win32-cross)
	    eah_win32_cross_build=yes
	    ;;
	--win32-cygwin)
	    eah_win32_cygwin_build=yes
	    ;;
	--release)
	    shift
	    case $1 in
		--* | "")
		    echo "ERROR: --release requires a version-argument"
		    fail
		    ;;
		*)
		    eah_release="$1"
		    ;;
	    esac
	    eah_recompile=yes
	    ;;
	--debug)
	    eah_debug=yes
	    CFLAGS="-g3 -ggdb -fno-common $CFLAGS"
	    CXXFLAGS="-g3 -ggdb -fno-common -fno-inline $CXXFLAGS"
	    CPPFLAGS="-D_DEBUG -DUSE_BOINC_DEBUG -Wall -W -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings $CPPFLAGS"
	    ;;
	--profile)
	    CFLAGS="$CFLAGS -pg"
	    CXXFLAGS="$CXXFLAGS -pg"
	    eah_recompile=yes
	    ;;
	--recompile)
	    eah_recompile=yes
	    ;;
	--new-lal)
	    eah_new_lal=yes;
	    ;;
	--no-cvs-lal)
	    eah_cvs_lal=no;
	    ;;
	--new-boinc)
	    eah_new_boinc=yes;
	    eah_boinc_rpdist=no;
	    ;;
	--no-graphics)
	    eah_graphics=no;
	    ;;
	--no-dlopen)
	    eah_dlopen=no;
	    ;;
	--with-ssl=*)
	    with_ssl="$1";
	    ssldir=`echo "$1" | sed 's/--with-ssl=//'`
	    CPPFLAGS="$CPPFLAGS -I$ssldir/include"
	    LIBS="$LIBS -L$ssldir/lib"
	    ;;
	*) 
	    echo "Unknown option or paramter: '$1'";
	    usage
	    ;;
    esac
    shift
done

# some sanity checks on options
if [ "$eah_cuda" = "yes" ] && [ "$eah_opencl" = "yes" ] ; then
    echo "ERROR: enabling both CUDA and OpenCL makes no sense - exiting"
    exit -1
fi
if [ "$eah_32bit" = "yes" ] && [ "$eah_64bit" = "yes" ] ; then
    echo "ERROR: enabling both 32bit and 64bit makes no sense - exiting"
    exit -1
fi


# check the CUDA version
if [ "$eah_cuda" = "yes" ]; then
    if [ "$eah_win32_cygwin_build" = "yes" ]; then
	export MSCDIR="`env | sed -n 's/VS[0-9]*COMNTOOLS=//p'`"
	export PATH="$PATH:`cygpath -u '$MSCDIR'`:/cygdrive/c/CUDA/bin"
	export CUDA_INSTALL_PATH=/cygdrive/c/CUDA
    fi
    cuda_vers=`nvcc --version | sed -n 's/.* release \([0-9]\.[0-9]*\).*/\1/p'`
    if echo $cuda_vers | awk -F. '{ major=2; minor=2;
                                    if ($1 < major) {exit -1};
                                    if ($1 > major) {exit 0};
                                    if ($2 < minor) {exit -1};
                                  }' 
    then
	:
    else
	echo "ERROR: There's something wrong with your CUDA setup."
	echo "detected CUDA version is $cuda_vers"
	echo "Einstein@home currently only compiles with CUDA 2.2 or newer"
	exit -1
    fi
fi

## some post-processing for release-building: 
## --release without --debug implies --optimize
if [ -n "${eah_release}" ] && [ "${eah_debug}" = no ]; then
    eah_optimize=yes
fi

## turn on optimization-switches
if [ "${eah_optimize}" = yes ]; then
    CPPFLAGS="$CPPFLAGS"
    CFLAGS="-O3 $CFLAGS"
    CXXFLAGS="-O3 $CXXFLAGS"
    if [ "${eah_os}" = darwin -a "${build_cpu}" = powerpc ]; then  ## some darwin-powerpc-specific optimizations
	darwin_opt=""
	CFLAGS="${darwin_opt} ${CFLAGS}"
	CXXFLAGS="${darwin_opt} ${CXXFLAGS}"
    fi
fi

## ---------------------------------
## check availability of xterm for log-output:
if [ "${eah_show_log}" = "X" ]; then
    echo $ECHO_N "Checking xterm... $ECHO_C"
    if (! type xterm >& /dev/null) || (! xterm -geometry 1x1+0+0 -e /bin/sh exit >& /dev/null); then 
	eah_show_log=no;
	echo
	echo "failed. No log-output will be shown!"
	echo
    else
	echo "ok."
    fi
elif [ "${eah_show_log}" = "X" ]; then
    echo $ECHO_N "Checking rxvt... $ECHO_C"
    if (! type rxvt >& /dev/null) || (! rxvt -geometry 1x1+0+0 -e /bin/sh exit >& /dev/null); then 
	eah_show_log=no;
	echo
	echo "failed. No log-output will be shown!"
	echo
    else
	echo "ok."
    fi
elif [ "${eah_show_log}" = "mac" ]; then
    rm -rf /tmp/follow.command
    ln -s /usr/bin/true /tmp/follow.command
    if [ `uname -r | sed 's/\..*//'` -gt 8 ]; then
	eah_mac_term="-a Terminal ${eah_here}/follow_log.terminal"
    else
	eah_mac_term="-a Terminal ${eah_here}/follow_log.term"
    fi
    if open ${eah_mac_term} 2>/dev/null; then
	sleep 1
	rm -rf /tmp/follow.command
    else
	eah_show_log=no;
	echo
	echo "failed. No log-output will be shown!"
	echo
    fi
fi


## -----------------------------------
## set source- and build-locations
if [ -z "$SOURCE_LOCATION" ]; then
    if [ "$eah_win32_cross_build" = "yes" ]; then
	if [ "$eah_64bit" = "yes" ]; then
	    SOURCE_LOCATION=${eah_here}/EaH_sources_win64
	else
	    SOURCE_LOCATION=${eah_here}/EaH_sources_win32
	fi
    elif [ "$eah_64bit" = "yes" ]; then
	SOURCE_LOCATION=${eah_here}/EaH_sources_64
    elif [ "$eah_32bit" = "yes" ]; then
	SOURCE_LOCATION=${eah_here}/EaH_sources_32
    else
	SOURCE_LOCATION=${eah_here}/EaH_sources
    fi
fi

if [ -z "$BUILD_LOCATION" ]; then
    if [ "$eah_win32_cross_build" = "yes" ]; then
	if [ "$eah_64bit" = "yes" ]; then
	    BUILD_LOCATION=${eah_here}/EaH_build_win64
	else
	    BUILD_LOCATION=${eah_here}/EaH_build_win32
	fi
    elif [ "$eah_64bit" = "yes" ]; then
	BUILD_LOCATION=${eah_here}/EaH_build_64
    elif [ "$eah_32bit" = "yes" ]; then
	BUILD_LOCATION=${eah_here}/EaH_build_32
    else
	BUILD_LOCATION=${eah_here}/EaH_build
    fi
fi

if [ "$eah_64bit" = "yes" ]; then
    CPPFLAGS="$CPPFLAGS -m64" # stupid, but necessary for BOINC
    CFLAGS="$CFLAGS -m64"
    CXXFLAGS="$CXXFLAGS -m64"
    LDFLAGS="$LDFLAGS -m64"
elif [ "$eah_32bit" = "yes" ]; then
    CPPFLAGS="$CPPFLAGS -m32" # stupid, but necessary for BOINC
    CFLAGS="$CFLAGS -m32"
    CXXFLAGS="$CXXFLAGS -m32"
    LDFLAGS="$LDFLAGS -m32"
fi


## ------------------------------
## any release gets dedicated build-directory name
if [ -n "${eah_release}" ]; then
    BUILD_LOCATION=${BUILD_LOCATION}_release_${APP_NAME}_${eah_release}
    if [ -d ${BUILD_LOCATION} ]; then
	echo "================================================================================"
	echo "WARNING: Building release ${eah_release} but the corresponding BUILD_LOCATION already exists:"
	echo 
	echo "${BUILD_LOCATION}"
	echo
	echo "==> Each release-binary *SHOULD REALLY* be built in a dedicated build-directory!"
	echo 
	echo "Only proceed if you know what you're doing, otherwise 'Ctrl-C' NOW and choose a"
	echo "new unique version-number for your release. You have been warned!!"
	echo "================================================================================"
	read eah_confirm
    fi
fi

## check that BUILD_LOCATION and SOURCE_LOCATION are absolute paths
tmppath1=${BUILD_LOCATION#/}
tmppath2=${SOURCE_LOCATION#/}
if [ "${tmppath1}" = "${BUILD_LOCATION}" ] || [ "${tmppath2}" = "${SOURCE_LOCATION}" ];  then
    echo 
    echo "ERROR: BUILD_LOCATION and SOURCE_LOCATION have to be absolute paths!"
    echo 
    fail
fi

## handle Mac OS X SDKs
if [ "$eah_cuda" = "yes" ] &&  uname -v | grep '^Darwin.*I386$' >/dev/null; then
    eah_macos_sdk_for_boinc_only=yes
elif [ "${eah_sdk}" = "10.4-i386" ]; then
    export LDFLAGS="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -arch i386 $LDFLAGS"
    export CPPFLAGS="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CPPFLAGS"
    export CFLAGS="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CFLAGS"
    export CXXFLAGS="-isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386 $CXXFLAGS"
    export SDKROOT="/Developer/SDKs/MacOSX10.4u.sdk"
    export MACOSX_DEPLOYMENT_TARGET=10.4
elif [ "${eah_sdk}" = "10.3-ppc" ]; then
    export LDFLAGS="-arch ppc -D_NONSTD_SOURCE -isystem /Developer/SDKs/MacOSX10.3.9.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.3.9.sdk $LDFLAGS"
    export CPPFLAGS="-arch ppc -D_NONSTD_SOURCE -isystem /Developer/SDKs/MacOSX10.3.9.sdk $CPPFLAGS"
    export CFLAGS="-arch ppc -D_NONSTD_SOURCE -isystem /Developer/SDKs/MacOSX10.3.9.sdk $CFLAGS"
    export CXXFLAGS="-arch ppc -D_NONSTD_SOURCE -isystem /Developer/SDKs/MacOSX10.3.9.sdk $CXXFLAGS"
    export SDKROOT="/Developer/SDKs/MacOSX10.3.9.sdk"
    export MACOSX_DEPLOYMENT_TARGET=10.3
elif false; then
    export LDFLAGS="-isysroot /Developer/SDKs/MacOSX10.3.9.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.3.9.sdk -arch ppc $LDFLAGS /Developer/SDKs/MacOSX10.3.9.sdk/usr/lib/libSystemStubs.a"
#   export LDFLAGS="-isysroot /Developer/SDKs/MacOSX10.3.9.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.3.9.sdk -arch ppc $LDFLAGS -L/Developer/SDKs/MacOSX10.3.9.sdk/usr/lib -lSystemStubs"
    export CPPFLAGS="-isysroot /Developer/SDKs/MacOSX10.3.9.sdk -arch ppc $CPPFLAGS"
    export CFLAGS="-isysroot /Developer/SDKs/MacOSX10.3.9.sdk -arch ppc $CFLAGS"
    export CXXFLAGS="-isysroot /Developer/SDKs/MacOSX10.3.9.sdk -arch ppc $CXXFLAGS"
    export SDKROOT="/Developer/SDKs/MacOSX10.3.9.sdk"
    export MACOSX_DEPLOYMENT_TARGET=10.3
fi


## ---------- compatibility-switch for compilation/build ----------
## We alwyas compile with __NO_CTYPE as default, because this
## avoids GLIBC_2.3 symbols being pulled in via ctype.h-macros,
## This is important for building portable binaries, and in any 
## case should have no negative effects, so we use it by default.

CPPFLAGS="$CPPFLAGS -D__NO_CTYPE"

## We now always compile for BOINC APIv6
CPPFLAGS="$CPPFLAGS -DBOINC_APIV6"

## Linux binaries compiled with this script should all feature extended statcktraces
CPPFLAGS="$CPPFLAGS -DEXT_STACKTRACE"

## --------------------------------------------------

BUILD_INSTALL=${BUILD_LOCATION}/extra_install

## locations for building/installing 'auxiliary' packages:
## these don't normally get rebuilt (unless --rebuild-aux is used)
AUX_SOURCES=${SOURCE_LOCATION}
AUX_INSTALL=${SOURCE_LOCATION}/aux_install


## --------------------------------------------------------------
## create required subdirectories for build
mkdir -p ${SOURCE_LOCATION} || fail
mkdir -p ${BUILD_LOCATION}  || fail
mkdir -p ${BUILD_LOCATION}/extra_sources || fail
mkdir -p ${AUX_INSTALL} || fail

# the following is only meant for Solaris:
# it assumes that a gcc is installed in /usr/local and thus the
# corresponding libs are in /usr/local/lib. Then we create a link to
# libstdc++.a, leaving out everything else, in particular the nasty
# libgcc_s.so
# for this to work we MUST NOT have -L/usr/local/lib in the LDFLAGS
# Bernd
if [ `uname` = SunOS -a -r /usr/local/lib/libstdc++.a -a \! -r ${AUX_INSTALL}/lib/libstdc++.a ] ; then
  mkdir -p ${AUX_INSTALL}/lib || fail
  ln -s /usr/local/lib/libstdc++.a ${AUX_INSTALL}/lib/libstdc++.a || fail
fi

## now adjust the PATH so that local installs will be found (autoconf,automake..):
PATH="${AUX_INSTALL}/bin:${PATH}"
echo export PKG_CONFIG_PATH="${BUILD_INSTALL}/lib/pkgconfig:${AUX_INSTALL}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export PKG_CONFIG_PATH="${BUILD_INSTALL}/lib/pkgconfig:${AUX_INSTALL}/lib/pkgconfig:${PKG_CONFIG_PATH}"
## ------------------------------------------------------------
## use pre-packaged LAL unless --new-lal was specified
LAL_TGZ="lal-4.0_current"
if [ "${eah_cvs_lal}" = yes ]; then
    lal_update_opts="$LAL_TAG"
    eah_LAL_PREFIX=${BUILD_LOCATION}/extra_install
elif [ "${eah_new_lal}" = yes ]; then
    eah_cvs_lal="yes"
    LAL_TAG="cvs-HEAD"
    lal_update_opts=-A
    eah_LAL_PREFIX=${BUILD_LOCATION}/extra_install
else
    LAL_TAG=${LAL_TGZ}.tar.gz
    eah_LAL_PREFIX=${AUX_INSTALL}
fi

## 'pin' BOINC-version to $BOINC_TAG unless --new-boinc was specified
if [ "${eah_new_boinc}" = yes ]; then
    BOINC_TAG="latest version (using HEAD)"
    BOINC_SVN_TAG="HEAD"
elif [ "${eah_boinc_rpdist}" = yes ]; then
    BOINC_TAG="using boinc-rpdist.tar.gz"
else
    boinc_update_opts="${BOINC_TAG}"
fi

if [ "$eah_win32_cross_build" = "yes" ]; then
    CPPFLAGS="-DMINGW_WIN32"
    CFLAGS="-gstabs3 $CFLAGS"
    CXXFLAGS="-gstabs3 $CXXFLAGS"
elif [ "$eah_win32_cygwin_build" = "yes" ]; then
    CPPFLAGS="-mno-cygwin -DMINGW_WIN32 $CPPFLAGS"
    CFLAGS="-gstabs3 -mno-cygwin $CFLAGS"
    CXXFLAGS="-gstabs3 -mno-cygwin $CXXFLAGS"
    LDFLAGS="-mno-cygwin $LDFLAGS"
else
    CFLAGS="-g $CFLAGS"
    CXXFLAGS="-g $CXXFLAGS"
fi

## ------------------------------------------------------------
## keep log-file of global build-options given to this script
## (allow appending to previous log-files) 
LOGFILE=${BUILD_LOCATION}/build.log
echo "----------------------------------------" >> "$LOGFILE"
echo "Date of build: "`date` >> "$LOGFILE"
echo "  as date tag: `date +%Y-%m-%dT%H:%M%z`"  >> "$LOGFILE"
echo >> "$LOGFILE"
## we have to strip the $'s from Id-string to avoid cvs-expansion
## in case the build-log gets checked into cvs!
## FIXME: the following code uses obsolete CVS ID tags.
## It should be modified to use git version information.
cvs_version=`echo '$Id$' | sed 's/[$]\(.*\)[$]/\1/g'`

echo "Version of the build-script: '${cvs_version}'" >> ${LOGFILE}
echo >> ${LOGFILE}
echo "Command-line: '${cmdline_log}'" >> ${LOGFILE}
echo >> ${LOGFILE}
echo "Built on: `hostname`" >> ${LOGFILE}
echo >> ${LOGFILE}

##----------------------------------------
## show user some summary of what we'll do

echo "System info:"				| tee -a "$LOGFILE"	
echo
echo "config_guess 	= ${eah_config_guess}" 	| tee -a "$LOGFILE"	
echo "build_cpu    	= ${build_cpu}" 	| tee -a "$LOGFILE"	
echo "build_vendor 	= ${build_vendor}" 	| tee -a "$LOGFILE"	
echo "build_os_full	= ${build_os_full}" 	| tee -a "$LOGFILE"	
echo "eah_system_name	= ${eah_system_name}"	| tee -a "$LOGFILE"	
echo | tee -a "$LOGFILE"


echo "Environment:" >> ${LOGFILE}
echo
echo "This script uses the following SHELL variables to regulate it's build"
echo
echo "CC                 =       ${CC}"			| tee -a ${LOGFILE}
echo "CXX                =       ${CXX}"		| tee -a ${LOGFILE}
echo "CFLAGS             =       ${CFLAGS}"		| tee -a ${LOGFILE}
echo "CXXFLAGS           =       ${CXXFLAGS}"		| tee -a ${LOGFILE}
echo "CPPFLAGS           =       ${CPPFLAGS}"		| tee -a ${LOGFILE}
echo "LDFLAGS            =       ${LDFLAGS}"		| tee -a ${LOGFILE}
echo "LIBS               =       ${LIBS}"		| tee -a ${LOGFILE}
echo "DEMODCC            =       ${DEMODCC}"		| tee -a ${LOGFILE}
echo "DEMODFLAGS         =       ${DEMODFLAGS}"		| tee -a ${LOGFILE}
echo "BOINCC             =       ${BOINCC}"		| tee -a ${LOGFILE}
echo "BUILD_LOCATION     =       ${BUILD_LOCATION}"	| tee -a ${LOGFILE}
echo "SOURCE_LOCATION    =       ${SOURCE_LOCATION}"	| tee -a ${LOGFILE}
echo "LAL_TAG            =       ${LAL_TAG}"		| tee -a ${LOGFILE}
echo "LAL_PATCH          =       ${LAL_PATCH}"		| tee -a ${LOGFILE}
if [ -n "$BOINC_SVN_TAG" ]; then
  echo "BOINC_SVN_TAG      =       ${BOINC_SVN_TAG}"	| tee -a ${LOGFILE}
else
  echo "BOINC_TAG          =       ${BOINC_TAG}"	| tee -a ${LOGFILE}
fi
echo "BOINC_PATCH        =       ${BOINC_PATCH}"	| tee -a ${LOGFILE}
echo "LALAPPS_TAG        =       ${LALAPPS_TAG}"	| tee -a ${LOGFILE}
echo "EAH_TAG            =       ${EAH_TAG}"		| tee -a ${LOGFILE}
echo "APP_NAME           =       ${APP_NAME}"		| tee -a ${LOGFILE}

if [ -n "${DEMODFLAGS}" -a -z  "${DEMODCC}" ] ; then
  echo
  echo "WARNING: DEMODFLAGS will be ignored if DEMODCC is not set!"
  echo "         Are you sure that this is what you want?"
fi

echo
echo "Step 1: autoconf, automake, libtool, pkg-config"
if [ "$eah_win32_cross_build" = "yes" ]; then
    echo "Step 1a: MinGW"
    echo "Step 1b: zlib, binutils"
fi
echo "Step 2: fftw3"
echo "Step 3: GSL"
echo "Step 4: BOINC"
echo "Step 5: LAL"
echo "Step 6: Einstein@Home"

echo
if [ -n "${eah_release}" ]; then
    echo "NOTE: Building a portable release-candidate! Version: ${eah_release}"
fi

echo
echo "If these choices are what you want, please type <ENTER>, else do <CONTROL-C>"
echo
echo "If you wish to start with a particular step N, just type N <ENTER>"
echo
read startstep

if test -z "${startstep}" ; then
    startstep=0
fi

echo "Selected start-step = ${startstep}"  >> ${LOGFILE}

## --------------------------------------------------------------
##
## prepare the general cmdline-arguments we'll use for running 
## the 'configure'-scripts
##
if [ -n "$CFLAGS" ]; then 
    eah_CFLAGS=`echo CFLAGS=\'${CFLAGS}\'`
fi
if [ -n "$CXXFLAGS" ]; then
    eah_CXXFLAGS=`echo CXXFLAGS=\'${CXXFLAGS}\'`
fi

## don't change the order of things below unless you REALLY know what you're doing
## e.g. LAL_PREFIX-paths need to come first, because they could be either in BUILD_INSTALL or AUX_INSTALL,
## depending on whether --new-lal was used or not. We have to make sure we use the right version if both should be available.
if [ "$eah_win32_cross_build" = "yes" ] || [ "$eah_win32_cygwin_build" = "yes" ]; then
    eah_cross_CPPFLAGS="-DEINSTEINATHOME_CROSS_BUILD"
fi
eah_CPPFLAGS=`echo CPPFLAGS=\'-I${BUILD_INSTALL}/include/BOINC -I${BUILD_INSTALL}/include/boinc -I${BUILD_INSTALL}/include -I${AUX_INSTALL}/include/bfd -I${AUX_INSTALL}/include ${CPPFLAGS} $eah_cross_CPPFLAGS\'`
eah_LDFLAGS=`echo LDFLAGS=\'-L${BUILD_INSTALL}/lib -L${AUX_INSTALL}/lib ${LDFLAGS}\'`
eah_LAL_FLAGS=`echo LAL_PREFIX=${eah_LAL_PREFIX}`

eah_configure_args="${eah_CPPFLAGS} ${eah_CFLAGS} ${eah_CXXFLAGS} ${eah_LDFLAGS} ${eah_LAL_FLAGS}"

## explicitely disable debugging when optimizing
if [ "${eah_optimize}" = yes ]; then
  eah_configure_args="$eah_configure_args --disable-debug"
fi

## BSD'sMacOSX needs some special flags
if [ "${eah_os}" = darwin ]; then
    eah_configure_args="${eah_configure_args} --with-apple-opengl-framework --with-x=no"
fi
## FreeBSD also seems to need some special flags...? (from Stein Sandbech)
if [ "${eah_os}" = freebsd ]; then
    eah_configure_args="${eah_configure_args} --with-x=no"
fi

if [ "$eah_win32_cross_build" = "yes" ]; then
    TARGET_HOST=i586-pc-mingw32
    BUILD_HOST="${eah_system_name}"
    CROSS_CONFIG_OPTS="--host=$TARGET_HOST --build=$BUILD_HOST"
    export TARGET_HOST BUILD_HOST CROSS_CONFIG_OPTS
elif [ "$eah_win32_cygwin_build" = "yes" ]; then
    TARGET_HOST=i586-pc-mingw32
    BUILD_HOST="${eah_system_name}"
    export TARGET_HOST BUILD_HOST
fi

if [ "$TARGET_HOST" = "i586-pc-mingw32" ] ; then
    eah_target_name=windows_intelx86
    eah_target_ext=.exe
else
    eah_target_name="${eah_system_name}"
fi

## ------------------------------
## OK. We are ready to start.

## ==============================
## run 6 main-steps
## ==============================
if [ "$eah_win32_cross_build" = "yes" ]; then
    if [ ${startstep} -le 1 ]; then step1; step1a; fi
    if [ ${startstep} -le 2 ]; then step1b; step2; fi
elif [ "$eah_win32_cygwin_build" = "yes" ]; then
    if [ ${startstep} -le 1 ]; then step1; fi
    if [ ${startstep} -le 2 ]; then step1b; step2; fi
elif [ "${eah_os}" = linux ]; then
    if [ ${startstep} -le 1 ]; then step1; fi
    if [ ${startstep} -le 2 ]; then step1b; step2; fi
else
    if [ ${startstep} -le 1 ]; then step1; fi
    if [ ${startstep} -le 2 ]; then step2; fi
fi
if [ ${startstep} -le 3 ]; then step3; fi
if [ ${startstep} -le 4 ]; then step4; fi
if [ ${startstep} -le 5 ]; then step5; fi
if [ ${startstep} -le 6 ]; then step6; fi

## ====================
## put final executables into build-directory
##====================

## ---------- copy binary into build-dir

f="${BUILD_LOCATION}/lalsuite/lalapps/src/pulsar/hough/src2/eah_HierarchicalSearch${eah_target_ext}"
test -f "$f" && cp -f "$f" "${BUILD_LOCATION}/cfsBOINC${eah_target_ext}"

## ---------- consistency-checks on final binary

## check consistency of debug/no-debug build between lal and lalapps
if [ "$eah_optimize" = "yes" -a -r "${BUILD_LOCATION}/cfsBOINC" ]; then
    laldebug=`nm "${BUILD_LOCATION}/cfsBOINC" | grep "[^X]LALFree"`;
    if [ -n "$laldebug" ]; then
	echo "laldebug = ${laldebug}"
	echo "================================================================================"
	echo "WARNING: it seems you linked a debug-LAL to a non-debug binary !!"
	echo "This will most likely lead to memory-checking problems if lalDebugLevel > 0"
	echo ""
	echo "It is *highly* recommended to rebuild LAL using '--rebuild-aux' in non-debug mode"
	echo "================================================================================"
	echo 
    fi
fi

## ==============================
## special treatment for release-binaries
##
## if --release: we also rename the binaries into 
## their canonical release-name: "einstein_$version_$eah_system_name[_nog]
## PLUS we do some portability checks under Linux
## ==============================
if [ -n "${eah_release}" ]; then
    ## rename into 'canonical' release-names
    eah_canonical_name=${APP_NAME}_${eah_release}_${eah_target_name}${eah_target_ext}
    eah_graphics_name=${APP_NAME}_${eah_release}_graphics_${eah_target_name}${eah_target_ext}
    cd ${BUILD_LOCATION}
    mv ./cfsBOINC${eah_target_ext} ./${eah_canonical_name}
    eah_release_binaries=${eah_canonical_name}
    if [ -r starsphere${eah_target_ext} ]; then
	mv ./starsphere${eah_target_ext} ./${eah_graphics_name}
	eah_release_binaries="${eah_release_binaries} ${eah_graphics_name}"
    fi
    if [ -r build.log ]; then
	eah_build_log="${eah_canonical_name}.log"
	cp ./build.log "./${eah_build_log}"
    fi
    ## now bundle the release-files into a tar.gz:
    eah_next="tar cf ${eah_canonical_name}.tar ${eah_release_binaries} ${eah_build_log}"
    echo "${eah_next}"
    eval ${eah_next} || fail
    eah_next="gzip ${eah_canonical_name}.tar"
    echo ${eah_next}
    eval ${eah_next} || fail 
    ## and move the final release-tar.gz into the 'pwd' where this script was started
    eah_next="mv ${eah_canonical_name}.tar.gz ${eah_here}"
    echo ${eah_next}
    eval ${eah_next}|| fail

    ## check portability of the linux-release
    if [ "${eah_os}" = linux ]; then
	glibc23_deps=`nm ${eah_release_binaries} | grep GLIBC_2.3`
	gcc33_deps=`nm ${eah_release_binaries} | grep GCC_3.3`
	if [ -n "${glibc23_deps}" ] || [ -n "${gcc33_deps}" ]; then
	    echo "================================================================================" 	| tee -a ${LOGFILE}
	    echo "WARNING: your linux release-binary contains GLIBC_2.3 and/or GCC_3.3 symbols:"	| tee -a ${LOGFILE}
	    echo ${glibc23_deps} 									| tee -a ${LOGFILE}
	    echo ${gcc33_deps}										| tee -a ${LOGFILE}	
	    echo
	    echo "==> This will seriously affect the portability of this binary."			| tee -a ${LOGFILE}
	    echo "It is not recommended to distribute this application under BOINC"			| tee -a ${LOGFILE}
	    echo "================================================================================"	| tee -a ${LOGFILE}
	fi ## glibc_2.3|gcc_3.3 symbols found
	ill_dynlibs=`ldd ${eah_release_binaries} | grep 'stdc++'`
	if [ -n "${ill_dynlibs}" ]; then
	    echo "================================================================================"	| tee -a ${LOGFILE}
	    echo "WARNING: your linux release-binary contains dynamic libstc++ dependencies"		| tee -a ${LOGFILE}
	    echo
	    echo "==> This will seriously affect the portability of this binary."			| tee -a ${LOGFILE}
	    echo "It is not recommended to distribute this application under BOINC"			| tee -a ${LOGFILE}
	    echo "================================================================================"	| tee -a ${LOGFILE}
	fi ## stdc++ dynamic dependencies found

    fi ## eah_os==linux

fi ## if eah_release

echo "**********************************************************************"
echo "*"
echo "* Ok, we're done. If everything went fine, you should find the"
echo "* Einstein@Home-executables in your BUILD_LOCATION (see above)"
echo "*"
if [ -n ${eah_release} ]; then
echo "*  You should also find the release-bundle in the local directory, "
echo "* under the name ${eah_canonical_name}.tar.gz"
echo "*"
fi
echo "**********************************************************************"
