#!/bin/sh
# script to build lalsuite (to be run by Metronome on the execute node)

# this script expects:
#
# - pwd is the directory under which lalsuite should be built
# - condor installed in PATH (and CONDOR_CONFIG defined if necessary)
# - $NMI_git_repo defined (by Metronome, from our submit file)
# - $NMI_git_id defined (by Metronome, from our submit file)

# wrap this whole script in a block combining stdout & stderr to ease debugging
{
    # make sure DEBUG is defined
    if [[ -z "$DEBUG" ]]; then
	DEBUG=0
    fi

    # exit immediately if any command exits with a non-zero status.
    set -e
    # treat unset variables as an error when performing parameter expansion.
    set -u
    # print (unexpanded) shell input lines as they are read
    set -v
    # print (expanded) commands before they're executed
    set -x

    if [ -d head ]; then
	echo Moving old install out of the way...
	mv -v head head.old.$$
    fi

    # most of what's below is adapted from
    # <https://www.lsc-group.phys.uwm.edu/daswg/docs/howto/lal-install.html>
    # (circa 2010-02-22)

    BASE_DIR=$PWD
    export LSCSOFT_SRCDIR=$BASE_DIR/src
    export LSCSOFT_ROOTDIR=$BASE_DIR/opt/lscsoft

    mkdir -p ${LSCSOFT_SRCDIR}
    cd ${LSCSOFT_SRCDIR}

    # if local repo doesn't already exist, create it
    if [[ ! -d ${LSCSOFT_SRCDIR}/lalsuite/.git ]]; then
	cd ${LSCSOFT_SRCDIR}
	# making a mirror and then a local clone ensures all refs are available
	# avoids "reference is not a tree" errors
	# (e.g., for a84c74b7169788f3879722943872298f15ec222f)
	git clone --mirror $NMI_git_repo
	git clone ./lalsuite
    fi

    cd ${LSCSOFT_SRCDIR}/lalsuite
#    git checkout -f $NMI_git_id -b test.$NMI_git_branch
    git checkout $NMI_git_id
#    git clean -dqfx

    # hack to work around metaio change
    # see <https://sympa.ligo.org/wws/arc/daswg/2013-01/msg00875.html>
    # 1359608400 is approximate (log sluething or git bisect is needed to find exact date of change)
    if [[ $(git log -1 --format='%ct' $NMI_git_id) -lt 1359608400 ]]; then
	set +u
	export PKG_CONFIG_PATH=/opt/metaio-8.2.1/lib/pkgconfig:$PKG_CONFIG_PATH
	set -u
    fi

    mkdir -p ${LSCSOFT_ROOTDIR}/etc
    echo "export LSCSOFT_LOCATION=${LSCSOFT_ROOTDIR}" > ${LSCSOFT_ROOTDIR}/etc/lscsoftrc

    # for debugging, but also to detect and fail if Condor isn't in the PATH
    echo "condor_config_val LIB =" $(condor_config_val LIB)

    for subsys in lal lalframe lalmetaio lalsimulation lalburst laldetchar lalinspiral lalpulsar lalxml lalinference lalstochastic lalapps pylal glue; do
# 
	if [[ ! -d ${LSCSOFT_SRCDIR}/lalsuite/$subsys ]]; then
            echo "Hmm, no $subsys source dir found; moving on..."
            continue;
	fi
	SUBSYS=$(echo $subsys | tr '[:lower:]' '[:upper:]')
	SUBSYS_PREFIX=$LSCSOFT_ROOTDIR/$subsys
	eval export ${SUBSYS}_PREFIX=$LSCSOFT_ROOTDIR/$subsys
	cd ${LSCSOFT_SRCDIR}/lalsuite/$subsys
        
	set +u
	. ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	set -u
	
	if [[ $subsys == "pylal" ]]; then
	    python setup.py build
            python setup.py install --prefix=$SUBSYS_PREFIX
	elif [[ $subsys == "glue" ]]; then
            # nasty kludge to work around glue install bug known to exist in
            # s6abc_lowmass tag 85f56e4a1555a60fe2ee98dde0b6e22afface3ad
            if [[ -L glue/misc/generate_vcs_info.py ]]; then
		rm glue/misc/generate_vcs_info.py
		cp lal/lib/generate_vcs_info.py glue/misc/
		touch glue/misc/__init__.py
            fi
            python setup.py install --prefix=$SUBSYS_PREFIX
	else
            ./00boot
            if [[ $subsys == "lalapps" ]]; then
    		./configure --prefix=$SUBSYS_PREFIX --enable-condor
	    else
#                ./configure --prefix=$SUBSYS_PREFIX --enable-swig-python
		./configure --prefix=$SUBSYS_PREFIX
            fi
	    if [[ "$DEBUG" != 0 ]]; then
  		make -j V=1 install
	    else
		make -j2 V=1 install
	    fi
	fi

	echo "# setup $SUBSYS for development:  " >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	echo "export ${SUBSYS}_LOCATION=\$LSCSOFT_LOCATION/$subsys" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	echo "if [ -f "\$${SUBSYS}_LOCATION/etc/${subsys}-user-env.sh" ]; then" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	echo "  source \$${SUBSYS}_LOCATION/etc/${subsys}-user-env.sh" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	echo "fi" >> ${LSCSOFT_ROOTDIR}/etc/lscsoftrc
	cd -
    done

    # useful info to have logged for later debugging
    env
    ${LSCSOFT_ROOTDIR}/lal/bin/lal-version

# end of stdout/stderr-combining block
} 2>&1
