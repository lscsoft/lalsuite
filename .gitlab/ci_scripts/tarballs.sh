# ----------------------------------------------------------------------
# LALSuite GitLab-CI: make tarballs for each subpackage
# ----------------------------------------------------------------------

tarball="${CI_JOB_GROUP_NAME##*:}"
echo "tarball=${tarball}"
case ${tarball} in

    # make the tarball for LAL only
    # - this job will run much faster than the `tarballs:pkg` job so we use
    #   it to release the LAL package-level build jobs as early as possible
    lal)
        pushd lal/
        ./00boot
        source ${LCI_SCRIPTS}/build_env.sh
        eval ./configure "${LCI_CONFIGURE_FLAGS}"
        make dist
        mv *.tar.* ${TARBALL_DIR}/
        ;;

    # make tarballs for all packages except LAL
    pkg)
        ./00boot
        source ${LCI_SCRIPTS}/build_env.sh
        eval ./configure "${LCI_CONFIGURE_FLAGS}"
        for subdir in ${LCI_PKGLIST_X_LAL}; do
            pushd ${subdir}
            make dist
            mv *.tar.* ${TARBALL_DIR}/
            popd
        done
        ;;

    top)
        ./00boot
        source ${LCI_SCRIPTS}/build_env.sh
        eval ./configure "${LCI_CONFIGURE_FLAGS}"
        make dist
        mv *.tar.* ${TARBALL_DIR}/
        ;;

    *)
        echo "ERROR: unknown tarball choice '${tarball}'"
        exit 1
        ;;

esac
