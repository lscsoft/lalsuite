# ----------------------------------------------------------------------
# LALSuite GitLab-CI: top-level build using plain make
# ----------------------------------------------------------------------

# install build dependencies
TOP_BUILDDEP_SUBDIRS="${LCI_PKGLIST}"
source ${LCI_SCRIPTS}/top_builddep_${LCI_PKG}.sh

# setup build environment
source ${LCI_SCRIPTS}/build_env.sh

# configure
eval ./configure "${LCI_CONFIGURE_FLAGS}"

# build targets
for target in ${MAKE_TARGETS:-distcheck}; do
    make -j${CPU_COUNT} ${target}
done
