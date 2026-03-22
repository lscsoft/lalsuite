# ----------------------------------------------------------------------
# LALSuite GitLab-CI: build package using plain make
# ----------------------------------------------------------------------

# -- setup -------------------------------------------------------------

# set up environment
for userenv in $(find ${PACKAGE_ROOT_DIR} -name 'lal*-user-env.sh'); do
    source ${userenv}
done

# install build dependencies
TOP_BUILDDEP_SUBDIRS="."
source ${LCI_SCRIPTS}/top_builddep_${LCI_PKG}.sh

# set ./configure flags
export LCI_CONFIGURE_FLAGS_PKG_BUILD_MAKE="--prefix=${PACKAGE_DIR}"

# -- build -------------------------------------------------------------

# setup build environment
source ${LCI_SCRIPTS}/build_env.sh

# configure
eval ./configure "${LCI_CONFIGURE_FLAGS}"

# build, test, and install package
make -j${CPU_COUNT} all
make -j${CPU_COUNT} check
make -j${CPU_COUNT} install
make -j${CPU_COUNT} installcheck
