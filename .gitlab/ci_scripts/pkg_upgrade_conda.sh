# ----------------------------------------------------------------------
# LALSuite GitLab-CI: upgrade Conda packages
# ----------------------------------------------------------------------

# ensure we are creating a brand new environment
export CONDA_ENVS_PATH="${LCI_TMPDIR}/conda/envs"
rm -rf ${CONDA_ENVS_PATH}

# install latest release
${LCI_SCRIPTS}/retry conda create -n upgrade-test "python=${LCI_PYTHON_VERSION}" ${LCI_PKGLIST}

# copy upstream Conda packages to a local channel
upstream_condas=$(find ${PACKAGE_ROOT_DIR} -name '*.conda')
echo "===== upstream Conda packages"
printf "%s\n" ${upstream_condas}
echo "====="
local_channel="${LCI_TMPDIR}/conda/local-channel"
rm -rf ${local_channel}/
mkdir -p ${local_channel}/
cp -rv ${PACKAGE_ROOT_DIR}/*/*/ ${local_channel}/
${LCI_SCRIPTS}/retry conda index "${local_channel}"
conda config --add channels "${local_channel}"
${LCI_SCRIPTS}/retry conda search "*lal*" --channel "${local_channel}" --override-channels

# upgrade all packages
${LCI_SCRIPTS}/retry conda update -n upgrade-test --use-local ${LCI_PKGLIST}
