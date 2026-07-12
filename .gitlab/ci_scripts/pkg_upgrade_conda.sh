# ----------------------------------------------------------------------
# LALSuite GitLab-CI: upgrade Conda packages
# ----------------------------------------------------------------------

# remove package pins
conda config --remove-key pinned_packages

# activate environment for testing package upgrade
conda activate lalsuite-ci-upgrade
conda list --name lalsuite-ci-upgrade

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
${LCI_SCRIPTS}/retry conda search --quiet --channel "${local_channel}" --override-channels "*lal*"

# upgrade all packages
pkg_spec=""
for pkg_name in ${LCI_PKGLIST}; do
    pkg_version=$(conda list --name lalsuite-ci-upgrade --full-name ${pkg_name} --json | yq -r '.[0].version')
    pkg_spec="${pkg_spec} ${pkg_name}>${pkg_version}"
done
echo "===== Conda package spec"
echo ${pkg_spec}
echo "====="
${LCI_SCRIPTS}/retry conda install --quiet --name lalsuite-ci-upgrade ${pkg_spec}
lalapps_version
