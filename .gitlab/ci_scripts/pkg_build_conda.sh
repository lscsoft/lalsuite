# ----------------------------------------------------------------------
# LALSuite GitLab-CI: build Conda package
# ----------------------------------------------------------------------

# copy upstream Conda packages to a local channel
upstream_condas=$(find ${PACKAGE_ROOT_DIR} -name '*.conda')
if [ "X${upstream_condas}" != X ]; then
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
fi

# where to build things
export CONDA_BLD_PATH="${PWD}/conda-bld"

# render YAML file to use our tarball
sha256=$(openssl dgst -r -sha256 ${TARBALL} | cut -d\  -f1)
tar -xf ${TARBALL} --wildcards ${PACKAGE}-*/conda/ --strip-components=1
sed "s|@TARBALL@|${TARBALL}|g;s|@SHA256@|${sha256}|g" conda/meta.yaml.in > conda/meta.yaml

# create a feedstock from the conda recipe
git config --global user.name "${GITLAB_USER_NAME}"
git config --global user.email "${GITLAB_USER_EMAIL}"
${LCI_SCRIPTS}/retry conda smithy init conda/ --feedstock-directory ${PACKAGE}-feedstock
pushd ${PACKAGE}-feedstock/

# handle migrations that are bundled with the tarball
mkdir -p .ci_support/migrations
find recipe/migrations -type f -name "*.yaml" -exec cp -nv {} .ci_support/migrations/ \;

# regenerate the feedstock
${LCI_SCRIPTS}/retry conda smithy regenerate --no-check-uptodate
git ls-files
conda smithy recipe-lint --conda-forge || true  # lint, but don't fail

# specify `--python` to build only the Python version of the Conda container
export CONDA_BUILD_ARGS="--python \"${LCI_PYTHON_VERSION} *_cp*\""


echo "======================== Conda configuration ========================="

# use package-specific configuration
conda_config_name="LCI_CONDA_CFG_${PACKAGE}"
if [ "X${!conda_config_name}" != X ]; then
    LCI_CONDA_CFG="${!conda_config_name}"
    echo "Using configuration ${LCI_CONDA_CFG} from ${conda_config_name}"
else
    echo "Using configuration ${LCI_CONDA_CFG} from LCI_CONDA_CFG"
fi

conda_config_file=".ci_support/${LCI_CONDA_CFG}.yaml"
if [ -f ${conda_config_file} ]; then
    echo "Found configuration file ${conda_config_file}"
else

    # try to find Python-specific configuration
    conda_config_file="$(echo .ci_support/${LCI_CONDA_CFG}*python${LCI_PYTHON_VERSION}cp*.yaml)"
    if [ -f ${conda_config_file} ]; then
        echo "Found configuration file ${conda_config_file}"
    else
        echo "ERROR: ${LCI_CONDA_CFG} does not match any configuration files"
        exit 1
    fi

fi

echo "----------------------------------------------------------------------"

# setup build environment
source ${LCI_SCRIPTS}/build_env.sh

# build packages
# - we use eval here because CONDA_BUILD_ARGS contains multiple spaces
eval ${LCI_SCRIPTS}/retry --max-try 2 conda build \
    recipe/ \
    --dirty \
    --error-overlinking \
    --keep-old-work \
    --no-anaconda-upload \
    --prefix-length 128 \
    --variant-config-files ${conda_config_file} \
    ${CONDA_BUILD_ARGS}

popd

# move Conda packages
package_dirs=$(find ${CONDA_BLD_PATH} -name '*.conda' -printf '%h\n' | sort | uniq)
mv -v ${package_dirs} ${PACKAGE_DIR}/

# save Conda-specific build products
mkdir conda-specific/
mv -v ${PACKAGE}-feedstock conda-specific/
