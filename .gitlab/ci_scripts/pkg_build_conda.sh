# ----------------------------------------------------------------------
# LALSuite GitLab-CI: build Conda package
# ----------------------------------------------------------------------

# where to build things
export CONDA_BLD_PATH="${PWD}/conda-bld"

# render YAML file to use our tarball
sha256=$(openssl dgst -r -sha256 ${TARBALL} | cut -d' ' -f1)
tar -xf ${TARBALL} --wildcards ${PACKAGE}-*/conda/ --strip-components=1
sed "s|@TARBALL@|${TARBALL}|g;s|@SHA256@|${sha256}|g" conda/meta.yaml.in > conda/meta.yaml

# create a feedstock from the conda recipe
git config --global user.name "${GITLAB_USER_NAME}"
git config --global user.email "${GITLAB_USER_EMAIL}"
${LCI_SCRIPTS}/retry conda smithy init conda/ --feedstock-directory ${PACKAGE}-feedstock
pushd ${PACKAGE}-feedstock/

# handle migrations that are bundled with the tarball
mkdir -p .ci_support/migrations
find recipe/migrations -type f -name "*.yaml" -exec cp --update=none --verbose {} .ci_support/migrations/ \;

# regenerate the feedstock
${LCI_SCRIPTS}/retry conda smithy regenerate --no-check-uptodate
git ls-files
conda smithy recipe-lint --conda-forge || true  # lint, but don't fail

# specify `--python` to build only the Python version of the Conda container
export CONDA_BUILD_ARGS="--python \"${LCI_PYTHON_VERSION} *_${LCI_CONDA_ABI}\""

echo "======================== Conda configuration ========================="

# try to find configuration file using just platform tag
plat_tag="${LCI_CONDA_PLAT_BASE}"
pkg_plat_tag_name="LCI_CONDA_PLAT_${PACKAGE}"
if [ "X${!pkg_plat_tag_name}" != X ]; then
    plat_tag="${plat_tag}_${!pkg_plat_tag_name}"
    echo "Platform tag: ${plat_tag} from LCI_CONDA_PLAT_BASE and ${pkg_plat_tag_name}"
else
    echo "Platform tag: ${plat_tag} from LCI_CONDA_PLAT_BASE"
fi
conda_config_file="$(echo .ci_support/${plat_tag}*.yaml)"
if [ -f "${conda_config_file}" ]; then
    echo "Found configuration file '${conda_config_file}'"
else
    echo "INFO: '${conda_config_file}' does not match a unique configuration file"

    # try to find configuration file also using Python and ABI tags
    echo "Python tag: python${LCI_PYTHON_VERSION}"
    echo "ABI tag: ${LCI_CONDA_ABI}"
    conda_config_file="$(echo .ci_support/${plat_tag}_python${LCI_PYTHON_VERSION}*_${LCI_CONDA_ABI}.yaml)"
    if [ -f "${conda_config_file}" ]; then
        echo "Found configuration file '${conda_config_file}'"
    else
        echo "INFO: '${conda_config_file}' does not match a unique configuration file"
        echo "ERROR: could not determine Conda configuration"
        exit 1
    fi

fi

echo "----------------------------------------------------------------------"

# setup build environment
source ${LCI_SCRIPTS}/build_env.sh

# build packages
# - we use eval here because CONDA_BUILD_ARGS contains multiple spaces
eval conda build \
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
