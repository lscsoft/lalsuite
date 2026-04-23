# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install upstream Conda packages
# ----------------------------------------------------------------------

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
