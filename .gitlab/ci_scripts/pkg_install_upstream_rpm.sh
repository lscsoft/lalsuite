# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install upstream RPM packages
# ----------------------------------------------------------------------

upstream_rpms=$(find ${PACKAGE_ROOT_DIR} -name '*.rpm')
if [ "X${upstream_rpms}" != X ]; then
    echo "===== upstream RPMs"
    printf "%s\n" ${upstream_rpms}
    echo "====="
    ${LCI_SCRIPTS}/retry dnf install -y ${upstream_rpms}
fi
