# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install upstream .deb packages
# ----------------------------------------------------------------------

upstream_debs=$(find ${PACKAGE_ROOT_DIR} -name '*.deb')
if [ "X${upstream_debs}" != X ]; then
    echo "===== upstream .debs"
    printf "%s\n" ${upstream_debs}
    echo "====="
    ${LCI_SCRIPTS}/retry apt-get -y -q install ${upstream_debs}
fi
