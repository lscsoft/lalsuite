# ----------------------------------------------------------------------
# LALSuite GitLab-CI: upgrade .deb packages
# ----------------------------------------------------------------------

# create local repo for upstream .debs
upstream_debs=$(find ${PACKAGE_ROOT_DIR} -name '*.deb')
echo "===== upstream .debs"
printf "%s\n" ${upstream_debs}
echo "====="
mkdir -pv /srv/local-apt-repository
cp -v ${upstream_debs} /srv/local-apt-repository
/usr/lib/local-apt-repository/rebuild
${LCI_SCRIPTS}/retry apt-get -y -q update

# upgrade all packages
${LCI_SCRIPTS}/retry apt-get -y -q install \
    $(printf "lib%s-dev " ${LCI_PKGLIST_X_LALAPPS}) \
    $(printf "python3-%s " ${LCI_PKGLIST_X_LALAPPS}) \
    ${LCI_PKGLIST}
lalapps_version
