# ----------------------------------------------------------------------
# LALSuite GitLab-CI: upgrade .deb packages
# ----------------------------------------------------------------------

# update APT cache
${LCI_SCRIPTS}/retry apt-get -y -q update

# upgrade distribution
${LCI_SCRIPTS}/retry apt-get -y -q upgrade

# install latest release
${LCI_SCRIPTS}/retry apt-get -y -q install lalsuite lalsuite-dev

# remove lalsuite meta-packages
dpkg -r lalsuite lalsuite-dev

# create local repo for upstream .debs
upstream_debs=$(find ${PACKAGE_ROOT_DIR} -name '*.deb')
echo "===== upstream .debs"
printf "%s\n" ${upstream_debs}
echo "====="
mkdir -pv /srv/local-apt-repository
cp -v ${upstream_debs} /srv/local-apt-repository
/usr/lib/local-apt-repository/rebuild

# upgrade all packages
${LCI_SCRIPTS}/retry apt-get -y -q update
${LCI_SCRIPTS}/retry apt-get -y dist-upgrade
