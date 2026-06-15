# ----------------------------------------------------------------------
# LALSuite GitLab-CI: upgrade RPM packages
# ----------------------------------------------------------------------

# create local repo for upstream RPMs
local_repo="${CI_PROJECT_DIR}/local_repo"
upstream_rpms=$(find ${PACKAGE_ROOT_DIR} -name '*.rpm')
echo "===== upstream RPMs"
printf "%s\n" ${upstream_rpms}
echo "====="
mkdir -p ${local_repo}
cp -v ${upstream_rpms} ${local_repo}
createrepo --quiet --workers "${CPU_COUNT}" "${local_repo}"
cat > /etc/yum.repos.d/local_repo.repo <<EOF
[local_repo]
name=Local repo
baseurl=file://${local_repo}
enabled=1
gpgcheck=0
EOF

# upgrade all packages
${LCI_SCRIPTS}/retry dnf -y upgrade
lalapps_version
