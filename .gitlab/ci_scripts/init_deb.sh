# ----------------------------------------------------------------------
# LALSuite GitLab-CI: initialise build for .deb packages
# ----------------------------------------------------------------------

# configure the APT cache
cat <<EOF >/etc/apt/apt.conf.d/99cache.conf
Dir::Cache "${APT_CACHE_DIR}";
EOF
mkdir -p "${APT_CACHE_DIR}/archives/partial"

# set packaging environment
export DEBFULLNAME="${GITLAB_USER_NAME}"
export DEBEMAIL="${GITLAB_USER_EMAIL}"
