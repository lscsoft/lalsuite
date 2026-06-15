# ----------------------------------------------------------------------
# LALSuite GitLab-CI: initialise build for RPM packages
# ----------------------------------------------------------------------

# parallelise the build (rpmbuild uses this variable by default)
export RPM_BUILD_NCPUS=${CPU_COUNT}

# configure dnf to cache packages locally
dnf config-manager --save --setopt="cachedir=${RPM_CACHE_DIR}" --setopt="keepcache=1"
