# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install .deb dependencies for top-level build
# ----------------------------------------------------------------------

# update APT cache
${LCI_SCRIPTS}/retry apt-get -y -q update

# install build dependencies
for subdir in ${TOP_BUILDDEP_SUBDIRS}; do
    pushd ${subdir}

    # remove LALSuite dependencies on itself
    awk '
         /^Build-Depends:/ { bd=1; print; next }
         bd==1 && ( / lal/ || / liblal/ || / python.*-lal/ ) { next }
         /^[^ ]/ { bd=0 }
         { print }
        ' debian/control > debian/control.tmp
    mv -f debian/control.tmp debian/control

    # install build dependencies
    ${LCI_SCRIPTS}/retry apt-get build-dep -y -q .

    popd

done
