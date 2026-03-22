# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install RPM dependencies for top-level build
# ----------------------------------------------------------------------

# collect build dependencies from RPM spec files
specfiles=
subdir=
for subdir in ${TOP_BUILDDEP_SUBDIRS}; do

    # find RPM spec file
    specfile=$(ls -1 ${subdir}/*.spec)

    # remove LALSuite dependencies on itself
    sed -i -e '/BuildRequires: *\(lal\|liblal\|python.*-lal\)/d' ${specfile}

    specfiles="${specfiles} ${specfile}"

done

# install build dependencies
${LCI_SCRIPTS}/retry dnf builddep -y -q ${specfiles}
