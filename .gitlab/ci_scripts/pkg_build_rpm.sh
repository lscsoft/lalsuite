# ----------------------------------------------------------------------
# LALSuite GitLab-CI: build RPM package
# ----------------------------------------------------------------------

# install upstream RPMs
upstream_rpms=$(find ${PACKAGE_ROOT_DIR} -name '*.rpm')
if [ "X${upstream_rpms}" != X ]; then
    echo "===== upstream RPMs"
    printf "%s\n" ${upstream_rpms}
    echo "====="
    ${LCI_SCRIPTS}/retry dnf install -y ${upstream_rpms}
fi

# where to build things
topdir="${PWD}/"

# build src.rpm
rpmbuild -ts --define "_topdir ${topdir}" ${TARBALL}
srpm=$(ls -1 ${topdir}/SRPMS/${PACKAGE}-*.src.rpm)

# install build dependencies
${LCI_SCRIPTS}/retry dnf builddep -y ${srpm}

# print installed packages
dnf list installed --quiet

# setup build environment
source ${LCI_SCRIPTS}/build_env.sh

# build binary RPMs
rpmbuild --rebuild --noclean \
    --define "_topdir ${topdir}" \
    --define '__spec_clean_body %{nil}' \
    ${srpm}
rpms=$(ls -1 ${topdir}/RPMS/*/*${PACKAGE}-*.rpm)

# print package info
for rpmf in ${rpms}; do
    echo "===== ${rpmf}"
    rpm -qlp "${rpmf}"
    echo "Files:"
    rpm -qip "${rpmf}"
    echo "Provides:"
    rpm -qp --provides "${rpmf}"
    echo "Requires:"
    rpm -qp --requires "${rpmf}"
done
echo "====="

# install the packages we built
${LCI_SCRIPTS}/retry dnf -y -q install ${rpms}

# lint packages (RPM files only)
echo "----- rpmlintrc"
grep . rpmlintrc
echo "----- rpmlintrc"
rpmlint -f rpmlintrc ${rpms} > rpmlint.out || true
cat rpmlint.out

# create a codeclimate report from the rpmlint report
python3 -m rpmlint_codeclimate rpmlint.out -s ${PACKAGE}.spec -o code-quality-rpmlint.json
sed -i -e "s|${PACKAGE}.spec|${PACKAGE}/${PACKAGE}.spec.in|g" code-quality-rpmlint.json

# move RPM packages
mv -v ${rpms} ${PACKAGE_DIR}/
