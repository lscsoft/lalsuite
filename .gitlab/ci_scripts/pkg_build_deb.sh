# ----------------------------------------------------------------------
# LALSuite GitLab-CI: build .deb package
# ----------------------------------------------------------------------

# install upstream .debs
upstream_debs=$(find ${PACKAGE_ROOT_DIR} -name '*.deb')
if [ "X${upstream_debs}" != X ]; then
    echo "===== upstream .debs"
    printf "%s\n" ${upstream_debs}
    echo "====="
    mkdir -pv /srv/local-apt-repository
    cp -v ${upstream_debs} /srv/local-apt-repository
    /usr/lib/local-apt-repository/rebuild
fi

# update APT cache
${LCI_SCRIPTS}/retry apt-get -y -q update

# copy tarball and rename to orig tarball
suffix=$(basename ${TARBALL} | sed 's/.*\.\(tar\..*\)/\1/')
version=$(basename ${TARBALL} | sed 's/[^-]*-\(.*\)\.tar\..*/\1/' | tr '-' '~')
cp -v ${TARBALL} ../${PACKAGE}_${version}.orig.${suffix}

# update changelog
dch -v "${version}-1+local" -b 'Rebuilt automatically on git.ligo.org CI'

# install build dependencies
${LCI_SCRIPTS}/retry mk-build-deps \
    --tool "apt-get -y -q -o Debug::pkgProblemResolver=yes --no-install-recommends" \
    --install \
    --remove
rm -rfv *-build-deps_*

# enable CUDA support
if [ "X${LCI_CUDA}" = Xyes ]; then
    LALSUITE_DEB_WITH_CUDA=/usr/local/cuda
else
    LALSUITE_DEB_WITH_CUDA=
fi

# build .deb packages
debuild \
    --prepend-path=/usr/lib/ccache \
    --set-envvar="CCACHE_DIR=${CCACHE_DIR}" \
    --set-envvar="LALSUITE_DEB_WITH_CUDA=${LALSUITE_DEB_WITH_CUDA}" \
    -us -uc -r \
    --lintian-opts --color=always --allow-root

# create a codeclimate report from the lintian report
echo "---- *lintian-overrides"
find debian/ -name '*lintian-overrides' | sort | xargs grep -H .
echo "---- *lintian-overrides"
lintian --allow-root --fail-on none --info --color never ../${PACKAGE}_*.changes > lintian.out
python3 -m lintian_codeclimate lintian.out -d . -o code-quality-lintian.json
sed -i \
    -e "s|debian\/control|debian/control.in|g" \
    -e "s|debian\/|$PACKAGE/debian/|g" \
    code-quality-lintian.json

# print package info
cd ..
for debf in *.deb; do
    echo "===== ${debf}";
    dpkg --info "${debf}";
    dpkg --contents "${debf}";
done
echo "====="

# move .deb packages
mv -v *.deb *.dsc *.changes ${PACKAGE_DIR}/
