cp -v ${LAL_TEST_SRCDIR}/H1-POWER_S5-830046410-50274.xml.gz bucluster_test.xml.gz && \
lalburst_cluster --verbose --cluster-algorithm excesspower --program power bucluster_test.xml.gz && \
rm -vf bucluster_test.xml.gz
