cp -v ${dir}H1-POWER_S5-830046410-50274.xml.gz bucluster_test.xml.gz && \
lalapps_bucluster --verbose --cluster-algorithm excesspower --program power bucluster_test.xml.gz && \
rm -vf bucluster_test.xml.gz
