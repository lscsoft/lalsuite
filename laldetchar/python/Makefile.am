BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = laldetchar

dist_bin_SCRIPTS = \
	laldetchar-idq-append-file-to-cache \
	laldetchar-idq-post-svm-train

if HAVE_PYTHON

pybin_scripts = \
	laldetchar-hveto \
	laldetchar-hveto-summary-page \
	laldetchar-idq-add-excluded-vars \
	laldetchar-idq-build-auxmvc-vectors \
	laldetchar-idq-convertkwtosb \
	laldetchar-idq-evaluate-times \
	laldetchar-idq-fap-calibration \
	laldetchar-idq-gdb-glitch-tables \
	laldetchar-idq-gdb-timeseries \
	laldetchar-idq-prepare-training-auxmvc-samples \
	laldetchar-idq-realtime \
	laldetchar-idq-select-channels-from-kw-config \
	laldetchar-idq-summary \
	laldetchar-idq-svm-evaluate \
	laldetchar-idq-svm-train \
	laldetchar-idq-train \
	laldetchar-segments-from-data \
	$(END_OF_LIST)

endif
