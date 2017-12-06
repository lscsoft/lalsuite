BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = laldetchar

if HAVE_PYTHON

pybin_scripts = \
	laldetchar-hveto \
	laldetchar-hveto-summary-page \
	laldetchar-idq-add-excluded-vars \
	laldetchar-idq-build-auxmvc-vectors \
	laldetchar-idq-convertkwtosb \
	laldetchar-idq-evaluate-times \
	laldetchar-idq-calibration \
	laldetchar-idq-gdb-alert-handler \
	laldetchar-idq-gdb-glitch-tables \
	laldetchar-idq-gdb-processor \
	laldetchar-idq-gdb-timeseries \
	laldetchar-idq-prepare-training-auxmvc-samples \
	laldetchar-idq-realtime \
	laldetchar-idq-select-channels-from-kw-config \
	laldetchar-idq-summary \
	laldetchar-idq-svm-evaluate \
	laldetchar-idq-svm-train \
	laldetchar-idq-train \
	laldetchar-segments-from-data \
	laldetchar-idq-ann-convert \
	laldetchar-idq-ann-train \
	laldetchar-idq-ann-evaluate \
	laldetchar-idq-frames-to-segments \
	laldetchar-idq-gst-sender \
	laldetchar-idq-gst-relay \
	laldetchar-idq-gst-listener \
	laldetchar-idq-gdb-local-performance \
	$(END_OF_LIST)

endif
