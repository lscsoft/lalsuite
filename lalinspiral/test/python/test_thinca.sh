#!/bin/sh

BUILDDIR=${LAL_TEST_BUILDDIR:-"."}
SRCDIR=${LAL_TEST_SRCDIR:-"."}

for i in 1 2 ; do
	echo
	echo "====== Testing --min-instruments=${i} ======"
	echo

	# reset input files
	for f in thinca_min${i}.xml.gz ; do echo "reset $f" ; "cp" ${SRCDIR}/test_thinca_input.xml.gz ${f} ; done
	for f in brute_min${i}.xml.gz ; do echo "reset $f" ; "cp" ${SRCDIR}/test_thinca_input.xml.gz ${f} ; done

	# run thinca and brute force algorithm
	lalinspiral_thinca --verbose --threshold 0.005 --min-instruments ${i} thinca_min${i}.xml.gz
	${BUILDDIR}/thinca_brute_force_coinc --verbose --delta-t 0.005 --min-instruments ${i} brute_min${i}.xml.gz

	# print total event counts
	for f in thinca_min${i}.xml.gz ${SRCDIR}/test_thinca_correct_min${i}.xml.gz ; do echo "coincs in ${f}:  $(ligolw_print -t coinc_event ${f} | wc -l)" ; done
	for f in brute_min${i}.xml.gz ; do echo "coincs in ${f}:  $(ligolw_print -t coinc_event ${f} | wc -l)" ; done

	# compare thinca and brute force outputs to correct result
	${BUILDDIR}/thinca_compare_coinc_lists --verbose thinca_min${i}.xml.gz ${SRCDIR}/test_thinca_correct_min${i}.xml.gz
	${BUILDDIR}/thinca_compare_coinc_lists --verbose brute_min${i}.xml.gz ${SRCDIR}/test_thinca_correct_min${i}.xml.gz
done
