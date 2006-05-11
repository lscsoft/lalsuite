## Process this file with automake to produce Makefile.in
TESTS = NoTest.sh
EXTRA_DIST = NoTest.sh RingTest.c
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = .dvi-dep *.out *.dat snr-*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
