check : glue_ligolw_ilwd_verify iterutils_verify lal_verify ligolw_test01 ligolw_test03a ligolw_test03b ligolw_test03c ligolw_test04 ligolw_test05 ligolw_sqlite_test offsetvector_verify segmentsUtils_verify segments_verify test_ligolw_array test_ligolw_lsctables test_ligolw_table test_ligolw_tokenizer test_ligolw_utils test_ligolw_utils_segments
	@echo "All Tests Passed"

define printpassfail
	{ echo "Pass" ; true ; } || { echo "Fail" ; false ; }
endef

SERVEFILESPORT := 16543

define servefilesstart
	{ python -c "import SocketServer, SimpleHTTPServer ; SocketServer.TCPServer(('127.0.0.1', $(SERVEFILESPORT)), SimpleHTTPServer.SimpleHTTPRequestHandler).serve_forever()" & sleep 1 ; }
endef
define servefilesstop
	{ killall python ; sleep 1 ; true ; }
endef

glue_ligolw_ilwd_verify :
	@echo "=== start glue_ligolw_ilwd_verify ===>"
	python glue_ligolw_ilwd_verify.py && $(printpassfail)
	@echo "<=== end glue_ligolw_ilwd_verify ==="

iterutils_verify :
	@echo "=== start iterutils_verify ===>"
	python iterutils_verify.py && $(printpassfail)
	@echo "<=== end iterutils_verify ==="

lal_verify :
	@echo "=== start lal_verify ===>"
	python lal_verify.py && $(printpassfail)
	@echo "<=== end lal_verify ==="

ligolw_test01 :
	@echo "=== start ligolw_test01 ===>"
	# test reading array data from a file
	python ligo_lw_test_01.py && $(printpassfail)
	@echo "<=== end ligolw_test01 ==="

ligolw_test03a :
	@echo "=== start ligolw_test03a ===>"
	# test reading compressed and non-compressed files from stdin
	{ cat inspiral_event_id_test_in1.xml.gz | ligolw_print >/dev/null ; } && $(printpassfail)
	{ cat inspiral_event_id_test_in2.xml | ligolw_print >/dev/null ; } && $(printpassfail)
	@echo "<=== end ligolw_test03a ==="

ligolw_test03b :
	@echo "=== start ligolw_test03b ===>"
	# test reading compressed and non-compressed files
	ligolw_print inspiral_event_id_test_in1.xml.gz >/dev/null && $(printpassfail)
	ligolw_print inspiral_event_id_test_in2.xml >/dev/null && $(printpassfail)
	@echo "<=== end ligolw_test03b ==="

ligolw_test03c :
	@echo "=== start ligolw_test03c ===>"
	# test reading compressed and non-compressed files from http:// URLs
	$(servefilesstart) && { ligolw_print http://127.0.0.1:$(SERVEFILESPORT)/inspiral_event_id_test_in1.xml.gz >/dev/null && $(printpassfail) ; } ; $(servefilesstop)
	$(servefilesstart) && { ligolw_print http://127.0.0.1:$(SERVEFILESPORT)/inspiral_event_id_test_in2.xml >/dev/null && $(printpassfail) ; } ; $(servefilesstop)
	@echo "<=== end ligolw_test03c ==="

ligolw_test04 :
	@echo "=== start ligolw_test04 ===>"
	{ ligolw_cut --delete-table sngl_inspiral <inspiral_event_id_test_in1.xml.gz | cmp ligolw_cut_proof.xml ; } && $(printpassfail)
	@echo "<=== end ligolw_test04 ==="

ligolw_test05 :
	@echo "=== start ligolw_test05 ===>"
	# make sure XML writing code compresses files
	FILENAME=$(shell mktemp --suffix .xml.gz) && { ligolw_add --output $${FILENAME} ligolw_sqlite_test_input.xml.gz && gunzip --test $${FILENAME} && $(printpassfail) ; } ; rm -f $${FILENAME}
	@echo "<=== end ligolw_test05 ==="

ligolw_sqlite_test :
	@echo "=== start ligolw_sqlite_test ===>"
	sh ligolw_sqlite_test.sh && $(printpassfail)
	@echo "<=== end ligolw_sqlite_test ==="

offsetvector_verify :
	@echo "=== start offsetvector_verify ===>"
	python offsetvector_verify.py && $(printpassfail)
	@echo "<=== end offsetvector_verify ==="

segmentsUtils_verify :
	@echo "=== start segmentsUtils_verify ===>"
	python segmentsUtils_verify.py && $(printpassfail)
	@echo "<=== end segmentsUtils_verify ==="

segments_verify :
	@echo "=== start segments_verify ===>"
	python segments_verify.py && $(printpassfail)
	@echo "<=== end segments_verify ==="

test_ligolw_array :
	@echo "=== start test_ligolw_array ===>"
	./test_ligolw_array.py && $(printpassfail)
	@echo "<=== end test_ligolw_array ==="

test_ligolw_lsctables :
	@echo "=== start test_ligolw_lsctables ===>"
	./test_ligolw_lsctables.py && $(printpassfail)
	@echo "<=== end test_ligolw_lsctables ==="

test_ligolw_table :
	@echo "=== start test_ligolw_table ===>"
	./test_ligolw_table.py && $(printpassfail)
	@echo "<=== end test_ligolw_table ==="

test_ligolw_tokenizer :
	@echo "=== start test_ligolw_tokenizer ===>"
	./test_ligolw_tokenizer.py && $(printpassfail)
	@echo "<=== end test_ligolw_tokenizer ==="

test_ligolw_utils :
	@echo "=== start test_ligolw_utils ===>"
	./test_ligolw_utils.py && $(printpassfail)
	@echo "<=== end test_ligolw_utils ==="

test_ligolw_utils_segments :
	@echo "=== start test_ligolw_utils_segments ===>"
	./test_ligolw_utils_segments.py && $(printpassfail)
	@echo "<=== end test_ligolw_utils_segments ==="

clean :
	rm -f ligo_lw_test_01*png
