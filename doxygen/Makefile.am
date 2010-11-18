all: html

.PHONY: getfigs html pdf clean

getfigs:
	echo -n "Collecting figures from LAL ... "
	./collect_figures.sh
	echo -n "done."

html: getfigs
	echo -n "Running doxygen to produce HTML ... "
	doxygen doxygen-html.cfg
	echo "done."

pdf: getfigs
	echo -n "Running doxygen to produce PDF ... "
	doxygen doxygen-pdf.cfg
	echo -n "done. Now compiling the pdf manual 'refman.pdf' ..."
	cd latex && pool_size=4000000 make refman.pdf
	ln -fs latex/refman.pdf .
	echo "done."

clean:
	rm -rf refman.pdf ./latex ./html ./figures
