all: doxygen

doxygen:
	echo -n "Collecting figures from LAL ... "
	./collect_figures.sh
	echo -n "done. Now running doxygen ... "
	doxygen lsd_doxy.cfg
	echo "done."


refman:
	cd latex && make refman.pdf
	ln -fs latex/refman.pdf .
	echo "done".

clean:
	rm -rf ./latex ./html ./figures
