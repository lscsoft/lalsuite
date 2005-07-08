all: doxygen refman

doxygen:
	echo -n "Collecting figures from LAL ... "
	./collect_figures.sh
	echo -n "done. Now running doxygen ... "
	doxygen lal-doxy.cfg
	echo "done."


refman: 
	cd latex && make refman.pdf  
	ln -fs latex/refman.pdf . >& /dev/null
	echo "done".

clean: 
	rm -f ./figures/*.eps figures/*.png;
	rm -rf latex/* html/*
