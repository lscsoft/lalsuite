## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = inject.tex \
  inject_phase_diagram.eps inject_phase_diagram.pdf \
  inject_ecliptic.eps inject_ecliptic.pdf \
  inject_galactic.eps inject_galactic.pdf \
  inject_geodetic.eps inject_geodetic.pdf \
  inject_geographic.eps inject_geographic.pdf \
  inject_horizon.eps inject_horizon.pdf \
  inject_lat_long.eps inject_lat_long.pdf
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
