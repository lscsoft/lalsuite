## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Thresholds.h ExcessPower.h TFTransform.h

all-local: $(top_builddir)/include/lal/Thresholds.h $(top_builddir)/include/lal/ExcessPower.h $(top_builddir)/include/lal/TFTransform.h

$(top_builddir)/include/lal/Thresholds.h: Thresholds.h
	cp -f $(srcdir)/Thresholds.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/TFTransform.h: TFTransform.h
	cp -f $(srcdir)/TFTransform.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/ExcessPower.h: ExcessPower.h
	cp -f $(srcdir)/ExcessPower.h $(top_builddir)/include/lal

