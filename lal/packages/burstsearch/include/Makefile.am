## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Thresholds.h ExcessPower.h TFTransform.h

all-local: $(top_builddir)/include/Thresholds.h $(top_builddir)/include/ExcessPower.h $(top_builddir)/include/TFTransform.h

$(top_builddir)/include/Thresholds.h: Thresholds.h
	cp -f $(srcdir)/Thresholds.h $(top_builddir)/include

$(top_builddir)/include/TFTransform.h: TFTransform.h
	cp -f $(srcdir)/TFTransform.h $(top_builddir)/include

$(top_builddir)/include/ExcessPower.h: ExcessPower.h
	cp -f $(srcdir)/ExcessPower.h $(top_builddir)/include

