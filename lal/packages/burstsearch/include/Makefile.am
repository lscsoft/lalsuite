## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Thresholds.h ExcessPower.h TFTransform.h

all-local: $(top_srcdir)/include/Thresholds.h $(top_srcdir)/include/ExcessPower.h $(top_srcdir)/include/TFTransform.h

$(top_srcdir)/include/Thresholds.h: Thresholds.h
	cp -f Thresholds.h $(top_srcdir)/include

$(top_srcdir)/include/TFTransform.h: TFTransform.h
	cp -f TFTransform.h $(top_srcdir)/include

$(top_srcdir)/include/ExcessPower.h: ExcessPower.h
	cp -f ExcessPower.h $(top_srcdir)/include

