## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FlatMesh.h PulsarTimes.h StackMetric.h Resample.h \
                     ComputeSky.h LALDemod.h

all-local: $(top_builddir)/include/lal/FlatMesh.h   \
           $(top_builddir)/include/lal/PulsarTimes.h \
           $(top_builddir)/include/lal/StackMetric.h \
           $(top_builddir)/include/lal/Resample.h \
           $(top_builddir)/include/lal/ComputeSky.h \
           $(top_builddir)/include/lal/LALDemod.h

$(top_builddir)/include/lal/FlatMesh.h: FlatMesh.h
	cp -f $(srcdir)/FlatMesh.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/PulsarTimes.h: PulsarTimes.h
	cp -f $(srcdir)/PulsarTimes.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/StackMetric.h: StackMetric.h
	cp -f $(srcdir)/StackMetric.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Resample.h: Resample.h
	cp -f $(srcdir)/Resample.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/ComputeSky.h: ComputeSky.h
	cp -f $(srcdir)/ComputeSky.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALDemod.h: LALDemod.h
	cp -f $(srcdir)/LALDemod.h $(top_builddir)/include/lal
