## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = StochasticCrossCorrelation.h \
	 	     Dirichlet.h

#	FalseAlarm.h FalseDismissal.h FDismissalVsFAlarm.h \


all-local: $(top_builddir)/include/lal/StochasticCrossCorrelation.h   \
	   $(top_builddir)/include/lal/Dirichlet.h

#          $(top_builddir)/include/lal/FalseAlarm.h  \
#	   $(top_builddir)/include/lal/FalseDismissal.h  \
#	   $(top_builddir)/include/lal/FDismissalVsFAlarm.h  \

$(top_builddir)/include/lal/StochasticCrossCorrelation.h: StochasticCrossCorrelation.h
	cp -f $(srcdir)/StochasticCrossCorrelation.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Dirichlet.h: Dirichlet.h
	cp -f $(srcdir)/Dirichlet.h $(top_builddir)/include/lal


# $(top_builddir)/include/lal/FalseAlarm.h: FalseAlarm.h
# 	cp -f $(srcdir)/FalseAlarm.h $(top_builddir)/include/lal

# $(top_builddir)/include/lal/FalseDismissal.h: FalseDismissal.h
# 	cp -f $(srcdir)/FalseDismissal.h $(top_builddir)/include/lal

# $(top_builddir)/include/lal/FDismissalVsFAlarm.h: FDismissalVsFAlarm.h
# 	cp -f $(srcdir)/FDismissalVsFAlarm.h $(top_builddir)/include/lal
