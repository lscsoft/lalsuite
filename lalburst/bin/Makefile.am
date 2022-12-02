BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- C programs -------------

bin_PROGRAMS = \
	lalburst_version \
	$(END_OF_LIST)

lalburst_version_SOURCES = version.c

TESTS += \
	lalburst_version \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	lalburst_cluster \
	lalburst_coinc \
	lalburst_cut \
	lalburst_gen_timeslides \
	lalburst_injfind \
	lalburst_inj_pic \
	lalburst_plot_tisi \
	lalburst_power_calc_likelihood \
	lalburst_power_final \
	lalburst_power_meas_likelihood \
	lalburst_power_plot_binj \
	lalburst_power_plot_binjtf \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif

# -- Manuals ----------------

man1_MANS = $(help2man_MANS)
