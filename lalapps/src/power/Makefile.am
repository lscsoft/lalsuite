include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

BURSTPROGS = lalapps_binj
lalapps_binj_SOURCES = binj.c
if LALINSPIRAL
INSPPROGS = lalapps_power
lalapps_power_SOURCES = power.c
endif

bin_PROGRAMS = $(BURSTPROGS) $(INSPPROGS)

if HAVE_PYTHON
pybin_scripts = \
	lalapps_binj_pic \
	lalapps_binjfind \
	lalapps_bucluster \
	lalapps_bucut \
	lalapps_burca \
	lalapps_burca_tailor \
	lalapps_cafe \
	lalapps_gen_timeslides \
	lalapps_plot_tisi \
	lalapps_power_calc_likelihood \
	lalapps_power_final \
	lalapps_power_likelihood_pipe \
	lalapps_power_online_pipe \
	lalapps_power_pipe \
	lalapps_power_plot_binjtf \
	lalapps_power_plot_burca2 \
	lalapps_power_plot_burst \
	lalapps_power_plot_burstrate \
	lalapps_power_plot_binj \
	lalapps_power_plot_burca \
	lalapps_power_plot_detresponse \
	lalapps_power_veto \
	lalapps_xml_plotlalseries
pkgpython_PYTHON = power.py
endif

pkgdata_DATA = power_pipe.ini power_online_pipe.ini

EXTRA_DIST += \
	H1_L1_dt125_df128.png \
	Z_histogram.png \
	clustering_example.png \
	conditioning.png \
	freqdomainfilter.png \
	h_histogram.png \
	power_pipeline.png \
	plotbinj_L1_4.png \
	plotbinj_L1_5.png \
	plotburca2-P-793154128-2525061.png \
	power.dox \
	power_online_pipe.ini \
	power_pipe.ini \
	psds.png \
	sk_histogram.png \
	sksk.png \
	tiles_histogram.png \
	tiles_histogram_adjusted.png \
	timedomainfilter_04hz.png \
	timedomainfilter_16hz.png \
	windows.png \
	$(END_OF_LIST)
