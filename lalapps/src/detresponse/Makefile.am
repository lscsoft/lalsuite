## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_detresponse lalapps_test_gridding lalapps_test_earth lalapps_bar_response lalapps_ifo_response lalapps_ifo_skyavg_ts
lalapps_detresponse_SOURCES = detresponse.c detresponse.h util.c util.h set_source_params.h set_source_params.c generate_timeseries_response.c compute_skygrid.c cmdline.c cmdline.h skygrid.c skygrid.h
lalapps_test_gridding_SOURCES = test_gridding.c util.c util.h skygrid.c skygrid.h make_gridding.c make_gridding.h
lalapps_test_earth_SOURCES = test_earth.c util.c util.h skygrid.c skygrid.h
lalapps_bar_response_SOURCES = bar_response.c util.c util.h set_source_params.h set_source_params.c skygrid.c skygrid.h make_gridding.h make_gridding.c
lalapps_ifo_response_SOURCES = ifo_response.c util.c util.h set_source_params.h set_source_params.c skygrid.c skygrid.h make_gridding.h make_gridding.c
lalapps_ifo_skyavg_ts_SOURCES = ifo_skyavg_ts.c util.c util.h set_source_params.h set_source_params.c skygrid.c skygrid.h make_gridding.h make_gridding.c
# lalapps_dopplergrid_SOURCES = util.c util.h dopplergrid.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I.
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

check: $(bin_PROGRAMS)
	./lalapps_detresponse -S -NTEST -r0. -d1.57 -o1.023 -Dtest -s2001 -u2048 -i600. -v4 -e7
	./lalapps_detresponse -W -NTEST -o1.023 -Dtest -s2001 -u2048 -i600. -v4 -e7
