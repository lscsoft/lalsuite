## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
bin_PROGRAMS = lalapps_make_job_input_data lalapps_make_sfts
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_make_job_input_data_SOURCES = make_job_input_data.c
lalapps_make_sfts_SOURCES = make_sfts.c
else
bin_PROGRAMS = 
endif
EXTRA_DIST =  lalapps_make_job_input_data \
lalapps_make_sfts \
make_cluster_data_list \
make_job_input_data.c \
make_locked_data_list \
make_sfts.c \
README \
S1-H1segments.txt \
S1-H2segments.txt \
S1-L1segments.txt
INCLUDES = -I$(top_srcdir)/src
