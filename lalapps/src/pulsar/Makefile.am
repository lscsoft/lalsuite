## Process this file with automake to produce Makefile.in
# $Id$

SUBDIRS = FDS_isolated Injections

if FRAME
FRAMEPROGS = lalapps_make_job_input_data lalapps_make_sfts
else
FRAMEPROGS = 
endif

bin_PROGRAMS = $(FRAMEPROGS) 

if FRAME
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_make_job_input_data_SOURCES = make_job_input_data.c
lalapps_make_sfts_SOURCES = make_sfts.c
endif

EXTRA_DIST = \
make_locked_data_list \
make_cluster_data_list \
make_job_input_data.c \
make_sfts.c \
README \
SFT.00000 \
SFT.00001 \
SFT.00002 \
SFT.00003 \
SFT.00004 \
SFT.00005 \
SFT.00006 \
SFT.00007 \
SFT.00008 \
SFT.00009 \
S1-H1segments.txt \
S1-H2segments.txt \
S1-L1segments.txt \
S2-H1segments.txt \
S2-H2segments.txt \
S2-L1segments.txt
INCLUDES = -I$(top_srcdir)/src
