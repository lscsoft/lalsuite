## Process this file with automake to produce Makefile.in
# $Id$

SUBDIRS = FDS_isolated Injections FDS_griphyn hough FDS_binary

if FRAME
FRAMEPROGS = lalapps_make_job_input_data lalapps_make_sfts lalapps_Calibrate60sSFTs lalapps_Calibrate1800SFTs lalapps_FreqAverager
else
FRAMEPROGS = lalapps_Calibrate60sSFTs lalapps_Calibrate1800SFTs lalapps_FreqAverager
endif

bin_PROGRAMS = $(FRAMEPROGS)
#bin_PROGRAMS = $(FRAMEPROGS) SFTfileIOTest

#SFTfileIOTest_SOURCES = SFTfileIOTest.c SFTfileIO.c SFTfileIO.h 

if FRAME
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_make_job_input_data_SOURCES = make_job_input_data.c
lalapps_make_sfts_SOURCES = make_sfts.c
endif

lalapps_Calibrate60sSFTs_SOURCES = Calibrate60sSFTs.c Calibrate60sSFTs.h
lalapps_Calibrate1800SFTs_SOURCES = Calibrate1800SFTs.c Calibrate1800SFTs.h
lalapps_FreqAverager_SOURCES = FreqAverager.c FreqAverager.h

EXTRA_DIST = \
make_locked_data_list \
make_cluster_data_list \
make_job_input_data.c \
make_sfts.c \
Calibrate60sSFTs.c \
Calibrate60sSFTs.h \
Calibrate1800SFTs.c \
Calibrate1800SFTs.h \
FreqAverager.c \
FreqAverager.h \
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
