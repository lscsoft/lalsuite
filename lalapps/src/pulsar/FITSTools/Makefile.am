include $(top_srcdir)/gnuscripts/lalapps.am

AM_CPPFLAGS += $(PAGER_CPPFLAGS)

if CFITSIO

bin_PROGRAMS = \
	lalapps_fitscopy \
	lalapps_imarith \
	lalapps_imcopy \
	lalapps_imlist \
	lalapps_imstat \
	lalapps_listhead \
	lalapps_liststruc \
	lalapps_modhead \
	lalapps_tabcalc \
	lalapps_tablist \
	lalapps_tabmerge \
	lalapps_tabselect \
        $(END_OF_LIST)

lalapps_fitscopy_SOURCES = fitscopy.c
lalapps_imarith_SOURCES = imarith.c
lalapps_imcopy_SOURCES = imcopy.c
lalapps_imlist_SOURCES = imlist.c
lalapps_imstat_SOURCES = imstat.c
lalapps_listhead_SOURCES = listhead.c
lalapps_liststruc_SOURCES = liststruc.c
lalapps_modhead_SOURCES = modhead.c
lalapps_tabcalc_SOURCES = tabcalc.c
lalapps_tablist_SOURCES = tablist.c
lalapps_tabmerge_SOURCES = tabmerge.c
lalapps_tabselect_SOURCES = tabselect.c

endif
