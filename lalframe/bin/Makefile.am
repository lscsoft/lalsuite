BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- C programs -------------

bin_PROGRAMS = \
	lalfr-cksum \
	lalfr-stat \
	lalfr-dump \
	lalfr-print \
	lalfr-fmt \
	lalfr-cat \
	lalfr-split \
	lalfr-cut \
	lalfr-paste \
	lalfr-stream \
	lalfr-vis \
	lalframe_version \
	$(END_OF_LIST)

lalfr_cksum_SOURCES = cksum.c
lalfr_stat_SOURCES = stat.c
lalfr_dump_SOURCES = dump.c
lalfr_print_SOURCES = print.c
lalfr_fmt_SOURCES = fmt.c
lalfr_cat_SOURCES = cat.c utils.c utils.h
lalfr_split_SOURCES = split.c utils.c utils.h
lalfr_cut_SOURCES = cut.c utils.c utils.h
lalfr_paste_SOURCES = paste.c utils.c utils.h
lalfr_stream_SOURCES = stream.c
lalfr_vis_SOURCES = vis.c
lalframe_version_SOURCES = version.c

TESTS += \
	lalframe_version \
	$(END_OF_LIST)

lalfr_MANS = \
	lalfr-cksum.1 \
	lalfr-stat.1 \
	lalfr-dump.1 \
	lalfr-print.1 \
	lalfr-fmt.1 \
	lalfr-cat.1 \
	lalfr-split.1 \
	lalfr-cut.1 \
	lalfr-paste.1 \
	lalfr-stream.1 \
	lalfr-vis.1 \
	$(END_OF_LIST)

man1_MANS = \
	$(lalfr_MANS) \
	$(END_OF_LIST)

if HELP2MAN
man1_MANS += \
	lalframe_version.man \
	$(END_OF_LIST)
endif

EXTRA_DIST += \
	$(lalfr_MANS) \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif
