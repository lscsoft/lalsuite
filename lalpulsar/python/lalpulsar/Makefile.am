BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	piecewise_model \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	git_version.py \
	lineFileParser.py \
	simulateCW.py \
	NstarTools.py \
	PulsarParametersWrapper.py \
	simulateHeterodynedCW.py \
	knope_utils.py \
	pulsarhtmlutils.py \
	pulsarpputils.py \
	public_sft_directory.py \
	metric_utils.py \
	$(END_OF_LIST)

pymodule_DATA = \
	public_sft_directory_README.md \
	$(END_OF_LIST)


if DOXYGEN

BUILT_SOURCES += public_sft_directory_README.dox
MOSTLYCLEANFILES += public_sft_directory_README.dox
EXTRA_DIST += public_sft_directory_README.dox
public_sft_directory_README.dox: public_sft_directory_README.md Makefile
	$(AM_V_GEN)echo "/**" > $@; \
	echo "\page public_sft_directory_readme Public SFT Directory README" >>$@; \
	echo "\$(MARKDOWN2HTML_DOXYBLOCK)" >>$@; \
	$(MARKDOWN2HTML) $< >>$@; \
	echo "\end$(MARKDOWN2HTML_DOXYBLOCK)" >>$@; \
	echo "*/" >>$@

endif

endif
