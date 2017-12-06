EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

# Add compiled test programs to this variable
test_programs += ConfigFileTest
test_programs += LALMath3DPlotTest
test_programs += LALMathNDPlotTest
test_programs += PrintFTSeriesTest
test_programs += PrintVectorTest
test_programs += ReadFTSeriesTest
test_programs += SegmentsIOTest
test_programs += StreamInputTest
test_programs += StreamSeriesInputTest
test_programs += UserInputTest

# Add shell, Python, etc. test scripts to this variable
test_scripts += GzipTest.sh

# Add any helper programs required by tests to this variable
test_helpers += GzipTest

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	ConfigFile.cfg \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	SegmentsOutput1.data \
	$(END_OF_LIST)

EXTRA_DIST += \
	ConfigFileSample.cfg \
	ConfigFileSample2.cfg \
	ConfigFileSample3.cfg.gz \
	SegmentsInput1.data \
	StreamInput.data \
	StreamSeriesInput.data
	$(END_OF_LIST)
