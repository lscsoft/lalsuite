include $(top_srcdir)/gnuscripts/lallibs_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

noinst_HEADERS = \
	ArrayFactoriesTest_source.c \
	ArraySequenceFactoriesTest_source.c \
	VectorFactoriesTest_source.c \
	VectorSequenceFactoriesTest_source.c

check_PROGRAMS = \
	ArrayFactoriesTest \
	ArraySequenceFactoriesTest \
	VectorFactoriesTest \
	StringVectorTest \
	VectorSequenceFactoriesTest

TESTS = \
	$(check_PROGRAMS)
