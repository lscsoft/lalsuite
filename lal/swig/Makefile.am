include $(top_srcdir)/gnuscripts/lallibs_swig.am

swig_function_prefixes = XLAL LAL

swig_module_lal_libs = \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la \
	$(top_builddir)/lib/lal/liblal.la \
	$(END_OF_LIST)

pkginclude_HEADERS += \
	SWIGCommon.i \
	SWIGOctave.i \
	SWIGPython.i \
	$(END_OF_LIST)
