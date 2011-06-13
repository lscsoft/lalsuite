# SWIG makefile for generating header lists
# Author: Karl Wette, 2011

# link this makefile to Makefile.in for automatic remaking
$(srcdir)/Makefile.in : $(top_srcdir)/swig/swig-header.mk

# name of SWIG interface header include file
swig_iface_header = $(top_builddir)/swig/swig-iface-header-$(subst /,-,$(subdir)).swg

# generate SWIG interface header include file
all install : $(swig_iface_header)
$(swig_iface_header) : $(swig_headers)
	@echo "// $(subdir)" > $@
	@for file in $(swig_headers); do echo 'swiglal_include(lal/'"$$file"')' >> $@; done

# cleanup
swig_clean = $(swig_iface_header)
