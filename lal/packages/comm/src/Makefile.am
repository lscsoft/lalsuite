## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a

if MPI
liblal_a_SOURCES = Comm.c
else
liblal_a_SOURCES = 
endif

EXTRA_DIST       = Comm.c

CC=@MPICC@

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

