## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a
liblal_a_SOURCES = \
  CreateTimeFreqParam.c \
  CreateTimeFreqRep.c \
  DestroyTimeFreqParam.c \
  DestroyTimeFreqRep.c \
  Dwindow.c \
  TfrPswv.c \
  TfrRsp.c \
  TfrSp.c \
  TfrWv.c

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)
