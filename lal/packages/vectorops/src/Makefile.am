## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a
liblal_a_SOURCES = \
        CCVectorDivide.c \
        CCVectorMultiply.c \
        CCVectorMultiplyConjugate.c \
        CVectorAbs.c \
        CVectorAngle.c \
        UnwrapREAL4Angle.c \
        ZZVectorDivide.c   \
        ZZVectorMultiply.c \
        ZZVectorMultiplyConjugate.c \
        ZVectorAbs.c \
        ZVectorAngle.c \
        UnwrapREAL8Angle.c \
        SCVectorMultiply.c \
        SSVectorMultiply.c \
        DZVectorMultiply.c \
        DDVectorMultiply.c

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

