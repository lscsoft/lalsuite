## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a

VCREATEC = \
        I2CreateVector.c I4CreateVector.c \
        I8CreateVector.c \
        U2CreateVector.c U4CreateVector.c \
        U8CreateVector.c \
        SCreateVector.c DCreateVector.c CreateVector.c \
        CCreateVector.c ZCreateVector.c CHARCreateVector.c 

VDESTROYC = \
        I2DestroyVector.c I4DestroyVector.c \
        I8DestroyVector.c \
        U2DestroyVector.c U4DestroyVector.c \
        U8DestroyVector.c \
        SDestroyVector.c DDestroyVector.c DestroyVector.c \
        CDestroyVector.c ZDestroyVector.c CHARDestroyVector.c 

ACREATEC = \
        I2CreateArray.c I4CreateArray.c \
        I8CreateArray.c \
        U2CreateArray.c U4CreateArray.c \
        U8CreateArray.c \
        SCreateArray.c DCreateArray.c CreateArray.c \
        SCreateArray.c DCreateArray.c CreateArray.c \
        CCreateArray.c ZCreateArray.c

ADESTROYC = \
        I2DestroyArray.c I4DestroyArray.c \
        I8DestroyArray.c \
        U2DestroyArray.c U4DestroyArray.c \
        U8DestroyArray.c \
        SDestroyArray.c DDestroyArray.c DestroyArray.c \
        CDestroyArray.c ZDestroyArray.c


VSCREATEC = \
        I2CreateVectorSequence.c I4CreateVectorSequence.c \
        I8CreateVectorSequence.c \
        U2CreateVectorSequence.c U4CreateVectorSequence.c \
        U8CreateVectorSequence.c \
        SCreateVectorSequence.c DCreateVectorSequence.c \
        CreateVectorSequence.c \
        CCreateVectorSequence.c ZCreateVectorSequence.c \
        CHARCreateVectorSequence.c

VSDESTROYC = \
        I2DestroyVectorSequence.c I4DestroyVectorSequence.c \
        I8DestroyVectorSequence.c \
        U2DestroyVectorSequence.c U4DestroyVectorSequence.c \
        U8DestroyVectorSequence.c \
        SDestroyVectorSequence.c DDestroyVectorSequence.c \
        DestroyVectorSequence.c \
        CDestroyVectorSequence.c ZDestroyVectorSequence.c \
        CHARDestroyVectorSequence.c

liblal_a_SOURCES = $(VCREATEC) $(VDESTROYC) $(ACREATEC) $(ADESTROYC) \
                   $(VSCREATEC) $(VSDESTROYC)
MAINTAINERCLEANFILES = $(VCREATEC) $(VDESTROYC) $(ACREATEC) $(ADESTROYC) \
                       $(VSCREATEC) $(VSDESTROYC)
EXTRA_DIST = \
	CreateVector.m4 DestroyVector.m4 \
        CreateArray.m4 DestroyArray.m4 \
	CreateVectorSequence.m4 DestroyVectorSequence.m4

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

VTYPECODES = Z C D S I2 I4 I8 U2 U4 U8 CHAR ""
ATYPECODES = Z C D S I2 I4 I8 U2 U4 U8 ""

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

srcs:
	for type in $(VTYPECODES) ; do \
	m4 -DTYPECODE=$${type} DestroyVector.m4 > $${type}DestroyVector.c ; \
	m4 -DTYPECODE=$${type} CreateVector.m4 > $${type}CreateVector.c ; \
	m4 -DTYPECODE=$${type} CreateVectorSequence.m4 \
		> $${type}CreateVectorSequence.c ; \
	m4 -DTYPECODE=$${type} DestroyVectorSequence.m4 \
		> $${type}DestroyVectorSequence.c ; \
	done
	for type in $(ATYPECODES) ; do \
	m4 -DTYPECODE=$${type} DestroyArray.m4 > $${type}DestroyArray.c ; \
	m4 -DTYPECODE=$${type} CreateArray.m4 > $${type}CreateArray.c ; \
	done
