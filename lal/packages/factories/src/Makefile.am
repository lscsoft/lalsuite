## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a

VCREATEC = \
        I2CreateVector.c I4CreateVector.c \
        I8CreateVector.c \
        SCreateVector.c DCreateVector.c CreateVector.c \
        CCreateVector.c ZCreateVector.c CHARCreateVector.c 

VDESTROYC = \
        I2DestroyVector.c I4DestroyVector.c \
        I8DestroyVector.c \
        SDestroyVector.c DDestroyVector.c DestroyVector.c \
        CDestroyVector.c ZDestroyVector.c CHARDestroyVector.c 

VSCREATEC = \
        I2CreateVectorSequence.c I4CreateVectorSequence.c \
        I8CreateVectorSequence.c \
        SCreateVectorSequence.c DCreateVectorSequence.c \
        CreateVectorSequence.c \
        CCreateVectorSequence.c ZCreateVectorSequence.c \
        CHARCreateVectorSequence.c

VSDESTROYC = \
        I2DestroyVectorSequence.c I4DestroyVectorSequence.c \
        I8DestroyVectorSequence.c \
        SDestroyVectorSequence.c DDestroyVectorSequence.c \
        DestroyVectorSequence.c \
        CDestroyVectorSequence.c ZDestroyVectorSequence.c \
        CHARDestroyVectorSequence.c

liblal_a_SOURCES = $(VCREATEC) $(VDESTROYC) $(VSCREATEC) $(VSDESTROYC)
MAINTAINERCLEANFILES = $(VCREATEC) $(VDESTROYC) $(VSCREATEC) $(VSDESTROYC)
EXTRA_DIST = \
	CreateVector.m4 DestroyVector.m4 \
	CreateVectorSequence.m4 DestroyVectorSequence.m4

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

TYPECODES = Z C D S I2 I4 I8 CHAR ""

.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

srcs:
	for type in $(TYPECODES) ; do \
	m4 -DTYPECODE=$${type} DestroyVector.m4 > $${type}DestroyVector.c ; \
	m4 -DTYPECODE=$${type} CreateVector.m4 > $${type}CreateVector.c ; \
	m4 -DTYPECODE=$${type} CreateVectorSequence.m4 \
		> $${type}CreateVectorSequence.c ; \
	m4 -DTYPECODE=$${type} DestroyVectorSequence.m4 \
		> $${type}DestroyVectorSequence.c ; \
	done
