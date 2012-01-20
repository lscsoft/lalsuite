## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libvectorops.la

libvectorops_la_SOURCES = \
	MatrixDivide.c \
	MatrixMultiply.c \
	MatrixPower.c \
	MiscMatlab.c \
	VectorMultiply.c \
	VectorPolar.c \
	VectorIndexRange.c

noinst_HEADERS = \
	MatrixDivide_source.c \
	MatrixMultiply_source.c \
	MatrixPower_source.c \
	MiscMatlab_source.c \
	VectorIndexRange_source.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
