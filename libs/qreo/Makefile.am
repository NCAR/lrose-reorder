###############################################
# makefile.am template for automake
###############################################

# compile flags

AM_CFLAGS = -I.
AM_CFLAGS += -I../includes
AM_CFLAGS += -I../../../include
AM_CFLAGS += -I$(prefix)/include
AM_CFLAGS += -m32
AM_CFLAGS += -DNETCDF

AM_FFLAGS = $(AM_CFLAGS)
AM_CXXFLAGS = $(AM_CFLAGS)

# target library

lib_LIBRARIES = libreorder_qreo.a

# sources

libreorder_qreo_a_SOURCES = \
	fofclib.c \
	foflib.F \
	qreo1.F \
	qreo2.F \
	qreo3.F \
	unixflib.F \
	unixqreo.F \
	xcedcio.c

