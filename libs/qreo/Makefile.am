###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../includes -DNETCDF -m32
AM_FFLAGS = -m32
AM_CXXFLAGS = $(AM_CFLAGS)
AM_LDFLAGS = -L.

# target libraries

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

