###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../publicIncludes -m32
AM_FFLAGS = -I. -I../publicIncludes -m32
AM_CXXFLAGS = $(AM_CFLAGS)

AM_LDFLAGS  = -L.
AM_LDFLAGS += -L../../Libraries/dap
AM_LDFLAGS += -L../../Libraries/qreo
AM_LDFLAGS += -L../../Libraries/s2

LDADD = -lreorder_qreo
LDADD += -lreorder_dap
LDADD += -lreorder_s2

# target libraries

bin_PROGRAMS = qreou

# sources

qreou_SOURCES = \
	fofclibPlus.c \
	qreoroot_uf.F \
	qufap.F \
	unixqufap.F


