###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../includes -m32
AM_FFLAGS = -I. -I../includes -m32
AM_CXXFLAGS = $(AM_CFLAGS)

AM_LDFLAGS  = -L.
AM_LDFLAGS += -L../../libs/dap
AM_LDFLAGS += -L../../libs/qreo
AM_LDFLAGS += -L../../libs/s2
AM_LDFLAGS += -L$(prefix)/lib

LDADD = -lreorder_qreo
LDADD += -lreorder_dap
LDADD += -lreorder_s2

LDADD += -lnetcdf
LDADD += -lhdf5_hl
LDADD += -lhdf5
LDADD += -lz
LDADD += -lbz2

# target libraries

bin_PROGRAMS = qreod

# sources

qreod_SOURCES = \
	qreoroot_dorade.F

