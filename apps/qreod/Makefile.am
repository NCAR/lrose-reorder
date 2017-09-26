###############################################
# makefile.am template for automake
###############################################

# compile flags

AM_CFLAGS = -I.
AM_CFLAGS += -I../includes
AM_CFLAGS += -I../../../include
AM_CFLAGS += -I$(prefix)/include
AM_CFLAGS += -m32

AM_FFLAGS = $(AM_CFLAGS)
AM_CXXFLAGS = $(AM_CFLAGS)

# linker flags

AM_LDFLAGS  = -L.
AM_LDFLAGS += -L../../libs/dap
AM_LDFLAGS += -L../../libs/qreo
AM_LDFLAGS += -L../../libs/s2
AM_LDFLAGS += -L../../../lib
AM_LDFLAGS += -L$(prefix)/lib
AM_LDFLAGS += -m32
AM_LDFLAGS += -Wl,-rpath=$(prefix)/lib

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

