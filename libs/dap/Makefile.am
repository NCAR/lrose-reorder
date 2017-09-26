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

# target library

lib_LIBRARIES = libreorder_dap.a

# sources

libreorder_dap_a_SOURCES = \
	dap_common.c \
	qdap.F

