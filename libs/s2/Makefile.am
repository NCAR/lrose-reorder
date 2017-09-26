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

lib_LIBRARIES = libreorder_s2.a

# sources

libreorder_s2_a_SOURCES = \
	by_products.c \
	cimm_dd.c \
	dd_catalog.c \
	dd_crackers.c \
	dd_der_flds.c \
	dd_files.c \
	dd_io_mgmt.c \
	dd_ncdf.c \
	dd_swp_files.c \
	dd_time.c \
	dd_uf.c \
	dda_common.c \
	ddb_common.c \
	ddin.c \
	ddout.c \
	dorade_ncdf.c \
	dorade_share.c \
	dorade_tape.c \
	dorade_uf.c \
	elda_dd.c \
	etl_dd.c \
	fof_dd.c \
	gecho.c \
	gneric_dd.c \
	gpro_data.c \
	ground_echo.c \
	hrd_dd.c \
	nc_dd.c \
	nssl_mrd.c \
	product_x.c \
	stdhrd.c \
	swp_file_acc.c \
	tdwr_dd.c \
	uf_dd.c \
	xwsrqc.c
