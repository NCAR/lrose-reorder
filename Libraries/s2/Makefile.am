###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../publicIncludes -m32
AM_FFLAGS = -I. -I../publicIncludes -m32
AM_CXXFLAGS = $(AM_CFLAGS)
AM_LDFLAGS = -L.

# target libraries

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
