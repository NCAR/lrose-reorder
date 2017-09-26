###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../includes -DNETCDF -m32
AM_FFLAGS = -m32
AM_CXXFLAGS = $(AM_CFLAGS)
AM_LDFLAGS = -L.

# target libraries

lib_LIBRARIES = libreorder_dap.a

# pkginclude_HEADERS = \
# 	../includes/CellSpacing.h \
# 	../includes/CellSpacingFP.h \
# 	../includes/CellVector.h \
# 	../includes/Comment.h \
# 	../includes/Correction.h \
# 	../includes/FieldLidar.h \
# 	../includes/FieldParam.h \
# 	../includes/FieldRadar.h \
# 	../includes/LidarDesc.h \
# 	../includes/LidarParam.h \
# 	../includes/Parameter.h \
# 	../includes/Pdata.h \
# 	../includes/Platform.h \
# 	../includes/Qdata.h \
# 	../includes/RadarDesc.h \
# 	../includes/Ray.h \
# 	../includes/SConscript \
# 	../includes/Sweep.h \
# 	../includes/TimeSeries.h \
# 	../includes/Volume.h \
# 	../includes/Waveform.h \
# 	../includes/Xtra_stuff.h \
# 	../includes/dd_defines.h \
# 	../includes/dd_files.h \
# 	../includes/dd_general_info.h \
# 	../includes/dd_math.h \
# 	../includes/dd_time.h \
# 	../includes/dd_types.h \
# 	../includes/dgi_func_decl.h \
# 	../includes/dorade_data_descriptors.h \
# 	../includes/dorade_headers.h \
# 	../includes/dorade_share.h \
# 	../includes/function_decl.h \
# 	../includes/input_sweepfiles.h \
# 	../includes/piraq.h \
# 	../includes/piraqx.h \
# 	../includes/point_in_space.h \
# 	../includes/solo_utils.h \
# 	../includes/super_SWIB.h

# sources

libreorder_dap_a_SOURCES = \
	dap_common.c \
	qdap.F

