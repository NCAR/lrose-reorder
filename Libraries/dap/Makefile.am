###############################################
# makefile.am template for automake
###############################################

# compile flags - include header subdirectory

AM_CFLAGS = -I. -I../publicIncludes -m32
AM_FFLAGS = -m32
AM_CXXFLAGS = $(AM_CFLAGS)
AM_LDFLAGS = -L.

# target libraries

lib_LIBRARIES = libreorder_dap.a

pkginclude_HEADERS = \
	../publicIncludes/CellSpacing.h \
	../publicIncludes/CellSpacingFP.h \
	../publicIncludes/CellVector.h \
	../publicIncludes/Comment.h \
	../publicIncludes/Correction.h \
	../publicIncludes/FieldLidar.h \
	../publicIncludes/FieldParam.h \
	../publicIncludes/FieldRadar.h \
	../publicIncludes/LidarDesc.h \
	../publicIncludes/LidarParam.h \
	../publicIncludes/Parameter.h \
	../publicIncludes/Pdata.h \
	../publicIncludes/Platform.h \
	../publicIncludes/Qdata.h \
	../publicIncludes/RadarDesc.h \
	../publicIncludes/Ray.h \
	../publicIncludes/SConscript \
	../publicIncludes/Sweep.h \
	../publicIncludes/TimeSeries.h \
	../publicIncludes/Volume.h \
	../publicIncludes/Waveform.h \
	../publicIncludes/Xtra_stuff.h \
	../publicIncludes/dd_defines.h \
	../publicIncludes/dd_files.h \
	../publicIncludes/dd_general_info.h \
	../publicIncludes/dd_math.h \
	../publicIncludes/dd_time.h \
	../publicIncludes/dd_types.h \
	../publicIncludes/dgi_func_decl.h \
	../publicIncludes/dorade_data_descriptors.h \
	../publicIncludes/dorade_headers.h \
	../publicIncludes/dorade_share.h \
	../publicIncludes/function_decl.h \
	../publicIncludes/input_sweepfiles.h \
	../publicIncludes/piraq.h \
	../publicIncludes/piraqx.h \
	../publicIncludes/point_in_space.h \
	../publicIncludes/solo_utils.h \
	../publicIncludes/super_SWIB.h

# sources

libreorder_dap_a_SOURCES = \
	dap_common.c \
	qdap.F

