#################################################################
# Autoconf template for reorder
# Process this file with autoconf to produce a configure script.
##################################################################

# init

AC_PREREQ([2.50])
AC_INIT([reorder], [1.40], [http://github.com/NCAR/lrose-reorder])

# Create the VERSION file, which contains the package version from
# AC_INIT.

echo -n AC_PACKAGE_VERSION>VERSION
AC_SUBST(PACKAGE_VERSION)
AC_MSG_NOTICE([netCDF AC_PACKAGE_VERSION])

# Keep libtool macros in an m4 directory.
#AC_CONFIG_MACRO_DIR([m4])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# We will output a config.h for the C library, and a fortran inc file,
# nfconfig.inc.
AC_CONFIG_HEADERS([config.h])

# initialize automake
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])

# Fix for automake 1.13
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Check for the existance of this file before proceeding.
AC_CONFIG_SRCDIR([libs/dap/qdap.F])

# user options
#MSG_NOTICE([checking user options])

# Did the user specify a location for the HDF5 library?
#AC_MSG_CHECKING([whether a location for the HDF5 library was specified])
#AC_ARG_WITH([hdf5],
#              [AS_HELP_STRING([--with-hdf5=<directory>],
#                              [Specify location of HDF5 library. Configure will \
#                              expect to find subdirs include and lib.])],
#            [HDF5DIR=$with_hdf5])
#AC_MSG_RESULT([$HDF5DIR])
#AC_SUBST(HDF5DIR, [$HDF5DIR])

# Did the user specify a location for the NETCDF library?
#AC_MSG_CHECKING([whether a location for the NETCDF library was specified])
#AC_ARG_WITH([netcdf],
#              [AS_HELP_STRING([--with-netcdf=<directory>],
#                              [Specify location of NETCDF library. Configure will \
#                              expect to find subdirs include and lib.])],
#            [NETCDFDIR=$with_netcdf])
#AC_MSG_RESULT([$NETCDFDIR])
#AC_SUBST(NETCDFDIR, [$NETCDFDIR])

# Add the HDF5 directory to includes and link flags
#if test ! "x$HDF5DIR" = x; then
#   LDFLAGS = "${LDFLAGS} -L${HDF5DIR}/lib"
#   CPPFLAGS = "${CPPFLAGS} -I${HDF5DIR}/include"
#fi    

# Add the NETCDF directory to includes and link flags
#if test ! "x$NETCDFDIR" = x; then
#   LDFLAGS = "${LDFLAGS} -L${NETCDFDIR}/lib"
#   CPPFLAGS = "${CPPFLAGS} -I${NETCDFDIR}/include"
#fi    

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77([gfortran g77])
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

# Checks for library functions.
AC_HEADER_STDC
AC_FUNC_STAT
AC_FUNC_UTIME_NULL
AC_CHECK_FUNCS([floor gettimeofday memset modf pow strerror])

# set makefile list

AC_CONFIG_FILES([
  Makefile
  libs/Makefile
  libs/dap/Makefile
  libs/qreo/Makefile
  libs/s2/Makefile
  apps/Makefile
  apps/qreod/Makefile
  apps/qreou/Makefile
])

AC_OUTPUT
