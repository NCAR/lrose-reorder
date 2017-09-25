#################################################################
# Autoconf template for reorder
# Process this file with autoconf to produce a configure script.
##################################################################

# init

AC_PREREQ([2.50])
AC_INIT([reorder], [1.38], [http://github.com/NCAR/lrose-reorder])

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
AC_CONFIG_SRCDIR([Libraries/dap/qdap.F])

# user options
#MSG_NOTICE([checking user options])

# Did the user specify a location for the HDF5 library?
AC_MSG_CHECKING([whether a location for the HDF5 library was specified])
AC_ARG_WITH([hdf5],
              [AS_HELP_STRING([--with-hdf5=<directory>],
                              [Specify location of HDF5 library. Configure will \
                              expect to find subdirs include and lib.])],
            [HDF5DIR=$with_hdf5])
AC_MSG_RESULT([$HDF5DIR])
AC_SUBST(HDF5DIR, [$HDF5DIR])

# Did the user specify a location for the NETCDF library?
AC_MSG_CHECKING([whether a location for the NETCDF library was specified])
AC_ARG_WITH([netcdf],
              [AS_HELP_STRING([--with-netcdf=<directory>],
                              [Specify location of NETCDF library. Configure will \
                              expect to find subdirs include and lib.])],
            [NETCDFDIR=$with_netcdf])
AC_MSG_RESULT([$NETCDFDIR])
AC_SUBST(NETCDFDIR, [$NETCDFDIR])

# 32-but build

CFLAGS = "-m32 -DNETCDF"
FFLAGS = "-m32 -DNETCDF"

# Add the HDF5 directory to includes and link flags
if test ! "x$HDF5DIR" = x; then
   LDFLAGS = "${LDFLAGS} -L${HDF5DIR}/lib"
   CPPFLAGS = "${CPPFLAGS} -I${HDF5DIR}/include"
fi    

# Add the NETCDF directory to includes and link flags
if test ! "x$NETCDFDIR" = x; then
   CFLAGS = " -DNETCDF"
   FFLAGS = " -DNETCDF"
   LDFLAGS = "${LDFLAGS} -L${NETCDFDIR}/lib"
   CPPFLAGS = "${CPPFLAGS} -I${NETCDFDIR}/include"
fi    

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
  Libraries/Makefile
  Libraries/dap/Makefile
  Libraries/qreo/Makefile
  Libraries/s2/Makefile
  Applications/Makefile
  Applications/qreod/Makefile
  Applications/qreou/Makefile
])

AC_OUTPUT



#AC_MSG_RESULT( NETCDFINCS: $NETCDFINCS )

#if test -z "$NETCDFLIBS" ; then
#if test -n "$NETCDFLIBDIR"; then
#   NETCDFLIBS="-L$NETCDFLIBDIR -lnetcdf -lnetcdff -ludunits2 -lhdf5 -lhdf5_hl -lz"
#else
#   NETCDFLIBS="-lnetcdf -lnetcdff -ludunits2 -lhdf5 -lhdf5_hl -lz"
#fi
#fi
#AC_MSG_RESULT( NETCDFLIBS: $NETCDFLIBS )

# We don't use any C checks at the moment, but without it autoconf warns
# about the existence of C source files without a C compiler check.  So
# setup the C compiler but without any of the usual checks.  (This were
# included by the initial autoscan.)
#AC_PROG_CC
# Checks for header files.
#AC_LANG([C])
#AC_CHECK_HEADERS([fcntl.h malloc.h memory.h netinet/in.h stdlib.h string.h strings.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_TYPE_OFF_T
#AC_TYPE_PID_T
#AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_MKTIME
#AC_FUNC_REALLOC
#AC_FUNC_WAIT3
#AC_CHECK_FUNCS([atexit bzero floor ftruncate gethostname gettimeofday memmove memset modf pow select strtoul])

#AM_CONDITIONAL([MAKE_WITH_NETCDF], [test "$with_netcdf" != no])

#AC_SUBST(NETCDFINCDIR)
#AC_SUBST(NETCDFLIBDIR)
#AC_SUBST(NETCDFLIBS)
#AC_SUBST(NETCDFINCS)

