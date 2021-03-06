README for REORDER VERSION 1.38-3.10 (LJM Jan 14, 2015)
=======================================================

Installing 32-bit system packages
=================================

Reorder is built as a 32-bit application, i.e. using the -m32 compile flag.

So you will need to install the 32-bit version of the supporting libraries.

For example on RedHat-based systems (including Centos, Scientific Linux)
run the following:

yum install -y ftp git svn cvs scons \
  gcc gcc-c++ gcc-gfortran \
  libgfortran.i686 libquadmath.i686 \
  glibc-devel.i686 libX11-devel.i686 libXext-devel.i686 \
  libtiff-devel.i686 libpng-devel.i686 \
  libstdc++-devel.i686 libtiff-devel.i686 \
  zlib-devel.i686 expat-devel.i686 flex-devel.i686 \
  fftw-devel.i686 bzip2-devel.i686

On Debian, you need to run the following:


  /usr/bin/dpkg --add-architecture i386
  apt-get update

and use apt-get to install the following:

  apt-get install libx11-6:i386 \
                  libstdc++-4.9-dev:i386 \
                  libquadmath0-dev:i386 \
                  libpng12-dev:i386 \
                  libx11-dev:i386 \
                  libxext-dev:i386 \
                  lib32stdc++-4.9-dev \
                  libstdc++5:i386 \
                  libstdc++6:i386

Quickstart install - in ~/reorder
=================================

  mkdir ~/reorder
  cd ~/reorder
  git clone https://github.com/ncar/lrose-netcdf
  git clone https://github.com/ncar/lrose-reorder
  cd lrose-netcdf/
  ./build_and_install_netcdf.m32 -x ~/reorder/
  cd ../lrose-reorder
  ./configure --prefix=${home}/reorder
  make install

Binaries will be in:

  ~/reorder/bin

Choose a location for your install
==================================

This is a 32-bit install, which is somewhat unusual these days.

Therefore it is a good idea to keep it separate from other installs.

Choose a separate location for this install.

Examples would be:

  ${home}/reorder
  /opt/local/reorder
  /usr/local/reorder

For /usr/local and /opt/local, you will need to be root.

For the remainder of this doc, we will assume the install location is:

  ${home}/reorder

Installing 32-bit versions of HDF5 and NetCDF libraries
=======================================================

32-bit versions for the HDF5 and NetCDF libraries are generally not
available in the form of packages.

Therefore, you will need to build them from the lrose distribution.

Check out the distribution:

  mkdir -p ~/git
  cd ~/git
  git clone https://github.com/ncar/lrose-netcdf

Build and install:

  cd lrose-netcdf
  ./build_and_install_netcdf.m32 -x ~/reorder

where the -x argument indicates where you want to perform your install.

Building with configure (recommended)
=====================================

Check out the distribution:

  mkdir -p ~/git
  cd ~/git
  git clone https://github.com/ncar/lrose-reorder

Build and install:

  cd lrose-reorder
  ./configure --prefix=${home}/reorder
  make -j 8
  make install

The binaries will be in:

  ${home}/reorder/bin

Building with scons
===================

If you build with scons, netcdf must be in a standard location.
This is tricky since this is a 32-bit build and will probably
conflict with the 64-bit install.

So probably it is better to use the configure build.

Note that the SCons construction tool is now required for building reorder.
Pre-built packages for SCons are available for most common operating systems;
check the common repositories for your OS and/or the download site at 
http://www.scons.org.

Unpack the reorder tar file and move to the resulting directory:

    $ tar xvzf Reorder_3.00.tar.gz
    $ cd Reorder/trunk

To build reorder without netCDF support:

    $ scons
    
Or, if you have netCDF libraries installed and want netCDF support in reorder:

    $ scons useNetcdf=True
    
The two reorder programs qreod and qreou will be created in the 'bin' 
directory below your current directory.

NOTE: The build method assumes that GNU compilers (gcc, gfortran, and/or g77)
will be used.  If you do not have GNU compilers, the build may still work,
but you need to explicitly disable use of a GNU-specific compiler flag by adding
'use_m32=False' to the appropriate scons command line above.


