# LROSE Reorder

Reorder is a package for gridding ray-based radar data.

Reorder was developed by Jay Miller at NCAR/MMM and NCAR/EOL.

## Version

The current version is 1.38-3.10 (LJM Jan 14, 2015)

See the ChangeLog for release notes and details.

## Introduction

As part of this project, we have started with Reorder that has been used for a number of years by several members of the scientific community.  All aspects of Reorder are currently being scrutinized in order to more fully understand how exactly the current code is performing.  This effort is expected to lead to a decision to either modernize and overhaul Reorder or to abandon it entirely.  If the latter choice is made, then the necessary requirements for a replacement gridder will be identified with input from the user community and incorporated into a newly designed package, keeping in mind that its lifetime is expected to be 5-10 years.

## Gridding Scheme

Reorder uses either a distance-weighting scheme (Cressman or exponential), equal (uniform) weighting, or the closest point to populate the output Cartesian grid.  The input samples are looped over in range (radial direction) one beam at a time to populate the output Cartesian grid with SUMS of (distance-weighting*field) and SUMS of (distance-weighting).  Distances between radar sample locations (RAE) and output grid points (XYZ) are calculated for all output grid points within a box centered on the radar sample location and having dimensions specified by the user.  Such an approach results in estimates at output grid points that are filtered versions of the input samples.  Strictly, then, this is not an interpolator and therefore has certain characteristics and limitations.

There are three important aspects to the Reorder gridding algorithm.  The first two are the dimensions and orientation of the box (region of influence), and the third is the distance-weighting formulation itself.  The user specifies three attributes, along with the values to be used (parameters).  The attributes are either XRADIUS, YRADIUS, and ZRADIUS or RGRADIUS, AZRADIUS, and ELRADIUS.

1. Dimensions of the box:  The X, Y, and Z radii are half-distances in the X, Y, and Z directions centered on the current range gate location in the Cartesian grid.  The R, A, and E radii were apparently intended to be half-distances in the R, A, and E directions with the idea of matching the resolution of the angular sampling.  In principle the numbers of samples used to populate the output grid would then be constant with range, leading to a constant reduction in error and preserving angular (not linear) resolution.  However, a careful look at the actual coded formulation has revealed this is not quite true.  This RAE formulation included an attempt to rotate the box such that it was aligned with the local range-azimuth-elevation directions.  The result of this failed attempt was to produce a pattern in the number of input measurements (count attribute) used that varied in the azimuthal direction.  This  pattern when displayed in the Cartesian output domain resembled a "clover-leaf" with minima near 0, 90, 180, and 270 degrees in azimuth and was not constant in the range direction.
2. Orientation of the box:  The box is implicitly oriented relative to the output grid axes because of the way that affected indices in the x, y, and z coordinate directions are calculated.  A range of x, y, and z indices from a minimum to a maximum value are calculated according to the integer in a minimum and maximum in output grid coordinate equal to range gate coordinate +/- radius value.
3. Distance-weighting formulation:

## Coding deficiencies

Some of the deficiencies in Reorder that have been identified so far are:

* This program has been around a number of years, certainly as far back as the early 1990s.  This has resulted in many "add-ons" that have not been well documented either internally in the source code or externally in any detail.  This makes it rather difficult after the fact for any person unfamiliar with the code to determine the real purpose of some functions and routines.  Such is the inevitable consequence of many legacy software packages that are not continuously maintained and updated, especially in keeping with more modern programming languages, standards, and compilers.
* Since the current source code is a mix of Fortran 77 and C, it requires that any single program maintainer be skilled in both programming languages.
* The original code is extremely modular, perhaps too much so, with many small modules (Fortran subroutines and functions), making it hard to trace the various paths through the code.
* This modularity further complicates any understanding of the code since oftentimes the variable names used in the parameter list in the call to a subroutine or function and those in the subroutine or function are different.  This programming style makes it very difficult to follow any particular variable through the code when either debugging existing code or adding any new functionality.
* The parts of code most relevant to the gridding algorithm from user-specified attributes and any associated parameters to the gridding (filtering) to writing the output file of gridded values represents only a small fraction, perhaps as little as 10%, of the total lines of code.  Much of the remaining code is used for parsing character strings and manipulating bytes using outdated Fortran methods.
* Selection of the radii of influence, either as XYZ or RAE radii, is problematic.  The linear density of input range-azimuth-elevation locations is not uniform in the output Cartesian space, making it impossible to have a single specification of radii of influence that will work at all radar slant ranges.

## Modifications

Modifications that have been made are:

* Commented out a call to a Fortran subroutine that was intended to reorient the RAE box to be along the radial, azimuth, and elevation directions.  The algorithm inside this routine led to a four-leafed pattern in the numbers of samples being used.  This pattern had minima near 0, 90, 180, and 270 deg azimuths and resembled a four-leaf clover.  The numbers of samples used should be nearly constant with azimuth at a fixed range.  Also, the numbers of samples used should be nearly constant in the range direction since the RAE radii were intended to match the original sampling resolutions.  These constancies in the expected numbers of samples is not at all present when the "count" field that can be generated during the gridding process is viewed.
* Added a new hybrid scheme for the radii of influence (box dimensions) centered on and surrounding the range gate sampling location.  The user can now specify both Cartesian (XYZ) and Spherical (RAE) radii that will be used in this hybrid scheme.  The idea behind an XYZ set of radii is to maintain constant box dimensions as range increases.  Large numbers of samples will be used near the radar because of the much smaller linear dimensions associated with typical angular sampling resolutions.  The use of an RAE set of radii on the other hand more nearly matches the original angular resolution, and the linear size of the box increases with range.  This approach leads to roughly equal numbers of samples as range increases to be used in populating the output grid.  The hybrid scheme uses the maxima of the XYZ and linear RAE radii so that near the radar, the dimensions of the box are more like XYZ radii and farther out in range they are more like RAE radii.
* Fixed a small error in the signage of longitude to follow the normal convention of having East (West) longitude as positive (negative).
* Added SAVE xmin, ymin, zmin inside qreo2.F:subroutine filter (my version 1.38) after Chris Burghart found and fixed a problem (version 1.39) with Reorder compiled with gfortran.  Michael Bell had found that, with a gfortran compilation, the grid was empty when compared with the g77 compiled version run on his hurricane Katrina test case using KLIX WSR-88D data.

## Sample script

Click on run-reorder.txt to see a sample script for gridding Dorade sweep files from S-POL during the STEPS-2000 program near Goodland KS.

## Installing 32-bit system packages

Reorder is built as a 32-bit application, i.e. using the -m32 compile flag.

So you will need to install the 32-bit version of the supporting libraries.

For example on RedHat-based systems (including Centos, Scientific Linux)
run the following:

```
yum install -y ftp git svn cvs scons \
  gcc gcc-c++ gcc-gfortran \
  libgfortran.i686 libquadmath.i686 \
  glibc-devel.i686 libX11-devel.i686 libXext-devel.i686 \
  libtiff-devel.i686 libpng-devel.i686 \
  libstdc++-devel.i686 libtiff-devel.i686 \
  zlib-devel.i686 expat-devel.i686 flex-devel.i686 \
  fftw-devel.i686 bzip2-devel.i686
```

On Debian, you need to run the following:

```
  /usr/bin/dpkg --add-architecture i386
  apt-get update
  apt-get install libx11-6:i386 \
                  libstdc++-4.9-dev:i386 \
                  libquadmath0-dev:i386 \
                  libpng12-dev:i386 \
                  libx11-dev:i386 \
                  libxext-dev:i386 \
                  lib32stdc++-4.9-dev \
                  libstdc++5:i386 \
                  libstdc++6:i386
```

## Quickstart install - in ~/reorder

```
  mkdir ~/reorder
  cd ~/reorder
  git clone https://github.com/ncar/lrose-netcdf
  git clone https://github.com/ncar/lrose-reorder
  cd lrose-netcdf/
  ./build_and_install_netcdf.m32 -x ~/reorder/
  cd ../lrose-reorder
  ./configure --prefix=${home}/reorder
  make install
```

Binaries will be in:

```
  ~/reorder/bin
```

## Choose a location for your install

This is a 32-bit install, which is somewhat unusual these days.

Therefore it is a good idea to keep it separate from other installs.

Choose a separate location for this install.

Examples would be:

```
  ${home}/reorder
  /opt/local/reorder
  /usr/local/reorder
```

For ```/usr/local``` and ```/opt/local```, you will need to be root.

For the remainder of this doc, we will assume the install location is:

```
  ${home}/reorder
```

# Installing 32-bit versions of HDF5 and NetCDF libraries

32-bit versions for the HDF5 and NetCDF libraries are generally not
available in the form of packages.

Therefore, you will need to build them from the lrose distribution.

Check out the distribution:

```
  mkdir -p ~/git
  cd ~/git
  git clone https://github.com/ncar/lrose-netcdf
```

Build and install:

```
  cd lrose-netcdf
  ./build_and_install_netcdf.m32 -x ~/reorder
```

where the -x argument indicates where you want to perform your install.

## Building with configure (recommended)

Check out the distribution:

```
  mkdir -p ~/git
  cd ~/git
  git clone https://github.com/ncar/lrose-reorder
```

Build and install:

```
  cd lrose-reorder
  ./configure --prefix=${home}/reorder
  make -j 8
  make install
```

The binaries will be in:

```
  ${home}/reorder/bin
```

# Building with scons

If you build with scons, netcdf must be in a standard location.
This is tricky since this is a 32-bit build and will probably
conflict with the 64-bit install.

So probably it is better to use the configure build.

Note that the SCons construction tool is now required for building reorder.
Pre-built packages for SCons are available for most common operating systems;
check the common repositories for your OS and/or the download site at 
http://www.scons.org.

Unpack the reorder tar file and move to the resulting directory:

```
    $ tar xvzf Reorder_3.00.tar.gz
    $ cd Reorder/trunk
```

To build reorder without netCDF support:

```
    $ scons
```

Or, if you have netCDF libraries installed and want netCDF support in reorder:

```
    $ scons useNetcdf=True
```

The two reorder programs qreod and qreou will be created in the 'bin' 
directory below your current directory.

NOTE: The build method assumes that GNU compilers (gcc, gfortran, and/or g77)
will be used.  If you do not have GNU compilers, the build may still work,
but you need to explicitly disable use of a GNU-specific compiler flag by adding
'use_m32=False' to the appropriate scons command line above.



