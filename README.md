SWIPE' pitch estimator, v. 1.5
==============================

Based on Camacho, Arturo. A sawtooth waveform inspired pitch estimator for
speech and music. Doctoral dissertation, University of Florida. 2007.

Implemented in C by Kyle Gorman <kylebgorman@gmail.com>

How to cite:
------------

Please cite this dissertation, and if possible include a URL to the source.

How to install:
---------------

For all platforms: `mkdir build && cd build` To compile as a library, type `cmake ..` at the terminal. To compile the with a program binary type `cmake -DBUILD_BINARY=TRUE ..`

Linux: All the large libraries should be available as packages if you're using a "modern" distro. For instance, on a Ubuntu system (Ubuntu 9.04, "Jaunty Jackalope", kernel 2.6.28-13-generic), I ran:

    sudo apt-get install libblas-dev liblapack-dev libfftw3-dev libsndfile1-dev swig

This installs the necessary libraries and all their dependencies. Similar
incantations are available for other Linux distributions.

Mac OS X: The linear algebra libraries ([C]LAPACK, BLAS) ship with Mac OS X. [fftw3](http://www.fftw.org/).

If you are superuser and wish to install globally the autoconf method should work fine:

    tar -xvzf downloadedPackage.tar.gz
    cd folderOfPackageCreatedByUnTARring/
    ./configure; make; make install;

These two libraries are also available via Fink and DarwinPorts.

Windows/CYGWIN: Unsupported. Send details of any successes, however.

Miscellany:
-----------

This library has now been incorporated into the excellent [Speech Signal Processing Toolkit](http://sp-tk.sourceforge.net/)

Somewhat of a hackjob when ported to C++, but includes optimzations to some loops using eigen/strict bounds. Uses fitpack to fit and evaluate the spline.

Reference
Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993
Developer
Paul Dierckx, Department of Computer Science, K.U. Leuven, Celestijnenlaan 200 A, B-3001, Heverlee, Belgium
Paul.Dierckx@cs.kuleuven.ac.be