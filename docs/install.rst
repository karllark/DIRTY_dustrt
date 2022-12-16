.. _install:

#############
Install DIRTY
#############

Code
====

Start with downloading the latest release of the code from `DIRTY releases
<https://github.com/karllark/DIRTY_dustrt/releases>`_.
The latest development version can be obtained from the `github repository
<https://github.com/karllark/DIRTY_dustrt>`_
instead, but note that your mileage may vary with this version.

DIRTY uses the cfitsio library for FITS I/O.  This may already be
installed on your system, but maybe not.  It can be installed locally,
but make sure to set up `pkg-config` or run `./configure --with-cfitsio[=PATH]`
appropriately so that the cfitsio library can be found during compilation.

Run make in the DIRTY directory to generate the `dirty` executable.

Dependencies
============

Conda
-----

.. code::

    conda create -n dirty -c conda-forge autoconf automake cfitsio compilers pkgconfig sphinx
    conda activate dirty

Debian/Ubuntu
-------------

.. code::

    apt install build-essential libcfitsio-dev pkg-config python-sphinx

RHEL/CentOS 7
-------------

.. code::

    yum install epel-release
    yum install autoconf automake cfitsio-devel gcc-g++ pkg-config python3-sphinx

RHEL/CentOS 8+/Fedora
---------------------

Appstream and EPEL lack the `python3-sphinx` package (broken?). Use a virtualenv to supplement the missing compile-time dependency.

.. code::

    dnf install epel-release
    dnf install autoconf automake cfitsio-devel gcc-g++ pkg-config python3
    python3 -mvenv venv
    source venv/bin/activate
    pip install sphinx

Installation (Release)
======================

.. code::

    git clone https://github.com/karllark/DIRTY_dustrt
    cd DIRTY_dustrt
    git tag -l -n
    git checkout [tag]
    ./bootstrap
    ./configure --prefix=$HOME/programs/dirty
    make
    make install
    export PATH=$HOME/programs/dirty/bin:$PATH

Installation (Development)
==========================

.. code::

    git clone https://github.com/karllark/DIRTY_dustrt
    cd DIRTY_dustrt
    ./bootstrap
    ./configure --enable-debug
    make


Dust Grain Info
===============

The associated dust grain files can be found at XXXXX.


