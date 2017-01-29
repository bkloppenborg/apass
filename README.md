APASS Data Reduction and Processing Software
=====

# Prerequisites

The APASS data reduction has the following dependencies:

* Fortran 77 compiler (gfortran is sufficent)
* SuperMongo 2.4.7 64-bit
* X11
* SLAlib (Fortran variant)


## SuperMongo

filtran requires [SuperMongo](http://www.astro.princeton.edu/~rhl/sm/) 2.4.27
(or possibly newer).

Note that recent versions of C will complain about the definition of `FREE` in
`src/bison/new.h`. As such, you can either replace the file:

    wget -O src/bison/new.h ftp://ftp.andrew.cmu.edu/pub/AUIS/src/overhead/bison/new.h

or if that link is broken, replace the definition of `FREE` in
`src/bison/new.h` with the following:

    #ifdef __STDC__
    #define	FREE(x)		(x ? (void) free((char *) (x)) : (void)0)
    #else
    #define FREE(x) 	((x) != 0 && (free ((char *) (x)), 0))
    #endif

After this, simply 

    cd sm_2_4_7
    set_opts    # take defaults, maybe move the install to ~/local
    vim makefile

Check that `CC = cc -Wall -Dlinux` is set. You probably don't need the swap
definition on modern systems.

    make
    make install

## SLAlib

filtran also requires the Fortran version of
[SLALib](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html). Compile it
however you wish. A simple CMake script to build it is as follows:

    cmake_minimum_required(VERSION 2.8)
    project(slalib Fortran)
    
    #set(CMAKE_Fortran_FLAGS -m32)
    
    file(GLOB SOURCE "*.f")
    
    add_library(sla_f ${SOURCE})

Which should be inserted into the same directory as the `.f` files.
After you have this, you can build as:

    mkdir build
    cd build
    cmake ..
    make

This will compile the `libsla_f.a` library and place it in the `build`
directory.

# Build instructions

    cd build
    cmake ..
    make
