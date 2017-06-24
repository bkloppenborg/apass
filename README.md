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

After this, configure SuperMongo

    cd sm_2_4_7
    set_opts    # take defaults, maybe move the install to ~/local
    
Edit the `src/options.h` file to remove the `#error` lines at the top of
the file and add the `FORTRAN_APPEND _` line. After this, verify that
the configuration in the makefile make sense:
    
    vim Makefile

Specifically, check that `CC = cc -Wall -Dlinux` is set. You probably 
don't need the swap definition on modern systems.

    make
    make install
    
After the installation step, check the external Fortran function naming
convention used by your compiler. Fortran functions will either be prefixed by
`f`, suffixed by one or two underscores. To determine the naming convention, run
something similar to the following in your shell:

    $ cd path_to_plotsub_install_dir
    
    $ nm -a libplotsub_d.a | grep sm_device
    0000000000000eed T sm_device
    000000000000071d T sm_device__
    
The first command locates your library. It should be in `/usr/local` or whatever
directory you specified in the `set_opts` step above. In this case the compiler
used the double-underscore convention. By default the build scripts here assume
the double-underscore convention. If your compiler is different, modify the  
`fortran/CmakeLists.txt` accordingly.

## SLAlib

filtran also requires the Fortran version of
[SLALib](http://star-www.rl.ac.uk/docs/sun67.htx/sun67.html). Compile it
however you wish. A simple CMake script to build it is as follows:

    cmake_minimum_required(VERSION 2.8)
    project(slalib Fortran)
    
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

# Build and installation instructions

To build and install the software to `/usr/local` (CMake's default), simply:

    cd build
    cmake ..
    make
    sudo make install

To change the output directory of the installation, specify
`-DCMAKE_INSTALL_PREFIX=/path/to/installation/directory` at cmake time.
If you simply want to use the built binaries, look in subdirectories of `build`
after running `make` (e.g. `build/fortran/filtran`)
