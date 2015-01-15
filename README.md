COUGAR README
=============


Installation:
-------------------------
To install COUGAR, download the source and ensure the following prerequisites are met

* c++11 compiler like a recent version of g++
* HDF5 (1.8.0) or newer
* boost_program_options
* boost_math_tr1
* boost_unit_test_framework
* python 2.7 (with packages numpy, matplotlib and h5py)

If HDF5 is unavailable, the code can be compiled without it. 

To compile without HDF5, type `make`. To compile with HDF5, type 'make gs_solver_hdf'. 

Testing is performed by typing `make runtests` 

Operation
-----------

    ./gs_solver

Try out the command line options!

    --version
and

    --help
    
The options listed in help can be specified on the command line or in grad-shafranov.cfg (the default config file).

Not all the options are implemented yet.

Of course, if you specify an illegal value for an option, or a file which does not exist or is not the right format, the program may terminate with an Error message, or with a crash.

To plot the output, type

    python python/viz.py cougar.out.tsv 
    
if the default output filename is used.  More generally, type

    python python/viz.py <filename> 
    
To access the data in python, import viz, then use

    F = viz.plot(filename, format)

F will be a dictionary of the grid and data arrays.
