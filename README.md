COUGAR README
=============

Nonstandard Requirements:
-------------------------
* c++11 compiler like a recent version of g++
* boost_program_options
* python 2.7 (with library numpy, matplotlib and h5py)

to compile, type `make`

to run tests, type `make runtests` 

to run type

    ./gs_solver

Try out the command line options!

    --version
and

    --help
    
The options listed in help can be specified on the command line or in grad-shafranov.cfg (the default config file).

Not all the options are implemented yet.

Of course, if you specify an illegal value for an option, or a file which does not exist or is not the right format, the program may terminate with an Error message, or with a crash.

To plot the output, type

    python python/viz.py cougar.out.tsv tsv
    
if the default output filename is used.  More generally, type

    python python/viz.py <filename> <tsv|hdf>
    
To access the data in python, import viz, then use

    F = viz.plot(filename, format)

F will be a dictionary of the grid and data arrays.
