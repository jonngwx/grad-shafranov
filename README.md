COUGAR User manual
==================


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

In order to run COUGAR, a number of physical parameters must be provided. These are listed here and can also be found by calling the solver with the `--help` command line option. The parameters can be provided using the command line or a configuration file using the `-c` flag (the default file is `grad-shafranov.cfg`). The exact format of the command line options and configuration file is provided in the documentation.

* Main program and solver configuration: The positions and values of the external currents must be provided in a `.tsv` file, as must be the positions of the physical limiters. Currents must be outside the computational boundary, while limiters must be at least two cells from the edge of the boundary.
* Grid and geometry: The computational domain should be specified.
* Saddle point initial geometry: The positions of initial guesses for two saddle points of the flux function should be specified. These guesses should be within the computational grid. The saddle points are used to calculate the plasma boundary.
* The form and magnitude of the pressure and current: The magnetic field on-axis, pressure on-axis and the exponents for the forms of the pressure and toroidal field distributions should be given. 
* The initial current distribution: The total current, radius of the initial distribution and location of the centre of current density should be provided. 
* Output format and frequency: The format of the output, non-standard fields to write and the file name of the output file should be provided. There are also options to write output during the iterative process. 

To run the code, type 

    ./gs_solver

with the appropriate command line options. If you specify an illegal value for an option, or a file which does not exist or is not the right format, the program may terminate with an Error message, or with a crash.

To plot the output, type

    python python/viz.py cougar.out.tsv 
    
if the default output filename is used.  More generally, type

    python python/viz.py <filename> 
    
To access the data in python, import viz, then use

    F = viz.plot(filename, format)

F will be a dictionary of the grid and data arrays.
