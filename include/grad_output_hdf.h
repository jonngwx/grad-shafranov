#ifndef GRAD_OUTPUT_HDF
#define GRAD_OUTPUT_HDF


/**
 * \file Header file for grad_output_hdf
 * */

#include <assert.h> 
#include <stdlib.h>
#include "grad_output.h"
#include "field.h"
#include <hdf5.h>

class Grad_Output_Hdf : public Grad_Output{
public:
    /**
     * Constructor for output class
     * @param f pointer to field containing flux function
     * @param p pointer to field containing pressure function
     * @param g pointer to field containing g function
     * @param outputs string of comma separated output options
     * */
    Grad_Output_Hdf(Field* f, Grid* grid, RHSfunc* p, RHSfunc* g, const char* outputs);
    ~Grad_Output_Hdf();

    /**
     * Writes output to file.
     * @param filename name of output file
     * */
    void write_output(const char* filename);
 private:
    /**
     * @brief Converts an array of pointers to a 1d array
     * @param f array of pointers to convert
     * @param x 1d array to write to
     * @param nx dimensions of 2d array in x
     * @param ny dimensions of 2d array in y
     * */
    void twod_to_oned(const double * const *f, double *x, int nx, int ny);

    inline void check(herr_t status){
        assert(status >= 0);
    };
    
};

#endif
