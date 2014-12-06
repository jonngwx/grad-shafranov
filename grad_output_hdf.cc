#include "grad_output_hdf.h"
#include "field.h"
#include "grid.h"
#include <stdio.h>
#include "hdf5.h"
#include "hdf5_hl.h"

Grad_Output_Hdf::Grad_Output_Hdf(Field* f0, Grid* grid0, RHSfunc * p0, RHSfunc * g0, const char* outputs) {
    f = f0;
    p = p0;
    g = g0;
    grid = grid0;
    Grad_Output::parse_outputs(outputs);
}

Grad_Output_Hdf::~Grad_Output_Hdf(){
}

void Grad_Output_Hdf::write_output(const char* filename){
    hid_t       file_id;
    hsize_t     dims[2],dimr[1],dimz[1];
    herr_t      status;
    const int nr = grid->nr_;
    const int nz = grid->nz_;
    
    /* Create a new file using default properties. */
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Create the data space for the dataset. */
    dims[0] = nr; 
    dims[1] = nz; 
    dimr[0] = nr;
    //    dimr[1] = 1;
    dimz[0] = nz;
    //dimz[1] = 1;
    /* Create and write R. */
    status = H5LTmake_dataset(file_id,"/R",1,dimr,H5T_NATIVE_DOUBLE, grid->R_);
    
    
    /* Close the file. */
    status = H5Fclose(file_id);

}

