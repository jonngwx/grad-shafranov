#include "grad_output_hdf.h"
#include "field.h"
#include "grid.h"
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>

Grad_Output_Hdf::Grad_Output_Hdf(Field* f0, Field* jphi0, Grid* grid0, Field* p0, Field* g0, const char* outputs) {
    f = f0;
    jphi = jphi0;
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
    check(status);
    status = H5LTmake_dataset(file_id,"/z",1,dimz,H5T_NATIVE_DOUBLE, grid->z_);
    check(status);
    
    double *x = new double[nr*nz];
    // hackery due to need for contiguous memory
    twod_to_oned(f->f_, x, nr, nz);
    status = H5LTmake_dataset(file_id,"/psi",2,dims,H5T_NATIVE_DOUBLE,x);
    check(status);

    hsize_t dimpsil = 1;
    status = H5LTset_attribute_double(file_id, "/psi", "psi_l", &(f->f_l), dimpsil);
    check(status);
    status = H5LTset_attribute_double(file_id, "/psi", "psi_0", &(f->f_0), dimpsil);
    check(status);
    for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nz; ++j){
            x[i*nz + j] = p->f_[i][j];
        }
    }
    status = H5LTmake_dataset(file_id,"/p",2,dims,H5T_NATIVE_DOUBLE,x);
    for (int i = 0; i < nr; ++i){
        for (int j = 0; j < nz; ++j){
            x[i*nz + j] = g->f_[i][j];
        }
    }
    status = H5LTmake_dataset(file_id,"/g",2,dims,H5T_NATIVE_DOUBLE,x);

    /* Do custom output */
    for (auto i : output_list){
        switch(i){
        case CURRENT:
            twod_to_oned(jphi->f_, x, nr, nz);
            status = H5LTmake_dataset(file_id,"/j",2,dims,H5T_NATIVE_DOUBLE,x);
            break;
        case TOROIDAL_FIELD:
            for (int i = 0; i < nr; i++){
                for (int j = 0; j < nz; j++){
                    x[i*nz + j] = g->f_[i][j]/grid->R_[i];
                }
            }
            status = H5LTmake_dataset(file_id,"/bt",2,dims,H5T_NATIVE_DOUBLE,x);
            break;
        }
    }

    /* Close the file. */
    status = H5Fclose(file_id);
    delete [] x;

}

void Grad_Output_Hdf::twod_to_oned(const double * const *f, double *x, int nx, int ny){
    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            x[i*ny + j] = f[i][j];
        }
    }
}
