#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Interpolate
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/grid.h"
#include "../include/field.h"
#include "../include/interpolate.h"
#include <math.h>

const int OutsideInterp = -2;

/*!
 * Tests Interpolate::F functions using paraboloid
 */
BOOST_AUTO_TEST_CASE (Paraboloid_interp) {
    int nr = 50;
    int nz = 50;
    double Rmin = -5;
    double Rmax = 5;
    double zmin = -3;
    double zmax = 3;

    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    
    // Initialize Psi as paraboloid
    Field *psi = new Field(*grid);
    for (int i=0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            psi->f_[i][j] = grid->R_[i]*grid->R_[i] + grid->z_[j]*grid->z_[j];
        }
    }
    
    Interpolate *inter = new Interpolate(*grid, *psi);
    inter->updateInterpolation(0,0);
    
//    double psi_interp;
//    try {
//        psi_interp = inter->bicubicInterpolate(0.05,0.05);
//    }
//    catch(int i) {
//        if (i == OutsideInterp) {
//            printf("Interpolation outside of current gridcell\n");
//        }
//    }
//    BOOST_CHECK_CLOSE(0.005, psi_interp, 10);
//    
//    double psi_interp_r;
//    try {
//        psi_interp_r = inter->bicubicInterpolate_r(0.05,0.05);
//    }
//    catch(int i) {
//        if (i == OutsideInterp) {
//            printf("Interpolation outside of current gridcell\n");
//        }
//    }
//    BOOST_CHECK_CLOSE(0.1, psi_interp_r, 10);
//    
//    double psi_interp_rr;
//    try {
//        psi_interp_rr = inter->bicubicInterpolate_rr(0.05,0.05);
//    }
//    catch(int i) {
//        if (i == OutsideInterp) {
//            printf("Interpolation outside of current gridcell\n");
//        }
//    }
//    BOOST_CHECK_CLOSE(2, psi_interp_rr, 10);
    
    double psi_interp;
    try {
        psi_interp = inter->F(0.05,0.05);
    }
    catch(int i) {
        if (i == OutsideInterp) {
            printf("Interpolation outside of current gridcell\n");
        }
    }
    BOOST_CHECK_CLOSE(0.005, psi_interp, .5);
    
    double psi2_interp;
    try {
        psi2_interp = inter->F(0.01,0.01);
    }
    catch(int i) {
        if (i == OutsideInterp) {
            printf("Interpolation outside of current gridcell\n");
        }
    }
    BOOST_CHECK_CLOSE(0.0002, psi2_interp, .5);
    
    double psir2_interp;
    try {
        psir2_interp = inter->F_r(0.01, 0.01);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(0.02, psir2_interp, .5);
    
    double psiz_interp;
    try {
        psiz_interp = inter->F_z(0.05, 0.05);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(0.1, psiz_interp, .5);
    
    double psiz2_interp;
    try {
        psiz2_interp = inter->F_z(0.01, 0.01);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(0.02, psiz2_interp, .5);


    double psizz_interp;
    try {
        psizz_interp = inter->F_zz(0.01, 0.01);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(2, psizz_interp, .005);

    double psirr_interp;
    try {
        psirr_interp = inter->F_rr(0.01, 0.01);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(2, psirr_interp, .005);


}


/*!
* Tests Interpolate::F functions using a cubic function
 */
BOOST_AUTO_TEST_CASE (cubic_interp) {
    int nr = 50;
    int nz = 50;
    double Rmin = -5;
    double Rmax = 5;
    double zmin = -3;
    double zmax = 3;


    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    
    // Initialize Psi as some polynomial function
    Field *psi = new Field(*grid);
    for (int i=0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            double R = grid->R_[i];
            double z = grid->z_[j];
            psi->f_[i][j] = R*R*R+2*z*z*z+ 3*R*z*z + 4*R*R*z*z;
        }
    }
    
    Interpolate *inter = new Interpolate(*grid, *psi);
    double R0 = 2.3;
    double z0 = 1.;

 
    inter->updateInterpolation(R0,z0);
    double psi_interp;
    try {
        psi_interp = inter->F(R0,z0);
    }
    catch(int i) {
        if (i == OutsideInterp) {
            printf("Interpolation outside of current gridcell\n");
        }
    }
    BOOST_CHECK_CLOSE(R0*R0*R0+2*z0*z0*z0+3*R0*z0*z0 + 4*R0*R0*z0*z0, psi_interp, .05);
    
    double psir2_interp;
    try {
        psir2_interp = inter->F_r(R0,z0);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(3*R0*R0 + 8*R0*z0*z0 + 3*z0*z0, psir2_interp, .05);
    
    double psiz_interp;
    try {
        psiz_interp = inter->F_z(R0,z0);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(6*z0*z0+ 8*R0*R0*z0 + 6*R0*z0, psiz_interp, .05);
    
    double psizz_interp;
    try {
        psizz_interp = inter->F_zz(R0,z0);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
    BOOST_CHECK_CLOSE(12*z0 + 8*R0*R0 + 6*R0, psizz_interp, .05);

    double psirz_interp;
    try {
        psirz_interp = inter->F_rz(R0,z0);
    }
    catch(int i) {
        if (i == OutsideInterp) printf("Interpolation outside of current gridcell\n");
    }
   BOOST_CHECK_CLOSE( 16*R0*z0 + 6*z0, psirz_interp, .05);


}
