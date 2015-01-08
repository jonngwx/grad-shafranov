#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Critical
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/grid.h"
#include "../include/field.h"
#include "../include/critical.h"
#include <math.h>
#include "../include/boundary.h"

/*!
 * Tests Critical::Psi_magnetic using paraboloid
 */
BOOST_AUTO_TEST_CASE (Paraboloid_min) {
    // set up Critical
    double z_limiter1 = 0.75;
    double z_limiter2 = -0.75;
    int max_iter_crit = 10;
    double error_tol_crit = 0.1;
    int nr = 50;
    int nz = 50;
    double Rmin = -5;
    double Rmax = 5;
    double zmin = -3;
    double zmax = 3;
    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    // Initialize magnetic axis
    double R0 = 0.1;
    double z0 = 0.1;
    
    // Initialize Psi as paraboloid
    Field *psi = new Field(*grid);
    for (int i=0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            psi->f_[i][j] = grid->R_[i]*grid->R_[i] + grid->z_[j]*grid->z_[j];
        }
    }
    
    Critical *crit = new Critical(*grid, *psi, max_iter_crit, error_tol_crit, z_limiter1, z_limiter2, R0, z0);
    double rcrit, zcrit, Psi_min;

    // Calculate Psi_o using previous coordinates
    crit->Psi_magnetic(R0, z0, &rcrit, &zcrit, &Psi_min);
    printf("rcrit = %f\n", rcrit);
    printf("zcrit = %f\n", zcrit);
}