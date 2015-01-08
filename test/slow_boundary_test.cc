#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SlowBoundary
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/boundary.h"
#include "../include/slow_boundary.h" 
#include "../include/grid.h"
#include "../include/field.h"
#include <math.h>
#include "gauss_seidel.h"

/*!
 * Tests Slow Boundary in the case of no external
 * coils. This test uses the Shafranov-Solovev Solution on
 * Jardin p. 104 to calculate Jphi, and then calcs 
 * psi on the boundary using the first term of 
 * eqn 4.65 (which is missing a factor of mu0). 
 * These boundary values should match the
 * boundary of the original psi-field (eqn 4.39),
 * since it was already a solution to the GSE.  
 */

BOOST_AUTO_TEST_CASE (Shaf_Sol_Slow) {

    /***********************************
     * Specify parameters for the grid
     **********************************/
    double Rmin = 4.0;
    double Rmax = 8.0;
    double zmin = -2.0;
    double zmax = 2.0;
    int nr = 20;
    int nz = 20;
  
    double max_iter = 1000000;
    double epsilon = 1e-8;
    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    Field *psi = new Field(*grid);
    Field *psi_old = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *psib = new SlowBoundary(grid);
    /*********************************
     * Specify parameters for SS model
     *********************************/
    double B0 = 1.0; // B field in Teslas at mag axis
    double R0 = 6.0; // radial position of mag axis
    double kap0 = 1.1; // ellipticity of plasma
    double q0 = 2.01; // safety factor at mag axis
    double eps = 0.2; // inverse aspect ratio
    const double mu0 = 4*M_PI*1e-7; // mag perm of free space
    double psi_pv = eps*eps*kap0*R0*R0*B0/(2*q0); // defines psi at plas/vac interf.
    double R,Z;
    double coeff = B0/(R0*R0*kap0*q0);
    for (int i = 0; i < grid->nr_; ++i) {
      R = grid->R_[i];
      for (int j = 0; j < grid->nz_; ++j) {
	Z = grid->z_[j];
	psi_old->f_[i][j] = (coeff/2.0)*(R*R*Z*Z + (kap0*kap0/4.0)*pow((R*R-R0*R0),2));
	psi->f_[i][j] = psi_old->f_[i][j]+1;
	if ((psi_old->f_[i][j] >= 0) && (psi_old->f_[i][j] <= psi_pv)) {
	  jphi->f_[i][j] = R*coeff*(kap0*kap0 + 1)/mu0;
	}
	else {
	  jphi->f_[i][j] = 0;
	} 
      }
    }

    /*********************************
     * Calculate new boundary for psi
     *********************************/
    psib->CalcB(psi, jphi); // updates the boundary grid points of psi

    
    // Solve for psi given the proper boundary conditions
    EllipticSolver *solver = new GaussSeidel(*grid, *psi);
    solver->coeff();
    for (int n = 0; n < max_iter; ++n) {
        if (n==0) solver->step_1(*jphi);
        else solver->step(*jphi);
//	printf("n = %d\n, norm = %f\n",n,solver->norm());
        if (solver->norm() < epsilon) {
	  break;

	}
    }
//    printf("\n FINAL \n");
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
          if (i==0) {printf("psi = %f \n", psi->f_[i][j]);}
	  if ((psi_old->f_[i][j] >= 0) && (psi_old->f_[i][j] <= psi_pv)) {
            //printf("%f\t", psi->f_[i][j]);
	    BOOST_CHECK_CLOSE(psi_old->f_[i][j], psi->f_[i][j], 1);
	    
	  }
	}
        //printf("\n");
    }
}


