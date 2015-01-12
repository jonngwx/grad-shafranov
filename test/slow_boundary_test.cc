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
#include "../include/green_fcn.h"

/*!
 * This roundabout test involves the consistency
 * of green_fcn.cc, slow_boundary.cc, and gauss_seidel.cc.
 *  
 */

BOOST_AUTO_TEST_CASE (Shaf_Sol_Very_Slow) {

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
  //  double perim = 2.0*(nr + nz -2.0);
    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    Field *psi = new Field(*grid);
    Field *psi_old = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *psib = new SlowBoundary(psi, grid);
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
    psib->CalcB(psi); // updates the boundary grid points of psi

    
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

    // Compare selected points between psi and 
    // psi_green (psi calculated using green
    // functions competely.
    


}

BOOST_AUTO_TEST_CASE(Green_Fcn_Test){

    /***********************************
     * Specify parameters for the grid
     **********************************/
    double Rmin = 5.0;
    double Rmax = 7.0;
    double zmin = -1.0;
    double zmax = 1.0;
    int nr = 250;
    int nz = 250;
    int perim = 2*(nr + nz - 2);  

    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    double dr = grid->dr_;
    double dz = grid->dz_;
    Field *psi_old = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *b = new Boundary(psi_old, grid);
   
    /*********************************
     * Specify parameters for SS model
     *********************************/
    double B0 = 1.0; // B field in Teslas at mag axis
    double R0 = 6.0; // radial position of mag axis
    double kap0 = 1.1; // ellipticity of plasma
    double q0 = 2.01; // safety factor at mag axis
    double eps = 0.3; // inverse aspect ratio
    const double mu0 = 4*M_PI*1e-7; // mag perm of free space
    double psi_pv = eps*eps*kap0*R0*R0*B0/(2*q0); // defines psi at plas/vac interf.
    printf("psi_pv = %f \n", psi_pv);
    double R,z;
    double coeff = B0/(R0*R0*kap0*q0);
    for (int i = 0; i < grid->nr_; ++i) {
      R = grid->R_[i];
      for (int j = 0; j < grid->nz_; ++j) {
	z = grid->z_[j];
	psi_old->f_[i][j] = (coeff/2.0)*(R*R*z*z + (kap0*kap0/4.0)*pow((R*R-R0*R0),2));
	if ((psi_old->f_[i][j] >= 0) && (psi_old->f_[i][j] <= psi_pv)) {
	  jphi->f_[i][j] = R*coeff*(kap0*kap0 + 1)/mu0;
	}
	else {
	  jphi->f_[i][j] = 0;
          printf("out of bounds \n");
	} 
      }
    }
    int R_index = (grid->nr_/5)*4;
    int z_index = (grid->nz_/5)*4;
    double Rp = grid->R_[R_index];
    double zp = grid->z_[z_index]; 
    double psi_calc = 0;
    for (int i = 0; i < grid->nr_; ++i) {
      R = grid->R_[i];
      for (int j = 0; j < grid->nz_; ++j) {
        z = grid->z_[j];
        psi_calc += mu0*green_fcn(R,z,Rp,zp)*jphi->f_[i][j]*dr*dz;
      }
    }

    double *g_factor1 = new double[perim]();
    double *g_factor2 = new double[perim]();

    /******************************************
     * Initialize g-factors for area integral(s)
     ******************************************/
    int I,J;
    for (int l = 0; l < perim; ++l) {
      I = b->LtoI(l);
      J = b->LtoJ(l);
      R = grid->R_[I];
      z = grid->z_[J];
      if (l>=0 && l<=(nr-2)) {
        g_factor1[l] = -(dr/dz)*(green_fcn(R,z+dz/2.0,Rp,zp) - green_fcn(R,z-dz/2.0,Rp,zp));
        g_factor2[l] = -(dr/dz)*(psi_old->f_[I][1] - psi_old->f_[I][0]); 
      }
      else if (l>=(nr-1) && l<=(nr+nz-3)) {
        g_factor1[l] = (dz/dr)*(green_fcn(R+dr/2.0,z,Rp,zp) - green_fcn(R-dr/2.0,z,Rp,zp));
        g_factor2[l] = (dz/dr)*(psi_old->f_[nr-1][J] - psi_old->f_[nr-2][J]);
      }
      else if (l>=(nr+nz-2) && l<=(2*nr+nz-4)) {
        g_factor1[l] = (dr/dz)*(green_fcn(R,z+dz/2.0,Rp,zp) - green_fcn(R,z-dz/2.0,Rp,zp));
        g_factor2[l] = (dr/dz)*(psi_old->f_[I][nz-1] - psi_old->f_[I][nz-2]);
      }
      else {
        g_factor1[l] = -(dz/dr)*(green_fcn(R+dr/2.0,z,Rp,zp) - green_fcn(R-dr/2.0,z,Rp,zp));
        g_factor2[l] = -(dz/dr)*(psi_old->f_[1][J] - psi_old->f_[0][J]);
      }

      psi_calc += psi_old->f_[I][J]*g_factor1[l]/R;
      psi_calc += -green_fcn(R,z,Rp,zp)*g_factor2[l]/R;
    }

    BOOST_CHECK_CLOSE(psi_calc, psi_old->f_[R_index][z_index],5);
}






