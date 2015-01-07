#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SlowBoundary
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/boundary.h"
#include "../include/slow_boundary.h" 
#include "../include/grid.h"
#include "../include/field.h"
#include <math.h>

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
  double Rmin = 5.0;
  double Rmax = 7.0;
  double zmin = -1.0;
  double zmax = 1.0;
  int nr = 20;
  int nz = 20;

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
  double eps = 0.3; // inverse aspect ratio
  const double mu0 = 4*M_PI*1e-7; // mag perm of free space
  double psi_pv = eps*eps*kap0*R0*R0*B0/(2*q0); // defines psi at plas/vac interf.
  
  /************************************
   * Initialize psi, psi_old, and jphi 
   ************************************/
  double R,Z;
  double coeff = B0/(R0*R0*kap0*q0);
  for (int i = 0; i < grid->nr_; ++i) {
    R = grid->R_[i];
    for (int j = 0; j < grid->nz_; ++j) {
      Z = grid->z_[j];
      psi->f_[i][j] = (coeff/2.0)*(R*R*Z*Z + (kap0*kap0/4.0)*pow((R*R-R0*R0),2));
      psi_old->f_[i][j] = psi->f_[i][j];
      if ((psi->f_[i][j] >= 0) && (psi->f_[i][j] <= psi_pv)) {
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

  /*****************************************
   * Compare new bdy pt to its old version
   *****************************************/
  int perim = 2*(nr + nz - 2);
  int I=0;
  int J=0;
  for (int l=0; l<perim; ++l) {
    I = psib->LtoI(l);
    J = psib->LtoJ(l);
    BOOST_CHECK_CLOSE(psi_old->f_[I][J], psi->f_[I][J],4);
  }


}
