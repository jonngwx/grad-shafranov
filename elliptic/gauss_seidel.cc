#include "gauss_seidel.h"
#include "field.h"
#include "grid.h"
#include <assert.h>
#include <vector>
#include "elliptic_solver.h"
#include <stdio.h>


/*!
 * @file gauss_seidel.cc
 * @brief Base class implementation of GaussSeidel
 * @section DESCRIPTION
 */

GaussSeidel::GaussSeidel(const Grid &GridS, Field &Psi) :
  EllipticSolver(GridS, Psi) { }

/*!
 * @brief Calculates coefficients for iteration
 */
void GaussSeidel::coeff() {
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    B = (dr*dr*dz*dz)/(2*dz*dz + 2*dr*dr);
    C = 1/(dz*dz);
    for (int i = 0; i < Grid_.nr_; ++i) {
        A[i] = 1/(2*Grid_.R_[i]*dr) + 1/(dr*dr);
    }
}

void GaussSeidel::step_1(const Field &jphi){
  step(jphi);
}

/*!
 * @brief For first iteration - use Gauss Seidel with blending
 * @param jphi current evaluated at current Psi
 */
void GaussSeidel::step(const Field &jphi){
  const double mu0 = 0.0000012566370614; // in SI units
  
// Save Psi_ to Psi_prev
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  for (int i = 0; i < nr; ++i) {
    for(int j = 0; j < nz; ++j) {
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }
// Copy over boundary values
//  boundary(Psi_prev_, Psi_);
  for (int i = 1; i < nr-1; ++i) {
    for (int j = 1;j < nz-1; ++j) {
      Psi_.f_[i][j] = B*(-jphi.f_[i][j]*mu0*Grid_.R_[i] - Psi_prev_.f_[i+1][j]*A[i] + Psi_.f_[i-1][j]*A[i] + Psi_.f_[i][j-1]*C - Psi_prev_.f_[i][j+1]*C);
    }
  }
  iter(0.5);
}
