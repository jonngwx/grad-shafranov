#include "sor.h"
#include "elliptic_solver.h"
#include "grid.h"
#include "field.h"
#include <math.h>
#include <stdio.h>
#include <vector>


/*!
 * @file sor.cc
 * @brief Base class implementation of Successive Over-Reduction implementation of EllipticSolver
 * @section DESCRIPTION
 */

SOR::SOR(const Grid &GridS, Field &Psi, double omega_init):
  EllipticSolver(GridS, Psi),
  Psi_prev_prev_(GridS),
  omega_init_(omega_init) {
}

SOR::~SOR() {}

/*!
 * @brief Calculates coefficients for iteration
 */
void SOR::coeff() {
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    B = (dr*dr*dz*dz)/(2*dz*dz + 2*dr*dr);
    C = 1/(dz*dz);
    for (int i = 0; i < Grid_.nr_; ++i) {
        A[i] = 1/(2*Grid_.R_[i]*dr) + 1/(dr*dr);
    }
}

/*!
 * @brief For first iteration - use Gauss Seidel with blending
 * @param jphi current evaluated at current Psi
 */
void SOR::step_1(const Field &jphi) {
  const double mu0 = 0.0000012566370614; // in SI units
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  // Save Psi_ to Psi_prev
  for (int i = 0; i < nr; ++i) {
    for(int j = 0; j < nz; ++j) {
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
      Psi_prev_prev_.f_[i][j] = Psi_prev_.f_[i][j];
    }
  }

  // Iterate over non-boundary
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f_[i][j] = B*(-jphi.f_[i][j]*mu0*Grid_.R_[i] - Psi_prev_.f_[i+1][j]*A[i+1] + Psi_.f_[i-1][j]*A[i-1] + Psi_.f_[i][j-1]*C - Psi_prev_.f_[i][j+1]*C);
    }
  }
  iter(omega_init_);
}

/*!
 * @brief Iterate with over-relaxation parameter omega
 * @param jphi current evaluated at current Psi
 */
void SOR::step(const Field &jphi) {
  const double mu0 = 0.0000012566370614; // in SI units
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_prev_.f_[i][j] = Psi_prev_.f_[i][j];
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }
  // Enforce boundary condition on Psi_
//  boundary(Psi_, Psi_prev_);
  double om = omega();
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f_[i][j] = (1-om)*Psi_prev_.f_[i][j] + om*B*(-jphi.f_[i][j]*mu0*Grid_.R_[i] - Psi_prev_.f_[i+1][j]*A[i+1] + Psi_.f_[i-1][j]*A[i-1] + Psi_.f_[i][j-1]*C - Psi_prev_.f_[i][j+1]*C);
    }
  }
}

/*!
 * @brief Calculate over-relaxation parameter omega for finite differencing at each iteration
 */
double SOR::omega() {
  double delta = norm_max(Psi_prev_, Psi_prev_prev_)/norm_max(Psi_, Psi_prev_);
  return 2/(1 + sqrt(1 - delta));
}

