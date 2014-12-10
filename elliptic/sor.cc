#include "sor.h"
#include "elliptic_solver.h"
#include "grid.h"
#include "field.h"
#include <math.h>

/*!
 * @file sor.cc
 * @brief Base class implementation of Successive Over-Reduction implementation of EllipticSolver
 * @section DESCRIPTION
 */

SOR::SOR(const Grid &GridS, Field &Psi, double omega_init):
  EllipticSolver(GridS, Psi),
  omega_init_(omega_init) {
  Psi_prev_prev_ = new Field(GridS);
}

SOR::~SOR() {
  delete Psi_prev_prev_;
}

/*!
 * @brief Calculates coefficients for iteration
 */
void SOR::coeff() {
  double dr = Grid_->dr_;
  double dz = Grid_->dz_;
  double nr = Grid_->nr_;
  double nz = Grid_->nz_;
  double *r = Grid_->R_;
  double e = 2*(1/(dr*dr) + 1/(dz*dz));
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      e = 2*(1/(dr*dr) + 1/(dz*dz));
      a[i][j] = (1/(dr*dr) + 1/(2*r[i]*dr))*(1/e);
      b[i][j] = (1/(dz*dz))*(1/e);
      c[i][j] = b[i][j];
      d[i][j] = (1/(dr*dr) - 1/(2*r[i]*dr))*(1/e);
      f[i][j] = 1;
    }
  }
}

/*!
 * @brief For first iteration - use Gauss Seidel with blending
 * @param jphi current evaluated at current Psi
 */
void SOR::step_1(const Field &jphi) {
  double nr = Grid_->nr_;
  double nz = Grid_->nz_;
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_->f_[i][j] = Psi_->f_[i][j];
    }
  }
  boundary(*Psi_, *Psi_prev_);
  // Iterate over non-boundary
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_->f_[i][j] = a[i][j]*Psi_prev_->f_[i+1][j] + b[i][j]*Psi_->f_[i-1][j] + c[i][j]*Psi_->f_[i][j+1] + d[i][j]*Psi_->f_[i][j-1] - f[i][j]*jphi.f_[i][j];
    }
  }
  iter(omega_init_);
}

/*!
 * @brief Iterate with over-relaxation parameter omega
 * @param jphi current evaluated at current Psi
 */
void SOR::step(const Field &jphi) {
  double nr = Grid_->nr_;
  double nz = Grid_->nz_;
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_prev_->f_[i][j] = Psi_prev_->f_[i][j];
      Psi_prev_->f_[i][j] = Psi_->f_[i][j];
    }
  }
  // Enforce boundary condition on Psi_
  boundary(*Psi_, *Psi_prev_);
  double om = omega();
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_->f_[i][j] = om*(a[i][j]*Psi_->f_[i-1][j] + b[i][j]*Psi_->f_[i][j-1] + c[i][j]*Psi_prev_->f_[i][j+1] + d[i][j]*Psi_prev_->f_[i+1][j] + f[i][j]*jphi.f_[i][j]) + (1-om)*Psi_prev_->f_[i][j];
    }
  }
}

/*!
 * @brief Calculate over-relaxation parameter omega for finite differencing at each iteration
 */
double SOR::omega() {
  double delta = norm_max(*Psi_prev_, *Psi_prev_prev_)/norm_max(*Psi_, *Psi_prev_);
  return 2/(1 + sqrt(1 - delta));
}

