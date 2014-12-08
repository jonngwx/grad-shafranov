#include "gauss-seidel.h"
#include "field.h"
#include "grid.h"


GaussSeidel::GaussSeidel(const Grid &GridS, Field &Psi) :
  EllipticSolver(GridS, Psi) { }

/*!
 * Calculate coefficients for iteration from grid parameters
 */
void coeff(const Grid &GridS) {
  
}

/*!
 * Perform one iteration
 */
void GaussSeidel::step(const Field &jphi){
// Save Psi_ to Psi_prev
  double nr = Grid_->nr_;
  double nz = Grid_->nz_;
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_->f_[i][j] = Psi_->f_[i][j];
    }
  }
  
  boundary(*Psi_prev_, *Psi_);

  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      Psi_->f_[i][j] = a[i][j]*Psi_prev_->f_[i+1][j] + b[i][j]*Psi_->f_[i-1][j] + c[i][j]*Psi_prev_->f_[i][j+1] + d[i][j]*Psi_->f_[i][j-1] + e[i][j]*Psi_prev_->f_[i][j] + f[i][j]*jphi.f_[i][j];
    }
  }
  iter();
}