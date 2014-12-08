#include "elliptic_solver.h"
#include <math.h>
#include "field.h"

EllipticSolver::EllipticSolver(const Grid &Grid, Field &Psi){
  Grid_ = &Grid;
  Psi_ = &Psi;
  Psi_prev_ = new Field(Grid.nr_,Grid.nz_);
  a = new double*[Grid.nr_];
  b = new double*[Grid.nr_];
  c = new double*[Grid.nr_];
  d = new double*[Grid.nr_];
  f = new double*[Grid.nr_];
  for(int i = 0; i < Grid.nr_; ++i) {
    a[i] = new double[Grid.nz_]();
    b[i] = new double[Grid.nz_]();
    c[i] = new double[Grid.nz_]();
    d[i] = new double[Grid.nz_]();
    f[i] = new double[Grid.nz_]();
  }
}

EllipticSolver::~EllipticSolver(){
  delete Psi_prev_;
  for (int i = 0; i < Grid_->nr_; ++i) {
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] f;
  }
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
  delete [] f;
}

/*!
 * Calculate maximum of |Psi - Psi_prev| over grid
 */
double EllipticSolver::norm_max(const Field &Psi, const Field &Psi_prev) {
  double max = 0;
  int nr = Grid_->nr_;
  int nz = Grid_->nz_;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0 ;j < nz; ++j) {
      if (abs(Psi.f_[i][j] - Psi_prev.f_[i][j]) > max)
        max = abs(Psi.f_[i][j] - Psi_prev.f_[i][j]);
    }
  }
  return max;
}

/*!
 * Calculates 2-norm of diffence between Psi and Psi_prev for convergence testing
 */
double EllipticSolver::norm() {
  double sum = 0;
  double nr = Grid_->nr_;
  double nz = Grid_->nz_;
  for (int i = 0; i < nr; ++i) {
    for(int j = 0; j < nz; ++j) {
      sum += (Psi_->f_[i][j]-Psi_prev_->f_[i][j])*(Psi_->f_[i][j]-Psi_prev_->f_[i][j]);
    }
  }
  return sqrt(sum);
}

/*!
 * Norm of residuals
 */
double EllipticSolver::residuals(const Field &Psi, const Field &Psi_prev) {
  return 0;
}

/*!
 * Blend with old solution
 */
void EllipticSolver::iter(double omega) {
  int nr = Grid_->nr_;
  int nz = Grid_->nz_;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0 ;j < nz; ++j) {
      Psi_->f_[i][j] = omega*Psi_->f_[i][j] + (1-omega)*Psi_prev_->f_[i][j];
    }
  }
}

/*!
 * Enforce boundary condition
 */
void EllipticSolver::boundary(Field &Psi, const Field &Psi_prev) {
  int nr = Grid_->nr_;
  int nz = Grid_->nz_;
  for (int i = 0; i < nr; ++i) {
    Psi.f_[i][0] = Psi_prev.f_[i][0];
    Psi.f_[i][nz-1] = Psi_prev.f_[i][nz-1];
  }
  for (int i = 0; i < nz; ++i) {
    Psi.f_[0][i] = Psi_prev.f_[i][0];
    Psi.f_[nr-1][i] = Psi_prev.f_[nr-1][0];
  }
}