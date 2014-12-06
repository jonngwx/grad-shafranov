#include "elliptic_solver.h"
#include <math.h>
#include "field.h"

// Initialize solver with current Psi
void init(const Field &Psi) {
  Psi_ = Psi;
}

// Calculate maximum of |Psi - Psi_prev| over grid
double EllipticSolver::norm_max(const Field &Psi, const Field &Psi_prev) {
  double max = 0;
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      if (abs(Psi.f_[i][j] - Psi_prev.f_[i][j]) > max)
        max = abs(Psi.f_[i][j] - Psi_prev.f_[i][j]);
    }
  }
  return max;
}

// Norm of residuals
double EllipticSolver::residuals(const Field &Psi, const Field &Psi_next) {
  
}

// Blend with old solution
void EllipticSolver::iter() {
  double alpha = 0.5;
  
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      Psi_.f_[i][j] = alpha*Psi_.f_[i][j] + (1-alpha)*Psif_prev.f[i][j];
    }
  }
}

// Enforce boundary condition
void EllipticSolver::boundary(const Field &Psi, const Field &Psi_next) {
  for (int i = 0; i < nr_; ++i) {
    Psi_next.f_[i][0] = Psi.f_[i][0];
    Psi_next.f_[i][nz_-1] = Psi.f_[i][nz_-1];
  }
  for (int i = 0; i < nz_; ++i) {
    Psi_next.f_[0][i] = Psi.f_[i][0];
    Psi_next.f_[nr_-1][i] = Psi.f_[nr_-1][0];
  }
}

double EllipticSolver::epsilon() {
  return epsilon_;
}