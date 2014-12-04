#include elliptic-solver.h
#include <math.h>

// Initialize solver with current Psi
void init(const Field &Psi) {
  Psi_ = Psi;
}

// Calculate maximum of |Psi - Psi_prev| over grid
double EllipticSolver::norm_max(const Field &Psi, const Field &Psi_prev) {
  double max = 0;
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      if (abs(Psi.f[i][j] - Psi_prev.f[i][j]) > max)
        max = abs(Psi.f[i][j] - Psi_prev.f[i][j]);
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
      Psi_.f[i][j] = alpha*Psi_.f[i][j] + (1-alpha)*Psi_prev.f[i][j];
    }
  }
}

// Enforce boundary condition
void EllipticSolver::boundary(const Field &Psi, const Field &Psi_next) {
  for (int i = 0; i < nr_; ++i) {
    Psi_next.f[i][0] = Psi.f[i][0];
    Psi_next.f[i][nz-1] = Psi.f[i][nz-1];
  }
  for (int i = 0; i < nz_; ++i) {
    Psi_next.f[0][i] = Psi.f[i][0];
    Psi_next.f[nr-1][i] = Psi.f[nr-1][0];
  }
}

double EllipticSolver::epsilon() {
  return epsilon_;
}