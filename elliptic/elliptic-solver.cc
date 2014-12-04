#include elliptic-solver.h
#include <math.h>

EllipticSolver::EllipticSolver() {
  
}

// Calculate maximum of |Psi_n - Psi_n-| over grid
double EllipticSolver::norm_max(const Field &Psi_n, const Field &Psi_n-) {
  double max = 0;
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      if (abs(Psi_n.f[i][j] - Psi_n-.f[i][j]) > max)
        max = abs(Psi_n.f[i][j] - Psi-.f[i][j]);
    }
  }
  return max;
}

// Enforce boundary condition
EllipticSolver::boundary(const Field &Psi_n, const Field &Psi_n+) {
  for (int i = 0; i < nr_; ++i) {
    Psi_n.f[i][0] = Psi_n+.f[i][0];
    Psi_n.f[i][nz-1] = Psi_n+.f[i][nz-1];
  }
  for (int i = 0; i < nz_; ++i) {
    Psi_n.f[0][i] = Psi_n+.f[i][0];
    Psi_n.f[nr-1][i] = Psi_n+.f[nr-1][0];
  }
}