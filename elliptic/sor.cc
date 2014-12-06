#include "sor.h"
#include "grid.h"
#include "field.h"
#include <math.h>

SOR::SOR(const Grid &GridS, double omega_init, double epsilon) :
  Grid_(GridS),
  nr_(GridS.nr_),
  nz_(GridS.nz_),
  omega_init_(omega_init),
  epsilon_(epsilon) {
    a = new double[nr_]();
    b = new double[nr_]();
    c = new double[nr_]();
    d = new double[nr_]();
    f = new double[nr_]();
    for(int i = 0; i< nr_; ++i) {
      a[i] = new double[nz_]();
      b[i] = new double[nz_]();
      c[i] = new double[nz_]();
      d[i] = new double[nz_]();
      f[i] = new double[nz_]();
    }
}

SOR::~SOR() {
  for (int i = 0; i < nr_; ++i) {
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

// Enforce boundary condition
void SOR::EllipticSolver::boundary(const Field &Psi, const Field &Psi_next) {
    for (int i = 0; i < nr_; ++i) {
        Psi_next.f_[i][0] = Psi.f_[i][0];
        Psi_next.f_[i][nz-1] = Psi.f_[i][nz-1];
    }
    for (int i = 0; i < nz_; ++i) {
        Psi_next.f_[0][i] = Psi.f_[i][0];
        Psi_next.f_[nr-1][i] = Psi.f_[nr-1][0];
    }
}

// For first iteration - use Gauss Seidel with blending
void SOR::SOR_1(const Field &jphi) {
  boundary(Psi, Psi_next);
  
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f_[i][j] = a[i][j]*Psi_prev_.f_[i+1][j] + b[i][j]*Psi_[i-1][j] + c[i][j]*Psi_.f_[i][j+1] + d[i][j]*Psi_[i][j-1] - e[i][j]*jphi.f_[i][j];
    }
  }
  iter();
}

// Calculate coefficients for iteration from grid parameters
void SOR::coeff(const Grid &GridS) {
    double dr = Grid.dr_;
    double dz = Grid.dz_;
    double e = 2*(1/(dr^2) + 1/(dz^2));
    for (int i = 1; i < nr-1; ++i) {
        for(int j = 1; j < nz-1; ++j) {
            e = 2*(1/(dr^2) + 1/(dz^2));
            a[i][j] = (1/(dr^2) + 1/(2*r[i]*dr))*(1/e);
            b[i][j] = (1/(dz^2))*(1/e);
            c[i][j] = b[i][j];
            d[i][j] = (1/(dr^2) - 1/(2*r[i]*dr))*(1/e);
            f[i][j] = 1;
        }
    }
}

// Iterate with over-relaxation parameter omega
void SOR::step(const Field &jphi) {
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_prev_.f_[i][j] = Psi_prev_.f_[i][j];
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }
  // Enforce boundary condition on Psi_
  // Field *J = RHS(Psi_prev_);
  boundary(Psi_, Psi_prev_);
  double omega = omega();
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f_[i][j] = omega*(a[i][j]*Psi_.f_[i-1][j] + b[i][j]*Psi_.f_[i][j-1] + c[i][j]*Psi_prev_.f_[i][j+1] + d[i][j]*Psi_prev_.f_[i+1][j] + f[i][j]*jphi.f_[i][j]) + (1-omega)*Psi_prev_.f_[i][j];
    }
  }
}

double SOR::norm() {
    double sum = 0;
    for (int i = 0; i < nr; ++i) {
        for(int j = 0; j < nz; ++j) {
            sum += (Psi_.f_[i][j]-Psi_prev_.f_[i][j])^2;
        }
    }
    return sqrt(sum);
}

// Calculate over-relaxation parameter
double SOR::omega() {
  double delta = norm_max(Psi_prev_, Psi_prev_prev_)/norm_max(Psi_, Psi_prev_);
  return 2/(1 + sqrt(1 - delta));
}

#endif