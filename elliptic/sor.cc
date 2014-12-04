#include sor.h
#include <math.h>

SOR::SOR(const Grid &GridS, double omega_init, double epsilon) :
  nr_(GridS.nr),
  nz_(GridS.nz),
  omega_init_(omega_init),
  max_iter_(max_iter),
  epsilon_(epsilon),
  {
    a = new double[nr]();
    b = new double[nr]();
    c = new double[nr]();
    d = new double[nr]();
    e = new double[nr]();
    f = new double[nr]();
    for(int i = 0; i< nr; ++i) {
      a[i] = new double[nz]();
      b[i] = new double[nz]();
      c[i] = new double[nz]();
      d[i] = new double[nz]();
      e[i] = new double[nz]();
      f[i] = new double[nz]();
    }
}

SOR::~SOR() {
  for (int i = 0; i < nr_; ++i) {
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] e;
    delete [] f;
  }
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
  delete [] e;
  delete [] f;
}

// Calculates coefficients for iteration
void coeff(const Field &Psi_n) {
  
}

// For first iteration - use Gauss Seidel with blending
void SOR_1() {
  // Field *J = RHS(Psi_);
  boundary(Psi, Psi_next);
  
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f[i][j] = a[i][j]*Psi_prev_.f[i+1][j] + b[i][j]*Psi_[i-1][j] + c[i][j]*Psi_.f[i][j+1] + d[i][j]*Psi_[i][j-1] - e[i][j]*J.f[i][j];
    }
  }
  iter();
}

// Iterate with over-relaxation parameter omega
void step() {
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_prev_.f[i][j] = Psi_prev_.f[i][j];
      Psi_prev_.f[i][j] = Psi_.f[i][j];
    }
  }
  // Enforce boundary condition on Psi_
  // Field *J = RHS(Psi_prev_);
  boundary(Psi_, Psi_prev_);
  double omega = omega();
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_.f[i][j] = omega*(a[i][j]*Psi_prev_.f[i+1][j] + b[i][j]*Psi_.f[i-1][j] + c[i][j]*Psi_prev_.f[i][j+1] + d[i][j]*Psi_.f[i][j-1] + e[i][j]*Psi_prev_[i][j] + f[i][j]*J.f[i][j]);
    }
  }
}

// Calculate over-relaxation parameter
double omega() {
  double delta = norm_max(Psi_prev_, Pxi_prev_prev_)/norm_max(Psi_, Psi_prev_);
  return 2/(1 + sqrt(1 - delta));
}

#endif