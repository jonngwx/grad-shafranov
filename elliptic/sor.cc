#include sor.h
#include <math.h>

SOR::SOR(const Field &Psi, double omega_init, int max_iter, double epsilon) :
  nr_(Psi.nr),
  nz_(Psi.nz),
  omega_init_(omega_init),
  max_iter_(max_iter),
  epsilon_(epsilon),
  {
    a = new double[nr]();
    b = new double[nr]();
    c = new double[nr]();
    d = new double[nr]();
    e = new double[nr]();
    for(int i = 0; i< nr; ++i) {
      a[i] = new double[nz]();
      b[i] = new double[nz]();
      c[i] = new double[nz]();
      d[i] = new double[nz]();
      e[i] = new double[nz]();
    }
}

SOR::~SOR() {
  for (int i = 0; i < nr_; ++i) {
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] e;
  }
  delete [] a;
  delete [] b;
  delete [] c;
  delete [] d;
  delete [] e;
}

// Calculates coefficients for iteration
void coeff(const Field &Psi_n) {
  
}

// For first iteration - use Gauss Seidel with blending
void SOR_1(const Field &Psi_n, const Field &Psi_n+) {
  boundary(Psi_n, Psi_n+);
  // Iterate over non-boundary region
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_n+.f[i][j] =
    }
  }
}

// Iterate with over-relaxation parameter omega
void step(const Field &Psi_n, const Field &Psi_n+) {
  boundary(Psi_n, Psi_n+);
  
}

// Calculate over-relaxation parameter
double omega(const Field &Psi_n-, const Field &Psi_n, const Field &Psi_n+) {
  double delta = norm_max(Psi_n, Pxi_n-)/norm_max(Psi_n+, Psi_n);
  return 2/(1 + sqrt(1 - delta));
}

#endif