#include "include/sor.h"
#include "include/elliptic_solver.h"
#include "grid.h"
#include "field.h"
#include <math.h>

SOR::SOR(const Grid &GridS, const Field &Psi, const Field &Psi_prev, const Field &Psi_prev_prev, double omega_init, double epsilon) :
  Grid_(GridS),
  Psi_(Psi),
  Psi_prev_(Psi_prev),
  Psi_prev_prev_(Psi_prev_prev),
  nr_(GridS.nr_),
  nz_(GridS.nz_),
  epsilon_(epsilon),
  omega_init_(omega_init)
   {
    *a = new double[nr_];
    *b = new double[nr_];
    *c = new double[nr_];
    *d = new double[nr_];
    *f = new double[nr_];
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

// Calculates coefficients for iteration
void SOR::coeff() {
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    double *r = Grid_.R_;
    double e = 2*(1/(dr*dr) + 1/(dz*dz));
    for (int i = 1; i < nr_-1; ++i) {
        for(int j = 1; j < nz_-1; ++j) {
            e = 2*(1/(dr*dr) + 1/(dz*dz));
            a[i][j] = (1/(dr*dr) + 1/(2*r[i]*dr))*(1/e);
            b[i][j] = (1/(dz*dz))*(1/e);
            c[i][j] = b[i][j];
            d[i][j] = (1/(dr*dr) - 1/(2*r[i]*dr))*(1/e);
            f[i][j] = 1;
        }
    }
}

// For first iteration - use Gauss Seidel with blending
void SOR::SOR_1(const Field &jphi) {
    boundary(Psi_, Psi_prev_);
    
    // Iterate over non-boundary region
    for (int i = 1; i < nr_-1; ++i) {
        for(int j = 1; j < nz_-1; ++j) {
            Psi_.f_[i][j] = a[i][j]*Psi_prev_.f_[i+1][j] + b[i][j]*Psi_.f_[i-1][j] + c[i][j]*Psi_.f_[i][j+1] + d[i][j]*Psi_.f_[i][j-1] - f[i][j]*jphi.f_[i][j];
        }
    }
    iter(omega_init_);
}

// Enforce boundary condition
void SOR::boundary(const Field &Psi, const Field &Psi_next) {
    for (int i = 0; i < nr_; ++i) {
        Psi_next.f_[i][0] = Psi.f_[i][0];
        Psi_next.f_[i][nz_-1] = Psi.f_[i][nz_-1];
    }
    for (int i = 0; i < nz_; ++i) {
        Psi_next.f_[0][i] = Psi.f_[i][0];
        Psi_next.f_[nr_-1][i] = Psi.f_[nr_-1][0];
    }
}

// Calculate maximum of |Psi - Psi_prev| over grid
double SOR::norm_max(const Field &Psi, const Field &Psi_prev) {
    double max = 0;
    for (int i = 0; i < nr_; ++i) {
        for (int j = 0 ;j < nz_; ++j) {
            if (abs(Psi.f_[i][j] - Psi_prev.f_[i][j]) > max)
                max = abs(Psi.f_[i][j] - Psi_prev.f_[i][j]);
        }
    }
    return max;
}

double SOR::epsilon() {
    return epsilon_;
}

// Iterate with over-relaxation parameter omega
void SOR::step(const Field &jphi) {
  // Save Psi_ to Psi_prev and Psi_prev to Psi_prev_prev
  for (int i = 1; i < nr_-1; ++i) {
    for(int j = 1; j < nz_-1; ++j) {
      Psi_prev_prev_.f_[i][j] = Psi_prev_.f_[i][j];
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }
  // Enforce boundary condition on Psi_
  // Field *J = RHS(Psi_prev_);
  boundary(Psi_, Psi_prev_);
  double om = omega();
  // Iterate over non-boundary region
  for (int i = 1; i < nr_-1; ++i) {
    for(int j = 1; j < nz_-1; ++j) {
      Psi_.f_[i][j] = om*(a[i][j]*Psi_.f_[i-1][j] + b[i][j]*Psi_.f_[i][j-1] + c[i][j]*Psi_prev_.f_[i][j+1] + d[i][j]*Psi_prev_.f_[i+1][j] + f[i][j]*jphi.f_[i][j]) + (1-om)*Psi_prev_.f_[i][j];
    }
  }
}

// Blend with old solution
void SOR::iter(double omega) {
    
    for (int i = 0; i < nr_; ++i) {
        for (int j = 0 ;j < nz_; ++j) {
            Psi_.f_[i][j] = omega*Psi_.f_[i][j] + (1-omega)*Psi_prev_.f_[i][j];
        }
    }
}

double SOR::norm() {
    double sum = 0;
    for (int i = 0; i < nr_; ++i) {
        for(int j = 0; j < nz_; ++j) {
            sum += (Psi_.f_[i][j]-Psi_prev_.f_[i][j])*(Psi_.f_[i][j]-Psi_prev_.f_[i][j]);
        }
    }
    return sqrt(sum);
}

// Calculate over-relaxation parameter
double SOR::omega() {
  double delta = norm_max(Psi_prev_, Psi_prev_prev_)/norm_max(Psi_, Psi_prev_);
  return 2/(1 + sqrt(1 - delta));
}

