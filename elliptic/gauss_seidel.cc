#include "gauss-seidel.h"
#include "field.h"
#include "grid.h"


GaussSeidel::GaussSeidel(const Grid &GridS, double epsilon) :
  nr_(GridS.nr_),
  nz_(GridS.nz_),
  epsilon_(epsilon) {
    a = new double[nr]();
    b = new double[nr]();
    c = new double[nr]();
    d = new double[nr]();
    e = new double[nr]();
    for(int i = 0; i < nr; ++i) {
      a[i] = new double[nz]();
      b[i] = new double[nz]();
      c[i] = new double[nz]();
      d[i] = new double[nz]();
      e[i] = new double[nz]();
    }
}

GaussSeidel::~GaussSeidel(){
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

void coeff(const Grid &GridS) {
  
}

void GaussSeidel::step(const Field &jphi){
// Save Psi_ to Psi_prev
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }
  
  boundary(Psi_prev_, Psi_);

  for (int i = 0; i < nr_; ++i) {
    for (int j = 0 ;j < nz_; ++j) {
      Psi_.f_[i][j] = a[i][j]*Psi_prev_.f_[i+1][j] + b[i][j]*Psi_.f_[i-1][j] + c[i][j]*Psi_prev_.f_[i][j+1] + d[i][j]*Psi_.f_[i][j-1] + e[i][j]*Psi_prev_[i][j] + f[i][j]*jphi.f_[i][j];
    }
  }
  iter();
}