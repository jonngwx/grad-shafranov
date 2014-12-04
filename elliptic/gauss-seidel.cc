#include gauss-seidel.h

GaussSeidel::GaussSeidel(const Field &Psi_n, int max_iter, double epsilon) :
  nr_(Psi_n.nr),
  nz_(Psi_n.nz),
  max_iter_(max_iter),
  epsilon_(epsilon) {
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

~GaussSeidel::GaussSeidel(){
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

void coeff(const Field &Psi_n) {
  
}

void GaussSeidel::step(const Field &Psi_n, const Field &Psi_n+){
    
}