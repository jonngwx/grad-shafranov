#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include gauss-seidel.h

class GaussSeidel : public EllipticSolver {
public:
  GaussSeidel(const Field &Psi_n, int max_iter, double epsilon);
  void step(const Field &Psi_n, const Field &Psi_n+);
  ~GaussSeidel();
private:
  // Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **e;
};

#endif