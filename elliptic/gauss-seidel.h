#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include elliptic-solver.h

class GaussSeidel : public EllipticSolver {
public:
  GaussSeidel(const Grid &GridS, double epsilon);
  void step(const Field &Psi, const Field &Psi_next);
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