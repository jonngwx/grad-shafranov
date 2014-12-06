#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "elliptic_solver.h"
#include "field.h"
#include "grid.h"

class GaussSeidel : public EllipticSolver {
public:
  GaussSeidel(const Grid &GridS, double epsilon);
  void step(const Field &jphi);
  ~GaussSeidel();
private:
  Grid &Grid_;
  Field &Psi_;
  Field &Psi_prev_;
  const int nr_;
  const int nz_;
  const double epsilon_;
  // Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **e;
};

#endif
