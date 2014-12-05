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
  // Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **e;
};

#endif
