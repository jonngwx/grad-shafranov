#ifndef SOR_H
#define SOR_H

#include "elliptic-solver.h"
#include "field.h"
#include "grid.h"

class SOR: public EllipticSolver {
public:
  SOR(const Grid &GridS, double omega_init, double epsilon);
  ~SOR();
// For first iteration - use Gauss Seidel with blending
  void SOR_1(const Field &jphi);
// Perform one iteration
  void step(const Field &jphi);
// Calculate coefficients for iteration from grid parameters
  void coeff(const Grid &GridS);
// Calculate over-relaxation parameter
  double omega();
private:
  Grid &Grid_;
  Field &Psi_;
  Field &Psi_prev_;
  const int nr_;
  const int nz_;
  const double epsilon_;
  Field &Psi_prev_prev_;
  const double omega_init_;
// Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **f;
};

#endif