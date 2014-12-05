#ifndef ELLIPTIC_SOLVER_H
#define ELLIPTIC_SOLVER_H

#include "field.h"
#include "grid.h"

class EllipticSolver {
public:
  virtual ~EllipticSolver();
// Initialize solver with current Psi
  void init(const Field &Psi, const Grid &Grid);
// Norm with last soluation
  double norm_max(const Field &Psi, const Field &Psi_next);
// Norm of residuals
  double residuals(const Field &Psi, const Field &Psi_next);
// Calculates coefficients for iteration
  virtual void coeff(const Grid &GridS) = 0;
// Blend with Psi_ with Psi_prev to iterate
  void iter();
// Enforce boundary condition for n+
  void boundary(const Field &Psi, const Field &Psi_next);
// Get epsilon
  void epsilon();
private:
  Grid &Grid_;
  Field &Psi_;
  Field &Psi_prev_;
  const int nr_;
  const int nz_;
  const double epsilon_;
};

#endif
