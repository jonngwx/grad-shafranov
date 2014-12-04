#ifndef ELLIPTIC_SOLVER_H
#define ELLIPTIC_SOLVER_H

#include elliptic-solver.h

class EllipticSolver {
public:
  virtual ~EllipticSolver();
// Norm with last soluation
  double norm(const Field &Psi_n, const Field &Psi_n-);
// Norm of residuals
  double residuals(const Field &Psi_n, const Field &Psi_n-);
// Calculates coefficients for iteration
  virtual void coeff(const Field &Psi_n) = 0
// Blend with old solution
  void iter(const Field &Psi_n, const Field &Psi_n-);
private:
  const int dimen_;
  const double epsilon_;
  const int max_iter_;
};

#endif