#ifndef ELLIPTIC_SOLVER_H
#define ELLIPTIC_SOLVER_H

#include "field.h"
#include "grid.h"
#include <vector>

class EllipticSolver {
public:
  EllipticSolver(const Grid &Grid, Field &Psi);
  virtual ~EllipticSolver();
/*!
 * Norm with last soluation
 */
  double norm_max(const Field &Psi, const Field &Psi_prev);
/*!
 * Blend Psi_ with Psi_prev_ to iterate with parameter omega
 */
  void iter(double omega);
/*!
 * Enforce boundary condition for Psi_ using Psi_prev_
 */
  void boundary(Field &Psi, const Field &Psi_prev);
/*!
 * Calculates 2-norm of diffence between Psi_ and Psi_prev_ for convergence testing
 */
  double norm();
/*!
 * Perform first iteration
 */
  virtual void step_1(const Field &jphi) = 0;
/*!
 * Perform one iteration
 */
  virtual void step(const Field &jphi) = 0;
/*!
 * Calculates coefficients for iteration
 */
  virtual void coeff() = 0;
/*!
 * Norm of residuals
 */
  double residuals(const Field &Psi, const Field &Psi_prev);
  
protected:
  const Grid &Grid_;
  Field &Psi_;
  Field Psi_prev_;
// Coefficients
  std::vector<double> A;
  double B, C;
};
#endif
