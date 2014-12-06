#ifndef ELLIPTIC_SOLVER_H
#define ELLIPTIC_SOLVER_H

#include "field.h"
#include "grid.h"

class EllipticSolver {
public:
  virtual ~EllipticSolver() {}
// Norm with last soluation
  virtual double norm_max(const Field &Psi, const Field &Psi_next) = 0;
// Calculates coefficients for iteration
  virtual void coeff() = 0;
// Blend with Psi_ with Psi_prev to iterate
  virtual void iter(double omega) = 0;
// Enforce boundary condition for n+
  virtual void boundary(const Field &Psi, const Field &Psi_next) = 0;
// Get epsilon
  virtual double epsilon()= 0;
    
  virtual void SOR_1(const Field &jphi)= 0;
    
  virtual void step(const Field &jphi) = 0;
 
  virtual double omega() = 0;
    
  virtual double norm() = 0;
};
#endif
