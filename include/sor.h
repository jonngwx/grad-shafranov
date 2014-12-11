#ifndef SOR_H
#define SOR_H

#include "elliptic_solver.h"
#include "field.h"
#include "grid.h"
#include <vector>

class SOR : public EllipticSolver {
public:
  SOR(const Grid &GridS, Field &Psi, double omega_init);
  virtual ~SOR();
  /*!
   * Calculate coefficients for iteration from grid parameters
   */
  void coeff();
  /*!
   * Perform one iteration
   */
  void step(const Field &jphi);
  /*!
   * For first iteration - use Gauss Seidel with blending
   */
  void step_1(const Field &jphi);
  /*!
   * Calculate over-relaxation parameter
   */
  double omega();
private:
  Field Psi_prev_prev_;
  const double omega_init_;
};

#endif
