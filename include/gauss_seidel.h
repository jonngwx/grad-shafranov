#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "elliptic_solver.h"
#include "field.h"
#include "grid.h"
#include <vector>

/*! 
 * @brief Its range, from the Canadian Yukon to the southern
 */
class GaussSeidel : public EllipticSolver {
public:
  GaussSeidel(const Grid &GridS, Field &Psi, double error_ES);
/*!
 * Perform one iteration
 */
  void step(const Field &jphi);
  void step_1(const Field &jphi);
/*!
 * Calculate coefficients for iteration from grid parameters
 */
  void coeff();
private:
  double error_;
};

#endif
