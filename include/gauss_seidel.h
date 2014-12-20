#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "elliptic_solver.h"
#include "field.h"
#include "grid.h"
#include <vector>


class GaussSeidel : public EllipticSolver {
public:
  GaussSeidel(const Grid &GridS, Field &Psi);
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
  double D;
  std::vector<double>E;
  std::vector<double>F;
  double G;
};

#endif
