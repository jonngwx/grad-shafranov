/*!
 * @file gauss_seidel.h
 * @brief Declarations for class GaussSeidel
 */

#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include <vector>
#include "field.h"
#include "grid.h"
#include "elliptic_solver.h"

/*! 
 * @brief An implementation of a Gauss-Seidel matrix solver method.
 *
 * Here is a reference to what it is and how it works.
 *
 */
class GaussSeidel : public EllipticSolver {
public:
  /*! 
   * @brief Constructor for GaussSeidel
   * How does it use these variables??
   */
  GaussSeidel(const Grid &GridS, Field &Psi, double error_ES);
  ~GaussSeidel();
  /*!
   * @brief A step of the Gauss-Seidel algorithm
   * @param jphi current evaluated at current Psi
   *
   * This is how it works.
   */
  void step(const Field &jphi);
  /*!
   * @brief [deprecated] identical to step 
   * @param jphi current evaluated at current Psi
   */
  void step_1(const Field &jphi);
  /*!
   * @brief Calculate coefficients for iteration from grid parameters
   */
  void coeff();
private:
  double error_; //!< MEANS SOMETHING
};

#endif
