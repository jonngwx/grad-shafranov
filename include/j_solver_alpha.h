/*!
 * @file j_solver_alpha.h
 * @author Peter J. Bolgert
 * @brief Header declarations for JSolverAlpha.
 */


#ifndef J_SOLVER_ALPHA_H_
#define J_SOLVER_ALPHA_H_

#include "grid.h"
#include "field.h"
#include "j_solver.h"
/*! 
 * @brief Calculates current according to simple
 * formulas. It conserves total current and the 
 * pressure at the magnetic axis.
 */
class JSolverAlpha : public JSolver{
 public:
  /**
   * @brief Constructor for JSolverAlpha
   * @param P0 Pressure on axis
   * @param g0 B_T*R on axis
   * @param n1 exponent of pressure function
   * @param n2 exponent of toroidal field function
   * @param Ip Total plasma current
   * @param grid pointer to grid with axis data
   */
  JSolverAlpha(double P0, double g0, double n1, double n2, double Ip,
               Grid *grid);
  ~JSolverAlpha();
  void update(Field *jphi, Field *psi, Field *p, Field *g);

 private:
  double n1_;
  double n2_;
};

#endif  // J_SOLVER_ALPHA_H_
