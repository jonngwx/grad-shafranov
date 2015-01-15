/*!
 * @file j_solver.h
 * @author Jonathan Ng
 * @brief Header declarations for JSolverNSTX.
 */

#ifndef J_SOLVER_NSTX_H_
#define J_SOLVER_NSTX_H_

#include "grid.h"
#include "field.h"
#include "j_solver.h"

/*! 
 * @brief Calculates current in NSTX. Assumes p' = Ax(1-x), gg' = c x 
 */
class JSolverNSTX : public JSolver {
 public:
  /**
   * @brief Constructor for JSolverAlpha
   * @param P0 Pressure on axis
   * @param g0 B_T*R on axis
   * @param n2 exponent of toroidal field function
   * @param Ip Total plasma current
   * @param grid pointer to grid with axis data
   */
  JSolverNSTX(double P0, double g0, double Ip, double n2, Grid *grid);
  ~JSolverNSTX();
  /*! which variables are 'in' and which are 'out' ? */
  void update(Field *jphi, Field *psi, Field *p, Field *g);
    
 private:
  double P1_; //!< coefficient of pressure function
  double n2_; //!< exponent in gg' function
};

#endif  // J_SOLVER_NSTX_H_
