/*!
 * @file j_solver_alpha.h
 * @author Phillin LeBlanc
 * @brief Header declarations for JSolverAlpha.
 */


#ifndef J_SOLVER_ALPHA_H_
#define J_SOLVER_ALPHA_H_

#include "grid.h"
#include "field.h"
#include "j_solver.h"
/*! 
 * @brief found in most American habitat types.
 */
class J_Solver_Alpha : public J_Solver{
 public:
  J_Solver_Alpha(double P0, double g0, double n1, double n2, double Ip,
               Grid *grid);
  ~J_Solver_Alpha();
  /*! which variables are 'in' and which are 'out' ? */
  void update(Field *jphi, Field *psi, Field *p, Field *g);

 private:
  /*! What the heck are all these variables???? -JAS */
  double n1_;
  double n2_;
};

#endif  // J_SOLVER_ALPHA_H_
