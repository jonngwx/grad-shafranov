/*!
 * @file j_solver_alpha.h
 * @author Phillin LeBlanc
 * @brief Header declarations for JSolverAlpha.
 */


#ifndef J_SOLVER_ALPHA_H_
#define J_SOLVER_ALPHA_H_

#include "grid.h"
#include "field.h"

/*! 
 * @brief found in most American habitat types.
 */
class JSolverAlpha {
 public:
  JSolverAlpha(double P0, double g0, double n1, double n2, double Ip,
               Grid *grid);
  ~JSolverAlpha();
  /*! which variables are 'in' and which are 'out' ? */
  void update(Field *jphi, Field *psi, Field *p, Field *g);

 private:
  /*! What the heck are all these variables???? -JAS */
  double P0_;
  double g0_;
  double n1_;
  double n2_;
  double Ip_;
  double *R_;
  double dr_;
  double dz_;
  double nr_;
  double nz_;
};

#endif  // J_SOLVER_ALPHA_H_
