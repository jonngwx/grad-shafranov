/*!
 * @file j_solver.h
 * @author Phillin LeBlanc
 * @brief Header declarations for JSolverAlpha.
 */


#ifndef J_SOLVER_H_
#define J_SOLVER_H_

#include "grid.h"
#include "field.h"

/*! 
 * @brief Interface for J_solver!
 */
class J_Solver {
 public:
    virtual ~J_Solver(){};
    /*! which variables are 'in' and which are 'out' ? */
    virtual void update(Field *jphi, Field *psi, Field *p, Field *g)=0;
    
 protected:
  /*! What the heck are all these variables???? -JAS */
  double P0_;
  double g0_;
  double Ip_;
  double *R_;
  double dr_;
  double dz_;
  double nr_;
  double nz_;
};

#endif  // J_SOLVER_H_
