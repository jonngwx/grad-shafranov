/*!
 * @file j_solver.h
 * @author Phillin LeBlanc
 * @brief Header declarations for JSolverAlpha.
 */


#ifndef J_SOLVER_NSTX_H_
#define J_SOLVER_NSTX_H_

#include "grid.h"
#include "field.h"

/*! 
 * @brief Calculates current in NSTX. Assumes p' = Ax(1-x), gg' = c x 
 */
class J_Solver_NSTX : public J_Solver {
 public:
    J_Solver_NSTX(double P0, double g0, double Ip, Grid *grid);
    ~J_Solver_NSTX();
    /*! which variables are 'in' and which are 'out' ? */
    void update(Field *jphi, Field *psi, Field *p, Field *g);
    
 private:
    double P1_;
    
};

#endif  // J_SOLVER_NSTX_H_
