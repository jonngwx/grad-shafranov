/*!
 * @file j_solver.h
 * @author Phillin LeBlanc
 * @brief Header declarations for JSolver.
 */


#ifndef J_SOLVER_H_
#define J_SOLVER_H_

#include "grid.h"
#include "field.h"

/*! 
 * @brief Interface for JSolver
 */
class JSolver {
 public:
    virtual ~JSolver(){};
    /*! which variables are 'in' and which are 'out' ? */
    virtual void update(Field *jphi, Field *psi, Field *p, Field *g)=0;
    
 protected:
  /*! What the heck are all these variables???? -JAS */
  double P0_; //!< Pressure at the magnetic axis
  double g0_; //!< Toroidal magnetic field * R at magnetic axis
  double Ip_; //!< Total plasma current (A)
  double *R_; //!< Array of radial grid points
  double dr_; //!< grid spacing in R
  double dz_; //!< grid spacing in z
  double nr_; //!< number of nodes in R
  double nz_; //!< number of nodes in z
};

#endif  // J_SOLVER_H_
