/*!
 * @file j_solver.h
 * @author Jonathan Ng
 * @brief Header declarations for JSolver.
 */


#ifndef J_SOLVER_H_
#define J_SOLVER_H_

#include "grid.h"
#include "field.h"

/*! 
 * @brief Base class for methods which calculate jphi as a function of phi
 */
class JSolver {
 public:
    virtual ~JSolver(){};
    /**
     * @brief Updates the current density, pressure and toroidal field given the magnetic flux
     * @param jphi pointer to field containing current density
     * @param psi pointer to field containing flux
     * @param p pointer to field containing pressure
     * @param g pointer to field containing toroidal field function
     */
    virtual void update(Field *jphi, Field *psi, Field *p, Field *g)=0;
    
 protected:
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
