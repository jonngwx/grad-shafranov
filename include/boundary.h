/*!
 * @file boundary.h
 * @author Peter J. Bolgert 
 * @brief Header file for the Boundary class.
 */
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "field.h"

/*!
 * @brief Contains a method to update the boundary values of Psi.
 *
 * This class is meant to be subclassed so that different types of boundary-value calculating methods can be implemented.
 */
class Boundary {
  public:
    Boundary(Grid *grid);
    virtual ~Boundary();

    /*!
     * @brief Updates boundary cells of psi based on the plasma current jphi.
     */
    virtual int CalcB(Field *psi, Field* jphi);
    /*! 
     * @brief Gets the horizontal index of the l'th boundary cell
     *
     * as counted ccw from the bottom left.
     * @param l Number of the boundary cell
     * @return The horizontal index of that cell. 
     */
    int LtoI(int l);

    /*! 
     * @brief Gets the vertical index of the l'th boundary cell
     * as counted ccw from the bottom left.
     * @param l Number of the boundary cell
     * @return The vertical index of that cell. 
     */
    int LtoJ(int l);
  protected:
    int nr_; //!< Number of grid cells in the r direction.
    int nz_; //!< Number of grid cells in the z direction.
};

#endif // BOUNDARY_H_
