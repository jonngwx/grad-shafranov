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
     * @param l Number of the boundary cell
     * @return The horizontal index of that cell. 
     * The index is counted ccw from the bottom left.
     * \verbatim
       L is numbered like:  So this function will return:
      
       6 5 4                0 1 2
       7   3                0   2
       0 1 2                0 1 2
       \endverbatim
     */
    int LtoI(int l) const;

    /*! 
     * @brief Gets the vertical index of the l'th boundary cell
     * @param l Number of the boundary cell
     * @return The vertical index of that cell. 
     * The index is counted ccw from the bottom left.
     * \verbatim
       L is numbered like:  So this function will return:
      
       6 5 4                2 2 2
       7   3                1   1
       0 1 2                0 0 0
       \endverbatim
     */
    int LtoJ(int l) const;
  protected:
    const int nr_; //!< Number of grid cells in the r direction.
    const int nz_; //!< Number of grid cells in the z direction.
};

#endif // BOUNDARY_H_
