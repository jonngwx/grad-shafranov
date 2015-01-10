/*!
 * @file grid.h
 * @author Felis concolor
 * @brief Header declarations for Grid.
 */

#ifndef GRID_H
#define GRID_H
#include <stdlib.h>

/*!
 * @brief Stores information about the solution grid.
 *
 * Stores the x and y axes of grid used in the solver and
 * methods to determine cell indices from coordinates.
 */
class Grid {
 public:
  /*!
   * Constructor of grid
   * @param R0 value of R at left boundary
   * @param Rend value of R at right boundary
   * @param z0 value of z at lower boundary
   * @param zend value of z at upper boundary
   * @param nr number of cells in the R direction
   * @param nz number of cells in the z direction
   */
  Grid(double R0, double Rend, double z0, double zend, int nr, int nz);
  ~Grid();
  
  /*!
   * @brief Converts radial position in meters to i-coordinate in grid-cell-space
   * which ranges from 0. to (nr_ - 1.0).
   * @param[in] r A radial position in meters.
   * @return The (floating-point-) number of cells in the radial direction 
   * out from the 0th cell that this location corresponds to. The returned value is clamped between 0.0 and (nr_ - 1.0).
   */
  double celli(double r) const;

  /*!
   * @brief Converts z position in meters to j-coordinate in grid-cell-space
   * which ranges from 0. to (nz_ - 1.0).
   * @param[in] z A z position in meters.
   * @return The (floating-point-) number of cells in the vertical direction 
   * up from the 0th cell that this location corresponds to. The returned value is clamped between 0.0 and (nz_ - 1.0).
   */
  double cellj(double z) const;

  const int nr_; //!< Number of cells in R direction.
  const int nz_; //!< Number of cells in z direction.
  double *R_;    //!< Array of grid cell radial locations in meters.
  double *z_;    //!< Array of grid cell vertical locations in meters.
  double dr_;    //!< Size of a grid cell in the r direction.
  double dz_;    //!< Size of a grid cell in the z direction.
};

#endif
