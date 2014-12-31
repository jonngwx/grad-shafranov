/*!
 * @file Header file for Field class
 */

#ifndef FIELD_H
#define FIELD_H

#include <stdlib.h>
#include "grid.h"

/*! 
 * @brief Container for 2d data and grid used in the solver.
 */
class Field {
  public:
  /*!
   * @brief Constructor of Field
   * @param grid The grid which this field lives in
   */
  Field(const Grid &grid);
  ~Field();

  const Grid * grid_; /** < pointer to grid data */
  double **f_; /** < pointer to field data */
  double f_l; /** < value of f at plasma edge */
  double f_0;  /** < value of f at mag axis */
};


#endif
