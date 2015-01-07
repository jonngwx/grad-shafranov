/*!
 * @file field.h
 * @author Cabrerae
 * @brief Header file for Field class
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
   * @param grid The grid which this field lives in.
   */
  Field(const Grid &grid);
  ~Field();

  const Grid * grid_; //!< Pointer to grid data: stores the physical dimension information.
  double **f_; //!< Pointer to field data: a 2D array of data.
  double f_l; //!< Value of the field at plasma edge.
  double f_0;  //!< Value of the field at the magnetic axis.
};


#endif
