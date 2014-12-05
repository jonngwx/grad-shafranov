#ifndef FIELD_H
#define FIELD_H

#include <stdlib.h>

/**
 * \file Header file for Field class
 * */

class Field {
  public:
  Field(int nr, int nz);
  ~Field();

  const int nr_; /** < number of points in R direction */
  const int nz_; /** < number of points in z direction */
  double **f_; /** < pointer to field data */
  /**
   * Returns a the value of the ith boundary cell
   * TODO: do we really need this?
   * 
   * @param i index on the boundary?
   * */
  double f_boundary(int i);
};


#endif
