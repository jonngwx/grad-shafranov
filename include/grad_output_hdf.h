/*!
 * @file grad_output_hdf.h
 * @author Jonathan Ng
 * @brief Header file for grad_output_hdf.
 */

#ifndef GRAD_OUTPUT_HDF
#define GRAD_OUTPUT_HDF

#include <assert.h>
#include <stdlib.h>
#include "grad_output.h"
#include "field.h"
#include <hdf5.h>

/*!
 * @brief Implementation of Grad_Output which writes data to a hdf5 file.
 */
class GradOutputHdf : public GradOutput {
 public:
  /**
   * Constructor for output class
   * @param f pointer to field containing flux function
   * @param jphi pointer to field containing current
   * @param p pointer to field containing pressure function
   * @param grid pointer to grid containing axis data
   * @param g pointer to field containing g function
   * @param outputs c string of comma separated output options
   * */
  GradOutputHdf(Field* f, Field* jphi, Grid* grid, Field* p, Field* g,
                  const char* outputs);
  ~GradOutputHdf();

  /**
   * Writes output to file.
   * @param filename name of output file
   * */
  void write_output(const char* filename);

 private:
  /**
   * @brief Converts an array of pointers to a 1d array
   * @param[in] f array of pointers to convert
   * @param[out] x 1d array to write to
   * @param[in] nx dimensions of 2d array in x
   * @param[in] ny dimensions of 2d array in y
   * */
  void twod_to_oned(const double* const* f, double* x, int nx, int ny);

  /**
   * @brief checks if hdf operations are successful
   * @param status return value of hdf function calls
   */
  inline void check(herr_t status) {
    assert(status >= 0);
  };
};

#endif
