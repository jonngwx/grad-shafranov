/*!
 * @file Header file for grad_output_txt
 * @author Mountain Lion
 */
#ifndef GRAD_OUTPUT_TXT
#define GRAD_OUTPUT_TXT
 
#include <stdlib.h>
#include "grad_output.h"
#include "field.h"

/*! 
 * @brief Implementation of Grad_Output which writes to a text file.
 */
class GradOutputTxt : public GradOutput{
public:
  /**
   * Constructor for output class
   * @param f pointer to field containing flux function
   * @param jphi pointer to field containing current
   * @param grid pointer to grid containing axis data
   * @param p pointer to field containing p
   * @param g pointer to field containing g
   * @param outputs string of comma separated output options
   * */
    GradOutputTxt(Field* f, Field* jphi, Grid* grid, Field* p, Field* g, const char* outputs);
  ~GradOutputTxt();

  /**
   * Writes output to file.
   * @param filename name of output file
   * */
  void write_output(const char* filename);
  
};

#endif
