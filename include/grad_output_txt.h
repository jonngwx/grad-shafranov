#ifndef GRAD_OUTPUT_TXT
#define GRAD_OUTPUT_TXT


/**
 * \file Header file for grad_output_txt
 * */
 
#include <stdlib.h>
#include "grad_output.h"
#include "field.h"

class Grad_Output_Txt : public Grad_Output{
public:
  /**
   * Constructor for output class
   * @param f pointer to field containing flux function
   * @param p pointer to field containing pressure function
   * @param g pointer to field containing g function
   * @param outputs string of comma separated output options
   * */
  Grad_Output_Txt(Field* f, Field* p, Field* g, const char* outputs);
  ~Grad_Output_Txt();

  /**
   * Writes output to file.
   * @param filename name of output file
   * */
  void write_output(const char* filename);
  
};

#endif
