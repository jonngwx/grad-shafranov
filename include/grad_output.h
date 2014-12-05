#ifndef GRAD_OUTPUT_H_
#define GRAD_OUTPUT_H_

#include <stdlib.h>
#include "field.h"
#include "grid.h"
#include "rhs_func.h"

class Grad_Output{
  public:
  virtual ~Grad_Output(){};
  virtual void write_output(const char* filename)=0;
  
protected:
  Field *f;
  RHSfunc *p;
  RHSfunc *g;
  Grid *grid;
  void parse_outputs(const char *outputs);
  // output_list
};

#endif
