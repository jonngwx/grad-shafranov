#ifndef GRAD_OUTPUT_H_
#define GRAD_OUTPUT_H_

#include <stdlib.h>
#include <vector>
#include "field.h"
#include "grid.h"
#include <algorithm>

class Grad_Output{
  public:
  virtual ~Grad_Output(){};
  virtual void write_output(const char* filename)=0;
  
protected:
  Field *f;
  Field *jphi;
  Field *p;
  Field *g;
  Grid *grid;
  /**
  * @brief Parses the string of outputs to write to file..
  * @param outputs comma separated string of outputs. Currently doesn't do anything.
  * 
  */
  void parse_outputs(const char *outputs);
  // output_list
  enum Output {CURRENT, TOROIDAL_FIELD};
  std::vector< Output > output_list;
  bool inline find(Output out){
      return std::find(output_list.begin(),output_list.end(),out)!=output_list.end();
  }

};

#endif
