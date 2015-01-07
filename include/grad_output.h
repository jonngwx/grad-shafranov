/*!
 * @file grad_output.h
 * @author Puma concolor
 * @brief Declarations for class Grad_Output
 */

#ifndef GRAD_OUTPUT_H_
#define GRAD_OUTPUT_H_

#include <stdlib.h>
#include <vector>
#include "field.h"
#include "grid.h"
#include <algorithm>

/*! 
 * @brief Interface for writing the output of the solver to file. 
 */
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
  * @brief Parses the string of field names to write to file.
  * Puts them in the member vector output_list.
  * @param[in] outputs A comma-separated string of outputs.
  */
  void parse_outputs(const char *outputs);
  // output_list
  enum Output {CURRENT, TOROIDAL_FIELD};
  std::vector< Output > output_list; /*!< A list of the different fields to be output. */
  bool inline find(Output out){
      return std::find(output_list.begin(),output_list.end(),out)!=output_list.end();
  }

};

#endif
