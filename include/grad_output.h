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
class GradOutput{
  public:
  virtual ~GradOutput(){};
  virtual void write_output(const char* filename)=0;
  
protected:
  Field *f; //!< pointer to field containing flux function
  Field *jphi; //!< pointer to field containing current density
  Field *p; //!< pointer to field containing pressure function
  Field *g; //!< pointer to field containing g function
  Grid *grid; //!< pointer to grid containing axis data
  /**
  * @brief Parses the string of field names to write to file.
  * Puts them in the member vector output_list.
  * @param[in] outputs A comma-separated string of outputs.
  */
  void parse_outputs(const char *outputs);
  // output_list
  enum Output {CURRENT, TOROIDAL_FIELD};
  std::vector< Output > output_list; //!< A list of the different fields to be output.

  /*!
   * @brief Determines whether a given Output is present in output_list.
   * @param[in] out The variety of output to test for.
   * @return If it is present in output_list, then true, else false. 
   */
  bool inline find(Output out){
      return std::find(output_list.begin(),output_list.end(),out)!=output_list.end();
  }

};

#endif
