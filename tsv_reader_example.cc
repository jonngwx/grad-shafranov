// a simple test of reading in some parameters to a text file.

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "tsv_reader.h"

int main (int argc, char * argv[]) {

  TsvData * td;
  std::string file_name(argv[1]);
  td = NewTsvDataFromFile(file_name);

  for (int i = 0; i < td->num_entries; ++i) {
    printf("R is %f and Z is %f and value is %f\n", td->r_locations[i], td->z_locations[i], td->values[i]);
  }
  
  return 0;
}
