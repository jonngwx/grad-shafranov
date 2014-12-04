// a simple test of reading in some parameters to a text file.

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main (int argc, char * argv[]) {

  std::string file_name(argv[1]);
  TsvData * td = NewTsvDataFromFile(file_name);

  for (int i = 0; i < td->num_entries; ++i) {
    printf("R is %f and Z is %f and value is %f\n", td->data[0][i], td->data[1][i], td->data[2][i]);
  }
  
  return 0;
}
