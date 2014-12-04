// a simple test of reading in some parameters to a text file.

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main (int argc, char * argv[]) {

  std::string file_name(argv[1]);
  CoilData * cd = NewCoilDataFromFile(file_name, 1);

  for (int i = 0; i < cd->num_entries; ++i) {
    printf("R is %f and Z is %f and I is %f\n", cd->r_locs[i], cd->z_locs[i], cd->currents[i]);
  }
  DeleteCoilData(cd);  
  return 0;
}
