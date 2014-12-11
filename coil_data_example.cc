// a simple test of reading in some parameters to a text file.

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main (int argc, char * argv[]) {

  std::string file_name(argv[1]);
  CoilData cd;
  int status = cd.load_from_tsv(file_name,1);
  if ( status != 0){
    std::cout << "load from tsv not successful. exiting.\n";
    exit(1);
  }

  for (int i = 0; i < cd.num_rows(); ++i) {
    for(int j=0; j < cd.num_columns(); ++j){
      std::cout << cd.data(i,j) << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
