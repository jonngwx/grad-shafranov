// a simple test of reading in some parameters to a text file.

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main (int argc, char * argv[]) {

  std::string file_name(argv[1]);
  Table td;
  int success = td.load_from_tsv(file_name,1);
  if ( success != 0){
    std::cout << "load from tsv not successful. exiting.\n";
    exit(1);
  }

  for (int i = 0; i < td.num_rows(); ++i) {
    for(int j=0; j < td.num_columns(); ++j){
      std::cout << td.data(i,j) << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
