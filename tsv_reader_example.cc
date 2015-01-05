/*!
 * @file tsv_reader_example.cc
 * @author Jacob Schwartz
 * @brief This is a short file that shows how to use Table:
 *
 * like how to create a Table and then load it up with data from a file
 * and then access the relevant variables from Table: the number of columns, the
 *number of rows, and the data itself.
 *
 * I'm not sure if we should put this and coil_data_example in an examples/
 *folder or get rid of them entirely.
 */
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main(int argc, char* argv[]) {

  // Create a table and load it with data
  std::string file_name(argv[1]);
  Table td;
  int success = td.load_from_tsv(file_name, 1);
  if (success != 0) {
    std::cout << "load from tsv not successful. exiting.\n";
    exit(1);
  }

  // Print out the data
  for (int i = 0; i < td.num_rows(); ++i) {
    for (int j = 0; j < td.num_columns(); ++j) {
      std::cout << td.data(i, j) << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
