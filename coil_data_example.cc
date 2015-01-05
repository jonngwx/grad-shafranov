/*!
 * @file coil_data_example.cc
 * @author Jacob Schwartz
 * @brief A simple program showing how to use CoilData.
 *
 * Demonstrates how to create a CoilData, load data from a file,
 * and then extract the useful parameters.
 *
 * Not sure if we should move this and tsv_reader_example.cc to a folder called
 * examples/ or if we should get rid of them entirely.
 */
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main(int argc, char* argv[]) {

  // Load data into a CoilData from file.
  std::string file_name(argv[1]);
  CoilData cd;
  int status = cd.load_from_tsv(file_name, 1);
  if (status != 0) {
    std::cout << "load from tsv not successful. exiting.\n";
    exit(1);
  }

  // Print out the data.
  for (int i = 0; i < cd.num_rows(); ++i) {
    for (int j = 0; j < cd.num_columns(); ++j) {
      std::cout << cd.data(i, j) << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
