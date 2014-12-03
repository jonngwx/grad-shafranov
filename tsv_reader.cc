// a simple test of reading in some parameters to a text file.

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>

#include "include/tsv_reader.h"

TsvData * NewTsvDataFromFile(const std::string tsv_file_name){

  printf("Will try reading from: %s\n", tsv_file_name.c_str());

  std::ifstream filein(tsv_file_name.c_str());
  std::string line;

  if (!filein.is_open()) {
    perror("Error while opening file");
  }
  std::vector< std::vector<double> > tsv_table;
  while (std::getline(filein, line)) {
    std::istringstream ss(line);
    std::istream_iterator<std::string> begin(ss), end;

    //putting all the tokens in the vector  
    std::vector<std::string> arrayTokens(begin, end);

    std::vector<double> doubleVector(arrayTokens.size());
    std::transform(arrayTokens.begin(), arrayTokens.end(), doubleVector.begin(), [](std::string const& val){ return std::stod(val);});

    //Should we now delete the arrayTokens vector?

    tsv_table.push_back(doubleVector);

  }
  if (filein.bad()) {
    perror("Error while reading file");
  }
  
  // convert vector of vectors into struct of 3 arrays.
  TsvData * td = new TsvData();
  assert (td != NULL);
  
  int num_tsv_lines = tsv_table.size();
  td->r_locations  = new double[num_tsv_lines];
  assert (td->r_locations != NULL);
  td->z_locations  = new double[num_tsv_lines];
  assert (td->z_locations != NULL);
  td->values = new double[num_tsv_lines];
  td->num_entries = num_tsv_lines;
  
  //To access values
  for (std::vector<double>::size_type i = 0; i != tsv_table.size(); ++i) {
    td->r_locations[i] = tsv_table[i][0];
    td->z_locations[i] = tsv_table[i][1];
    td->values[i] = tsv_table[i][2];
  }
  
  return td;
}
