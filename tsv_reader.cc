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

TsvData * NewTsvDataFromFile(const std::string tsv_file_name) {
  printf("Not skipping any header lines\n");
  return NewTsvDataFromFile(tsv_file_name, 0);
}

TsvData * NewTsvDataFromFile(const std::string tsv_file_name, const int header_lines) {

  printf("Will try reading from: %s\n", tsv_file_name.c_str());

  std::ifstream filein(tsv_file_name.c_str());
  std::string line;

  if (!filein.is_open()) {
    perror("Error while opening file");
  }
  std::vector< std::vector<double> > tsv_table;
  
  //skip 'header_lines' lines
  for (int i = 0; i < header_lines; ++i) {
    std::getline(filein,line);
  }
  std::vector<double> doubleVector; 
  while (std::getline(filein, line)) {
    std::istringstream ss(line);
    std::istream_iterator<std::string> begin(ss), end;

    //putting all the tokens in the vector  
    std::vector<std::string> arrayTokens(begin, end);

    doubleVector.resize(arrayTokens.size());
    
    //convert the vector of strings to a vector of doubles
    std::transform(arrayTokens.begin(), arrayTokens.end(), doubleVector.begin(), [](std::string const& val){ return std::stod(val);});

    //Should we now delete the arrayTokens vector?

    tsv_table.push_back(doubleVector);

  }
  if (filein.bad()) {
    perror("Error while reading file");
  }



  int num_tsv_lines = tsv_table.size();
  int num_cols = doubleVector.size(); 
  
  TsvData * td = new TsvData();
  assert (td != NULL);
  td->num_entries = num_tsv_lines;
  td->num_columns = num_cols;
  
  td->data = new double *[num_cols];
  for(int i = 0; i < num_cols; ++i) {
    td->data[i] = new double[num_tsv_lines];
    assert(td->data[i] != NULL);
  }
  
  for (int i = 0; i < num_cols; ++i) {
    for (int j = 0; j < num_tsv_lines; ++j){
      td->data[i][j] = tsv_table[j][i];
    }
  }
  
  return td;
}

void DeleteTsvData(TsvData * td){
  assert(td != NULL);
  printf("Deleting a tsv struct\n");
  for(int i = 0; i < td->num_columns; ++i){
    assert(td->data[i] != NULL);
    delete td->data[i];
  }
  delete td;
}

CoilData *  NewCoilDataFromFile(const std::string tsv_file_name){
  return NewCoilDataFromFile(tsv_file_name, 0);
}

CoilData * NewCoilDataFromFile(const std::string tsv_file_name, const int header_lines){
  TsvData * td = NewTsvDataFromFile(tsv_file_name, header_lines);
  assert(td->num_columns==3);
  CoilData * cd = new CoilData();
  cd->data=td->data;
  cd->num_columns=td->num_columns;
  cd->num_entries=td->num_entries;
  cd->r_locs = cd->data[0];
  cd->z_locs = cd->data[1];
  cd->currents = cd->data[2];

  return cd;
}

void DeleteCoilData(CoilData * cd){
  DeleteTsvData(cd);
}
