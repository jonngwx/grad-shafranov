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

struct TsvData {
  double * r_locations;
  double * z_locations;
  double * currents;
  int num_entries;
};

typedef struct TsvData TsvData;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    printf("USAGE: %s <name of file to read>\n", argv[1]);
    exit(1);
  }

  const std::string input_file_name = argv[1];
  printf("Will try reading from: %s\n", input_file_name.c_str());

  std::ifstream filein(input_file_name.c_str());
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
  td->currents = new double[num_tsv_lines];
  td->num_entries = num_tsv_lines;
  
  	//To access values
  for (std::vector<double>::size_type i = 0; i != tsv_table.size(); ++i) {
    td->r_locations[i] = tsv_table[i][0];
    td->z_locations[i] = tsv_table[i][1];
    td->currents[i] = tsv_table[i][2];
  }
  
  for (int i = 0; i < td->num_entries; ++i) {
    printf("R is %f and Z is %f and I is %f\n", td->r_locations[i], td->z_locations[i], td->currents[i]);
  }
  
  return 0;
}
