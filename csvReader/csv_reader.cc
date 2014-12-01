// a simple test of reading in some parameters to a text file.

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

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
  while (std::getline(filein, line)) {
    std::istringstream ss(line);
    std::istream_iterator<std::string> begin(ss), end;
    
    //putting all the tokens in the vector  
    std::vector<std::string> arrayTokens(begin, end);
    
    std::vector<double> doubleVector(arrayTokens.size());
    std::transform(arrayTokens.begin(), arrayTokens.end(), doubleVector.begin(), [](std::string const& val){ return std::stod(val);});
    
    //Should we now delete the arrayTokens vector?
    
    //print 'em all out.
    for ( int i = 0; i < doubleVector.size(); i++) {
      std::cout << doubleVector[i] << " ";
    }
    std::cout << "\n";
    
  }
  if (filein.bad()) {
    perror("Error while reading file");
  }
  return 0;
}
