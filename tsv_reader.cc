/*!
 * @file tsv_reader.cc
 * @brief implementation for Table and CoilData and PGData classes.
 * @author Jacob Schwartz
 */
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "include/tsv_reader.h"

int Table::load_from_tsv(const std::string tsv_file_name, int header_lines) {

  std::ifstream filein(tsv_file_name.c_str());
  std::string line;

  if (!filein.is_open()) {
    perror("Error while opening file");
    return kOpenFileError;
  }

  int lines_read_in = 0;
  // skip 'header_lines' lines
  for (int i = 0; i < header_lines; ++i) {
    std::getline(filein, line);
    lines_read_in++;
  }
  std::vector<double> doubleVector;
  while (std::getline(filein, line)) {
    std::istringstream ss(line);

    /*if the first character is a # ignore the line and read in another. */
    if(ss.peek() == '#'){
      continue;
    } else {
      lines_read_in++;
    }

    std::istream_iterator<std::string> begin(ss), end;
    // putting all the tokens in the vector
    std::vector<std::string> arrayTokens(begin, end);

    if (lines_read_in == header_lines + 1) {
      num_columns_ = arrayTokens.size();
    } else if (arrayTokens.size() != num_columns_) {
      std::cout << "Inconsistent number of columns in line " << lines_read_in
                << "\n";
      return kInconsistentColumnNumberError;
    }
    doubleVector.resize(arrayTokens.size());

    // convert the vector of strings to a vector of doubles
    std::transform(arrayTokens.begin(), arrayTokens.end(), doubleVector.begin(),
                   [](std::string const& val) { return std::stod(val); });

    data_.push_back(doubleVector);
  }

  if (filein.bad()) {
    perror("Error while reading file");
    return kFileReadError;
  }

  num_rows_ = data_.size();

  return 0;
}

int CoilData::load_from_tsv(const std::string tsv_file_name, int header_lines) {
  int status = Table::load_from_tsv(tsv_file_name, header_lines);
  if (status != 0) {
    return status;
  } else {
    if (num_columns_ != 3) {
      std::cout << "Error: CoilData must have three columns.\n";
      return kNotThreeColumnsError;
    } else {
      return 0;
    }
  }
}

int PGData::load_from_tsv(const std::string tsv_file_name, int header_lines) {
  int status = Table::load_from_tsv(tsv_file_name, header_lines);
  if (status != 0) {
    return status;
  } else {
    if (num_columns_ != 2) {
      std::cout << "Error: PGData must have three columns.\n";
      return kNotTwoColumnsError;
    } else {
      return 0;
    }
  }
}
