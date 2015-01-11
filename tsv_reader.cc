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
  num_columns_ = 0;
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

  if(num_rows_ <= 0) { 
    return kNoDataFoundError;
  } else {
    return 0;
  }
}

int CoilData::load_from_tsv(const std::string tsv_file_name, int header_lines) {
  int status = Table::load_from_tsv(tsv_file_name, header_lines);
  num_coil_subregions_ = 0;
  if (status != 0) {
    return status;
  } else {
    if (num_columns_ == 3) { //simple coil_data.tsv format
      coil_data_ = data_;
      return 0;
    } else if (num_columns_ == 14) { //compressed coil_data.tsv format
      return GenerateCoilData();
    } else {
      std::cout << "Error: CoilData must have three columns for the simple format or fourteen columns for the 'compressed' format.\n";
      return kNotCorrectNumColumnsError;
    }
  }
}

int CoilData::GenerateCoilData(){
  std::vector<double> one_coil_subregion;
  double coil_dr, coil_dz; // The sub-coil-spacing in r and z.
  double temp_r, temp_z; // One sub-coil center location.
  num_coil_subregions_ = 0;
  for (size_t i = 0; i < num_rows_; i++){

    coil_dr = CoilRegionW(i)/CoilRegionNR(i);
    coil_dz = CoilRegionH(i)/CoilRegionNZ(i);

    for (int nr = 0; nr < CoilRegionNR(i); nr++) {
      temp_r = CoilRegionR(i) - CoilRegionW(i)/2 + (nr + 0.5) * coil_dr;
      for (int nz = 0; nz < CoilRegionNZ(i); nz++) {
        temp_z = CoilRegionZ(i) - CoilRegionH(i)/2 + (nz + 0.5) * coil_dz;
        one_coil_subregion = {temp_r, temp_z, CoilRegionCurrent(i)};
        coil_data_.push_back(one_coil_subregion);
        num_coil_subregions_++;
      }
    }
  }
  return 0; // May have nonzero returns in the future (to indicate an error) in case of bad input. 
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
