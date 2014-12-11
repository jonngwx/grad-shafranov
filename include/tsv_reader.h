/**
 * @file tsv_data.h
 *
 * Contains declarations for the struct Table
 * and methods pertaining to it.
 * Also describes its sub-structs CoilData and PGData
 * and methods for creating and destroying them.
 *
 * @author Jacob Schwartz
 * @bug Should we change this from a struct to a class? Implement proper
 *constructors and destructors?
 */

#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <vector>
#include <string>

/**
 * A basic container for a 2D table or array of doubles.
 * data can be set and accessed like
 * data(row, col) */
class Table {
 protected:
  const static int kOpenFileError = 1;
  const static int kInconsistentColumnNumberError = 2;
  const static int kFileReadError     = 3;
  size_t num_columns_; /**< For keeping track of the array's size. */
  size_t num_rows_;    /**< For keeping track of the array's size. */
  std::vector< std::vector<double> > data_;
 public:
  /**
   * @brief Reads in a tsv file and creates a Table.
   * Detects the number of rows & columns in the file programmatically.
   * Behavior if some rows are shorter than others has not been tested.
   */
  virtual int load_from_tsv(const std::string filename, int header_lines = 0);
  int num_columns() const { return num_columns_; }
  int num_rows() const { return num_rows_; }
  double data(int row, int column) const { return data_[row][column]; }
};

class CoilData : public Table {
 private:
   const static int kNotThreeColumnsError = 4;
 public:
  int load_from_tsv(const std::string filename, int header_lines = 0) override;
  double r(int i) const { return data_[i][0]; }
  double z(int i) const { return data_[i][1]; }
  double current(int i) const { return data_[i][2]; }
};

class PGData : public Table {
 public:
  double psi(int i) const { return data_[i][0]; }
  double value(int i) const { return data_[i][1]; }
};

#endif
