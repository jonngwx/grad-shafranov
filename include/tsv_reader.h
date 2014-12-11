/**
 * @file tsv_data.h
 *
 * Contains declarations for the class Table
 * and methods pertaining to it.
 * Also describes its sub-classes CoilData and PGData.
 *
 * @author Jacob Schwartz
 */

#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <vector>
#include <string>

/*!
 * @brief A basic container for a 2D table or array of doubles.
 * data can be set and accessed like
 * data(row, col)
 */
class Table {
 protected:
  const static int kOpenFileError = 1;
  const static int kInconsistentColumnNumberError = 2;
  const static int kFileReadError = 3;
  size_t num_columns_; /**< For keeping track of the array's size. */
  size_t num_rows_;    /**< For keeping track of the array's size. */
  std::vector<std::vector<double> > data_;

 public:
  /*!
   * @brief Reads in a tsv file and creates a Table.
   * Detects the number of rows & columns in the file programmatically.
   * Behavior if some rows are shorter than others has not been tested.
   * @param [in] filename the filename, as a string
   * @param [in] header_lines how many inital lines to skip
   * @return 0 if successful, nonzero if error.
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
  /*!
   * @brief Reads in a tsv file and creates a CoilData,
   * Inherits from Table, then checks that there is the proper number of
   * columns.
   * @param [in] filename the filename, as a string
   * @param [in] header_lines how many inital lines to skip
   * @return 0 if successful, nonzero if error.
   */
  int load_from_tsv(const std::string filename, int header_lines = 0) override;
  double r(int i) const { return data_[i][0]; }
  double z(int i) const { return data_[i][1]; }
  double current(int i) const { return data_[i][2]; }
};

class PGData : public Table {
 private:
  const static int kNotTwoColumnsError = 4;

 public:
  /*!
  * @brief Reads in a tsv file and creates a PGData,
  * Inherits from Table, then checks that there is the proper number of columns.
  * @param [in] filename the filename, as a string
  * @param [in] header_lines how many inital lines to skip
  * @return 0 if successful, nonzero if error.
  */
  int load_from_tsv(const std::string filename, int header_lines = 0) override;
  double psi(int i) const { return data_[i][0]; }
  double value(int i) const { return data_[i][1]; }
};

#endif
