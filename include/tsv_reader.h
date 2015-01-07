/*!
 * @file tsv_reader.h
 * @author Jacob Schwartz
 * @brief Header declarations for Table, CoilData, and PGData.
 */
#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <vector>
#include <string>

/*!
 * @brief A basic container for a 2D table or array of doubles.
 *
 * data can be accessed like
 * data(row, col) and cannot currently be set.
 * There is also no standard constructor right now.
 * To use, do Table tb; tb.load_from_tsv(<filename>);
 */
class Table {
 protected:
  const static int kOpenFileError = 1;
  const static int kInconsistentColumnNumberError = 2;
  const static int kFileReadError = 3;
  const static int kNoDataFoundError = 4;
  size_t num_rows_;    //!< Number of rows of data.
  size_t num_columns_; //!< Number of columns of data.
  std::vector<std::vector<double> > data_; //!< A 2D-vector in which the data is stored, row-major.

 public:
  /*! 
   * @brief Constructor for Table.
   */
  Table() : num_columns_(0), num_rows_(0) {};
  ~Table() {};
  /*!
   * @brief Reads in a tsv file and populates a Table.
   *
   * Detects the number of rows & columns in the file programmatically.
   * Checks that there are the same number of columns for each row or stops and
   * returns an error.
   *
   * Skips header_lines lines to start. Then, while reading the body of the tsv,
   * will skip any line for which the first character is a #.
   * Inline or end-of-line comments are not allowed.
   *
   * @param [in] filename the filename, as a string
   * @param [in] header_lines how many inital lines to skip
   * @return 0 if successful, nonzero if error.
   */
  virtual int load_from_tsv(const std::string filename, int header_lines = 0);
  int num_columns() const { return num_columns_; }
  int num_rows() const { return num_rows_; }

  /*!
   * @brief get data from the table
   * @param [in] row integer of row (starting from 0)
   * @param [in] column integer of column (starting from 0)
   * @return the data point (which is a double)
   *
   * Does not do any bounds checking.
   */
  double data(int row, int column) const { return data_[row][column]; }
};

/*!
 * @brief A container for the r and z location (m) and currents (A) of external
 *coils.
 *
 * Load from a (3-column) tsv of doubles like Table does. Can access the data in
 * the same way or use getter functions r(i) z(i) and current(i)
 */
class CoilData : public Table {
 private:
  const static int kNotThreeColumnsError = 5;

 public:
  /*!
   * @brief Reads in a tsv file and populates a CoilData.
   *
   * Inherits from Table, then checks that there is the proper number of
   * columns.
   * @param [in] filename the filename, as a string
   * @param [in] header_lines how many inital lines to skip
   * @return 0 if successful, nonzero if error.
   */
  int load_from_tsv(const std::string filename, int header_lines = 0) override;

  /*!
   * @brief Getter for r location (meters) of i'th coil
   * @param [in] i the number of the coil
   * @return r location of the coil from the axis (meters)
   */
  double r(int i) const { return data_[i][0]; }

  /*!
   * @brief Getter for z location (meters) of i'th coil
   * @param [in] i the number of the coil
   * @return z location of the coil above the midplane (meters)
  */
  double z(int i) const { return data_[i][1]; }

  /*!
   * @brief Getter for current (Amps) of i'th coil
   * @param [in] i the number of the coil
   * @return current in amps of the coil
   */
  double current(int i) const { return data_[i][2]; }
};

/*!
 * @brief (THIS CLASS IS NOT USED)  A container for one sort of data describing initial P and G
 *
 * Useful for the 'array' type of P and G input.
 * Contains either P or G as a function of psi.
 * Reads from a 2-column text file: each row has a value of psi
 * and a value of p-or-g.
 *
 */
class PGData : public Table {
 private:
  const static int kNotTwoColumnsError = 4;

 public:
  /*!
  * @brief Reads in a tsv file and populates a PGData.
  *
  * Inherits from Table, then checks that there is the proper number of columns.
  * @param [in] filename the filename, as a string
  * @param [in] header_lines how many inital lines to skip
  * @return 0 if successful, nonzero if error.
  */
  int load_from_tsv(const std::string filename, int header_lines = 0) override;

  /*!
   * @brief Get the i'th value of psi
   * @param [in] i which value of psi (starting from 0)
   * @return the i'th value of psi
   */
  double psi(int i) const { return data_[i][0]; }

  /*!
   * @brief Get the i'th value of (p or g)
   * @param [in] i which value of psi (starting from 0)
   * @return the i'th value of (p or g)
   *
   * Since this PGTable stores either p or g it can only return the
   * corresponding value of course...
   */
  double value(int i) const { return data_[i][1]; }
};

#endif
