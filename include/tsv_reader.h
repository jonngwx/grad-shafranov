/**
 * @file tsv_data.h
 *
 * Contains declarations for the struct Table 
 * and methods pertaining to it.
 * Also describes its sub-structs CoilData and PGData
 * and methods for creating and destroying them.
 *
 * @author Jacob Schwartz
 * @bug Should we change this from a struct to a class? Implement proper constructors and destructors? 
 */

#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <string>

/** 
 * A basic container for a 2D table or array of doubles.
 * data can be set and accessed like
 * data[i][j].
 */
struct Table {
  double ** data; /**< The two-d array of actual data. */
  int num_columns; /**< For keeping track of the array's size. */
  int num_entries; /**< For keeping track of the array's size. */
};
typedef struct Table Table;

struct CoilData : Table {
  double * r_locs;
  double * z_locs;
  double * currents;
};
typedef struct CoilData CoilData;

struct PGData : Table {
  double * psis;
  double * values;
};
typedef struct PGData PGData;

/**
 * @brief Reads in a tsv file and creates a Table.
 * Detects the number of rows & columns in the file programmatically.
 * Behavior if some rows are shorter than others has not been tested.
 */
Table * NewTableFromFile(const std::string tsv_file_name, const int header_lines);
Table * NewTableFromFile(const std::string tsv_file_name);

/**
 * @brief Deletes a table.
 * May be memory leaks...
 */
void DeleteTable(Table * td);

/**
 * @brief Reads in a tsv file and creates a CoilData.
 * Calls NewTableFromFile and checks that there are 3 columns.
 * @see NewTableFromFile
 */
CoilData * NewCoilDataFromFile(const std::string tsv_file_name, const int header_lines);
CoilData * NewCoilDataFromFile(const std::string tsv_file_name);

/** @brief Deletes a CoilData
 *  Calls the superstruct's DeleteTable.
 *  @see DeleteTable.
 */ 
void DeleteCoilData(CoilData * cd);

/**
 * @brief Reads in a tsv file and creates a PGData.
 * Calls NewTableFromFile and checks that there are 3 columns.
 * @see NewTableFromFile
 */
PGData * NewPGDataFromFile(const std::string tsv_file_name, const int header_lines);
PGData * NewPGDataFromFile(const std::string tsv_file_name);

/** @brief Deletes a PGData
 *  Calls the superstruct's DeleteTable.
 *  @see DeleteTable.
 */ 
void DeletePGData(PGData * pgd);


#endif
