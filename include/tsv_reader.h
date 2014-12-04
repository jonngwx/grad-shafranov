#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <string>

struct TsvData {
  double ** data;
  int num_columns;
  int num_entries;
};

typedef struct TsvData TsvData;

struct CoilData : TsvData {
  double * r_locs;
  double * z_locs;
  double * currents;
};

typedef struct CoilData CoilData;

TsvData * NewTsvDataFromFile(const std::string tsv_file_name);
TsvData * NewTsvDataFromFile(const std::string tsv_file_name, const int header_lines);
void DeleteTsvData(TsvData * td);

CoilData * NewCoilDataFromFile(const std::string tsv_file_name);
CoilData * NewCoilDataFromFile(const std::string tsv_file_name, const int header_lines);
void DeleteCoilData(CoilData * cd);

#endif
