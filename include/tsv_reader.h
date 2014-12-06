#ifndef TSV_DATA_H_
#define TSV_DATA_H_

#include <string>

struct Table {
  double ** data;
  int num_columns;
  int num_entries;
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



Table * NewTableFromFile(const std::string tsv_file_name);
Table * NewTableFromFile(const std::string tsv_file_name, const int header_lines);
void DeleteTable(Table * td);

CoilData * NewCoilDataFromFile(const std::string tsv_file_name);
CoilData * NewCoilDataFromFile(const std::string tsv_file_name, const int header_lines);
void DeleteCoilData(CoilData * cd);

PGData * NewPGDataFromFile(const std::string tsv_file_name);
PGData * NewPGDataFromFile(const std::string tsv_file_name, const int header_lines);
void DeletePGData(PGData * pgd);


#endif
