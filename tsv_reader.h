#ifndef TSV_READER_H_
#define TSV_READER_H_

#include <string>

struct TsvData {
  double * r_locations;
  double * z_locations;
  double * values;
  int num_entries;
};

typedef struct TsvData TsvData;

TsvData * NewTsvDataFromFile(const std::string tsv_file_name);

#endif
