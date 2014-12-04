// a simple test of reading in some parameters to a text file.

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "include/tsv_reader.h"

int main (int argc, char * argv[]) {

  std::string file_name(argv[1]);
  PGData * pgd = NewPGDataFromFile(file_name, 1);

  for (int i = 0; i < pgd->num_entries; ++i) {
    printf("Psi is %f and value is %f\n", pgd->psis[i], pgd->values[i]);
  }
  DeletePGData(pgd);  
  return 0;
}
