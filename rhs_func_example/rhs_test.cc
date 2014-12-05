// a test of RHSfunc and other things

#include "tsv_reader.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "rhs_func.h"

using namespace std;

int main (int argc, char *argv[]){

  string file_name(argv[1]);  
  PGData * pgd = NewPGDataFromFile(file_name, 1);

  for (int i = 0; i < pgd->num_entries; ++i) {
    printf("Psi is %f and value is %f\n", pgd->psis[i], pgd->values[i]);
  }


  string pgtype = "type1";
  RHSfunc *p = new RHSfunc(pgtype, pgd);

  double value = p->eval(1.0);
  printf("The value of p(1.0) is %f\n",value);
  
  double prime = p->eval_prime(1.0);
  printf("The value of p'(1.0) is %f\n",prime);

  DeletePGData(pgd);

  return 0;
}
