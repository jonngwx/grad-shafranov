#include "include/rhs_func.h"
#include <string>
#include "include/tsv_reader.h"
using namespace std;

RHSfunc::RHSfunc(string pgtype, PGData *data)
  : pgtype_(pgtype),
    data_(data)
{}

RHSfunc::~RHSfunc()
{}

double RHSfunc::eval(double var) {

  return 3; //for now, lol 
}

double RHSfunc::eval_prime(double var) {

  return 4; //for now, lol
}
