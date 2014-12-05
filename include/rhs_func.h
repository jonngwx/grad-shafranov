#ifndef RHSFUNC_H_
#define RHSFUNC_H_

#include "tsv_reader.h"
#include <string>

using namespace std;

class RHSfunc {
  public:
    RHSfunc(string pgtype, PGData *data); 
    ~RHSfunc();
    double eval(double var);
    double eval_prime(double var);
  private:
    string pgtype_;
    PGData *data_;

};

#endif // RHSFUNC_H_

