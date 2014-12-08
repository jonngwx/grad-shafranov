#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "field.h"

class Boundary {
  public:
    virtual ~Boundary() {}
    virtual int CalcB(Field &psi, Field &jphi) = 0;
};

#endif // BOUNDARY_H_
