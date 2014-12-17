#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "field.h"

class Boundary {
  public:
    Boundary(Grid *grid);
    virtual ~Boundary();
    virtual int CalcB(Field *psi, Field* jphi);
    int LtoI(int l);
    int LtoJ(int l);
  protected:
    int nr_;
    int nz_;
};

#endif // BOUNDARY_H_
