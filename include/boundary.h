#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "field.h"

class Boundary {
  public:
    Boundary(Grid *grid);
    virtual ~Boundary();
    virtual int CalcB(Field *psi, Field* jphi);
    /*! What does this method do? What does it return? */
    /*! What are L? What are I and J? */
    int LtoI(int l);
    int LtoJ(int l);
  protected:
    /*! What are these things? */
    int nr_;
    int nz_;
};

#endif // BOUNDARY_H_
