#ifndef RHSARRAY_H_
#define RHSARRAY_H_

class RhsArray {
  public:
    RhsArray(Grid &grid, double **p, double **g); 
    ~RhsArray();
    int calcRHS(Field &psi, Field &rhsField);
  private:
    /* to be determined */

};

#endif // RHSARRAY_H_

