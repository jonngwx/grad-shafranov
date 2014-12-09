#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "field.h"
#include "grid.h"
#include "tsv_reader.h"
#include "boundary.h"

class SlowBoundary : public Boundary {
  public:
    SlowBoundary(Grid* grid, CoilData* cond_data); 
    ~SlowBoundary();
    int CalcB(Field* psi, Field* jphi);
  private:
    int nr_;
    int nz_;
    CoilData* cond_data_;
    int perim_;
    double ***g_plasma_;
    void LtoIJ(int ar[], int l);
};

#endif // SLOWBOUNDARY_H_

