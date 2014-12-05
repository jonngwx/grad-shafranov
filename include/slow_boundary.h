#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "field.h"
#include "grid.h"
#include "tsv_reader.h"
#include "boundary.h"

class SlowBoundary : public Boundary {
  public:
    SlowBoundary(Grid &grid, CoilData &cond_data); 
    ~SlowBoundary();
    int CalcB(Field &psi, Field &jphi);
  private:
    Grid &grid_;
    CoilData &cond_data_;

};

#endif // SLOWBOUNDARY_H_

