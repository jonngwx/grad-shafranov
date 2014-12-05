#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include "include/field.h"

SlowBoundary::SlowBoundary(Grid &grid, CoilData &cond_data)
  : grid_(grid),
    cond_data_(cond_data) {

}

SlowBoundary::~SlowBoundary()
{}

int SlowBoundary::CalcB(Field &psi, Field &jphi) {

  // lots to add here.  

  return 0;
}
