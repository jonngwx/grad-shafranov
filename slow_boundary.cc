#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include "include/field.h"

SlowBoundary::SlowBoundary(Grid* grid, CoilData* cond_data)
  : grid_(grid),
    cond_data_(cond_data) {
  
  int perim_ = 2*(grid_->nr_ + grid_->nz_ - 2);
  printf("%d\n", perim_);  
  // Initialize Green's Function Array
  g_plasma_ = new double**[grid_->nr_];
}

SlowBoundary::~SlowBoundary()
{}

int SlowBoundary::CalcB(Field &psi, Field &jphi) {

  // lots to add here.  

  return 0;
}
