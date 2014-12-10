#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include "include/field.h"
#include <stdio.h>

SlowBoundary::SlowBoundary(Grid* grid, CoilData* cond_data)
  : nr_(grid->nr_),
    nz_(grid->nz_),
    cond_data_(cond_data) {
  perim_ = 2*(nr_ + nz_ - 2);

  // Initialize Green's Function Array
  g_plasma_ = new double**[nr_];
  for (int i=0; i < nr_; ++i){
    g_plasma_[i] = new double*[nz_];
    for (int j=0; j < nz_; ++j){
      g_plasma_[i][j] = new double[perim_]();
      for (int l=0; l < perim_; ++l){
        g_plasma_[i][j][l] = 0.2;
      }
    }
  }

}

SlowBoundary::~SlowBoundary()
{}

int SlowBoundary::CalcB(Field* psi, Field* jphi) {
  int IJ[2];
  printf("perim_ is %d.\n",perim_); 
    for (int l=0; l < perim_; ++l){
    LtoIJ(IJ, l);
    printf("For l = %d, i = %d, j = %d \n", l, IJ[0],IJ[1]);
    psi->f_[IJ[0]][IJ[1]] =2000; 
  } 

  return 0;

}



void SlowBoundary::LtoIJ(int ar[],int l) { 
  int i,j; 

  if (l>=0 && l<=(nr_-2)) {
    i = l;
    j = 0;
  }
  else if (l>=(nr_-1) && l<=(nr_+nz_-3)) {
    i = nr_ - 1;
    j = l - (nr_-1);
  }
  else if (l>=(nr_+nz_-2) && l<=(2*nr_+nz_-4)) {
    i = (2*nr_+nz_-3) - l;
    j = nz_ - 1;
  }
  else {
    i = 0;
    j = (2*nr_+2*nz_-4)  - l;
  }
  
  ar[0] = i;
  ar[1] = j;
}

