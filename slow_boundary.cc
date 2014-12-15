#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include "include/field.h"
#include <stdio.h>
#include "include/green_fcn.h"

SlowBoundary::SlowBoundary(Grid* grid, CoilData* cond_data)
  : nr_(grid->nr_),
    nz_(grid->nz_),
    R_(grid->R_),
    z_(grid->z_),
    dr_(grid->dr_),
    dz_(grid->dz_),
    cond_data_(cond_data) {
  perim_ = 2*(nr_ + nz_ - 2);

  // Initialize Green's Function Array
  g_plasma_ = new double**[nr_];
  for (int i=0; i < nr_; ++i){
    g_plasma_[i] = new double*[nz_];
    for (int j=0; j < nz_; ++j){
      g_plasma_[i][j] = new double[perim_]();
      for (int l=0; l < perim_; ++l){
         g_plasma_[i][j][l] = green_fcn(R_[i],z_[j],R_[LtoI(l)],z_[LtoJ(l)]);
      }
    }
  }

}

SlowBoundary::~SlowBoundary()
{}

int SlowBoundary::CalcB(Field* psi, Field* jphi) {
  //printf("perim_ is %d.\n",perim_); 
    for (int l=0; l < perim_; ++l){
      printf("For l = %d, i = %d, j = %d \n", l, LtoI(l),LtoJ(l));
      psi->f_[LtoI(l)][LtoJ(l)] = 0; 
      for (int i=0; i < nr_; ++i) {
          for (int j=0; j < nz_; ++j) {
              psi->f_[LtoI(l)][LtoJ(l)] += g_plasma_[i][j][l]*(jphi->f_[i][j]);
          }
      }
    psi->f_[LtoI(l)][LtoJ(l)] *= (dr_*dz_); 
  } 

  return 0;

}



int SlowBoundary::LtoI(int l) { 
  int i; 
 
  if (l>=0 && l<=(nr_-2)) {
    i = l;
  }
  else if (l>=(nr_-1) && l<=(nr_+nz_-3)) {
    i = nr_ - 1;
  }
  else if (l>=(nr_+nz_-2) && l<=(2*nr_+nz_-4)) {
    i = (2*nr_+nz_-3) - l;
  }
  else {
    i = 0;
  }
  return i;
}

int SlowBoundary::LtoJ(int l) {
  int j;
 
  if (l>=0 && l<=(nr_-2)) {
    j = 0;
  } 
  else if (l>=(nr_-1) && l<=(nr_+nz_-3)) {
    j = l - (nr_-1);
  }
  else if (l>=(nr_+nz_-2) && l<=(2*nr_+nz_-4)) {
    j = nz_ - 1;
  }
  else {
    j = (2*nr_ + 2*nz_ - 4) - l;
  }
  return j;
}

