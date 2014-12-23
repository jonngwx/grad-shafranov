#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include "include/field.h"
#include <stdio.h>
#include "include/green_fcn.h"

SlowBoundary::SlowBoundary(Grid* grid, CoilData* cond_data)
  : Boundary(grid),
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
  
  /*What is this magic number? Provide some documentation please.*/
  double mu0 = 0.0000012566370614;
//    printf("perim_ is %d.\n",perim_); 
    for (int l=0; l < perim_; ++l){
//      printf("For l = %d, i = %d, j = %d \n", l, LtoI(l),LtoJ(l));
      psi->f_[LtoI(l)][LtoJ(l)] = 0; 
      for (int i=0; i < nr_; ++i) {
          for (int j=0; j < nz_; ++j) {
              psi->f_[LtoI(l)][LtoJ(l)] += mu0*g_plasma_[i][j][l]*(jphi->f_[i][j]);
          }
      }
    psi->f_[LtoI(l)][LtoJ(l)] *= (dr_*dz_);
//    printf("boundary value at i = %i and j = %i is %f \n", LtoI(l), LtoJ(l), psi->f_[LtoI(l)][LtoJ(l)]); 
  } 

  return 0;

}

