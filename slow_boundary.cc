/*!
 * @file slow_boundary.cc
 * @author Peter J. Bolgert
 * @brief Implementation for the SlowBoundary class.
 */
#include <stdio.h>
#include <math.h>

#include "include/boundary.h"
#include "include/tsv_reader.h"
#include "include/grid.h"
#include "include/field.h"
#include "include/green_fcn.h"
#include "include/slow_boundary.h"

SlowBoundary::SlowBoundary(Field* psi, Grid* grid)
    : Boundary(psi, grid),
      R_(grid->R_),
      z_(grid->z_),
      dr_(grid->dr_),
      dz_(grid->dz_) {

  cond_data_ = new CoilData();
//  printf("num rows = %i \n", cond_data_->num_coil_subregions());

  // Initialize Green's Function Array for Plasma Response
  g_plasma_ = new double** [nr_];
  for (int i = 0; i < nr_; ++i) {
    g_plasma_[i] = new double* [nz_];
    for (int j = 0; j < nz_; ++j) {
      g_plasma_[i][j] = new double[perim_]();
      for (int l = 0; l < perim_; ++l) {
        g_plasma_[i][j][l] = green_fcn(R_[i], z_[j], R_[LtoI(l)], z_[LtoJ(l)]);
      }
    }
  }
}

SlowBoundary::SlowBoundary(Field* psi, Grid* grid, CoilData* cond_data)
    : SlowBoundary(psi, grid) {

  *cond_data_ = *cond_data;

  // Initialize Green's Fcn Array for Coils
  g_coils_ = new double* [cond_data_->num_coil_subregions()]; //num_coil_subregions = num of coils
  for (int c = 0; c < cond_data_->num_coil_subregions(); ++c) {
    g_coils_[c] = new double[perim_]();
    for (int l = 0; l < perim_; ++l) {
      g_coils_[c][l] = green_fcn(cond_data_->r(c),cond_data_->z(c), R_[LtoI(l)], z_[LtoJ(l)]);
    }
  }
  
}

SlowBoundary::~SlowBoundary() {
 
  // Delete g_plasma_
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0; j < nz_; ++j) {
      delete[] g_plasma_[i][j];
    }
    delete[] g_plasma_[i];
  }
  delete[] g_plasma_;

  // Delete g_coils_
  for (int c = 0; c < cond_data_->num_coil_subregions(); ++c) {
    delete[] g_coils_[c];
  }
  delete[] g_coils_;
  delete cond_data_;
}

int SlowBoundary::CalcB(Field* jphi) {
  double mu0 = 4 * M_PI * 1e-7; /* magnetic permeability of free space */
  // printf("perim_ is %d.\n",perim_);
  int n_rows_ = cond_data_->num_coil_subregions();

  for (int l = 0; l < perim_; ++l) {
    // printf("For l = %d, i = %d, j = %d \n", l, LtoI(l),LtoJ(l));
    psib_old[l] = psi_->f_[LtoI(l)][LtoJ(l)];
    psi_->f_[LtoI(l)][LtoJ(l)] = 0;
    for (int i = 0; i < nr_; ++i) {
      for (int j = 0; j < nz_; ++j) {
        psi_->f_[LtoI(l)][LtoJ(l)] +=
            mu0 * g_plasma_[i][j][l] * (jphi->f_[i][j]);
      }
    }
    psi_->f_[LtoI(l)][LtoJ(l)] *= (dr_ * dz_);
    
    // add psi at bdy due to coils
    for (int c = 0; c < n_rows_; ++c) {
      //printf("l = %i, c = %i \n", l, c);
      psi_->f_[LtoI(l)][LtoJ(l)] += mu0 * g_coils_[c][l] * cond_data_->current(c);
    }
  }
  return 0;
}
