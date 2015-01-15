/*!
 * @file boundary.cc
 * @author Peter J. Bolgert 
 * @brief Implementation for the Boundary class.
 */
#include "include/boundary.h"
#include <math.h>
#include <stdio.h>
#include <grid.h>
#include "include/field.h"

Boundary::Boundary(Field *psi, Grid* grid) : psi_(psi), nr_(grid->nr_), nz_(grid->nz_) {
  perim_ = 2*(nr_ + nz_ - 2);
  psib_old = new double[perim_];
}

Boundary::~Boundary() { delete [] psib_old;}

int Boundary::CalcB(Field* jphi) { return 0; }

int Boundary::LtoI(int l) const {
  int i = 0;

  if (l >= 0 && l <= (nr_ - 2)) {
    i = l;
  } else if (l >= (nr_ - 1) && l <= (nr_ + nz_ - 3)) {
    i = nr_ - 1;
  } else if (l >= (nr_ + nz_ - 2) && l <= (2 * nr_ + nz_ - 4)) {
    i = (2 * nr_ + nz_ - 3) - l;
  } else {
    i = 0;
  }
  return i;
}

int Boundary::LtoJ(int l) const {
  int j = 0;

  if (l >= 0 && l <= (nr_ - 2)) {
    j = 0;
  } else if (l >= (nr_ - 1) && l <= (nr_ + nz_ - 3)) {
    j = l - (nr_ - 1);
  } else if (l >= (nr_ + nz_ - 2) && l <= (2 * nr_ + nz_ - 4)) {
    j = nz_ - 1;
  } else {
    j = (2 * nr_ + 2 * nz_ - 4) - l;
  }
  return j;
}

double Boundary::norm() {
  double sum = 0;
  
  for (int l = 0; l<perim_; ++l) {
    sum += pow((psi_->f_[LtoI(l)][LtoJ(l)]-psib_old[l]),2);
  }
  
  return sqrt(sum);


}
