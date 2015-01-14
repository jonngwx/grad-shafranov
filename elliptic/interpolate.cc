/*
 * @file interpolate.cc
 * @brief Implementation for class Interpolate.
 */
#include "interpolate.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>

const int OutsideInterp = -2;
const int OutsideGrid = -1;

Interpolate::Interpolate(Grid &GridS, Field &F) :
F_(F),
Grid_(GridS),
dr_(Grid_.dr_),
dz_(Grid_.dz_),
P_(4, std::vector<double>(4)) {}

Interpolate::~Interpolate() {}

void Interpolate::CorrectCellBoundsCheck(double r, double z) const {
  if ((r - r_curr_) > dr_ || (z - z_curr_) > dz_ || (r - r_curr_) < 0 || (z - z_curr_) < 0) {
    /* 
    printf("dr_ = %f\n", dr_);
    printf("r_curr_ = %f\n", r_curr_);
    printf("r = %f\n", r);
    printf("dz_ = %f\n", dz_);
    printf("z_curr_ = %f\n", z_curr_);
    printf("z = %f\n", z);
    */
    throw OutsideInterp;
  }
}

double Interpolate::F(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_; //radial location relative to the cell that the target point is in: goes from 0 to 1.
  double r2 = rc*rc;
  double r3 = r2*rc;
  double zc = (z - z_curr_)/dz_; //vertical location relative to the cell that the target point is in: goes from 0 to 1.
  double z2 = zc*zc;
  double z3 = z2*zc;
  return (a00 + a10*rc + a20*r2 + a30*r3) +
         (a01 + a11*rc + a21*r2 + a31*r3)*zc +
         (a02 + a12*rc + a22*r2 + a32*r3)*z2 +
         (a03 + a13*rc + a23*r2 + a33*r3)*z3;
}

double Interpolate::F_r(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_;
  double r2 = rc*rc;
  double zc = (z - z_curr_)/dz_;
  double z2 = zc*zc;
  double z3 = z2*zc;
  return ((a10 + 2*a20*rc + 3*a30*r2) +
          (a11 + 2*a21*rc + 3*a31*r2)*zc +
          (a12 + 2*a22*rc + 3*a32*r2)*z2 +
          (a13 + 2*a23*rc + 3*a33*r2)*z3)/dr_;
}

double Interpolate::F_rr(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_;
  double zc = (z - z_curr_)/dz_;
  double z2 = zc*zc;
  double z3 = z2*zc;
  return ((2*a20 + 6*a30*rc) +
          (2*a21 + 6*a31*rc)*zc +
          (2*a22 + 6*a32*rc)*z2 +
          (2*a23 + 6*a33*rc)*z3)/(dr_*dr_);
}

double Interpolate::F_rz(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_;
  double r2 = rc*rc;
  double zc = (z - z_curr_)/dz_;
  double z2 = zc*zc;
  return ((a11 + 2*a21*rc + 3*a31*r2) +
        2*(a12 + 2*a22*rc + 3*a32*r2)*zc +
        3*(a13 + 2*a23*rc + 3*a33*r2)*z2)/(dr_*dz_);
}

double Interpolate::F_z(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_;
  double r2 = rc*rc;
  double r3 = r2*rc;
  double zc = (z - z_curr_)/dz_;
  double z2 = zc*zc;
  return ((a01 + a11*rc + a21*r2 + a31*r3) +
        2*(a02 + a12*rc + a22*r2 + a32*r3)*zc +
        3*(a03 + a13*rc + a23*r2 + a33*r3)*z2)/dz_;
}

double Interpolate::F_zz(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double rc = (r - r_curr_)/dr_;
  double r2 = rc*rc;
  double r3 = r2*rc;
  double zc = (z - z_curr_)/dz_;
  return (2*(a02 + a12*rc + a22*r2 + a32*r3) +
          6*(a03 + a13*rc + a23*r2 + a33*r3)*zc)/(dz_*dz_);
}

//These coefficients come from www.paulinternet.nl/?page=bicubic.
void Interpolate::updateCoefficients() {
  a00 = P_[1][1];
  a01 = -.5*P_[1][0] + .5*P_[1][2];
  a02 = P_[1][0] - 2.5*P_[1][1] + 2*P_[1][2] - .5*P_[1][3];
  a03 = -.5*P_[1][0] + 1.5*P_[1][1] - 1.5*P_[1][2] + .5*P_[1][3];
  a10 = -.5*P_[0][1] + .5*P_[2][1];
  a11 = .25*P_[0][0] - .25*P_[0][2] - .25*P_[2][0] + .25*P_[2][2];
  a12 = -.5*P_[0][0] + 1.25*P_[0][1] - P_[0][2] + .25*P_[0][3] + .5*P_[2][0] - 1.25*P_[2][1] + P_[2][2] - .25*P_[2][3];
  a13 = .25*P_[0][0] - .75*P_[0][1] + .75*P_[0][2] - .25*P_[0][3] - .25*P_[2][0] + .75*P_[2][1] - .75*P_[2][2] + .25*P_[2][3];
  a20 = P_[0][1] - 2.5*P_[1][1] + 2*P_[2][1] - .5*P_[3][1];
  a21 = -.5*P_[0][0] + .5*P_[0][2] + 1.25*P_[1][0] - 1.25*P_[1][2] - P_[2][0] + P_[2][2] + .25*P_[3][0] - .25*P_[3][2];
  a22 = P_[0][0] - 2.5*P_[0][1] + 2*P_[0][2] - .5*P_[0][3] - 2.5*P_[1][0] + 6.25*P_[1][1] - 5*P_[1][2] + 1.25*P_[1][3] + 2*P_[2][0] - 5*P_[2][1] + 4*P_[2][2] - P_[2][3] - .5*P_[3][0] + 1.25*P_[3][1] - P_[3][2] + .25*P_[3][3];
  a23 = -.5*P_[0][0] + 1.5*P_[0][1] - 1.5*P_[0][2] + .5*P_[0][3] + 1.25*P_[1][0] - 3.75*P_[1][1] + 3.75*P_[1][2] - 1.25*P_[1][3] - P_[2][0] + 3*P_[2][1] - 3*P_[2][2] + P_[2][3] + .25*P_[3][0] - .75*P_[3][1] + .75*P_[3][2] - .25*P_[3][3];
  a30 = -.5*P_[0][1] + 1.5*P_[1][1] - 1.5*P_[2][1] + .5*P_[3][1];
  a31 = .25*P_[0][0] - .25*P_[0][2] - .75*P_[1][0] + .75*P_[1][2] + .75*P_[2][0] - .75*P_[2][2] - .25*P_[3][0] + .25*P_[3][2];
  a32 = -.5*P_[0][0] + 1.25*P_[0][1] - P_[0][2] + .25*P_[0][3] + 1.5*P_[1][0] - 3.75*P_[1][1] + 3*P_[1][2] - .75*P_[1][3] - 1.5*P_[2][0] + 3.75*P_[2][1] - 3*P_[2][2] + .75*P_[2][3] + .5*P_[3][0] - 1.25*P_[3][1] + P_[3][2] - .25*P_[3][3];
  a33 = .25*P_[0][0] - .75*P_[0][1] + .75*P_[0][2] - .25*P_[0][3] - .75*P_[1][0] + 2.25*P_[1][1] - 2.25*P_[1][2] + .75*P_[1][3] + .75*P_[2][0] - 2.25*P_[2][1] + 2.25*P_[2][2] - .75*P_[2][3] - .25*P_[3][0] + .75*P_[3][1] - .75*P_[3][2] + .25*P_[3][3];
}

void Interpolate::updateP(double r, double z) {
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  int is = (int)(Grid_.celli(r));
  int js = (int)(Grid_.cellj(z));

  if (is-1 < 0 || is+2 >= nr || js-1 < 0 || js+2 >=nz) {
      throw OutsideGrid;
  }

  r_curr_ = (Grid_.R_[is]);
  z_curr_ = (Grid_.z_[js]);
  // Fill in P_
  for (int i = 0; i < 4 ; i++) {
    for (int j = 0; j < 4 ; j++) {
      P_[i][j] = F_.f_[is + i - 1][js + j - 1];
    }
  }
}

void Interpolate::updateInterpolation(double r, double z){
  updateP(r,z);
  updateCoefficients();
}

void Interpolate::PrintAmnCoefficients(){
  printf("\n%f    + %fz    + %fz^2    + %fz^3\n", a00, a01, a02, a03);
  printf(  "%fr   + %frz   + %frz^2   + %frz^3\n", a10, a11, a12, a13);
  printf(  "%fr^2 + %fr^2z + %fr^2z^2 + %fr^2z^3\n", a20, a21, a22, a23);
  printf(  "%fr^3 + %fr^3z + %fr^3z^2 + %fr^3z^3\n\n", a30, a31, a32, a33);
}
