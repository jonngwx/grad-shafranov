#include "interpolate.h"
#include <math.h>
#include "field.h"
#include "grid.h"
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>

const int OutsideInterp = -2;
const int OutsideGrid = -1;

Interpolate::Interpolate(Grid &GridS, Field &Psi) :
Psi_(Psi),
Grid_(GridS),
dr_(Grid_.dr_),
dz_(Grid_.dz_),
P(4, std::vector<double>(4)) {}

Interpolate::~Interpolate() {}

void Interpolate::CorrectCellBoundsCheck(double r, double z) const {
  if ((r - r_curr) > dr_ || (z - z_curr) > dz_ || (r - r_curr) < 0 || (z - z_curr) < 0) {
    /* 
    printf("dr_ = %f\n", dr_);
    printf("r_curr = %f\n", r_curr);
    printf("r = %f\n", r);
    printf("dz_ = %f\n", dz_);
    printf("z_curr = %f\n", z_curr);
    printf("z = %f\n", z);
    */
    throw OutsideInterp;
  }
}

double Interpolate::Psi_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double r2 = r_*r_;
  double r3 = r2*r_;
  double z_ = (z - z_curr)/dz_;
  double z2 = z_*z_;
  double z3 = z2*z_;
  return (a00 + a10*r_ + a20*r2 + a30*r3) +
         (a01 + a11*r_ + a21*r2 + a31*r3)*z_ +
         (a02 + a12*r_ + a22*r2 + a32*r3)*z2 +
         (a03 + a13*r_ + a23*r2 + a33*r3)*z3;
}

double Interpolate::Psir_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double r2 = r_*r_;
  double z_ = (z - z_curr)/dz_;
  double z2 = z_*z_;
  double z3 = z2*z_;
  return ((a10 + 2*a20*r_ + 3*a30*r2) +
          (a11 + 2*a21*r_ + 3*a31*r2)*z_ +
          (a12 + 2*a22*r_ + 3*a32*r2)*z2 +
          (a13 + 2*a23*r_ + 3*a33*r2)*z3)/dr_;
}

double Interpolate::Psirr_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double z_ = (z - z_curr)/dz_;
  double z2 = z_*z_;
  double z3 = z2*z_;
  return ((2*a20 + 6*a30*r_) +
          (2*a21 + 6*a31*r_)*z_ +
          (2*a22 + 6*a32*r_)*z2 +
          (2*a23 + 6*a33*r_)*z3)/(dr_*dr_);
}

double Interpolate::Psirz_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double r2 = r_*r_;
  double z_ = (z - z_curr)/dz_;
  double z2 = z_*z_;
  return ((a11 + 2*a21*r_ + 3*a31*r2) +
        2*(a12 + 2*a22*r_ + 3*a32*r2)*z_ +
        3*(a13 + 2*a23*r_ + 3*a33*r2)*z2)/(dr_*dz_);
}

double Interpolate::Psiz_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double r2 = r_*r_;
  double r3 = r2*r_;
  double z_ = (z - z_curr)/dz_;
  double z2 = z_*z_;
  return ((a01 + a11*r_ + a21*r2 + a31*r3) +
        2*(a02 + a12*r_ + a22*r2 + a32*r3)*z_ +
        3*(a03 + a13*r_ + a23*r2 + a33*r3)*z2)/dz_;
}

double Interpolate::Psizz_interp(double r, double z) const {
  CorrectCellBoundsCheck(r,z);
  double r_ = (r - r_curr)/dr_;
  double r2 = r_*r_;
  double r3 = r2*r_;
  double z_ = (z - z_curr)/dz_;
  return (2*(a02 + a12*r_ + a22*r2 + a32*r3) +
          6*(a03 + a13*r_ + a23*r2 + a33*r3)*z_)/(dz_*dz_);
}

void Interpolate::updateCoefficients() {
  a00 = P[1][1];
  a01 = -.5*P[1][0] + .5*P[1][2];
  a02 = P[1][0] - 2.5*P[1][1] + 2*P[1][2] - .5*P[1][3];
  a03 = -.5*P[1][0] + 1.5*P[1][1] - 1.5*P[1][2] + .5*P[1][3];
  a10 = -.5*P[0][1] + .5*P[2][1];
  a11 = .25*P[0][0] - .25*P[0][2] - .25*P[2][0] + .25*P[2][2];
  a12 = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + .5*P[2][0] - 1.25*P[2][1] + P[2][2] - .25*P[2][3];
  a13 = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .25*P[2][0] + .75*P[2][1] - .75*P[2][2] + .25*P[2][3];
  a20 = P[0][1] - 2.5*P[1][1] + 2*P[2][1] - .5*P[3][1];
  a21 = -.5*P[0][0] + .5*P[0][2] + 1.25*P[1][0] - 1.25*P[1][2] - P[2][0] + P[2][2] + .25*P[3][0] - .25*P[3][2];
  a22 = P[0][0] - 2.5*P[0][1] + 2*P[0][2] - .5*P[0][3] - 2.5*P[1][0] + 6.25*P[1][1] - 5*P[1][2] + 1.25*P[1][3] + 2*P[2][0] - 5*P[2][1] + 4*P[2][2] - P[2][3] - .5*P[3][0] + 1.25*P[3][1] - P[3][2] + .25*P[3][3];
  a23 = -.5*P[0][0] + 1.5*P[0][1] - 1.5*P[0][2] + .5*P[0][3] + 1.25*P[1][0] - 3.75*P[1][1] + 3.75*P[1][2] - 1.25*P[1][3] - P[2][0] + 3*P[2][1] - 3*P[2][2] + P[2][3] + .25*P[3][0] - .75*P[3][1] + .75*P[3][2] - .25*P[3][3];
  a30 = -.5*P[0][1] + 1.5*P[1][1] - 1.5*P[2][1] + .5*P[3][1];
  a31 = .25*P[0][0] - .25*P[0][2] - .75*P[1][0] + .75*P[1][2] + .75*P[2][0] - .75*P[2][2] - .25*P[3][0] + .25*P[3][2];
  a32 = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + 1.5*P[1][0] - 3.75*P[1][1] + 3*P[1][2] - .75*P[1][3] - 1.5*P[2][0] + 3.75*P[2][1] - 3*P[2][2] + .75*P[2][3] + .5*P[3][0] - 1.25*P[3][1] + P[3][2] - .25*P[3][3];
  a33 = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .75*P[1][0] + 2.25*P[1][1] - 2.25*P[1][2] + .75*P[1][3] + .75*P[2][0] - 2.25*P[2][1] + 2.25*P[2][2] - .75*P[2][3] - .25*P[3][0] + .75*P[3][1] - .75*P[3][2] + .25*P[3][3];
}

void Interpolate::updateP(double r, double z) {
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  int is = (int)(Grid_.celli(r));
  int js = (int)(Grid_.cellj(z));

  if (is-1 < 0 || is+2 >= nr || js-1 < 0 || js+2 >=nz) {
      throw OutsideGrid;
  }

  r_curr = (Grid_.R_[is]);
  //printf("R0 = %f\n", Grid_.R_[0]);
  //printf("r_curr = %f\n", r_curr);
  z_curr = (Grid_.z_[js]);
  //printf("z0 = %f\n", Grid_.z_[0]);
  //printf("z_curr = %f\n", z_curr);
  // Fill in P
  for (int i = 0; i < 4 ; i++) {
    for (int j = 0; j < 4 ; j++) {
      P[i][j] = Psi_.f_[is + i - 1][js + j - 1];
    }
  }
}


void Interpolate::updateInterpolation(double r, double z){
    updateP(r,z);
    updateCoefficients();
}
