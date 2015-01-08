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
P(4, std::vector<double>(4)) {}

Interpolate::~Interpolate() {}

/*!
 * @brief Bivariate interpolation of Psi
 *
 * preserves smoothness as described in Akima, 1974
 *
 * Do you mean "A Method of Bivariate Interpolation and Smooth Surface Fitting Based on Local Procedures", Communications of the ACM, Volume 17 Issue 1, Jan. 1974 ? I will assume that this is the paper you were referring to. I will call that paper [[Akima 1974]] with the double [[]].
 * I'm not quite sure why you chose variable names here which are different from those in the paper... I will do my best to annotate this method. -JAS
 *
 * Calculates Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
 * Determines Psi = sum(a_ij*r^i z^j) inside rectangle defined by (r[i], r[i+1])x(z[i], z[i+1])
 */
double Interpolate::Psi_interp(double r, double z) {
    double r_, z_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    double r2 = r_*r_;
    double r3 = r2*r_;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    double z3 = z2*z_;
    return (a00 + a10*r_ + a20*r2 + a30*r3) +
    (a01 + a11*r_ + a21*r2 + a31*r3)*z_ +
    (a02 + a12*r_ + a22*r2 + a32*r3)*z2 +
    (a03 + a13*r_ + a23*r2 + a33*r3)*z3;
}

double Interpolate::Psir_interp(double r, double z) {
    double r_, z_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    double r2 = r_*r_;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    double z3 = z2*z_;
    return (a10 + 2*a20*r_ + 3*a30*r2) +
    (a11 + 2*a21*r_ + 3*a31*r2)*z_ +
    (a12 + 2*a22*r_ + 3*a32*r2)*z2 +
    (a13 + 2*a23*r_ + 3*a33*r2)*z3;
}

double Interpolate::Psirr_interp(double r, double z) {
    double z_, r_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    double z3 = z2*z_;
    return (2*a20 + 6*a30*r) +
    (2*a21 + 6*a31*r)*z_ +
    (2*a22 + 6*a32*r)*z2 +
    (2*a23 + 6*a33*r)*z3;
}

double Interpolate::Psirz_interp(double r, double z) {
    double r_, z_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    double r2 = r_*r_;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    return (a11 + 2*a21*r_ + 3*a31*r2) +
    2*(a12 + 2*a22*r_ + 3*a32*r2)*z_ +
    3*(a13 + 2*a23*r_ + 3*a33*r2)*z2;
}

double Interpolate::Psiz_interp(double r, double z) {
    double z_, r_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    double r2 = r_*r_;
    double r3 = r2*r_;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    return (a01 + a11*r_ + a21*r2 + a31*r3) +
    2*(a02 + a12*r_ + a22*r2 + a32*r3)*z_ +
    3*(a03 + a13*r_ + a23*r2 + a33*r3)*z2;
}

double Interpolate::Psizz_interp(double r, double z) {
    double r_, z_;
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    if ((r - r_curr) == 0) r_ = 0;
    if ((z - z_curr) == 0) z_ = 0;
    if ((r - r_curr) > dr || (z - z_curr) > dz || (r - r_curr) < 0 || (z - z_curr) < 0) {
        printf("dr = %f\n", dr);
        printf("r_curr = %f\n", r_curr);
        printf("r = %f\n", r);
        printf("dz = %f\n", dz);
        printf("z_curr = %f\n", z_curr);
        printf("z = %f\n", z);
        throw OutsideInterp;
    }
    r_ = (r - r_curr)/dr;
    double r2 = r_*r_;
    double r3 = r2*r_;
    z_ = (z - z_curr)/dz;
    double z2 = z_*z_;
    return 2*(a02 + a12*r_ + a22*r2 + a32*r3) +
    6*(a03 + a13*r_ + a23*r2 + a33*r3)*z2;
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
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    double nr = Grid_.nr_;
    double nz = Grid_.nz_;
    int is = (int)(Grid_.celli(r));
    int js = (int)(Grid_.cellj(z));
    
    if (is-1 < 0 || is+2 >= nr || js-1 < 0 || js+2 >=nz) {
        throw OutsideGrid;

    }

    r_curr = is*dr + (Grid_.R_[0]);
    printf("R0 = %f\n", Grid_.R_[0]);
    printf("r_curr = %f\n", r_curr);
    z_curr = js*dz + (Grid_.z_[0]);
    printf("z0 = %f\n", Grid_.z_[0]);
    printf("z_curr = %f\n", z_curr);
    // Fill in P
    P[0][0] = Psi_.f_[is-1][js-1];
    P[0][1] = Psi_.f_[is-1][js];
    P[0][2] = Psi_.f_[is-1][js+1];
    P[0][3] = Psi_.f_[is-1][js+2];
    P[1][0] = Psi_.f_[is][js-1];
    P[1][1] = Psi_.f_[is][js];
    P[1][2] = Psi_.f_[is][js+1];
    P[1][3] = Psi_.f_[is][js+2];
    P[2][0] = Psi_.f_[is+1][js-1];
    P[2][1] = Psi_.f_[is+1][js];
    P[2][2] = Psi_.f_[is+1][js+1];
    P[2][3] = Psi_.f_[is+1][js+2];
    P[3][0] = Psi_.f_[is+2][js-1];
    P[3][1] = Psi_.f_[is+2][js];
    P[3][2] = Psi_.f_[is+2][js+1];
    P[3][3] = Psi_.f_[is+2][js+2];
}
