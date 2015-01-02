#include "critical.h"
#include <math.h>
#include "field.h"
#include "grid.h"
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

const int OutsideGrid = -1;

Critical::Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0) :
Psi_(Psi),
max_iter(max_iter),
Grid_(GridS),
epsilon(epsilon),
z_limiter1(z_limiter1),
z_limiter2(z_limiter2),
R0(R0),
z0(z0),
P(4, std::vector<double>(4)) {
  Rl = R0;
  zl = z_limiter1;
  // Interpolate Psi_ at (R0, z_limiter1)
  try {
    updateP(R0,z_limiter1);
  }
  catch(int i) {
    printf("Error: limiter1 outside grid\n");
  }
  updateCoefficients();
  Psi_lim1 = Psi_interp(R0, z_limiter1);
  // Interpolate Psi_ at (R0, z_limiter2)
  try {
    updateP(R0, z_limiter2);
  }
  catch(int i) {
    printf("Error: limiter2 outside grid\n");
  }
  updateCoefficients();
  Psi_lim1 = Psi_interp(R0, z_limiter2);
  printf("INITIAL\n");
  printf("Rl = %f\n", Rl);
  printf("zl = %f\n", zl);
  printf("R0 = %f\n", R0);
  printf("z0 = %f\n", z0);
}

Critical::~Critical() {}

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
double Critical::Psi_interp(double r, double z) {
  double r2 = r*r;
  double r3 = r2*r;
  double z2 = z*z;
  double z3 = z2*z;
  return (a00 + a01*r + a02*r2 + a03*r3) +
  (a10 + a11*r + a12*r2 + a13*r3)*z +
  (a20 + a21*r + a22*r2 + a23*r3)*z2 +
  (a30 + a31*r + a32*r2 + a33*r3)*z3;
}

double Critical::Psir_interp(double r, double z) {
  double r2 = r*r;
  double z2 = z*z;
  double z3 = z2*z;
  return (a01 + 2*a02*r + 3*a03*r2) +
         (a11 + 2*a12*r + 3*a13*r2)*z +
         (a21 + 2*a22*r + 3*a23*r2)*z2 +
         (a31 + 2*a32*r + 3*a33*r2)*z3;
}

double Critical::Psirr_interp(double r, double z) {
  double z2 = z*z;
  double z3 = z2*z;
  return (2*a02 + 6*a03*r) +
  (2*a12 + 6*a13*r)*z +
  (2*a22 + 6*a23*r)*z2 +
  (2*a32 + 6*a33*r)*z3;
}

double Critical::Psirz_interp(double r, double z) {
  double r2 = r*r;
  double z2 = z*z;
  return (a11 + 2*a12*r + 3*a13*r2) +
       2*(a21 + 2*a22*r + 3*a23*r2)*z +
       3*(a31 + 2*a32*r + 3*a33*r2)*z2;
}

double Critical::Psiz_interp(double r, double z) {
  double r2 = r*r;
  double r3 = r2*r;
  double z2 = z*z;
  return (a10 + a11*r + a12*r2 + a13*r3) +
  2*(a20 + a21*r + a22*r2 + a23*r3)*z +
  3*(a30 + a31*r + a32*r2 + a33*r3)*z2;
}

double Critical::Psizz_interp(double r, double z) {
  double r2 = r*r;
  double r3 = r2*r;
  return 2*(a20 + a21*r + a22*r2 + a23*r3) +
  6*(a30 + a31*r + a32*r2 + a33*r3)*z;
}

/*!
 * @brief returns dr, dz to progress toward critical point in Psi
 */
void Critical::Psi_search(double r, double z, double *dr, double *dz) {
  double D;
  double Psi_rr = Psirr_interp(r, z);
  double Psi_rz = Psirz_interp(r, z);
  double Psi_zz = Psizz_interp(r, z);
  double Psi_z = Psiz_interp(r, z);
  double Psi_r = Psir_interp(r, z);
  // Update Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
  D = Psi_rr*Psi_zz - pow(Psi_rz,2);
  assert(D != 0);
  *dr = (-Psi_zz*Psi_r + Psi_rz*Psi_z)*(1.0/D);
  *dz = (Psi_rz*Psi_r - Psi_rr*Psi_z)*(1.0/D);
}

void Critical::updateCoefficients() {
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

void Critical::updateP(double r, double z) {
  int is, js;
  int i = (int)Grid_.celli(r);
  int j = (int)Grid_.cellj(z);
  if (Grid_.celli(r) - i >= 0) is = i - 1;
  else is = i - 2;
  if (Grid_.cellj(r) - j >= 0) js = j - 1;
  else js = j - 2;
  // Error checking - within grid bounds?
  if (is < 0 || (is+3) > (Grid_.nr_-1) || js < 0 || (js+3) > (Grid_.nz_-1)) {
    throw OutsideGrid;
  }
  printf("is = %d\n", is);
  printf("js = %d\n", js);
  // Fill in P
  P[0][0] = Psi_.f_[is][js];
  P[0][1] = Psi_.f_[is][js+1];
  P[0][2] = Psi_.f_[is][js+2];
  P[0][3] = Psi_.f_[is][js+3];
  P[1][0] = Psi_.f_[is+1][js];
  P[1][1] = Psi_.f_[is+1][js+1];
  P[1][2] = Psi_.f_[is+1][js+2];
  P[1][3] = Psi_.f_[is+1][js+3];
  P[2][0] = Psi_.f_[is+2][js];
  P[2][1] = Psi_.f_[is+2][js+1];
  P[2][2] = Psi_.f_[is+2][js+2];
  P[2][3] = Psi_.f_[is+2][js+3];
  P[3][0] = Psi_.f_[is+3][js];
  P[3][1] = Psi_.f_[is+3][js+1];
  P[3][2] = Psi_.f_[is+3][js+2];
  P[3][3] = Psi_.f_[is+3][js+3];
}

/*!
 * @brief Perform search for critical points beginning with initial
 * guess r, z
 */
void Critical::Psi_magnetic(double r, double z, double *rcrit, double *zcrit, double *Psi_min) {
  double Psi_min_;
  double dr, dz;
  double Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz, D;
  for (int i = 0; i < max_iter; ++i) {
    try {
      updateP(r,z);
    }
    // If r or z outside grid, use original coordinates
    catch(int i) {
      if (i == OutsideGrid) break;
    }
    updateCoefficients();
    Psi_r = Psir_interp(r,z);
    Psi_z = Psiz_interp(r,z);
    Psi_rz = Psirz_interp(r,z);
    Psi_zz = Psizz_interp(r,z);
    Psi_rr = Psirr_interp(r,z);
    // Calculate |del Psi(r,z)| - if within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      // Second derivative test
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a minimum
      if (D > 0) {
        double Psi_min_ = Psi_interp(r,z);
        *Psi_min = Psi_min_;
        *rcrit = r;
        *zcrit = z;
        Psi_min_ = (Psi_interp(*rcrit, *zcrit));
        *Psi_min = Psi_min_;
        return;
      }
      Psi_search(r, z, &dr, &dz);
      r += dr;
      z += dz;
      // Check if outside limiters
      if (z >= z_limiter1 || z <= z_limiter2) break;
      // Check if outside grid boundaries
      if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
    }
  }
  // If search failed, use original coordinates of magnetic axis
  Psi_min_ = Psi_interp(R0, z0);
  *Psi_min = Psi_min_;
    
  *rcrit = R0;
  *zcrit = z0;
  Psi_min_ = (Psi_interp(*rcrit, *zcrit));
  *Psi_min = Psi_min_;
  return;
}

/*!
 * @brief Perform search for critical points beginning with initial
 * guess r, z
 */
void Critical::Psi_limiter(double r, double z, double *rcrit, double *zcrit, double *Psi_min) {
  double dr, dz;
  double Psi_lim1, Psi_lim2;
  double Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz, D, Psi_crit;
  
  // Calculate minimum over limiters
  Psi_lim1 = Psi_interp(Rl, z_limiter1);
  Psi_lim2 = Psi_interp(Rl, z_limiter2);
  if(Psi_lim1 < Psi_lim2) {
    *rcrit = Rl;
    *zcrit = z_limiter1;
    *Psi_min = Psi_lim1;
  }
  else {
    *rcrit = R0;
    *zcrit = z_limiter2;
    *Psi_min = Psi_lim2;
  }
  for (int i = 0; i < max_iter; ++i) {
    assert(!isnan(r));
    assert(!isnan(z));
    try {
      updateP(r,z);
    }
    // If r or z outside grid, use limiters
    catch(int i) {
      if (i == OutsideGrid) break;
    }
    // Calculate |del Psi(r,z)|
    // Update Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
    updateCoefficients();
    Psi_r = Psir_interp(r,z);
    Psi_z = Psiz_interp(r,z);
    Psi_rz = Psirz_interp(r,z);
    Psi_zz = Psizz_interp(r,z);
    Psi_rr = Psirr_interp(r,z);
    // If within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a saddle point
      // compare with limiters and return
      if (D < 0) {
        printf("D < 0\n");
        Psi_crit = Psi_interp(r, z);
        if(Psi_crit < *Psi_min) {
          *rcrit = r;
          *zcrit = z;
          *Psi_min = Psi_crit;
        }
        return;
      }
      // If not a saddle point, use limiters
      break;
    }
    Psi_search(r, z, &dr, &dz);
    assert(!isnan(dr));
    assert(!isnan(dz));
    r += dr;
    z += dz;
    // Check if outside limiters
    if (z >= z_limiter1 || z <= z_limiter2) break;
    // Check if outside grid boundaries
    if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
  }
  // Alternate case - use previous value of limiter
  *rcrit = Rl;
  *zcrit = zl;
  *Psi_min = Psi_.f_l;
  return;
}

/*!
 * @brief Performs critical point search; updates Psi_l and Psi_o
 */
void Critical::update() {
  double rcrit, zcrit, Psi_min;
  // Calculate Psi_l using previous coordinates for limiter
  assert(!isnan(Rl));
  assert(!isnan(zl));
  if (Psi_lim1 < Psi_lim2) Psi_min = Psi_lim1;
  else Psi_min = Psi_lim2;
  Psi_limiter(Rl, zl, &rcrit, &zcrit, &Psi_min);
  // Update Psi_l, r_l, and z_l
  Psi_.f_l = Psi_min;
  Rl = rcrit;
  zl = zcrit;
  // Calculate Psi_o using previous coordinates
  Psi_magnetic(R0, z0, &rcrit, &zcrit, &Psi_min);
  Psi_.f_0 = Psi_min;
  R0 = rcrit;
  z0 = zcrit;
  printf("Rl = %f\n", Rl);
  printf("zl = %f\n", zl);
  printf("R0 = %f\n", R0);
  printf("z0 = %f\n", z0);
}



