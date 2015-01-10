#include "critical.h"
#include "interpolate.h"
#include "field.h"
#include "grid.h"
#include <math.h>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

const int OutsideGrid = -1;
const int OutsideInterp = -2;

Critical::Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0) :
Psi_(Psi),
max_iter(max_iter),
Grid_(GridS),
Inter_(Grid_, Psi_),
epsilon(epsilon),
z_limiter1(z_limiter1),
z_limiter2(z_limiter2),
R0(R0),
z0(z0) {
  Rl = R0;
  zl = z_limiter1;
  // Interpolate Psi_ at (R0, z_limiter1)
  try {
    Inter_.updateP(R0, z_limiter1);
  }
  catch(int i) {
//    if (i == OutsideGrid) printf("Error: limiter1 outside grid\n");
  }
  Inter_.updateCoefficients();
  Psi_lim1 = Inter_.Psi_interp(R0, z_limiter1);
  // Interpolate Psi_ at (R0, z_limiter2)
  try {
    Inter_.updateP(R0, z_limiter2);
  }
  catch(int i) {
    if (i == OutsideGrid) printf("Error: limiter2 outside grid\n");
  }
  Inter_.updateCoefficients();
  Psi_lim2 = Inter_.Psi_interp(R0, z_limiter2);

//  printf("INITIAL\n");
//  printf("Rl = %f\n", Rl);
//  printf("zl = %f\n", zl);
//  printf("R0 = %f\n", R0);
//  printf("z0 = %f\n", z0);
}

Critical::~Critical() {}

/*!
 * @brief returns dr, dz to progress toward critical point in Psi
 */
void Critical::Psi_search(double r, double z, double *dr, double *dz) {
  double D;
  double Psi_rr, Psi_rz, Psi_zz, Psi_z, Psi_r;
  Psi_rr = Inter_.Psirr_interp(r, z);

  Psi_rz = Inter_.Psirz_interp(r, z);
  
  Psi_zz = Inter_.Psizz_interp(r, z);
  Psi_z = Inter_.Psiz_interp(r, z);
  Psi_r = Inter_.Psir_interp(r, z);
  D = Psi_rr*Psi_zz - pow(Psi_rz,2);
  assert(D != 0);
  *dr = (-Psi_zz*Psi_r + Psi_rz*Psi_z)*(1.0/D);
  *dz = (Psi_rz*Psi_r - Psi_rr*Psi_z)*(1.0/D);
  //  printf("Psi_rr = %15.8f, Psi_rz = %15.8f, D = %15.8f, Psi_zz = %15.8f, Psi_z = %15.8f, Psi_r = %15.8f, dr = %15.8f, dz = %15.8f\n", Psi_rr, Psi_rz, D, Psi_zz, Psi_z, Psi_r ,*dr, *dz);
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
//      printf("PSI_MAGNETIC\n");
//      printf("r = %f\n", r);
//      printf("z = %f\n", z);
      Inter_.updateP(r,z);
    }
    // If r or z outside grid, use original coordinates
    catch(int i) {
      if (i == OutsideGrid) {
	printf("Interpolation outside grid\n");
	break;
      }
    }
    Inter_.updateCoefficients();
      Psi_r = Inter_.Psir_interp(r,z);
      Psi_z = Inter_.Psiz_interp(r,z);
      Psi_rz = Inter_.Psirz_interp(r,z);
      Psi_zz = Inter_.Psizz_interp(r,z);
      Psi_rr = Inter_.Psirr_interp(r,z);
    // Calculate |del Psi(r,z)| - if within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      // Second derivative test
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a minimum
      if (D > 0 && Psi_rr > 0) {
	Psi_min_ = Inter_.Psi_interp(r,z);
        *Psi_min = Psi_min_;
        *rcrit = r;
        *zcrit = z;
        return;
      }
    }
    Psi_search(r, z, &dr, &dz);
//    printf("dr = %f\n", dr);
//    printf("dz = %f\n", dz);
    r += dr;
    z += dz;
    // Check if outside limiters
    if (z >= z_limiter1 || z <= z_limiter2) break;
    // Check if outside grid boundaries
    if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
  }
  // If search failed, use original coordinates of magnetic axis
  Psi_min_ = Inter_.Psi_interp(R0, z0);

  *Psi_min = Psi_min_;
  *rcrit = R0;
  *zcrit = z0;
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
  Psi_lim1 = Inter_.Psi_interp(Rl, z_limiter1);
  Psi_lim2 = Inter_.Psi_interp(Rl, z_limiter2);
  
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
      Inter_.updateP(r,z);
    }
    // If r or z outside grid, use limiters
    catch(int i) {
      if (i == OutsideGrid) break;
    }
    // Calculate |del Psi(r,z)|
    // Update Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
    Inter_.updateCoefficients();
      Psi_r = Inter_.Psir_interp(r,z);
      Psi_z = Inter_.Psiz_interp(r,z);
      Psi_rz = Inter_.Psirz_interp(r,z);

      Psi_zz = Inter_.Psizz_interp(r,z);
      Psi_rr = Inter_.Psirr_interp(r,z);
    // If within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a saddle point
      // compare with limiters and return
      if (D < 0) {
          Psi_crit = Inter_.Psi_interp(r, z);

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
}
