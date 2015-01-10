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
R0(R0),
z0(z0),
phys_lim_R(R0),
phys_lim_zup(z_limiter1),
phys_lim_zdown(z_limiter2) {
  R_stag_up = R0;
  z_stag_up = z_limiter1;
  R_stag_down = R0;
  z_stag_down = z_limiter2;
  assert(phys_lim_zup > phys_lim_zdown);
  // Interpolate Psi_ at (R0, z_limiter1)
  try {
    Inter_.updateInterpolation(R0, z_limiter1);
  }
  catch(int i) {
    if (i == OutsideGrid) printf("Error: limiter1 outside grid\n");
  }
  Psi_stag_up = Inter_.Psi_interp(R0, z_limiter1);
  // Interpolate Psi_ at (R0, z_limiter2)
  try {
    Inter_.updateInterpolation(R0, z_limiter2);
  }
  catch(int i) {
    if (i == OutsideGrid) printf("Error: limiter2 outside grid\n");
  }
  Psi_stag_down = Inter_.Psi_interp(R0, z_limiter2);

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
      Inter_.updateInterpolation(r,z);
    }
    // If r or z outside grid, use original coordinates
    catch(int i) {
      if (i == OutsideGrid) {
	printf("Interpolation outside grid\n");
	break;
      }
    }
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
    if (z >= phys_lim_zup || z <= phys_lim_zdown) break;
    // Check if outside grid boundaries
    if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
  }
  // If search failed, use original coordinates of magnetic axis
  Inter_.updateInterpolation(R0,z0);
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
double Critical::Psi_limiter() {
  double Psi_phys_up, Psi_phys_down, Psi_phys;

  Inter_.updateInterpolation(phys_lim_R,phys_lim_zup);
  Psi_phys_up = Inter_.Psi_interp(phys_lim_R, phys_lim_zup);
  Inter_.updateInterpolation(phys_lim_R,phys_lim_zdown);
  Psi_phys_down = Inter_.Psi_interp(phys_lim_R, phys_lim_zdown);  

  Psi_phys = fmin(Psi_phys_up, Psi_phys_down);
  
  double r = R_stag_down; double z = z_stag_down;
  if (find_saddle(r,z)){
      R_stag_down = r;
      z_stag_down = z;
      Inter_.updateInterpolation(R_stag_down,z_stag_down);
      Psi_stag_down = Inter_.Psi_interp(R_stag_down,z_stag_down);
  } else {
      Psi_stag_down = Psi_phys + 1;
  }
  r = R_stag_up; z = z_stag_up;
  if (find_saddle(r,z)){
      R_stag_up = r;
      z_stag_up = z;
      Inter_.updateInterpolation(R_stag_up,z_stag_up);
      Psi_stag_up = Inter_.Psi_interp(R_stag_up,z_stag_up);
  } else {
      Psi_stag_up = Psi_phys + 1;
  }

  if (Psi_phys < Psi_stag_up && Psi_phys < Psi_stag_down) {
      // physical limiter is innermost
      printf("physical!, psi=%f\n", Psi_phys);
      return Psi_phys;
  } else {
      printf("stagnation\n");
      return fmin(Psi_stag_up, Psi_stag_down);
  }
}

/*!
 * @brief Performs critical point search; updates Psi_l and Psi_o
 */
void Critical::update() {
  double rcrit, zcrit, Psi_min;
  // Calculate Psi_l using previous coordinates for limiter

  // Update Psi_l, r_l, and z_l
  Psi_.f_l = Psi_limiter();
  // Calculate Psi_o using previous coordinates
  Psi_magnetic(R0, z0, &rcrit, &zcrit, &Psi_min);
  Psi_.f_0 = Psi_min;
  R0 = rcrit;
  z0 = zcrit;
}

bool Critical::find_saddle(double &r, double &z){
    double Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz, D;
    double dr =0;
    double dz=0;
  for (int i = 0; i < max_iter; ++i) {
    assert(!isnan(r));
    assert(!isnan(z));
    try {
      Inter_.updateInterpolation(r,z);
    }
    // If r or z outside grid, use limiters
    catch(int i) {
      if (i == OutsideGrid) break;
    }
    // Calculate |del Psi(r,z)|
    // Update Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
      Psi_r = Inter_.Psir_interp(r,z);
      Psi_z = Inter_.Psiz_interp(r,z);
      Psi_rz = Inter_.Psirz_interp(r,z);

      Psi_zz = Inter_.Psizz_interp(r,z);
      Psi_rr = Inter_.Psirr_interp(r,z);
      // If within tolerence
      if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
          D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
          // If critical point corresponds to a saddle point
          if (D < 0) {
              return true;
          }
          break;
      }
      Psi_search(r, z, &dr, &dz);
      assert(!isnan(dr));
      assert(!isnan(dz));
      r += dr;
      z += dz;
      // Check if outside limiters
      if (z >= Grid_.z_[Grid_.nz_-1] || z <= Grid_.z_[0]) break;
      // Check if outside grid boundaries
      if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
  }
  return false;
}
