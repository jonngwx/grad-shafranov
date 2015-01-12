#include "critical.h"
#include "interpolate.h"
#include "field.h"
#include "grid.h"
#include <math.h>
#include <vector>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "tsv_reader.h"
#include <algorithm>

const int OutsideGrid = -1;
const int OutsideInterp = -2;

Critical::Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, Table &limiters, double R_stag_up, double z_stag_up, double R_stag_down, double z_stag_down, double R0, double z0) :
Psi_(Psi),
max_iter(max_iter),
Grid_(GridS),
Inter_(Grid_, Psi_),
epsilon(epsilon),
R0(R0),
z0(z0),
R_stag_up_orig(R_stag_up),
z_stag_up_orig(z_stag_up),
R_stag_up(R_stag_up),
z_stag_up(z_stag_up),
R_stag_down_orig(R_stag_down),
z_stag_down_orig(z_stag_down),
R_stag_down(R_stag_down),
z_stag_down(z_stag_down),
limiters_(limiters) {
    // interpolate at stag point
    assert(limiters.num_rows()>0);
    assert(limiters.num_columns() == 2);
    for (int i = 0; i < limiters.num_rows();++i){
        Psi_phys_lim.push_back(0);
    }
  try {
      Inter_.updateInterpolation(R_stag_up, z_stag_up);
      Psi_stag_up = Inter_.Psi_interp(R_stag_up, z_stag_up);
  }
  catch(int i) {
      if (i == OutsideGrid) printf("Error: limiter1 outside grid\n");
      if (i == OutsideInterp) printf("Error: outside Interp\n");
      throw i;
  }

  // Interpolate Psi_ at (R0, z_limiter2)
  try {
      Inter_.updateInterpolation(R_stag_down, z_stag_down);
      Psi_stag_down = Inter_.Psi_interp(R_stag_down, z_stag_down);
  }
  catch(int i) {
      if (i == OutsideGrid) printf("Error: limiter2 outside grid\n");
      if (i == OutsideInterp) printf("Error: outside Interp\n");
      throw i;
  }

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
        //    printf("z = %f\n", z);
        Inter_.updateInterpolation(r,z);
        Psi_r = Inter_.Psir_interp(r,z);
        Psi_z = Inter_.Psiz_interp(r,z);
        Psi_rz = Inter_.Psirz_interp(r,z);
        Psi_zz = Inter_.Psizz_interp(r,z);
        Psi_rr = Inter_.Psirr_interp(r,z);

    }
    // If r or z outside grid, use original coordinates
    catch(int i) {
      if (i == OutsideGrid) {
          printf("Interpolation outside grid\n");
          break;
      }
      if (i == OutsideInterp) {
          printf("Interpolation failed\n");
          break;
      }
    }
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
    // Check if outside boundaries
    if (z >= Grid_.z_[Grid_.nz_-1]|| z <= Grid_.z_[0]) break;
    // Check if outside grid boundaries
    if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) break;
  }
  // If search failed, do this manually
  //  R0 = Grid_.R_[0]/2 + Grid_.R_[Grid_.nr_-1]/2;
  //  z0 = 0;
  Psi_min_ = Psi_.f_[0][0];
  //  Inter_.updateInterpolation(R0,z0);
  for (int i = 0; i < Grid_.nr_; ++i){
      for (int j = 0; j < Grid_.nz_; ++j){
          if (Psi_.f_[i][j] < Psi_min_){
              Psi_min_ = Psi_.f_[i][j];
              R0 = Grid_.R_[i];
              z0 = Grid_.z_[j];
          }
      }
  }
  //  Psi_min_ = Inter_.Psi_interp(R0, z0);
  printf("o point critical search failed, new min at %f, %f\n",R0,z0);
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
  double Psi_phys;

  for (int i = 0; i < limiters_.num_rows(); ++i){
      Inter_.updateInterpolation(limiters_.data(i,0), limiters_.data(i,1));
      Psi_phys_lim[i] = Inter_.Psi_interp(limiters_.data(i,0), limiters_.data(i,1));
      //printf("limiter %d at R = %f, z = %f \n", i+1, limiters_.data(i,0), limiters_.data(i,1));
  }

  Psi_phys = *std::min_element(Psi_phys_lim.begin(), Psi_phys_lim.end());

  double r = R_stag_down; double z = z_stag_down;
  if (find_saddle(r,z)){
      R_stag_down = r;
      z_stag_down = z;
      Inter_.updateInterpolation(R_stag_down,z_stag_down);
      Psi_stag_down = Inter_.Psi_interp(R_stag_down,z_stag_down);
  } else {
      //reset the search
      R_stag_down = R_stag_down_orig;
      z_stag_down = z_stag_down_orig;
      Psi_stag_down = Psi_phys + 1;
  }
  r = R_stag_up; z = z_stag_up;
  if (find_saddle(r,z)){
      R_stag_up = r;
      z_stag_up = z;
      Inter_.updateInterpolation(R_stag_up,z_stag_up);
      Psi_stag_up = Inter_.Psi_interp(R_stag_up,z_stag_up);
  } else {
      //reset the search
      R_stag_up = R_stag_up_orig;
      z_stag_up = z_stag_up_orig;
      Psi_stag_up = Psi_phys + 1;
  }

  if (Psi_phys < Psi_stag_up && Psi_phys < Psi_stag_down) {
      // physical limiter is innermost
      // printf("physical!, psi=%f\n", Psi_phys);
      return Psi_phys;
  } else {
      // printf("stagnation\n");
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
        break;
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
