/*!
 * @file critical.cc
 * @author Elizabeth J. Paul
 * @brief Implementation for the Critical class.
 */
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
      try {
	Inter_.updateInterpolation(limiters_.data(i,0), limiters_.data(i,1));
      } catch (int i){
	if (i == OutsideGrid) printf("Error: physical limiter outside/too close to computational boundary\n");
	if (i == OutsideInterp) printf("Error: outside Interp\n");
	throw i;
      }
        Psi_phys_lim.push_back(0);
    }
  try {
      Inter_.updateInterpolation(R_stag_up, z_stag_up);
      Psi_stag_up = Inter_.F(R_stag_up, z_stag_up);
  }
  catch(int i) {
      if (i == OutsideGrid) printf("Error: upper saddle point outside grid\n");
      if (i == OutsideInterp) printf("Error: outside Interp\n");
      throw i;
  }

  // Interpolate Psi_ at (R0, z_limiter2)
  try {
      Inter_.updateInterpolation(R_stag_down, z_stag_down);
      Psi_stag_down = Inter_.F(R_stag_down, z_stag_down);
  }
  catch(int i) {
      if (i == OutsideGrid) printf("Error: lower saddle point outside grid\n");
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
  Psi_rr = Inter_.F_rr(r, z);

  Psi_rz = Inter_.F_rz(r, z);
  
  Psi_zz = Inter_.F_zz(r, z);
  Psi_z = Inter_.F_z(r, z);
  Psi_r = Inter_.F_r(r, z);
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
        Inter_.updateInterpolation(r,z);
        Psi_r = Inter_.F_r(r,z);
        Psi_z = Inter_.F_z(r,z);
        Psi_rz = Inter_.F_rz(r,z);
        Psi_zz = Inter_.F_zz(r,z);
        Psi_rr = Inter_.F_rr(r,z);

    }
    catch(int i) {
    // If r or z outside grid, use original coordinates

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
	Psi_min_ = Inter_.F(r,z);
        *Psi_min = Psi_min_;
        *rcrit = r;
        *zcrit = z;
        return;
      }
      //      printf("not a minimum at R = %f, z = %f, D = %f, Psirr = %f, Psizz = %f, Psi = %f\n", r, z, D, Psi_rr, Psi_zz, Inter_.F(r,z));
      // if this fails, we are at a stationary point, but it's not a minimum so searching will keep bringing us back here and we break.
      break;
    }
    Psi_search(r, z, &dr, &dz);
    //    printf("dr = %f\n", dr);
    //    printf("dz = %f\n", dz);
    r += dr;
    z += dz;
    // Check if outside boundaries
    if (z >= Grid_.z_[Grid_.nz_-1]|| z <= Grid_.z_[0]) {
      printf("z outside boundaries %f\n",z);
      break;
    }
    // Check if outside grid boundaries
    if (r <= Grid_.R_[0] || r >= Grid_.R_[Grid_.nr_-1]) {
      printf("R outside boundaries %f\n",r);
      break;
    }
  }

  // If search fails, brute force it.
  Psi_min_ = Psi_.f_[0][0];
  for (int i = 0; i < Grid_.nr_; ++i){
      for (int j = 0; j < Grid_.nz_; ++j){
          if (Psi_.f_[i][j] < Psi_min_){
              Psi_min_ = Psi_.f_[i][j];
              R0 = Grid_.R_[i]; 
              z0 = Grid_.z_[j];
          }
      }
  }
  //  Psi_min_ = Inter_.F(R0, z0);
  printf("o point critical search failed, new min at %f, %f, psi = %f\n",R0,z0, Psi_min_);
  *Psi_min = Psi_min_;
  *rcrit = R0;
  *zcrit = z0;
  return;
}

/*!
 * @brief Perform search for limiter points
 */
double Critical::Psi_limiter() {
  double Psi_phys;

  for (int i = 0; i < limiters_.num_rows(); ++i){
      Inter_.updateInterpolation(limiters_.data(i,0), limiters_.data(i,1));
      Psi_phys_lim[i] = Inter_.F(limiters_.data(i,0), limiters_.data(i,1));
      //printf("limiter %d at R = %f, z = %f \n", i+1, limiters_.data(i,0), limiters_.data(i,1));
  }

  Psi_phys = *std::min_element(Psi_phys_lim.begin(), Psi_phys_lim.end());

  double r = R_stag_down; double z = z_stag_down;
  if (find_saddle(r,z)){
      R_stag_down = r;
      z_stag_down = z;
      Inter_.updateInterpolation(R_stag_down,z_stag_down);
      Psi_stag_down = Inter_.F(R_stag_down,z_stag_down);
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
      Psi_stag_up = Inter_.F(R_stag_up,z_stag_up);
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

/*!
 * @brief Finds saddle point using critical point search starting with r and z. 
 *
 * Returns false if search lands outside of horizontal limiters or grid 
 * boundaries.
 */

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
      Psi_r = Inter_.F_r(r,z);
      Psi_z = Inter_.F_z(r,z);
      Psi_rz = Inter_.F_rz(r,z);

      Psi_zz = Inter_.F_zz(r,z);
      Psi_rr = Inter_.F_rr(r,z);
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
