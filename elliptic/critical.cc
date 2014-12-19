#include "critical.h"
#include <math.h>
#include "field.h"
#include "grid.h"
#include <vector>
#include <stdio.h>
#include <assert.h>

Critical::Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0) :
  Psi_(Psi),
  max_iter(max_iter),
  Grid_(GridS),
  epsilon(epsilon),
  z_limiter1(z_limiter1),
  z_limiter2(z_limiter2),
  R0(R0),
  z0(z0) {
    alpha = new double*[Grid_.nr_];
    beta = new double*[Grid_.nr_];
    for(int i = 0; i < Grid_.nr_; ++i) {
      alpha[i] = new double[Grid_.nz_];
      beta[i] = new double[Grid_.nz_];
    }
    Rl = R0;
    zl = z_limiter1;
}

Critical::~Critical() {
  for (int i = 0; i < Grid_.nr_; ++i) {
    delete [] alpha[i];
    delete [] beta[i];
  }
  delete [] alpha;
  delete [] beta;
}

/*!
 * @brief Bivariate interpolation of Psi
 *
 * preserves smoothness as described in Akima, 1974
 * Calculates Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
 * Determines Psi = r^alpha z^beta inside rectangle defined by (r[i], r[i+1])x(z[i], z[i+1])
 */
void Critical::interpolate() {
  double c1, c2, c3, c4, c5, c6, c7, c8, d1, d2, d3, d4, d5, w_r1, w_r2, w_z1, w_z2, e1, e2, e3, e4;
  double *r = Grid_.R_;
  double *z = Grid_.z_;
  double **Psi = Psi_.f_;
  for (int i = 3; i < Grid_.nr_-3; ++i) {
    for(int j = 3; j < Grid_.nz_-3; ++j) {
      // C[i-1][j]
      c1 = (Psi[i][j] - Psi[i-1][j])/(r[i]-r[i-1]);
      // C[i][j]
      c2 = (Psi[i+1][j] - Psi[i][j])/(r[i+1]-r[i]);
      // C[i][j+1]
      c3 = (Psi[i+1][j+1] - Psi[i][j+1])/(r[i+1]-r[i]);
      // C[i-1][j-1]
      c4 = (Psi[i][j-1] - Psi[i-1][j-1])/(r[i]-r[i-1]);
      // c[i][j-1]
      c5 = (Psi[i+1][j-1] - Psi[i][j-1])/(r[i+1]-r[i]);
      // c[i-1][j+1]
      c6 = (Psi[i][j+1] - Psi[i-1][j+1])/(r[i]-r[i-1]);
      // c[i+1][j]
      c7 = (Psi[i+2][j] - Psi[i+1][j])/(r[i+2]-r[i+1]);
      // c[i-2][j]
      c8 = (Psi[i-2][j] - Psi[i-3][j])/(r[i-2]-r[i-3]);
      // D[i][j-1]
      d1 = (Psi[i][j] - Psi[i][j-1])/(z[j]-z[j-1]);
      // D[i][j]
      d2 = (Psi[i][j+1] - Psi[i][j])/(z[j+1]-z[j]);
      // D[i+1][j]
      d3 = (Psi[i+1][j+1] - Psi[i+1][j])/(z[j+1]-z[j]);
      // D[i][j+1]
      d4 = (Psi[i][j+2] - Psi[i][j+1])/(z[j+2]-z[j+1]);
      // D[i][j-2]
      d5 = (Psi[i][j-1] - Psi[i][j-2])/(z[j-1]-z[j-2]);
      w_r1 = abs(c7 - c2);
      w_r2 = abs(c1 - c8);
      if (w_r1 == w_r2 && w_r1 == 0) {
        w_r1 = 1;
        w_r2 = 1;
      }
      w_z1 = abs(d4 - d2);
      w_z2 = abs(d1 - d5);
      if (w_z1 == w_z2 && w_z1 == 0) {
        w_z1 = 1;
        w_z2 = 1;
      }
      double Psi_r = (w_r1*c1 + w_r2*c2)/(w_r1 + w_r2);
      double Psi_z = (w_z1*d1 + w_z2*d2)/(w_z1 + w_z2);
      // e[i][j]
      e1 = (c3 - c2)/(z[j+1]-z[j]);
      // e[i-1][j-1] = (c[i-1][j] - c[i-1][j-1])/(z[j]-z[j-1])
      e2 = (c1 - c4)/(z[j]-z[j-1]);
      // e[i][j-1] = (c[i][j] - c[i][j-1])/(z[j-1]-z[j-2])
      e3 = (c2 - c5)/(z[j-1] - z[j-2]);
      // e[i-1][j] = (c[i-1][j+1] - c[i-1][j])/(z[j+1] - z[j])
      e4 = (c6 - c1)/(z[j+1] - z[j]);
      double Psi_rz = (w_r1*(w_z1*e2 + w_z2*e4) + w_r2*(w_z1*e3 + w_z2*e1))/((w_r1 + w_r2)*(w_z1 + w_z2));
      assert(Psi_z != 0);
      assert(Psi_r != 0);
      alpha[i][j] = Psi_rz*r[i]/Psi_z;
      beta[i][j] = Psi_rz*z[j]/Psi_r;
      if (isnan(beta[i][j])) {
        printf("BETA NAN\n");
        printf("i = %d\n", i);
        printf("j = %d\n", j);
      }
      assert(!isnan(beta[i][j]));
      assert(!isnan(alpha[i][j]));
    }
  }
}

/*!
 * @brief Find alpha for given position
 */
double Critical::cell_alpha(double r, double z) {
  // Find containing cell's indices
  double i = Grid_.celli(r);
  double j = Grid_.cellj(z);
  // Determine where bivariate polynomial is defined
  if ((int)Grid_.celli(r) - i < 0) --i;
  if ((int)Grid_.cellj(z) - j < 0) --j;
  return alpha[(int)i][(int)j];
}

/*!
 * @brief Find beta for given position
 */
double Critical::cell_beta(double r, double z) {
  // Find containing cell's indices
  double i_ = Grid_.celli(r);
  double j_ = Grid_.cellj(z);
  int i = (int)i_;
  int j = (int)j_;
  if (i - i_ < 0) --i;
  if (j - j_ < 0) --j;
  assert(i >= 0);
  assert(i < Grid_.nr_);
  assert(j >= 0);
  assert(j < Grid_.nz_);
  if (isnan(beta[i][j])) {
    printf("IS NAN\n");
    printf("i = %d\n", i);
    printf("j = %d\n", j);
  }
  assert(!isnan(beta[i][j]));
  return beta[i][j];
}

/*!
 * @brief Interpolated Psi defined for all r, z
 */
double Critical::Psi_interp(double r, double z) {
  double alpha = cell_alpha(r,z);
  double beta = cell_beta(r,z);
  return pow(r,alpha)*pow(z,beta);
}

/*!
 * @brief returns dr, dz to progress toward critical point in Psi
 */
void Critical::Psi_search(double r, double z, double *dr, double *dz) {
  double beta, alpha, Psi_zz, Psi_rr, Psi_rz, Psi_r, Psi_z, D;
  // Psi = r^alpha z^beta
  // D = Psi_rr*Psi_zz - Psi_rz^2
  beta = cell_beta(r,z);
  alpha = cell_alpha(r,z);
  Psi_zz = beta*(beta-1)*pow(z,beta-2)*pow(r,alpha);
  Psi_rr = alpha*(alpha-1)*pow(r,alpha-2)*pow(z,beta);
  Psi_rz = alpha*beta*pow(r,alpha-1)*pow(z,beta-1);
  Psi_r = (alpha)*pow(r,alpha-1)*pow(z,beta);
  Psi_z = (beta)*pow(z,beta-1)*pow(r,alpha);
  D = Psi_rr*Psi_zz - pow(Psi_rz,2);
  assert(D != 0);
  assert(!isnan(beta));
  assert(!isnan(alpha));
  *dr = (-Psi_zz*Psi_r + Psi_rz*Psi_z)*(1.0/D);
  *dz = (Psi_rz*Psi_r - Psi_rr*Psi_z)*(1.0/D);
}

/*!
 * @brief Perform search for critical points beginning with initial
 * guess r, z
 */
void Critical::Psi_magnetic(double r, double z, double *rcrit, double *zcrit, double *Psi_min) {
  double Psi_min_;
  double dr, dz, alpha, beta;
  double Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz, D;
  for (int i = 0; i < max_iter; ++i) {
    // Calculate |del Psi(r,z)|
    beta = cell_beta(r,z);
    alpha = cell_alpha(r,z);
    Psi_r = (alpha)*pow(r,(alpha-1))*pow(z, beta);
    Psi_z = (beta)*pow(z,(beta-1))*pow(r, alpha);
    // If within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      // Second derivative test
      Psi_rr = (alpha)*(alpha-1)*pow(r, (alpha-2))*pow(z,beta);
      Psi_zz = (beta)*(beta-1)*pow(z, (beta-2))*pow(r, alpha);
      Psi_rz = (alpha)*pow(r,(alpha-1))*(beta)*pow(z, (beta-1));
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a minimum
      if (D > 0) {
        rcrit = &r;
        zcrit = &z;
        Psi_min_ = (Psi_interp(*rcrit, *zcrit));
        Psi_min = &Psi_min_;
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
  rcrit = &R0;
  zcrit = &z0;
  Psi_min_ = (Psi_interp(*rcrit, *zcrit));
  Psi_min = &Psi_min_;
  return;
}

/*!
 * @brief Perform search for critical points beginning with initial
 * guess r, z
 */
void Critical::Psi_limiter(double r, double z, double *rcrit, double *zcrit, double *Psi_min) {
  double dr, dz, alpha, beta;
  double Psi_lim1, Psi_lim2;
  double Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz, D, Psi_crit;
  
  // Calculate minimum over limiters
  Psi_lim1 = Psi_interp(R0, z_limiter1);
  Psi_lim2 = Psi_interp(R0, z_limiter2);
  if(Psi_lim1 < Psi_lim2) {
    rcrit = &R0;
    zcrit = &z_limiter1;
    Psi_min = &Psi_lim1;
  }
  else {
    rcrit = &R0;
    zcrit = &z_limiter2;
    Psi_min = &Psi_lim2;
  }
  
  for (int i = 0; i < max_iter; ++i) {
    assert(!isnan(r));
    assert(!isnan(z));
    // Calculate |del Psi(r,z)|
    beta = cell_beta(r,z);
    alpha = cell_alpha(r,z);
    Psi_r = (alpha)*pow(r,(alpha-1))*pow(z, beta);
    Psi_z = (beta)*pow(z,(beta-1))*pow(r, alpha);
    // If within tolerence
    if (sqrt(Psi_r*Psi_r + Psi_z*Psi_z) < epsilon){
      // Second derivative test
      Psi_rr = (alpha)*(alpha-1)*pow(r, (alpha-2))*pow(z,beta);
      Psi_zz = (beta)*(beta-1)*pow(z, (beta-2))*pow(r, alpha);
      Psi_rz = (alpha)*pow(r,(alpha-1))*(beta)*pow(z, (beta-1));
      D = Psi_rr*Psi_zz - Psi_rz*Psi_rz;
      // If critical point corresponds to a saddle point
      // compare with limiters and return
      if (D < 0) {
        Psi_crit = Psi_interp(r, z);
        if(Psi_crit < *Psi_min) {
          rcrit = &r;
          zcrit = &z;
          Psi_min = &Psi_crit;
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


  
