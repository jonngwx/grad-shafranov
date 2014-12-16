#include "critical.h"
#include <math.h>
#include "field.h"
#include "grid.h"
#include <vector>
#include <stdio.h>

Critical::Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2) :
  Psi_(Psi),
  Grid_(GridS),
  max_iter(max_iter),
  epsilon(epsilon),
  z_limiter1(z_limiter1),
  z_limiter2(z_limiter2) {
    alpha = new double*[Grid_.nr_];
    beta = new double*[Grid_.nr_];
    for(int i = 0; i < Grid_.nr_; ++i) {
      alpha[i] = new double[Grid_.nz_];
      beta[i] = new double[Grid_.nz_];
    }
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
  for (int i = 1; i < Grid_.nr_-1; ++i) {
    for(int j = 1; j < Grid_.nz_-1; ++j) {
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
      alpha[i][j] = Psi_rz*r[i]/Psi_z;
      beta[i][j] = Psi_rz*z[j]/Psi_r;
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
  double i = Grid_.celli(r);
  double j = Grid_.cellj(z);
  // Determine where bivariate polynomial is defined
  if ((int)Grid_.celli(r) - i < 0) --i;
  if ((int)Grid_.cellj(z) - j < 0) --j;
  return beta[(int)i][(int)j];
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
  *dr = (-Psi_zz*Psi_r + Psi_rz*Psi_z)*(1.0/D);
  *dz = (Psi_rz*Psi_r - Psi_rr*Psi_z)*(1.0/D);
}

/*!
 * @brief Perform search for critical points beginning with initial guess r, z
 */
void Critical::Psi_critical(double r, double z, double *rcrit, double *zcrit) {
  double dr, dz, alpha, beta;
  for (int i = 0; i < max_iter; ++i) {
    // Calculate |del Psi(r,z)|
    beta = cell_beta(r,z);
    alpha = cell_alpha(r,z);
    double Psi_r = (alpha)*pow(r,alpha-1)*pow(z,beta);
    double Psi_z = (beta)*pow(z,beta-1)*pow(r,alpha);
    // If within tolerence, return
    if (sqrt(pow(Psi_r,2) + pow(Psi_z,2)) < epsilon){
      *rcrit = r;
      *zcrit = z;
      return;
    }
    Psi_search(r, z, &dr, &dz);
    r += dr;
    z += dz;
  }
  printf("Critical point not found within %d iterations \n", max_iter);
  exit(1);
}

/*!
 * @brief Performs critical point search; updates Psi_l and Psi_o
 */
void Critical::update() {
    
}


  
