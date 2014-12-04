#include "critical.h"
#include <math.h>

Critical::Critical(const Grid &GridS, const Field &Psi) :
  Psi(Psi),
  GridS(GridS) {
    Psi_r = new double[GridS->nr];
    Psi_z = new double[GridS->nr];
    Psi_rr = new double[GridS->nr];
    Psi_zz = new double[GridS->nr];
    Psi_rz = new double[GridS->nr];
    alpha = new double[Grid->nr];
    beta = new double[Grid->nr];
    for(int i = 0; i < GridS->nr; ++i) {
      Psi_r[i] = new double[GridS->nz];
      Psi_z[i] = new double[GridS->nz];
      Psi_rr[i] = new double[GridS->nz];
      Psi_zz[i] = new double[GridS->nz];
      Psi_rz[i] = new double[GridS->nz];
      alpha[i] = new double[GridS->nz];
      beta[i] = new double[GridS->nz];
    }
}

void Critical::~Critical() {
  for (int i = 0; i < GridS->nr; ++i) {
    delete [] Psi_r[i];
    delete [] Psi_z[i];
    delete [] Psi_rr[i];
    delete [] Psi_zz[i];
    delete [] Psi_rz[i];
    delete [] alpha[i];
    delete [] beta[i];
  }
  delete [] Psi_r;
  delete [] Psi_z;
  delete [] Psi_rr;
  delete [] Psi_zz;
  delete [] Psi_rz;
  delete [] alpha;
  delete [] beta;
}

// Bivariate interpolation of Psi s.t. preserves smoothness
// as described in Akima, 1974
// Calculates Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz
// Determines Psi = r^alpha z^beta inside rectangle defined by
// (r[i], r[i+1])x(z[i], z[i+1])
void Critical::interpolate() {
  double c1, c2, c3, d1, d2, d3, w_r1, w_r2, w_z1, w_z2, e1, e2;
  
  for (int i = 1; i < nr-1; ++i) {
    for(int j = 1; j < nz-1; ++j) {
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
      // D[i][j-1]
      d1 = (Psi[i][j] - Psi[i][j-1])/(z[j]-z[j-1]);
      // D[i][j]
      d2 = (Psi[i][j+1] - Psi[i][j])/(z[j+1]-z[j]);
      // D[i+1][j]
      d3 = (Psi[i+1][j+1] - Psi[i+1][j])/(z[j+1]-z[j]);
      w_r1 = abs(c[i+1][j] - c[i][j]);
      w_r2 = abs(c[i-1][j] - c[i-2][j]);
      if (w_r1 == w_r2 && w_r1 == 0) {
        w_r1 = 1;
        w_r2 = 1;
      }
      w_z1 = abs(d[i][j+1] - d[i][j]);
      w_z2 = abs(d[i][j-1] - d[i][j-2]);
      if (w_z1 == w_z2 && w_z1 == 0) {
        w_z1 = 1;
        w_z2 = 1;
      }
      Psi_r[i][j] = (w_r1*c1 + w_r2*c2)/(w_r1 + w_r2);
      Psi_z[i][j] = (w_z1*d1 + w_z2*d2)/(w_z1 + w_z2);
      // e[i][j]
      e1 = (c3 - c2)/(z[j+1]-z[j]);
      // e[i-1][j-1] = (c[i-1][j] - c[i-1][j-1])/(z[j]-z[j-1])
      e2 = (c1 - c4)/(z[j]-z[j-1]);
      // e[i][j-1] = (c[i][j] - c[i][j-1])/(z[j-1]-z[j-2])
      e3 = (c2 - c5)/(z[j-1] - z[j-2]);
      // e[i-1][j] = (c[i-1][j+1] - c[i-1][j])/(z[j+1] - z[j])
      e4 = (c6 - c1)/(z[j+1] - z[j]);
      Psi_rz[i][j] = (w_r1*(w_z1*e2 + w_z2*e4) + w_r2*(w_z1*e3 + w_z2*e1))/((w_r1 + w_r2)(w_z1 + w_z2));
      alpha[i][j] = Psi_rz[i][j]*r[i][j]/Psi_z[i][j];
      beta[i][j] = Psi_rz[i][j]*z[i][j]/Psi_r[i][j];
    }
  }
}

// Interpolated Psi defined for all r, z
void Psi_interp(double r, double z) {
  double a, b;
  // Find containing cell's indices
  int i = GridS->celli(r, &a);
  int j = GridS->cellj(z, &b);
  // Determine where bivariate polynomial is defined
  if (a < i) --i;
  if (b < j) --j;
  double alpha = alpha[i][j];
  double beta = beta[i][j];
  // Psi = r^alpha z^beta
  return (r^alpha)*(z^beta);
}
  
  
