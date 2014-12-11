#include "elliptic_solver.h"
#include <math.h>
#include "field.h"
#include <stdio.h>
#include <vector>

/*!
 * @file elliptic_solver.cc
 * @brief Base class implementation of EllipticSolver
 * @section DESCRIPTION
 */

EllipticSolver::EllipticSolver(const Grid &Grid, Field &Psi) :
  Grid_(Grid),
  Psi_(Psi),
  Psi_prev_(Grid),
  A(Grid.nr_){}

EllipticSolver::~EllipticSolver(){}

/*!
 * @brief Calculate maximum of |Psi - Psi_prev| over grid
 * @param Psi current value of Psi
 * @param Psi_prev previous value of Psi
 */
double EllipticSolver::norm_max(const Field &Psi, const Field &Psi_prev) {
  double max = 0;
  int nr = Grid_.nr_;
  int nz = Grid_.nz_;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0 ;j < nz; ++j) {
      if (abs(Psi.f_[i][j] - Psi_prev.f_[i][j]) > max)
        max = abs(Psi.f_[i][j] - Psi_prev.f_[i][j]);
    }
  }
  return max;
}

/*!
 * @brief Calculates 2-norm of diffence between Psi and Psi_prev for convergence testing
 */
double EllipticSolver::norm() {
  double sum = 0;
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  for (int i = 0; i < nr; ++i) {
    for(int j = 0; j < nz; ++j) {
      sum += (Psi_.f_[i][j]-Psi_prev_.f_[i][j])*(Psi_.f_[i][j]-Psi_prev_.f_[i][j]);
    }
  }
  return sqrt(sum);
}

/*!
 * @brief Returns norm of residuals for convergence testing
 * @param Psi current value of Psi
 * @param Psi_prev previous value of Psi
 */
double EllipticSolver::residuals(const Field &Psi, const Field &Psi_prev) {
  return 0;
}

/*!
 * @brief Blend with old solution
 * @param omega blending parameter used for iteration
 */
void EllipticSolver::iter(double omega) {
  int nr = Grid_.nr_;
  int nz = Grid_.nz_;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0 ;j < nz; ++j) {
      Psi_.f_[i][j] = omega*Psi_.f_[i][j] + (1-omega)*Psi_prev_.f_[i][j];
    }
  }
}

/*!
 * @brief Enforce boundary condition by copying boundary values from
 * Psi to Psi_prev
 * @param Psi current value of Psi
 * @param Psi_prev previous value of Psi
 */
void EllipticSolver::boundary(Field &Psi, const Field &Psi_prev) {
  int nr = Grid_.nr_;
  int nz = Grid_.nz_;
  for (int i = 0; i < nr; ++i) {
    Psi.f_[i][0] = Psi_prev.f_[i][0];
    Psi.f_[i][nz-1] = Psi_prev.f_[i][nz-1];
  }
  for (int i = 0; i < nz; ++i) {
    Psi.f_[0][i] = Psi_prev.f_[0][i];
    Psi.f_[nr-1][i] = Psi_prev.f_[nr-1][0];
  }
}
