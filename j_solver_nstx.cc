/*!
 * @file j_solver_nstx.cc
 * @brief Implementation of JSolverNSTX
 */
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "j_solver_nstx.h"

JSolverNSTX::JSolverNSTX(double P0, double g0, double Ip, double n2, Grid *grid)
    : n2_(n2) {  // enforce boundary conditions on P
  P0_ = P0;
  g0_ = g0;
  Ip_ = Ip;
  R_ = grid->R_;
  dr_ = grid->dr_;
  dz_ = grid->dz_;
  nr_ = grid->nr_;
  nz_ = grid->nz_;
  P1_ = 6 * P0;
}

JSolverNSTX::~JSolverNSTX() {}

void JSolverNSTX::update(Field *jphi, Field *psi, Field *p, Field *g) {
  const double mu0 = 4 * M_PI * 1e-7;  // the permeability of free space
  double temp1 = 0;
  double temp2 = 0;
  double psi_s = 0; /* "psi-squiggle": a normalized psi equal to
  0 at the plasma-vacuum boundary and 1 on the magnetic axis.*/

  double delta_psi = psi->f_l - psi->f_0;
  double jtot_old = 0;
  // calc temps, then alpha g
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0; j < nz_; ++j) {
      jtot_old += jphi->f_[i][j] * dr_ * dz_;
      psi_s = (psi->f_l - psi->f_[i][j]) / delta_psi;
      /* If this is a point inside the plasma */
      if (psi_s > 0 && psi_s <= 2) {
        temp1 += R_[i] * P1_ * psi_s * (1 - psi_s);
        temp2 += pow(psi_s, n2_) / R_[i];
        p->f_[i][j] = P1_ * (pow(psi_s, 2) / 2 -
                             pow(psi_s, 3) / 3);  // update pressure field
      } else {
        p->f_[i][j] = 0;
      }
    }
  }
  double alpha_g = mu0 * (-temp1 + Ip_ * delta_psi / (dr_ * dz_)) / (temp2);
  // ggprime = -0.0140 (1-psi)^2 + .0258 (1-psi) +- .01114;
  // update fields g, and jphi
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0; j < nz_; ++j) {
      psi_s = (psi->f_l - psi->f_[i][j]) / delta_psi;
      // double ompsi = 1 - psi_s;
      if (psi_s > 0 && psi_s <= 1) {

        g->f_[i][j] = sqrt(g0_ * g0_ + 2. / n2_ * alpha_g *
                                           pow(psi_s, n2_));  // update g field
        //  double ggprime = -.0140*pow(ompsi,2) + .00258*ompsi -.01114;
        double ggprime = alpha_g * pow(psi_s, n2_);
        jphi->f_[i][j] =
            (P1_ * psi_s * (1 - psi_s) * R_[i] + 1 / mu0 / R_[i] * ggprime) /
            delta_psi;

      } else {
        g->f_[i][j] = g0_;
        jphi->f_[i][j] = 0;
      }
    }
  }

  double jtot = 0;
  for (int i = 0; i < nr_; ++i) {
    for (int j = 0; j < nz_; ++j) {
      jtot += jphi->f_[i][j];
    }
  }
  jtot *= (dr_ * dz_);
  // if (abs(Ip_ - jtot) >1) {
  //   printf("Ip = %f . summed jphi = %f\n, a_g = %f, dpsi = %f, psil = %f, psi0
  //     = %f, temp1 = %f\n", Ip_, jtot, alpha_g, delta_psi, psi->f_l, psi->f_0,
  //     temp1);
  // }
}
