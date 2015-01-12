#include "j_solver_nstx.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

J_Solver_NSTX::J_Solver_NSTX(double P0, double g0, double Ip, Grid *grid)
{ // enforce boundary conditions on P
  P0_ = P0;
  g0_ = g0;
  Ip_ = Ip;
  R_ = grid->R_;
  dr_ = grid->dr_;
  dz_ = grid->dz_;
  nr_ = grid->nr_;
  nz_ = grid->nz_;
  P1_ = -6*P0;
}

J_Solver_NSTX::~J_Solver_NSTX() {}

void J_Solver_NSTX::update(Field *jphi, Field *psi, Field *p, Field *g) {
  const double mu0 = 4 * M_PI * 1e-7; // the permeability of free space
  double temp1 = 0;
  double temp2 = 0;
  double psi_s = 0; /* "psi-squiggle": a normalized psi equal to
  0 at the plasma-vacuum boundary and 1 on the magnetic axis.*/

  double delta_psi = psi->f_l - psi->f_0;
 
  double jtot_old=0;
  // calc temps, then alpha g
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j < nz_; ++j) {
      jtot_old += jphi->f_[i][j] * dr_ * dz_;
      psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
      /* If this is a point inside the plasma */ 
      if (psi_s > 0 && psi_s <= 1) {
        temp1 += R_[i]*P1_*psi_s*(1-psi_s);
        temp2 += psi_s*psi_s / R_[i];
        p->f_[i][j] = P0_ + P1_*(pow(1-psi_s,2)/2 - pow(1-psi_s,3)/3); //update pressure field
      } else {
          p->f_[i][j] = 0;
      }
    }
  }
  double alpha_g = mu0*(-temp1 + Ip_*delta_psi/(dr_*dz_))/(temp2);

  // update fields g, and jphi
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j < nz_; ++j) {
      psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
      if (psi_s > 0 && psi_s <=1) {
        
        g->f_[i][j] = sqrt(g0_*g0_ + 2./3.*alpha_g*pow(psi_s,3)); //update g field
        jphi->f_[i][j] = (P1_*psi_s*(1-psi_s)*R_[i] + 1/mu0/R_[i]*alpha_g*psi_s*psi_s)/delta_psi;

        //P0_*temp1/delta_psi + g0_*g0_*alpha_g*temp2/(2*mu0*delta_psi); //update jphi
      } else {
          g->f_[i][j] = g0_;
          jphi->f_[i][j] = 0;
      }
    }
  }

  double jtot=0;
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j< nz_; ++j) {
      jtot += jphi->f_[i][j];  
    }
  }
  jtot *= (dr_*dz_);
  // if (abs(Ip_ - jtot) >1) {
  printf("Ip = %f . summed jphi = %f\n, a_g = %f, dpsi = %f, psil = %f, psi0 = %f, temp2 = %f\n", Ip_, jtot, alpha_g, delta_psi, psi->f_l, psi->f_0, temp2);
  // }
}
