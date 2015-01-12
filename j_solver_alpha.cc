#include "include/j_solver_alpha.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

J_Solver_Alpha::J_Solver_Alpha(double P0, double g0, double n1, double n2, double Ip, Grid *grid)
    : n1_(n1),
      n2_(n2)
{   
    P0_ = P0;
    g0_ = g0;
    Ip_ = Ip;
    R_ = grid->R_;
    dr_ = grid->dr_;
    dz_ = grid->dz_;
    nr_ = grid->nr_;
    nz_ = grid->nz_;
}

J_Solver_Alpha::~J_Solver_Alpha() {}

void J_Solver_Alpha::update(Field *jphi, Field *psi, Field *p, Field *g) {
  const double mu0 = 4 * M_PI * 1e-7; // the permeability of free space
  double temp1 = 0;
  double temp2 = 0;
  double psi_s = 0; /* "psi-squiggle": a normalized psi equal to
  0 at the plasma-vacuum boundary and 1 on the magnetic axis.*/

  //    printf("psi_l = %f, psi_0 = %f \n", psi->f_l, psi->f_0);
  /* delta_psi: The difference between the value of psi at the
   * plasma-vacuum boundary and the value on axis. */
  double delta_psi = psi->f_l - psi->f_0;
 
  double jtot_old=0;
  // calc temps, then alpha g
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j < nz_; ++j) {
      jtot_old += jphi->f_[i][j] * dr_ * dz_;
      psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
      /* If this is a point inside the plasma */ 
      if (psi_s > 0) {
          temp1 += R_[i]*n1_*pow(psi_s, n1_ - 1.0);
          temp2 += n2_*pow(psi_s, n2_ - 1.0) / R_[i];
          p->f_[i][j] = P0_*pow(psi_s,n1_); //update pressure field
      } else {
          p->f_[i][j] = 0;
      }
      // printf("psi = %f .... p = %f \n", psi->f_[i][j], p->f_[i][j]);
    }
  }
  double alpha_g = mu0*(-P0_*temp1 + Ip_*delta_psi/(dr_*dz_))/(0.5*g0_*g0_*temp2);
  // printf("alpha_g = %f \n", alpha_g);
  // printf("summed jphi before update = %f\n", jtot_old);

  // update fields g, and jphi
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j < nz_; ++j) {
      psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
      if (psi_s > 0) {
        
        g->f_[i][j] = g0_*sqrt(1 + alpha_g*pow(psi_s,n2_)); //update g field
        // printf("g = %f \n", g->f_[i][j]);
        // jphi = - p'*R - gg'/mu0/R
        // if gg' = .5 g0^2 alphag g_s'
        jphi->f_[i][j] = (n1_*P0_*pow(psi_s,n1_-1)*R_[i] + .5*g0_*g0_ / (mu0 * R_[i]) * alpha_g*pow(psi_s,n2_-1)*n2_)/delta_psi;

        //P0_*temp1/delta_psi + g0_*g0_*alpha_g*temp2/(2*mu0*delta_psi); //update jphi
      } else {
          g->f_[i][j] = g0_;
          jphi->f_[i][j] = 0;
      }
    }
  }

  // check that the sum of jphi is specified Ip (temporary)
  double jtot=0;
  for (int i=0; i < nr_; ++i) {
    for (int j=0; j< nz_; ++j) {
      jtot += jphi->f_[i][j];  
    }
  }
  jtot *= (dr_*dz_);
  //    if (abs(Ip_ - jtot) >1) {
  //  printf("Ip = %f . summed jphi = %f\n, a_g = %f, dpsi = %f, psil = %f, psi0 = %f, temp2 = %f\n", Ip_, jtot, alpha_g, delta_psi, psi->f_l, psi->f_0, temp2);
  //    }

  //    assert(abs(Ip_-jtot) < 1);

}
