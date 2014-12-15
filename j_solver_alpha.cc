#include "include/j_solver_alpha.h"
#include <stdio.h>
#include <math.h>

JSolverAlpha::JSolverAlpha(double P0, double g0, double n1, double n2, double Ip, Grid *grid)
    : P0_(P0),
      g0_(g0),
      n1_(n1),
      n2_(n2),
      Ip_(Ip),
      R_(grid->R_),
      dr_(grid->dr_),
      dz_(grid->dz_),
      nr_(grid->nr_),
      nz_(grid->nz_)
{}

JSolverAlpha::~JSolverAlpha()
{}

void JSolverAlpha::update(Field *jphi, Field *psi, Field *p, Field *g) {
    const double mu0 = 0.0000012566370614;
    double temp1 = 0;
    double temp2 = 0;
    double psi_s = 0; // "psi-squiggle"
    double delta_psi = psi->f_l - psi->f_0;

    // calc temps, then alpha g
    for (int i=0; i < nr_; ++i) {
        for (int j=0; j < nz_; j++) {
            psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
            temp1 += R_[i]*n1_*pow(psi_s, n1_ - 1);
            temp2 += n2_*pow(psi_s, n2_ - 1) / R_[i];
            p->f_[i][j] = P0_*pow(psi_s,n1_); //update pressure field
        }
    }
    double alpha_g = mu0*(-P0_*temp1 + Ip_*delta_psi/(dr_*dz_))/(0.5*g0_*g0_*temp2);

    // update fields g, and jphi
    for (int i=0; i < nr_; ++i) {
        for (int j=0; j < nz_; ++j) {
            psi_s = (psi->f_l - psi->f_[i][j])/delta_psi;
            g->f_[i][j] = g0_*sqrt(1 + alpha_g*pow(psi_s,n2_)); //update g field
            jphi->f_[i][j] = P0_*temp1/delta_psi + g0_*g0_*alpha_g*temp2/(2*mu0*delta_psi); //update jphi
        }
    }

}

 
