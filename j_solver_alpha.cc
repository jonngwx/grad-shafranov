#include "include/j_solver_alpha.h"

JSolverAlpha::JSolverAlpha(double P0, double g0, double n1, double n2, double Ip, Grid *grid)
    : P0_(P0),
      g0_(g0),
      n1_(n1),
      n2_(n2),
      Ip_(Ip),
      R_(grid->R_),
      dr_(grid->dr_),
      dz_(grid->dz_)

{}

JSolverAlpha::~JSolverAlpha()
{}

void JSolverAlpha::update(Field *jphi, Field *psi) {



}

 
