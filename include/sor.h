#ifndef SOR_H
#define SOR_H

#include "elliptic_solver.h"
#include "field.h"
#include "grid.h"

class SOR : public EllipticSolver {
public:
  SOR(const Grid &GridS, const Field &Psi, const Field &Psi_prev, const Field &Psi_prev_prev, double omega_init, double epsilon);
  ~SOR();
// For first iteration - use Gauss Seidel with blending
  void SOR_1(const Field &jphi);
// Calculate coefficients for iteration from grid parameters
  void coeff();
// Iterate
  void step(const Field &jphi);
  void boundary(const Field &Psi, const Field &Psi_next);
// Calculate over-relaxation parameter
  double omega();
  double epsilon();
  double norm();
  void iter(double omega);
  double norm_max(const Field &Psi, const Field &Psi_prev);
private:
    const Grid &Grid_;
    const Field &Psi_;
    const Field &Psi_prev_;
    const Field &Psi_prev_prev_;
    const int nr_;
    const int nz_;
    const double epsilon_;
  const double omega_init_;
// Coefficient arrays
  double **a;
  double **b;
  double **c;
  double **d;
  double **f;
};

#endif
