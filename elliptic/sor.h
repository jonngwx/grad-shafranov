#ifndef SOR_H
#define SOR_H

#include elliptic-solver.h

class SOR: public EllipticSolver {
public:
  SOR(const Grid &GridS, double omega_init, double epsilon);
  ~SOR();
// For first iteration - use Gauss Seidel with blending
  void SOR_1();
// Iterate
  void step();
// Calculate over-relaxation parameter
  double omega();
private:
  Field &Psi_prev_prev_;
  const double omega_init_;
// Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **e;
  const double **f;
};

#endif