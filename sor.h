#ifndef SOR_H
#define SOR_H

#include sor.h

class SOR: public EllipticSolver {
public:
  SOR(const Field &Psi, double omega_init, int max_iter, double epsilon);
  ~SOR();
// For first iteration - use Gauss Seidel with blending
  void SOR_1();
// Iterate
  void step(const Field &Psi_n, const Field &Psi_n+);
// Calculate over-relaxation parameter
  double omega(const Field &Psi_n-, const Field &Psi_n, const Field &Psi_n+);
private:
  const double omega_init_;
// Coefficient arrays
  const double **a;
  const double **b;
  const double **c;
  const double **d;
  const double **e;
};

#endif