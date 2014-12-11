#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Elliptic
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/grid.h"
#include "../include/field.h"
#include "../include/elliptic_solver.h"
#include "../include/sor.h"
#include "../include/gauss_seidel.h"
#include <math.h>

/*!
 * Tests SOR elliptic solver by comparison to vacuum 
 * solution constrained to be even about the Z = 0
 * midplane
 */
BOOST_AUTO_TEST_CASE (SOR_vacuum)
{
  double R0 = 0;
  double Rend = 10;
  double Z0 = -5;
  double Zend = 5;
  double nr = 100;
  double nz = 100;
  double omega_init = 0.5;
  double max_iter = 20;
  double epsilon = 0.1;
  Grid *grid = new Grid(R0, Rend, Z0, Zend, nr, nz);
  Field *psi = new Field(*grid);
  Field *jphi = new Field(*grid);
  
// Initialize psi and jphi
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      psi->f_[i][j] = 1;
      jphi->f_[i][j] = 0;
    }
  }
  
  EllipticSolver *solver = new SOR(*grid, *psi, omega_init);
  solver->coeff();
  for (int n = 0; n < max_iter; ++n) {
    if (n == 0) solver->step_1(*jphi);
    else solver->step(*jphi);
    if (solver->norm() < epsilon) break;
  }

// Vacuum solution multipole expansion
  Field *psi_sol = new Field(*grid);
  double Rc = 5;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      double R = grid->R_[i];
      double Z = grid->z_[j];
      psi_sol->f_[i][j] = Rc*Rc + 0.5*(R*R - Rc*Rc);
      psi_sol->f_[i][j] += 1/(8*Rc*Rc)*((R*R-Rc*Rc)*(R*R-Rc*Rc) - 4*R*R*Z*Z);
      psi_sol->f_[i][j] += 1/(24*pow(Rc,4))*(pow((R*R-Rc*Rc),3) - 12*R*R*Z*Z*(R*R - Rc*Rc) + 8*R*R*pow(Z,4));
    }
  }
  BOOST_CHECK_CLOSE(psi_sol->f_[50][50], psi->f_[50][50], 0.1);
}

/*!
 * Tests Gauss-Seidel elliptic solver by comparison to vacuum
 * solution constrained to be even about the Z = 0
 * midplane
 */
BOOST_AUTO_TEST_CASE (GS_vacuum)
{
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 100;
    double nz = 100;
    double max_iter = 20;
    double epsilon = 0.1;
    Grid *grid = new Grid(R0, Rend, Z0, Zend, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);
    
    // Initialize psi and jphi
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            psi->f_[i][j] = 1;
            jphi->f_[i][j] = 0;
        }
    }
    
    EllipticSolver *solver = new GaussSeidel(*grid, *psi);
    solver->coeff();
    for (int n = 0; n < max_iter; ++n) {
        if (n == 0) solver->step_1(*jphi);
        else solver->step(*jphi);
        if (solver->norm() < epsilon) break;
    }
    
    // Vacuum solution multipole expansion
    Field *psi_sol = new Field(*grid);
    double Rc = 5;
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            double R = grid->R_[i];
            double Z = grid->z_[j];
            psi_sol->f_[i][j] = Rc*Rc + 0.5*(R*R - Rc*Rc);
            psi_sol->f_[i][j] += 1/(8*Rc*Rc)*((R*R-Rc*Rc)*(R*R-Rc*Rc) - 4*R*R*Z*Z);
            psi_sol->f_[i][j] += 1/(24*pow(Rc,4))*(pow((R*R-Rc*Rc),3) - 12*R*R*Z*Z*(R*R - Rc*Rc) + 8*R*R*pow(Z,4));
        }
    }
    BOOST_CHECK_CLOSE(psi_sol->f_[50][50], psi->f_[50][50], 0.1);
}