/*!
 * @file gauss_seidel.cc
 * @brief Implementation of class GaussSeidel
 */
#include <assert.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "field.h"
#include "grid.h"
#include "elliptic_solver.h"
#include "gauss_seidel.h"

GaussSeidel::GaussSeidel(const Grid &GridS, Field &Psi, double error_ES) :
  EllipticSolver(GridS, Psi),
  error_(error_ES) { }

GaussSeidel::~GaussSeidel(){}

void GaussSeidel::coeff() {
    double dr = Grid_.dr_;
    double dz = Grid_.dz_;
    B = -(dr*dr*dz*dz)/(2*dz*dz + 2*dr*dr);
    D = -1/(dz*dz);
    for (int i = 0; i < Grid_.nr_; ++i) {
        A[i] = 1/(2*Grid_.R_[i]*dr) - 1/(dr*dr);
        C[i] = -1/(2*Grid_.R_[i]*dr) - 1/(dr*dr);
    }
}

void GaussSeidel::step_1(const Field &jphi){
  step(jphi);
}

void GaussSeidel::step(const Field &jphi){
  const double mu0 = 4 * M_PI * 1e-7; // in SI units
  double nr = Grid_.nr_;
  double nz = Grid_.nz_;
  
  // Copy psi to psi_prev
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      Psi_prev_.f_[i][j] = Psi_.f_[i][j];
    }
  }

  // Gauss seidel algorithm
  for (int k = 0; k < 1000; ++k) {
    // Save Psi_ to Psi_temp
    for (int i = 0; i < nr; ++i) {
      for(int j = 0; j < nz; ++j) {
        Psi_temp_.f_[i][j] = Psi_.f_[i][j];
      }
    }
    // Copy over boundary values
    //  boundary(Psi_prev_, Psi_);
    for (int i = 1; i < nr-1; ++i) {
      for (int j = 1;j < nz-1; ++j) {
        Psi_.f_[i][j] = B*(jphi.f_[i][j]*mu0*Grid_.R_[i] + A[i]*Psi_temp_.f_[i+1][j] + C[i]*Psi_.f_[i-1][j] + D*Psi_.f_[i][j-1] + D*Psi_temp_.f_[i][j+1]);
      }
    } 
    // Check convergence
    double sum = 0.0;
    for(int i = 0; i < nr; ++i) {
      for(int j = 0; j < nz; ++j) {
        sum += (Psi_.f_[i][j]-Psi_temp_.f_[i][j])*(Psi_.f_[i][j]-Psi_temp_.f_[i][j]);
      }
    }
    if (sqrt(sum) < error_) {break;}
    if (k == 999) {printf("Gauss-Seidel algorithm reached end without converging\n");}  
  }  // end Gauss seidel algorithm
  //iter(0.5);
}

