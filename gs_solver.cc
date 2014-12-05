#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "include/grid.h"
#include "include/field.h"
#include "include/rhs_func.h"
#include "include/grad_output.h"
#include "include/grad_output_txt.h"
#include "include/sor.h"
#include "include/slow_boundary.h"

#include "include/create_options.h"

using namespace std;

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g);

int main(int argc, char *argv[]){
    // JACOB, I'm going to manually set some constants right now so we can test everything else  :)
    int nr = 11;
    int nz = 11;
    double R0 = 5.0;
    double Rend = 10.0;
    double z0 = -3.0;
    double zend = 3.0;
    int maxIterM = 100;
    int maxIterN = 100;   
 
    PGData *pgd = NewPGDataFromFile("pg_data.tsv",1);
    CoilData *cd = NewCoilDataFromFile("coil_data.tsv",1);
    string pgtype = "array";
    
    Grid *grid = new Grid(R0, Rend, z0, zend, nr, nz);
    Field *psi = new Field(nr,nz);
    Field *psi_prev = new Field(nr,nz);
    Field *psi_next = new Field(nr,nz);
    Field *jphi = new Field(nr,nz);
    
    RHSfunc *p = new RHSfunc(pgtype, pgd);
    RHSfunc *g = new RHSfunc(pgtype, pgd);

    for (int i = 0; i < grid->nr_; ++i) {
      for (int j = 0; j < grid->nz_; ++j) {
        psi->f_[i][j] = i*j;
      }
    }


    // Elliptic solver for inner loop
//    EllipticSolver *solver = new SOR(grid, omega_init, epsilon);
    Boundary *psib = new SlowBoundary(*grid, *cd);    

    /** determine which output type */
    Grad_Output *grad_output = new Grad_Output_Txt(psi,grid,p,g,"this,is,a,test");

    // solve stuff
    for (int m = 0; m < maxIterM; ++m){
        calc_jphi(*grid, *jphi, *psi, *p, *g);
        psib->CalcB(*psi, *jphi); // PETER this should come after as the initial guess already has a self consistent boundary?
        // test convergence
        
//        solver->init(psi);
        for (int n = 0; n < maxIterN; ++n) {
            if (n != 0) calc_jphi(*grid, *jphi, *psi, *p, *g);
//            psi = solver->step(*jphi);
//            if (solver->norm() < solver->epsilon()) break;
        }
    }

    // output stuff
    
    grad_output->write_output("whatever");

    delete grad_output;
    delete grid;
    delete psi;
    delete psib;
//    delete solver;
    delete jphi;
}

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g){

  for (int i = 0; i < grid.nr_; ++i) {
      for (int j = 0; j < grid.nz_; ++j)  {
          jphi.f_[i][j] = p.eval(psi.f_[i][j]) + g.eval_prime(psi.f_[i][j]); // not the right formula, but you get the idea
      }
  }  

  return 0;
}
