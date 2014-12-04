#include <stdlib.h>
#include <stdio.h>
#include "field.h"
#include "grad_output.h"
#include "tsv_data.h"
#include "elliptic/sor.h"

int main(int argc, char *argv[]){
    // read inputs TO BE FILLED IN BY JACOB
    std::string filename = "";
    TsvData coils = NewTsvDataFromFile(filename); // is there an implicit alloc here?
    // make field data
    int nr,nz;
    double R0, Rend, z0, zend;
    Grid *grid = new Grid(R0, Rend, z0, zend, nr, nz);
    Field *Psi = new Field(nr,nz);
    Field *Psi_prev = new Field(nr,nz);
    Field *Psi_next = new Field(nr,nz);
    //RHS rhs = something

    // Elliptic solver for inner loop
    EllipticSolver *solver = new SOR(grid, omega_init, epsilon);
    Boundary *Psib = new SlowBoundary(nr,nz,Psi->dr, Psi->dz, coils);    

    /** determine which output type */
    Grad_Output grad_output = new Grad_Output_Txt(Psi,Psi,Psi,"what to write");
    // solve stuff
    for (int i = 0; i < maxIterM, ++i){
        Psib->CalcB(Psi); // PETER this should come after as the initial guess already has a self consistent boundary?
        // test convergence
        
        solver->init(*Psi);
        for (int j = 0; j < maxIterN, ++j) {
            Psi = solver->step();
            if (solver->norm() < solver->epsilon()) break;
        }
    }

    // output stuff

    grad_output.write_output("whatever");

    delete grad_output;
    delete grid;
    delete Psi;
    delete Psib;
    delete solver;
}
