#include <stdlib.h>
#include <stdio.h>
#include "field.h"
#include "grad_output.h"
#include "tsv_data.h"

int main(int argc, char *argv[]){
    // read inputs TO BE FILLED IN BY JACOB
    std::string filename = "";
    TsvData coils = NewTsvDataFromFile(filename); // is there an implicit alloc here?
    // make field data
    Field *Psi = new Field(R0, Rend, z0, zend, nr,nz);
    //RHS rhs = something

    // make solvers
    // EllipticSolver () = ???
    Boundary *Psib = new SlowBoundary(nr,nz,Psi->dr, Psi->dz, coils);    

    /** determine which output type */
    Grad_Output grad_output = new Grad_Output_Txt(Psi,Psi,Psi,"what to write");
    // solve stuff
    for (int i = 0; i < maxIterM,++i){
        Psib->CalcB(Psi); // PETER this should come after as the initial guess already has a self consistent boundary?
        // test convergence
        for (int j = 0; j < maxIterN,++j){
            // Elliptic_solver->solve(Psi)
            // test convergence
        }
    }

    // output stuff

    grad_output.write_output("whatever");

    delete grad_output;
    delete Psi;
    delete Psib;
}
