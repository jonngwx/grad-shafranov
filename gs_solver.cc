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
#include "include/tsv_reader.h"
#include "include/create_options.h"

using namespace std;

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g);

int main(int argc, char *argv[]){

  po::options_description visible_options("Allowed options");
  po::variables_map vm;
  int return_value = CreateOptions(argc, argv, visible_options, vm);
  if( return_value != 0 ) {
    exit(return_value);
  }

  if (vm.count("help")) {
    std::cout << visible_options << "\n";
    return 0;
  }
  
  if (vm.count("version")) {
      std::cout << "Grad Shafranov solver, version alpha\n";
      return 0;
  }

  int nr = vm["grid-elems-r"].as<int>();
  int nz = vm["grid-elems-z"].as<int>();
  double R0 =   vm["r-min"].as<double>();
  double Rend = vm["r-max"].as<double>();
  double z0 =   vm["z-min"].as<double>();
  double zend = vm["z-max"].as<double>();
  int maxIterM = vm["max-iter-M"].as<int>();
  int maxIterN = vm["max-iter-N"].as<int>();

  std::string pgtype;
  if (!vm.count("pgtype")){
    std::cout << "Must specify pgtype: (array, choice2, choice 3)  \n";
    exit(1);
  }
  else {
    pgtype = vm["pgtype"].as<string>();
    if (pgtype == "array"){
    }
    else if (pgtype == "choice2"){}
    else if (pgtype == "choice3"){}
    else {
      std::cout << "Error: pgtype unrecognized; exiting\n";
      exit(2);
    }
  }

  PGData *gd = NewPGDataFromFile(vm["p-filename"].as<string>(),1);
  PGData *pd = NewPGDataFromFile(vm["g-filename"].as<string>(),1);
  CoilData *cd = NewCoilDataFromFile(vm["coil-data-name"].as<string>(),1);

    Grid *grid = new Grid(R0, Rend, z0, zend, nr, nz);
    Field *psi = new Field(nr,nz);
    Field *psi_prev = new Field(nr,nz);
    Field *psi_prev_prev = new Field(nr,nz);
    Field *jphi = new Field(nr,nz);
    
    RHSfunc *p = new RHSfunc(pgtype, pd);
    RHSfunc *g = new RHSfunc(pgtype, gd);

    for (int i = 0; i < grid->nr_; ++i) {
      for (int j = 0; j < grid->nz_; ++j) {
        psi->f_[i][j] = i*j;
      }
    }


    // Elliptic solver for inner loop
    double omega_init = 0.5;
    double epsilon = 0.1;
    EllipticSolver *solver = new SOR(*grid, *psi, *psi_prev, *psi_prev_prev,  omega_init, epsilon);
    Boundary *psib = new SlowBoundary(*grid, *cd);

    /** determine which output type */
    Grad_Output *grad_output = new Grad_Output_Txt(psi,grid,p,g,"this,is,a,test");

    // solve stuff
    for (int m = 0; m < maxIterM; ++m){
        calc_jphi(*grid, *jphi, *psi, *p, *g);
        psib->CalcB(*psi, *jphi); // PETER this should come after as the initial guess already has a self consistent boundary?
        // test convergence
        
        solver->coeff();
        solver->SOR_1(*jphi);
        for (int n = 1; n < maxIterN; ++n) {
            calc_jphi(*grid, *jphi, *psi, *p, *g);
            solver->step(*jphi);
            if (solver->norm() < solver->epsilon()) break;
        }
    }

    // output stuff
    std::string full_output_name = vm["output-name"].as<string>()+".txt";
    grad_output->write_output(full_output_name.c_str());

    delete grad_output;
    delete grid;
    delete psi;
    delete psib;
    delete psi_prev;
    delete psi_prev_prev;
    delete jphi;
    DeletePGData(pd);
    DeletePGData(gd);
    DeleteCoilData(cd);
}

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g){
  const double mu0 = 0.0000012566370614; // in SI units
  for (int i = 0; i < grid.nr_; ++i) {
      for (int j = 0; j < grid.nz_; ++j)  {
          jphi.f_[i][j] = -grid.R_[i]*p.eval(psi.f_[i][j]) - g.eval(psi.f_[i][j])*g.eval_prime(psi.f_[i][j]/(mu0*grid.R_[i]));
      }
  }  

  return 0;
}
