#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include "include/grid.h"
#include "include/field.h"
#include "include/rhs_func.h"
#include "include/grad_output.h"
#include "include/grad_output_txt.h"
#include "include/sor.h"
#include "include/gauss_seidel.h"
#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/create_options.h"

#ifdef HDF_MODE
#include "grad_output_hdf.h"
#endif 

using namespace std;

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g);

int main(int argc, char *argv[])
{

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
    double Rmin =   vm["r-min"].as<double>();
    double Rmax = vm["r-max"].as<double>();
    double zmin =   vm["z-min"].as<double>();
    double zmax = vm["z-max"].as<double>();
    int maxIterM = vm["max-iter-M"].as<int>();
    int maxIterN = vm["max-iter-N"].as<int>();
    double epsilon = vm["error-tol-N"].as<double>();    

    std::string pgtype;
    if (!vm.count("pgtype")) {
        std::cout << "Must specify pgtype: (array, choice2, choice 3)  \n";
        exit(1);
    } else {
        pgtype = vm["pgtype"].as<string>();
        if (pgtype == "array") {
        } else if (pgtype == "choice2") {}
        else if (pgtype == "choice3") {}
        else {
            std::cout << "Error: pgtype unrecognized; exiting\n";
            exit(2);
        }
    }
    std::string output_type = vm["output-type"].as<string>();
  
    PGData gd; gd.load_from_tsv(vm["p-filename"].as<string>(),1);
    PGData pd; gd.load_from_tsv(vm["g-filename"].as<string>(),1);
    CoilData cd; cd.load_from_tsv(vm["coil-data-name"].as<string>(),1);
  
    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);

    RHSfunc *p = new RHSfunc(pgtype, &pd);
    RHSfunc *g = new RHSfunc(pgtype, &gd);
    
    double r_squared;
    double R0 = Rmin + (Rmax - Rmin)/2.0; // not necessarily true.  just for now
    double z0 = zmin + (zmax - zmin)/2.0; // see above comment
    double D = vm["j-phi-D"].as<double>();
    double Ip = vm["j-phi-Ip"].as<double>();

    // calc constant c to be consistent with specified Ip
    double a = pow(D/R0,2); // just for convenience.  should be less than 1
    double c = Ip*a/(4.188790205*R0*(pow(1-a,1.5) - (1 - 1.5*a)));
    printf("c = %f \n", c);

    for (int i = 0; i < grid->nr_; ++i) {
        for (int j = 0; j < grid->nz_; ++j) {
            r_squared = (grid->R_[i] - R0)*(grid->R_[i] - R0) + (grid->z_[j] - z0)*(grid->z_[j] - z0);
            if (r_squared < D*D) {
                jphi->f_[i][j] = (c/grid->R_[i])*(1 - r_squared/(D*D)) ;
            }
            else {
                jphi->f_[i][j] = 0;
            }  
        }
    }

    // Elliptic solver for inner loop
    double omega_init = 0.5;
    EllipticSolver *solver = new GaussSeidel(*grid, *psi);
    Boundary *psib = new SlowBoundary(grid, &cd);

    /** determine which output type */
    Grad_Output *grad_output;
    if (output_type == "tsv") {
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,"j,bt,bt,this,is,a,test");
    } else if (output_type == "hdf5") {
        #ifdef HDF_MODE
        grad_output = new Grad_Output_Hdf(psi,jphi,grid,p,g,"j,bt,this,is,a,test");
        #else 
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,"this,is,a,test");
        printf("Output type hdf not supported. Recompile with hdf libraries to enable. Defaulting to tsv. \n");
        #endif 
    } else {
        printf("Output type %s is not supported, use tsv or hdf5. Defaulting to tsv. \n", output_type.c_str());
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,"this,is,a,test");
    }

    // solve stuff
    solver->coeff();
    for (int m = 0; m < maxIterM; ++m) {
        psib->CalcB(psi, jphi);
        // test convergence

        // Iterate through elliptic solver
        for (int n = 0; n < maxIterN; ++n) {
            printf("n = %i \n", n);
            if (n == 0) solver->step_1(*jphi);
            solver->step(*jphi);
            calc_jphi(*grid, *jphi, *psi, *p, *g);
            printf("error norm = %f \n", solver->norm());
            if (solver->norm() < epsilon) break;
            if (n == maxIterN-1) {
                printf(" Elliptic solver reached maxIterN without convergence\n");
            }
        }
    }

  // output stuff
  std::string full_output_name = vm["output-name"].as<string>()+"."+output_type;
  grad_output->write_output(full_output_name.c_str());

  delete grad_output;
  delete grid;
  delete psi;
  delete psib;
  delete jphi;
  delete solver;
}

int calc_jphi(Grid &grid, Field &jphi, Field &psi, RHSfunc &p, RHSfunc &g)
{
    const double mu0 = 0.0000012566370614; // in SI units
    for (int i = 0; i < grid.nr_; ++i) {
        for (int j = 0; j < grid.nz_; ++j)  {
            jphi.f_[i][j] = -grid.R_[i]*p.eval(psi.f_[i][j]) - g.eval(psi.f_[i][j])*g.eval_prime(psi.f_[i][j]/(mu0*grid.R_[i]));
        }
    }

    return 0;
}
