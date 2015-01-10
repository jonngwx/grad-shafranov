/*! 
 * @file gs_solver.cc
 * @author The COUGAR Team
 * @brief The main program for COUGAR
 */

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <time.h>
#include "include/grid.h"
#include "include/field.h"
#include "include/grad_output.h"
#include "include/grad_output_txt.h"
#include "include/sor.h"
#include "include/gauss_seidel.h"
#include "include/slow_boundary.h"
#include "include/tsv_reader.h"
#include "include/create_options.h"
#include "include/j_solver_alpha.h"
#include "include/critical.h"

#ifdef HDF_MODE
#include "grad_output_hdf.h"
#endif 

using namespace std;

int main(int argc, char *argv[])
{
  
    clock_t time0 = clock();
    /*********************************************************
     * Load options from config file and command line into vm
     ********************************************************/
    po::options_description visible_options("Allowed options");
    po::variables_map vm;
    int return_value = CreateOptions(argc, argv, visible_options, vm);
    if( return_value != 0 ) {
        exit(return_value);
    }

    /************************************************
     * Process options that would print and then end
     * the program before any computation
     ***********************************************/
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
    int outputEveryN = vm["output-every-n"].as<int>();
    int outputEveryM = vm["output-every-m"].as<int>();
    double error_epsilon = vm["error-tol-N"].as<double>();

    /* commented out because pgtype is not used anywhere else in main
     * so this block is useless. If we're not planning
     * to use it we should delete it.
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
    }*/
  
    CoilData cd; cd.load_from_tsv(vm["coil-data-name"].as<string>(),1);

    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);
    Field *p = new Field(*grid);
    Field *g = new Field(*grid);

    /****************************************************************
     * Initialize the plasma current (jphi) and psi.
     * The initialization for current looks the same as Johnson 1979 Step A.
     * I'm not sure where this psi initialization comes from. -JAS
     ***************************************************************/
    double r_squared;
    double R0;
    if (vm.count("j-phi-R0")) {
      R0 = vm["j-phi-R0"].as<double>();
    } else {
      // not necessarily true. Just for now...
      R0 = (Rmax + Rmin)/2.0;
    }
    double z0;
    if (vm.count("j-phi-z0")) {
      z0 = vm["j-phi-z0"].as<double>();
    } else {
      // not necessarily true. Just for now...
      z0 = (zmax + zmin)/2.0; 
    }
    double D = vm["j-phi-D"].as<double>();
    double Ip = vm["j-phi-Ip"].as<double>();

    /* The initial plasma current profile is this parabola-like 
     *
     * J_phi = c/(R0 + r Cos[theta]) ( 1 - (r/D)^2 ) 
     *
     * where R0 is the radius of the magnetic axis,
     * r is the minor radius,
     * D is the half-diameter (seriously why have we named it D??)
     * and c is a normalization factor so that the sum of the plasma current
     * is Ip. c has units of Amps/meter.
     *
     * In order to figure out this normalization factor, integrate J_phi:
     *
     * In Mathematica, 
     *
     * Integrate[
      c/(R0 + r Cos[\[Theta]]) (1 - r^2/D0^2) r, {r, 0, D0}, {\[Theta], 0, 
        2 \[Pi]}, Assumptions -> {R0 > 0, r > 0, r < R0, 0 < D0, D0 < R0}]
     *
     * should work, then solve for (the output) = Ip.
     */
    double a = pow(D/R0,2); // just for convenience.  should be less than 1
    double c = Ip*a/((4.0*M_PI/3.0) * R0 * (pow(1-a,1.5) - (1 - 1.5*a)));
    printf("c = %f \n", c);

    for (int i = 0; i < grid->nr_; ++i) {
        for (int j = 0; j < grid->nz_; ++j) {
            r_squared = pow(grid->R_[i] - R0,2) + pow(grid->z_[j] - z0,2);
            if (r_squared < D*D) {
                jphi->f_[i][j] = (c/grid->R_[i])*(1 - r_squared/(D*D)) ;
            }
            else {
                jphi->f_[i][j] = 0;
            } 
            //What does the number 10 mean?
            //            psi->f_[i][j] = 10*exp(-pow(grid->R_[i] - R0,2))*exp(-pow(grid->z_[j],2)/10.);
        }
    }

    // set up JSolverAlpha class
    double P0 = vm["pgta-p0"].as<double>();
    double g0 = vm["pgta-g0"].as<double>();
    double n1 = vm["pgta-n1"].as<double>();
    double n2 = vm["pgta-n2"].as<double>();
    JSolverAlpha *jsa = new JSolverAlpha(P0,g0,n1,n2,Ip,grid);
    
    // Elliptic solver for inner loop
    EllipticSolver *solver = new GaussSeidel(*grid, *psi);
    Boundary *psib = new SlowBoundary(grid, &cd);

    // set up Critical
    double phys_lim_z_up = vm["phys_lim_z_up"].as<double>();
    double phys_lim_R_up = vm["phys_lim_R_up"].as<double>();
    double phys_lim_z_down = vm["phys_lim_z_down"].as<double>();
    double phys_lim_R_down = vm["phys_lim_R_down"].as<double>();

    double R_stag_up = vm["R_stag_up"].as<double>();
    double z_stag_up = vm["z_stag_up"].as<double>();
    double R_stag_down = vm["R_stag_down"].as<double>();
    double z_stag_down = vm["z_stag_down"].as<double>();

    int max_iter_crit = vm["max-iter-crit"].as<int>();
    double error_tol_crit = vm["error-tol-crit"].as<double>();
    Critical *crit = new Critical(*grid, *psi, max_iter_crit, error_tol_crit, phys_lim_R_up, phys_lim_z_up, phys_lim_R_down, phys_lim_z_down, R_stag_up, z_stag_up, R_stag_down,  z_stag_down, R0, z0);

    /** determine which output type: tsv or hdf5 */
    Grad_Output *grad_output;
    std::string output_type = vm["output-type"].as<string>();
    std::string output_list = vm["output-fields"].as<string>();
    if (output_type == "tsv") {
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,output_list.c_str());
    } else if (output_type == "hdf5") {
        #ifdef HDF_MODE
        grad_output = new Grad_Output_Hdf(psi,jphi,grid,p,g,output_list.c_str());
        #else 
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,output_list.c_str());
	output_type = "tsv";
        printf("Output type hdf not supported. Recompile with hdf libraries to enable. Defaulting to tsv. \n");
        #endif 
    } else {
        printf("Output type %s is not supported, use tsv or hdf5. Defaulting to tsv. \n", output_type.c_str());
        output_type = "tsv";
        grad_output = new Grad_Output_Txt(psi,jphi,grid,p,g,output_list.c_str());
    }
    clock_t time1 = clock();
    std::string output_filename_base = vm["output-name"].as<string>();

    // solve stuff   "What???" -JAS
    solver->coeff();
    /************************************************
     * Main program loop 
     ***********************************************/
    for (int m = 0; m < maxIterM; ++m) {
        psib->CalcB(psi, jphi);
        // output during calculation 
        if (outputEveryM > 0 && ((m % outputEveryM) == 0)){
            std::string partial_output_name = output_filename_base + ".m" + std::to_string(m) + "." + output_type;
            grad_output->write_output(partial_output_name.c_str());
            printf("Writing output for m = %d\n",m);
        }
        // test convergence
        // Iterate through elliptic solver
        for (int n = 0; n < maxIterN; ++n) {
            //printf("n = %i \n", n);
            if (n == 0) {
              solver->step_1(*jphi);
            } else {
              solver->step(*jphi);
            }
            crit->update();
            jsa->update(jphi, psi, p, g);
            //printf("error norm = %f \n", solver->norm());
            //printf("iteration # n = %d, m = %d\n", n, m);
            if (solver->norm() < error_epsilon){ 
              break;
            }

            // output during calculation
            if (outputEveryN > 0 && ((n % outputEveryN) == 0)){
                std::string partial_output_name = output_filename_base + std::to_string(n) + ".m" + std::to_string(m) + "." + output_type;
                grad_output->write_output(partial_output_name.c_str());
                printf("Writing output for n = %d, m = %d\n",n,m);
            }
            if (n == maxIterN-1) {
              printf(" Elliptic solver reached maxIterN without convergence\n");
            }
        }//end inner loop
        if (m == maxIterM-1){
            printf("boundary reached max iter M without convergence\n");
        }
    }//end outer loop
    clock_t time2 = clock();
  /************************************************
   * Write final output and close 
   ***********************************************/
  std::string full_output_name = output_filename_base + "." + output_type;
  grad_output->write_output(full_output_name.c_str());
  
  float t1 = ((float)(time1-time0))/CLOCKS_PER_SEC;
  float t2 = ((float)(time2-time1))/CLOCKS_PER_SEC;
  printf("init time = %f, solver time = %f\n",t1,t2);
  delete grad_output;
  delete crit;
  delete psi;
  delete psib;
  delete jphi;
  delete solver;
  delete p;
  delete g;
  delete jsa;
  delete grid;
}

