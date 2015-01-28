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
#include "j_solver.h"
#include "j_solver_nstx.h"

#ifdef HDF_MODE
#include "grad_output_hdf.h"
#endif 

using namespace std;

/* This function writes a paraboloid-like initial current profile to jphi. 
 *
 * J_phi = c/(R0 + r Cos[theta]) ( 1 - (r/D)^2 ) 
 *
 * where R0 is the radial location of the magnetic axis,
 * r is the minor radius (distance from (R0, z0)),
 * theta the angle as measured from the magnetic axis location (R0, z0),
 * D is the half-diameter,
 * and c is a normalization factor so that the total plasma current
 * is Ip. c has units of Amps/meter.
 *
 * In order to figure out this normalization factor, integrate J_phi:
 * In Mathematica, 
 * Integrate[
  c/(R0 + r Cos[\[Theta]]) (1 - r^2/D0^2) r, {r, 0, D0}, {\[Theta], 0, 
    2 \[Pi]}, Assumptions -> {R0 > 0, r > 0, r < R0, 0 < D0, D0 < R0}]
 *
 * should work, then solve for (the output) = Ip.
 */
void InitializeCurrentDensity(double Ip, double R0, double z0, double D, Field * jphi){
    const Grid * grid = jphi->grid_;
    double a = pow(D/R0,2); // just for convenience.  should be less than 1
    double c = Ip*a/((4.0*M_PI/3.0) * R0 * (pow(1-a,1.5) - (1 - 1.5*a)));
    double r_squared;
    for (int i = 0; i < grid->nr_; ++i) {
        for (int j = 0; j < grid->nz_; ++j) {
            r_squared = pow(grid->R_[i] - R0,2) + pow(grid->z_[j] - z0,2);
            if (r_squared < D*D) {
                jphi->f_[i][j] = (c/grid->R_[i])*(1 - r_squared/(D*D)) ;
            }
            else {
                jphi->f_[i][j] = 0;
            } 
        }
    }
}

/* 
 * This function creates a GradOutput of the correct type based on whether HDF is used and what the output_type is specified as.
 */
GradOutput * DetermineGradOutput(Field * psi, Field * jphi, Grid * grid, Field * p, Field * g, std::string & output_type, std::string & output_list){
  GradOutput * grad_output;
  if (output_type == "tsv") {
    grad_output = new GradOutputTxt(psi,jphi,grid,p,g,output_list.c_str());
  } else if (output_type == "hdf5") {
    #ifdef HDF_MODE
    grad_output = new GradOutputHdf(psi,jphi,grid,p,g,output_list.c_str());
    #else 
    grad_output = new GradOutputTxt(psi,jphi,grid,p,g,output_list.c_str());
    output_type = "tsv";
    printf("Output type hdf not supported. Recompile with hdf libraries to enable. Defaulting to tsv. \n");
    #endif 
  } else {
    printf("Output type %s is not supported, use tsv or hdf5. Defaulting to tsv. \n", output_type.c_str());
    output_type = "tsv";
    grad_output = new GradOutputTxt(psi,jphi,grid,p,g,output_list.c_str());
  }
  return grad_output;
}

/************************************************************************/

int main(int argc, char *argv[])
{
  clock_t time0 = clock();
  /*********************************************************
   * Load options from config file and command line into vm
   ********************************************************/
  po::options_description visible_options("Allowed options");
  po::variables_map vm;
  int error_code = CreateOptions(argc, argv, visible_options, vm);
  if( error_code!= 0 ) {
      exit(error_code);
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
      std::cout << "Grad Shafranov solver, version 0.9\n";
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
  double error_epsilon_N = vm["error-tol-N"].as<double>();
  double error_epsilon_M = vm["error-tol-M"].as<double>();

  CoilData cd; cd.load_from_tsv(vm["coil-data-name"].as<string>(),1);
  Table limiters; limiters.load_from_tsv(vm["limiter-name"].as<string>(),1);
  
  Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
  Field *psi = new Field(*grid);
  Field *jphi = new Field(*grid);
  Field *p = new Field(*grid);
  Field *g = new Field(*grid);

  /****************************************************************
   * Initialize the plasma current jphi.
   * The initialization for current looks the same as Johnson 1979 Step A.
   ***************************************************************/
  // Set up the magnetic axis location.
  double R0;
  if (vm.count("j-phi-R0")) {
    R0 = vm["j-phi-R0"].as<double>();
  } else {
    // Guess that the magnetic axis is at the middle of the grid. 
    R0 = (Rmax + Rmin)/2.0;
  }
  double z0;
  if (vm.count("j-phi-z0")) {
    z0 = vm["j-phi-z0"].as<double>();
  } else {
    z0 = (zmax + zmin)/2.0; 
  }

  double Ip = vm["j-phi-Ip"].as<double>();
  InitializeCurrentDensity(Ip,                          // Total plasma current
                           R0,                          // Plasma current r-center
                           z0,                          // Plasma current z-center
                           vm["j-phi-D"].as<double>(),  //The plasma current half-diameter
                           jphi);

  // Set up J_Solver class
  double P0 = vm["pgta-p0"].as<double>();
  double g0 = vm["pgta-g0"].as<double>();
  JSolver *js;
  std::string jname = vm["J-solver-type"].as<string>();
  if (jname == "nstx"){
    double n2 = vm["pgta-n2"].as<double>();
    js= new JSolverNSTX(P0,g0,Ip,n2,grid);
  } else {
    double n1 = vm["pgta-n1"].as<double>();
    double n2 = vm["pgta-n2"].as<double>();
    js = new JSolverAlpha(P0,g0,n1,n2,Ip,grid);
  }
  
  // Set up Elliptic solver for inner loop
  double error_ES = vm["error-tol-ES"].as<double>();
  EllipticSolver *solver = new GaussSeidel(*grid, *psi, error_ES);
  Boundary *psib = new SlowBoundary(psi, grid, &cd);

  // Set up Critical
  Critical *crit;
  try{
    crit = new Critical(*grid, *psi, vm["max-iter-crit"].as<int>(),
                                     vm["error-tol-crit"].as<double>(), 
                                     limiters,
                                     vm["R-stag-up"].as<double>(),
                                     vm["z-stag-up"].as<double>(),
                                     vm["R-stag-down"].as<double>(),
                                     vm["z-stag-down"].as<double>(),
                                     R0, z0);
  } catch (int i) {
    printf("Failed to setup critical! Abort!\n");
    return 1;
  }

  /** Determine which output type: tsv or hdf5 */
  std::string output_type = vm["output-type"].as<string>();
  std::string output_list = vm["output-fields"].as<string>();
  GradOutput *grad_output = DetermineGradOutput(psi, jphi, grid, p, g, output_type, output_list);
  std::string output_filename_base = vm["output-name"].as<string>();
  clock_t time1 = clock();

  /************************************************
   * Main program loop 
   ***********************************************/
  solver->coeff();
  for (int m = 0; m < maxIterM; ++m) {
    psib->CalcB(jphi);
    // output during calculation 
    if (m > 0 && outputEveryM > 0 && ((m % outputEveryM) == 0)){
      std::string partial_output_name = output_filename_base + ".m" + std::to_string(m) + "." + output_type;
      grad_output->write_output(partial_output_name.c_str());
      printf("Writing output for m = %d\n",m);
    }
    // test convergence
    if (m>0) {
      printf("Error for outer loop is %f \n", psib->norm());
      if (psib->norm() < error_epsilon_M){
        break;
      }
    }
    // Iterate through elliptic solver
    for (int n = 0; n < maxIterN; ++n) {
      //printf("n = %i \n", n);
      if (n == 0) {
        solver->step_1(*jphi);
      } else {
        solver->step(*jphi);
      }
      crit->update();
      js->update(jphi, psi, p, g);
      //printf("error norm = %f \n", solver->norm());
      //printf("iteration # n = %d, m = %d\n", n, m);
      if (solver->norm() < error_epsilon_N){ 
        printf("Inner loop converged at N = %d\n", n);
        break;
      }

      // output during calculation
      if (outputEveryN > 0 && ((n % outputEveryN) == 0)){
        std::string partial_output_name = output_filename_base + ".n" + std::to_string(n) + ".m" + std::to_string(m) + "." + output_type;
        grad_output->write_output(partial_output_name.c_str());
        printf("Writing output for n = %d, m = %d\n",n,m);
      }
      if (n == maxIterN-1) { printf(" Elliptic solver reached maxIterN without convergence"); }
    }//end inner loop
    if (m == maxIterM-1){ printf("boundary reached max iter M without convergence\n"); }
  }//end outer loop
  clock_t time2 = clock();

  /************************************************
   * Write final output and close 
   ***********************************************/
  std::string full_output_name = output_filename_base + "." + output_type;
  grad_output->write_output(full_output_name.c_str());
  
  delete grad_output;
  delete crit;
  delete psi;
  delete psib;
  delete jphi;
  delete solver;
  delete p;
  delete g;
  delete js;
  delete grid;

  float t1 = ((float)(time1-time0))/CLOCKS_PER_SEC;
  float t2 = ((float)(time2-time1))/CLOCKS_PER_SEC;
  printf("init time = %f, solver time = %f\n",t1,t2);
}

