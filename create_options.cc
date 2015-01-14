/*!
 * @file create_options.cc
 * @brief defines the CreateOptions function
 * @author Jacob Schwartz
 */

#include <iostream>
#include <fstream>
#include "create_options.h"

/*!
 * @brief loads options from command line and config file.
 * Called once from main. Encapsulated here to make main shorter.
 * @param[in] ac main's argc
 * @param[in] av main's argv
 * @param[out] visible a boost::program_options::options_description which will be loaded up with options marked 'visible', to be printed if the user says to print the help.
 * @param[out] vm a boost::program_options::variables_map which will be loaded up with keys and values of the options specified on the command line and in the config file. If an option is unspecified it will default to the value specified below in this function.
 * @return 0 to indicate success (in the future nonzero returns could indicate an error, but none have been implemented)
 */
int CreateOptions(int ac, char * av[], po::options_description &visible, po::variables_map &vm) {

  // Declare a group of options that will be 
  // allowed only on command line
  po::options_description generic("Generic options");
  generic.add_options()
      ("version,v", "print version std::string")
      ("help", "produce help message")
      ("config,c", po::value<std::string>()->default_value("grad-shafranov.cfg"),
            "Main configuration file name")
      ;

  // Declare a group of options that will be 
  // allowed both on command line and in
  // config file
  po::options_description config("Main program and solver configuration");
  config.add_options()
      ("coil-data-name", po::value< std::string>()->default_value("coil_data.tsv"), "Coil data file name")
      ("limiter-name", po::value< std::string>()->default_value("limiters.tsv"), "Limiter data file name")
      ("J-solver-type", po::value< std::string>()->default_value("alpha"), "Which solver to use to determine current density. [nstx | alpha]")
      //("pgtype", po::value<std::string>(), "p and g specification type (alpha, ..., ...)")
      //("p-filename", po::value<std::string>()->default_value("pg_data.tsv"), "p file name")
      //("g-filename", po::value<std::string>()->default_value("pg_data.tsv"), "g file name")
      ("max-iter-N,N", po::value<int>()->default_value(20), "Maximum inner loop iterations")
      ("error-tol-N", po::value<double>()->default_value(1.0e-4), "Error tolerance for inner loop")
      ("max-iter-M,M", po::value<int>()->default_value(1), "Maximum outer loop iterations")
      ("error-tol-M", po::value<double>()->default_value(1.0e-4), "Error tolerance for outer loop")
      ("max-iter-crit", po::value<int>()->default_value(100), "Max iterations for solving for the critical points")
      ("error-tol-crit",po::value<double>()->default_value(1.0e-3), "Convergence criterion for finding critical points")
      ("error-tol-ES", po::value<double>()->default_value(1.0e-5), "Error tolerance for internal iter. in GaussSeidel")
      ;

  po::options_description geometry_opts("Grid and geometry");
  geometry_opts.add_options()
      ("grid-elems-r", po::value<int>()->default_value(10), "Number of grid elements in r dimension")
      ("grid-elems-z", po::value<int>()->default_value(10), "Number of grid elements in z dimension")
      ("r-min", po::value<double>()->default_value(5), "Grid's minimum r location (meters)")
      ("r-max", po::value<double>()->default_value(10), "Grid's maximum r location (meters)")
      ("z-min", po::value<double>()->default_value(-3), "Grid's minimum z location (meters)")
      ("z-max", po::value<double>()->default_value(3), "Grid's maximum z location (meters)")
      ;

  po::options_description stagnation_point_opts("Stagnation point initial geometry");
  stagnation_point_opts.add_options()
      ("R-stag-down",po::value<double>()->default_value(.5), "R location of stagnation point in bottom half of tokamak")
      ("z-stag-down",po::value<double>()->default_value(-1.5), "z location of stagnation point in bottom half of tokamak")
      ("R-stag-up",po::value<double>()->default_value(.5), "R location of stagnation point in top half of tokamak")
      ("z-stag-up",po::value<double>()->default_value(1.5), "z location of stagnation point in top half of tokamak")
      ;
  /* 
   * n1, n2, p0 (pressure at magnetic axis), and g0 (R0*B0)
   */
  po::options_description pgta_opts("The form and magnitude of the pressure and current");
  pgta_opts.add_options()
      ("pgta-n1", po::value<double>()->default_value(1.0), "n1: An exponent for the form of the pressure distribution.")
      ("pgta-n2",  po::value<double>()->default_value(1.0), "n2: An exponent for the form of the current distribution.")
      ("pgta-p0",  po::value<double>()->default_value(12000.0), "p0: Pressure on axis, measured in Pascals")
      ("pgta-g0",  po::value<double>()->default_value(7.5), "g0: R0*B0 (rad of mag axis times field at that point)")
      ;
   
  /* for inputting j-phi
   * Ip, D, R0, z0, 
   */
  po::options_description j_phi("The initial current distribution");
  j_phi.add_options()
      ("j-phi-Ip", po::value<double>()->default_value(1000000.0), "Total plasma current")
      ("j-phi-D", po::value<double>()->default_value(1.0), "Radius of initial curr distribution")
      ("j-phi-R0", po::value<double>(), "Initial R value of magnetic axis")
      ("j-phi-z0", po::value<double>(), "Initial z value of magnetic axis")
      ;

  po::options_description outputs("Output format and frequency");
  outputs.add_options()
      ("output-fields", po::value<std::string>()->default_value("none"), "Fields to output: comma-separated list of J,Bt,...")
      ("output-type", po::value<std::string>()->default_value("tsv"), "tsv or hdf5 (won't work unless you've compiled with hdf5)")
      ("output-name", po::value<std::string>()->default_value("cougar.out"), "Prefix for output filename")
      ("output-every-n", po::value<int>()->default_value(-1), "Write the state to file every n inner loop iterations. -1 is off.")
      ("output-every-m", po::value<int>()->default_value(-1), "Write the state to file every m outer loop iterations. -1 is off.")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(geometry_opts).add(stagnation_point_opts).add(pgta_opts).add(j_phi).add(outputs);

  po::options_description config_file_options;
  config_file_options.add(config).add(geometry_opts).add(stagnation_point_opts).add(pgta_opts).add(j_phi).add(outputs);

  // Add to the visible options list (which will be returned
  // for possible printing later)
  visible.add(generic).add(config).add(geometry_opts).add(stagnation_point_opts).add(pgta_opts).add(j_phi).add(outputs);
  
  //store in vm
  store(po::command_line_parser(ac, av).
        options(cmdline_options).run(), vm);
  notify(vm);
  
  std::string config_file = vm["config"].as<std::string>();
  
  std::ifstream ifs(config_file.c_str());
  if (!ifs)
  {
    // We require a config file. Its absence will be dealt with in main.
  }
  else
  {
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  }

  return 0;
}
