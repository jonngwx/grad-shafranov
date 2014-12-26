#include <iostream>
#include <fstream>
#include "create_options.h"

int CreateOptions(int ac, char * av[], po::options_description &visible, po::variables_map &vm) {

  // Declare a group of options that will be 
  // allowed only on command line
  po::options_description generic("Generic options");
  generic.add_options()
      ("version,v", "print version std::string")
      ("help", "produce help message")
      ("config,c", po::value<std::string>()->default_value("grad-shafranov.cfg"),
            "name of a file of a configuration.")
      ;

  // Declare a group of options that will be 
  // allowed both on command line and in
  // config file
  po::options_description config("Configuration");
  config.add_options()
      ("coil-data-name,C", po::value< std::string>()->default_value("coil_data.tsv"), "coil data file name")
      ("pgtype", po::value<std::string>(), "p and g specification type (alpha, ..., ...)")
      ("p-filename", po::value<std::string>()->default_value("pg_data.tsv"), "p file name")
      ("g-filename", po::value<std::string>()->default_value("pg_data.tsv"), "g file name")
      ("grid-elems-r", po::value<int>()->default_value(10), "number of grid elements in r dimension")
      ("grid-elems-z", po::value<int>()->default_value(10), "number of grid elements in z dimension")
      ("r-min", po::value<double>()->default_value(5.0), "grid's minimum r location (meters)")
      ("r-max", po::value<double>()->default_value(10.0), "grid's maximum r location (meters)")
      ("z-min", po::value<double>()->default_value(-3.0), "grid's minimum z location (meters)")
      ("z-max", po::value<double>()->default_value(3.0), "grid's maximum z location (meters)")
      ("max-iter-N,N", po::value<int>()->default_value(20), "maximum inner iterations")
      ("error-tol-N", po::value<double>()->default_value(1.0e-4), "error tolerance inner loop")
      ("max-iter-M,M", po::value<int>()->default_value(1), "maximum outer iterations")
      ("error-tol-M", po::value<double>()->default_value(1.0e-4), "error tolerance outer loop")
      ("max-iter-crit", po::value<int>()->default_value(100), "max iterations for solving for the critical points")
      ("error-tol-crit",po::value<double>()->default_value(1.0e-3), "converge criterion for finding critical points")
      ("z_limiter1",po::value<double>()->default_value(2), "Horizontal limiter point for calculating critical points")
      ("z_limiter2",po::value<double>()->default_value(-2), "Horizontal limiter point for calculating critical points")
     ;

  /* For pgtype = "alpha"
   * n1, n2, p0 (pressure at magnetic axis), and g0 (R0*B0)
   */
  po::options_description pgtype_alpha("For pgtype=alpha");
  pgtype_alpha.add_options()
      ("pgta-n1", po::value<double>()->default_value(1.0), "n1")
      ("pgta-n2",  po::value<double>()->default_value(1.0), "n2")
      ("pgta-p0",  po::value<double>()->default_value(12000.0), "p0: pressure on axis")
      ("pgta-g0",  po::value<double>()->default_value(7.5), "g0: R0*B0 (loc of mag axis times field at that point)")
      ;
   
  /* for inputting j-phi
   * Ip, D, R0, z0, 
   */
  po::options_description j_phi("For specifying initial current distribution");
  j_phi.add_options()
      ("j-phi-Ip", po::value<double>()->default_value(1000000.0), "total plasma current")
      ("j-phi-D", po::value<double>()->default_value(2.0), "radius of initial curr dist")
      ("j-phi-R0", po::value<double>(), "initial R value of mag axis")
      ("j-phi-z0", po::value<double>(), "initial z value of mag axis")      ;

  po::options_description outputs("Output format");
  outputs.add_options()
      ("output-fields", po::value<std::string>()->default_value("none"), "comma-separated list of J,Bt,...")
      ("output-type", po::value<std::string>()->default_value("tsv"), "tsv or hdf5 (won't work unless you've compiled with hdf5)")
      ("output-name", po::value<std::string>()->default_value("cougar.out"), "prefix for output filename")
      ("output-every-n", po::value<int>()->default_value(-1), "when to write output for convergence testing of inner loop")
      ("output-every-m", po::value<int>()->default_value(-1), "when to write output for convergence testing of outer loop")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(pgtype_alpha).add(j_phi).add(outputs);

  po::options_description config_file_options;
  config_file_options.add(config).add(pgtype_alpha).add(j_phi).add(outputs);

  // Add to the visible options list (which will be returned
  // for possible printing later)
  visible.add(generic).add(config).add(pgtype_alpha).add(j_phi).add(outputs);
  
  //store in vm
  store(po::command_line_parser(ac, av).
        options(cmdline_options).run(), vm);
  notify(vm);
  
  std::string config_file = vm["config"].as<std::string>();
  
  std::ifstream ifs(config_file.c_str());
  if (!ifs)
  {
    // we require a config file. its abscence will be dealt with in main
  }
  else
  {
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  }

  return 0;
}
