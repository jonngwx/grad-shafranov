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
      ("coil-data-name,c", 
           po::value< std::string>()->default_value("coil_data.tsv"), 
           "coil data file name")
      ("pgtype", po::value<std::string>(), "p and p specification type")
      ("g-filename", po::value<std::string>(), "g file name")
      ("p-filename", po::value<std::string>(), "p file name")
      ("max-iter-N,N", po::value<int>()->default_value(10), 
            "maximum inner iterations")
      ("error-tol-N,eN", po::value<double>()->default_value(1e-5),
            "error tolerance inner loop")
      ("max-iter-M,M", po::value<int>()->default_value(10),
            "maximum outer iterations")
      ("error-tol-M,eM", po::value<double>()->default_value(1e-5),
            "error tolerance outer loop")
      ("output-type", po::value<std::string>()->default_value("tsv"), "tsv or hdf5")
      ("output name", po::value<std::string>()->default_value("cougar.out"), "prefix for output filename")
      ;

  // Hidden options, will be allowed both on command line and
  // in config file, but will not be shown to the user.
  po::options_description hidden("Hidden options");
  hidden.add_options()
      ("unused", po::value< std::vector<std::string> >(), "placeholder")
      ;

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(hidden);

  po::options_description config_file_options;
  config_file_options.add(config).add(hidden);

  // Add to the visible options list (which will be returned
  // for possible printing later)
  visible.add(generic).add(config);
  
  po::positional_options_description p;
  p.add("unused", -1);
  
  //store in vm
  store(po::command_line_parser(ac, av).
        options(cmdline_options).positional(p).run(), vm);
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
