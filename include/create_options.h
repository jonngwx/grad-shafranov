/**
 * @file create_options.h
 * @author Jacob Schwartz
 * @brief Header for the create_options method and a helper printer method.
 */
#ifndef CREATE_OPTIONS_H_
#define CREATE_OPTIONS_H_

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// @brief Prints a vector to an ostream, I think. 
/// From the boost library examples in boost_program_options file multiple_sources.cpp
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " ")); 
    return os;
}

/**
 * @brief Loads options from command line and config file.
 * @param[in] ac	main's argc: the number of command line tokens(?)
 * @param[in] av	main's argv: the string of the command line 
 * @param[out] visible	Loaded up with the list of Visible options (so it can be printed in main if help is called)
 * @param[out] vm	A key-value pair map of options.
 * @return 0 if successful, some other integer if a problem.
 */
int CreateOptions(int ac, char * av[], po::options_description &visible, po::variables_map &vm);

#endif
