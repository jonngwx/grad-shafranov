#ifndef CREATE_OPTIONS_H_
#define CREATE_OPTIONS_H_

#include <iostream>
#include <iterator>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Prints a vector to an ostream, I think. 
// From the boost library examples in boost_program_options
//  file multiple_sources.cpp
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " ")); 
    return os;
}

int CreateOptions(int ac, char * av[], po::options_description &visible, po::variables_map &vm);


#endif
