/*!
 * @file util_test.cc
 * @author Jacob Schwartz
 * @brief Test that util's linspace function works as advertised.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE util_test
#include <boost/test/unit_test.hpp>
#include <stdlib.h>
#include "util.h"

struct linsptest {
  linsptest() {
    BOOST_TEST_MESSAGE( "Setup fixture linsptest"); 
    num_nodes = 20;
    low_bound = 10.0;
    high_bound = 20.0;
    arr = new double[num_nodes];
    linspace(low_bound, high_bound, num_nodes, arr);
  }
  ~linsptest() {
    BOOST_TEST_MESSAGE( "Teardown fixture linsptest"); 
    delete arr; 
  }
  int num_nodes;
  double low_bound;
  double high_bound;
  double * arr;
};


BOOST_FIXTURE_TEST_SUITE( suite1, linsptest)

  BOOST_AUTO_TEST_CASE(test_low_bound){
    BOOST_CHECK_EQUAL(arr[0],low_bound);
  }

  BOOST_AUTO_TEST_CASE(test_high_bound){
    BOOST_CHECK_EQUAL(arr[num_nodes-1],high_bound);
  }

  BOOST_AUTO_TEST_CASE(test_one_in_the_middle){
    BOOST_CHECK_EQUAL(arr[1], low_bound + (high_bound - low_bound)/(num_nodes - 1.));
  }

BOOST_AUTO_TEST_SUITE_END()
