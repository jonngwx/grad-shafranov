/*! 
 * @file boundary_test.cc
 * @author Jacob Schwartz
 * @brief Some tests for the Boundary class using the boost_test library.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE boundary_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "../include/boundary.h"

struct boundarytest {
  boundarytest() { 
    BOOST_TEST_MESSAGE( "Setup fixture boundarytest");
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 3;
    double nz = 3;
    
    gr = new Grid(R0, Rend, Z0, Zend, nr, nz);
    f = new Field(*gr);
    b = new Boundary(f, gr); 
  }
  ~boundarytest(){ 
    BOOST_TEST_MESSAGE( "Teardown fixture boundarytest"); 
    delete gr;
    delete f;
    delete b;
  }

  Grid * gr;
  Field * f;
  Boundary * b;

};

BOOST_FIXTURE_TEST_SUITE( bdy_suite, boundarytest)

BOOST_AUTO_TEST_CASE(LToI) {
  const std::vector<int> correct_output = { 0, 1, 2, 2, 2, 1, 0, 0 };
  for (int l = 0; l < 8; l++){
    printf("hellp! l = %i \n", l);
    BOOST_CHECK_EQUAL(b->LtoI(l), correct_output[l]);
  }
}

BOOST_AUTO_TEST_CASE(LToJ) {
  const std::vector<int> correct_output = { 0, 0, 0, 1, 2, 2, 2, 1 };
  for (int l = 0; l < 8; l++){
    BOOST_CHECK_EQUAL(b->LtoJ(l), correct_output[l]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
