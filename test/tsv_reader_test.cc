/*! 
 * @file tsv_reader_test.cc
 * @author Jacob Schwartz
 * @brief Some tests for the Table and CoilData classes using the boost_test library.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tsv_reader_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "../include/tsv_reader.h"

struct tbtest {
  tbtest() { BOOST_TEST_MESSAGE( "Setup fixture tbtest"); }
  ~tbtest(){ BOOST_TEST_MESSAGE( "Teardown fixture tbtest"); }
  Table td;  
};

struct cdtest {
  cdtest() { BOOST_TEST_MESSAGE ("Setup fixture cdtest"); }
  ~cdtest(){ BOOST_TEST_MESSAGE( "Teardown fixture cdtest"); }
  CoilData cd;  
};

BOOST_FIXTURE_TEST_SUITE( suite1, tbtest)

BOOST_AUTO_TEST_CASE(two_good_cols) {
  int success = td.load_from_tsv("test/tsvReaderExamples/good_example_2_cols.tsv",1);
  BOOST_REQUIRE_EQUAL(success, 0);
  BOOST_CHECK_EQUAL(td.num_columns(), 2);
  BOOST_CHECK_EQUAL(td.num_rows(), 3);
  BOOST_CHECK_EQUAL(td.data(1,0), 3.4);
}

BOOST_AUTO_TEST_CASE(three_good_cols) {
  int success = td.load_from_tsv("test/tsvReaderExamples/good_example_3_cols.tsv",1);
  BOOST_REQUIRE_EQUAL(success, 0);
  BOOST_CHECK_EQUAL(td.num_columns(), 3);
}

BOOST_AUTO_TEST_CASE(good_example_comments){
  int success = td.load_from_tsv("test/tsvReaderExamples/good_example_comments.tsv",1);
  BOOST_REQUIRE_EQUAL(success, 0);
}

BOOST_AUTO_TEST_CASE(good_example_comments_2){
  int success = td.load_from_tsv("test/tsvReaderExamples/good_example_comments_2.tsv");
  BOOST_REQUIRE_EQUAL(success, 0);
}

BOOST_AUTO_TEST_CASE(bad_example_no_file){
  std::cout << "\nNext line should show an error:\n";
  int error_code = td.load_from_tsv("test/tsvReaderExamples/doesnt_exist.tsv",1);
  BOOST_REQUIRE_EQUAL(error_code,1);
}

BOOST_AUTO_TEST_CASE(bad_example_inconsistent_columns){
  std::cout << "\nNext line should show an error:\n";
  int error_code = td.load_from_tsv("test/tsvReaderExamples/bad_example_inconsistent_column_number.tsv",1);
  BOOST_CHECK_EQUAL(error_code, 2);
}

BOOST_AUTO_TEST_SUITE_END()
/***************
 * Test CoilData
 **************/
BOOST_FIXTURE_TEST_SUITE( suite2, cdtest)

BOOST_AUTO_TEST_CASE(bad_example_not_three_columns){
  std::cout << "\nNext line should show an error:\n";
  int error_code = cd.load_from_tsv("test/tsvReaderExamples/good_example_2_cols.tsv",1);
  BOOST_REQUIRE_EQUAL(error_code, 4);
}

BOOST_AUTO_TEST_SUITE_END()
