#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Output
#define protected public
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <math.h>
#include "field.h"
#include "grid.h"
#include "grad_output_txt.h"
#include <stdlib.h>
#include <stdio.h>

/*!
 * Tests the parsing of parse_output
 */
BOOST_AUTO_TEST_CASE (test_parse_output) {
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 10;
    double nz = 10;
    Grid *grid = new Grid(R0, Rend, Z0, Zend, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);
    Field * p =  new Field(*grid);
    Field * g = new Field(*grid);

    Grad_Output_Txt grad_output(psi,jphi,grid,p,g,"j");
    BOOST_CHECK_EQUAL(grad_output.output_list.size(),1);

    Grad_Output_Txt grad_output1(psi,jphi,grid,p,g,"gibberish");
    BOOST_CHECK_EQUAL(grad_output1.output_list.size(),0);

    Grad_Output_Txt grad_output2(psi,jphi,grid,p,g,"j,j,j,j,j");
    BOOST_CHECK_EQUAL(grad_output2.output_list.size(),1);

    Grad_Output_Txt grad_output3(psi,jphi,grid,p,g,"j,bt");
    BOOST_CHECK_EQUAL(grad_output3.output_list.size(),2);
    BOOST_CHECK_EQUAL(grad_output3.output_list[1],Grad_Output::TOROIDAL_FIELD);
    BOOST_CHECK_EQUAL(grad_output3.output_list[0],Grad_Output::CURRENT);
    delete p;
    delete g;
    delete psi;
    delete jphi;
    delete grid;

}


