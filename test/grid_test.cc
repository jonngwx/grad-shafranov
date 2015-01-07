#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Grid
#define protected public
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "grid.h"
#include <stdlib.h>

/*!
 * Tests grid celli
 */
BOOST_AUTO_TEST_CASE (test_celli) {
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 10;
    double nz = 10;

    Grid grid(R0, Rend, Z0, Zend, nr, nz);
    
    // check that lower bound is working
    double a = grid.celli(-10);
    BOOST_CHECK_CLOSE(a,0,6);

    // check that upper bound is working
    double b = grid.celli(50);
    BOOST_CHECK_CLOSE(b,nr-1,6);

    // check interpolation
    double c = grid.celli(2.5);
    BOOST_CHECK_CLOSE(c,2.25,.5);
}


/*!
 * Tests grid cellj
 */
BOOST_AUTO_TEST_CASE (test_cellj) {
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 10;
    double nz = 10;

    Grid grid(R0, Rend, Z0, Zend, nr, nz);
    
    // check that lower bound is working
    double a = grid.cellj(-10);
    BOOST_CHECK_CLOSE(a,0,6);

    // check that upper bound is working
    double b = grid.cellj(50);
    BOOST_CHECK_CLOSE(b,nr-1,6);

    // check interpolation
    double c = grid.cellj(2.5);
    BOOST_CHECK_CLOSE(c,6.7568,.5);
}

