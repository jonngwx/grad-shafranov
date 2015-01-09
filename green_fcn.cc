/*! 
 * @file green_fcn.cc
 * @author Peter J. Bolgert
* @brief implements green_fcn()
 */
#include <math.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace boost::math;

double green_fcn(double R1, double Z1, double R2, double Z2) {
  if ((R1 == R2) && (Z1 == Z2)) {
    return 0;
  } else {
    double k =
        sqrt(4.0 * R1 * R2 / ((R1 + R2) * (R1 + R2) + (Z1 - Z2) * (Z1 - Z2)));
    /*Numerical coefficient below is 1 over 2 pi */
    double g = -(1.0 / (2 * M_PI)) * (sqrt(R1 * R2) / k) *
               ((2 - k * k) * ellint_1(k) - 2 * ellint_2(k));
    return g;
  }
}
