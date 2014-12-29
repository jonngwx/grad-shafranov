/*! 
 * @file green_fcn.cc
 * @author ???
 * @brief defines cc
 */
#include <math.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace boost::math;

/*!
 * @brief get Psi at R1, Z1, due to current at R2, Z2
 *
 * Evaluates the Green's function that gives Psi at point R1, Z1,
 * from a toroidal current at radius R2 from axis and height Z2.
 *
 * @param[in] R1 radius of the field point
 * @param[in] Z1 height of the field point
 * @param[in] R2 radius of the source point
 * @param[in] Z2 height of the source point
 *
 * @return The value of psi from this one neat toroidal current element.
 */
double green_fcn(double R1, double Z1, double R2, double Z2) {
  if ((R1 == R2) && (Z1 == Z2)) {
    return 0;
  } else {
    double k =
        sqrt(4.0 * R1 * R2 / ((R1 + R2) * (R1 + R2) + (Z1 - Z2) * (Z1 - Z2)));
    /*Oh gosh what is that number?*/
    double g = (1.0 / 6.28318530718) * (sqrt(R1 * R2) / k) *
               ((2 - k * k) * ellint_1(k) - 2 * ellint_2(k));
    return g;
  }
}
