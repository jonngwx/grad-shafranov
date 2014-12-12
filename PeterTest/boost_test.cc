#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <stdio.h>
#include <math.h>

using namespace boost::math;

int main() {
    

    
    double K = ellint_1(0.5); // K(k)
    double E = ellint_2(0.5); // E(k)
    printf("K(0.5) = %f and E(0.5) = %f \n", K, E);
   







}
