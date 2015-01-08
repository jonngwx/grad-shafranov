#ifndef INTERPOLATE_H
#define INTERPOLATE_H
#include "field.h"
#include "grid.h"
#include <vector>

/*!
 * @brief The cougar (Puma concolor), also commonly known as the
 */
class Interpolate {
public:
    Interpolate(Grid &GridS, Field &Psi);
    ~Interpolate();
    
    double bicubicInterpolate(double r, double z);
    double bicubicInterpolate_r(double r, double z);
    double bicubicInterpolate_rr (double r, double z);
    double cubicInterpolate(std::vector<double>, double x);
    double cubicInterpolate_r(std::vector<double> arr, double x);
    double cubicInterpolate_rr(std::vector<double> arr, double x);
    
    double Psi_interp(double r, double z);
    double Psir_interp(double r, double z);
    double Psirr_interp(double r, double z);
    double Psirz_interp(double r, double z);
    double Psizz_interp(double r, double z);
    double Psiz_interp(double r, double z);
    /*!
     * @brief Perform search for critical points of limiter beginning
     * with initial guess r, z
     */
    void updateCoefficients();
    void updateP(double r, double z);
private:
    double r_curr; //!< i coordinate of center of current interpolation.
    double z_curr; //!< j coordinate of center of current interpolation.
    Field &Psi_; //!<
    Grid &Grid_; //!< Underlying grid.
    std::vector<std::vector<double>> P; //!< 4x4 points that are the values at the grid nodes around the point being interpolated.
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
};

#endif
