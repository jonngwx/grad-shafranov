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
    Interpolate(Grid &GridS, Field &F);
    ~Interpolate();
    
    /*
     * @brief Interpolate the value of the field F at the given point.
     * @param r The radial location to interpolate F at.
     * @param z The vertical location to interpolate F at.
     * @return The interpolated value of the field F at that point.
     */
    double F(double r, double z) const;
    /*
     * @brief Interpolate the value of dF/dr at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the derivative of the field dF/dr at that point.
     */
    double F_r(double r, double z) const;
    /*
     * @brief Interpolate the value of d^2F/dr^2 at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/dr^2 at that point.
     */
    double F_rr(double r, double z) const;
    /*
     * @brief Interpolate the value of d^2F/drdz at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/drdz at that point.
     */
    double F_rz(double r, double z) const; 
    /*
     * @brief Interpolate the value of d^2F/dz^2 at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/dz^2 at that point.
     */
    double F_zz(double r, double z) const;
    /*
     * @brief Interpolate the value of dF/dz at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the derivative of the field dF/dz at that point.
     */
    double F_z(double r, double z) const;
    
    /*!
     *@brief calls updateP and updateCoefficients(). This prevents errors as both need to be updated each step.
     */
    void updateInterpolation(double r, double z);
private:

    /*!
     * @brief Checks that the target point is within the region
     * that has currently-loaded coefficients, or throws an error.
     * @param[in] r The target point's radial location in meters.
     * @param[in] z The target point's vertical location in meters.
     */
    void CorrectCellBoundsCheck(double r, double z) const;

    /*!
     * @brief Loads the 16 main array grid points surrounding the target point into P.
     * @param[in] r The target point's radial location in meters.
     * @param[in] z The target point's vertical location in meters.
     *
     * P is a 4x4 vector that contains the grid values surrounding the
     * target point (which is typically not on a grid point).
     * It is a small slice of the main array.
     * This function loads those 16 surrounding points from the array
     * given the target point's location. 
     */
    void updateP(double r, double z);

    /*!
     * @brief Update the bicubic polynomial coefficients for interpolation in the target point's cell. 
     * 
     * This interpolator uses bicubic interpolation, which requires the 16 = 4x4 grid points that 
     * surround a target non-grid point. Based on the function values on the grid points we generate coefficients
     * \f$ a_{m,n}\f$ for the bicubic polynomial where m and n are indicies ranging from 0 to 3.
     * The bicubic polynomial is \f$ \Sigma a_{m,n} x^m y^n\f$.
     * Coefficient formulas are from:
     * www.paulinternet.nl/?page=bicubic
     */
    void updateCoefficients();
    double r_curr_; //!< i coordinate of center of current interpolation.
    double z_curr_; //!< j coordinate of center of current interpolation.
    Field &F_;   //!< The Field being interpolated over.
    Grid &Grid_; //!< Underlying grid.
    const double dr_; //!< Radial grid spacing. From Grid.
    const double dz_; //!< Vertical grid spacing. From Grid.
    std::vector<std::vector<double>> P_; //!< 4x4 points that are the values at the grid nodes around the point being interpolated.
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
};

#endif
