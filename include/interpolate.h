/*!
 * @file interpolate.h
 * @brief Declaration for class Interpolate.
 */
#ifndef INTERPOLATE_H
#define INTERPOLATE_H
#include "field.h"
#include "grid.h"
#include <vector>

/*!
 * @brief Interpolates to find the values and derivatives of a 2D data set between grid points.
 *
 * Uses bicubic interpolation as described on http://paulinternet.nl/?page=bicubic.
 * Given a Field (with its associated Grid) and then given a point within[1]
 * that field, even not on one of the grid nodes but in a 'cell' between four
 * nodes, we can use the values of the Field on the surrounding 16 = 4x4 nodes
 * to construct a polynomial \f$ \Sigma_{m=0}^3 \Sigma_{n=0}^3 a_{m,n} x^m y^n\f$
 * that represents the function over the whole grid cell.
 *
 * We can then calculate approximations to the function and to its first and second
 * derivatives over the whole cell.
 */
class Interpolate {
public:
    Interpolate(Grid &GridS, Field &F);
    ~Interpolate();
    
    /*!
     * @brief Interpolate the value of the field F at the given point.
     * @param r The radial location to interpolate F at.
     * @param z The vertical location to interpolate F at.
     * @return The interpolated value of the field F at that point.
     */
    double F(double r, double z) const;
    /*!
     * @brief Interpolate the value of dF/dr at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the derivative of the field dF/dr at that point.
     */
    double F_r(double r, double z) const;
    /*!
     * @brief Interpolate the value of d^2F/dr^2 at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/dr^2 at that point.
     */
    double F_rr(double r, double z) const;
    /*!
     * @brief Interpolate the value of d^2F/drdz at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/drdz at that point.
     */
    double F_rz(double r, double z) const; 
    /*!
     * @brief Interpolate the value of d^2F/dz^2 at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the second derivative of the field d^2F/dz^2 at that point.
     */
    double F_zz(double r, double z) const;
    /*!
     * @brief Interpolate the value of dF/dz at the given point.
     * @param r The radial location to interpolate at.
     * @param z The vertical location to interpolate at.
     * @return The interpolated value of the derivative of the field dF/dz at that point.
     */
    double F_z(double r, double z) const;
    
    /*!
     * @brief Calls updateP and updateCoefficients(). This prevents errors as both need to be updated each step.
     * @param[in] r a radial location inside the cell for which the bicubic polynomial will be generated. 
     * @param[in] z a vertical location inside the cell for which the bicubic polynomial will be generated. 
     */
    void updateInterpolation(double r, double z);

    void PrintAmnCoefficients();
private:

    /*!
     * @brief Checks that the target point is within the region
     * that has currently-loaded coefficients, or throws an error.
     * @param[in] r The target point's radial location in meters.
     * @param[in] z The target point's vertical location in meters.
     */
    void CorrectCellBoundsCheck(double r, double z) const;

    /*!
     * @brief Loads the 16 main array grid points surrounding the target point into P_.
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
     * The bicubic polynomial is \f$ \Sigma_{m=0}^3 \Sigma_{n=0}^3 a_{m,n} x^m y^n\f$.
     * Coefficient formulas are from:
     * http://www.paulinternet.nl/?page=bicubic
     */
    void updateCoefficients();
    double r_curr_; //!< i coordinate of cell of current interpolation.
    double z_curr_; //!< j coordinate of cell of current interpolation.
    Field &F_;   //!< The Field being interpolated over.
    Grid &Grid_; //!< Underlying grid.
    const double dr_; //!< Radial grid spacing. From Grid.
    const double dz_; //!< Vertical grid spacing. From Grid.
    std::vector<std::vector<double>> P_; //!< 4x4 points that are the values at the grid nodes around the point being interpolated.
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
};

#endif
