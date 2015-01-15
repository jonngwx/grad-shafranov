/*! 
 * @file green_fcn.h
 * @author Peter J. Bolgert
 * @brief delcaration for green_fcn()
 */

#ifndef GREEN_FCN_H
#define GREEN_FCN_H

/*!
 * @brief get Psi at R1, Z1, due to current at R2, Z2
 *
 * Evaluates the Green's function that gives Psi at point R1, Z1,
 * from a toroidal current at radius R2 from axis and height Z2.
 * This function is found in Johnson 1978, equation 9.
 *
 * @param[in] R1 radius of the field point
 * @param[in] Z1 height of the field point
 * @param[in] R2 radius of the source point
 * @param[in] Z2 height of the source point
 *
 * @return The value of psi from this one neat toroidal current element.
 */
double green_fcn(double R1, double Z1, double R2, double Z2);

#endif 
