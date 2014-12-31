#ifndef GRID_H
#define GRID_H
#include <stdlib.h>

/*!
 * @brief Stores information about the solution grid.
 *
 * Stores the x and y axes of grid used in the solver and 
 * methods to determine cell indices from coordinates.
 */
class Grid{
public:
    /*!
     * Constructor of grid 
     * @param R0 - value of R at left boundary
     * @param Rend - value of R at right boundary
     * @param z0 - value of z at lower boundary
     * @param zend - value of z at upper boundary
     * @param nr - number of grid points in the R direction
     * @param nz - number of points in the z direction
     *
     */
    Grid(double R0, double Rend, double z0, double zend, int nr, int nz);
    ~Grid();
    double celli(double r);
    double cellj(double z);
	
    const int nr_; /** < number of points in R direction */
    const int nz_; /** < number of points in z direction */
    double *R_; /** < pointer to array of radial grid points */
    double *z_; /** < pointer to array of vertical grid points */
    double dr_; /** < grid spacing in radial direction */
    double dz_; /** < grid spacing in vertical direction */  

};

#endif
