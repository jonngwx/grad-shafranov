#ifndef GRID_H
#define GRID_H
#include <stdlib.h>

class Grid{
public:
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
