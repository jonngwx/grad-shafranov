#ifndef GRID_H
#define GRID_H
#include <stdlib.h>

class Grid{
public:
    Grid(double R0, double Rend, double z0, double zend, int nr, int nz);
    ~Grid();
	
    const int nr; /** < number of points in R direction */
    const int nz; /** < number of points in z direction */
    double *R; /** < pointer to array of radial grid points */
    double *z; /** < pointer to array of vertical grid points */
    double dr; /** < grid spacing in radial direction */
    double dz; /** < grid spacing in vertical direction */  

}

#endif