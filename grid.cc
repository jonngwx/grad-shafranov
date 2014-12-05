#include "grid.h"

Grid::Grid(double R0, double Rend, double z0, double zend, int nr, int nz) : nr(nr), nz(nz), z0(z0), R0(R0) {
    R = new double[nr]();
    z = new double[nz]();	
    dr = (Rend-R0)/(nr-1);
    dz = (zend-z0)/(nz-1);
    // FIXME make the grid
}

Grid::~Grid(){
    delete [] z;
    delete [] R;
}

// Given r returns i coordinate of containing grid cell
int Grid::celli(double r) {
    return (int) (r - R0)*(1.0/dr);
}

// Given z returns j coordinate of containing grid cell
int Grid::cellj(double z) {
    return (int) (z - z0)*(1.0/dz);
}