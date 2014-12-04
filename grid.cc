#include "grid.h"

Grid::Grid(double R0, double Rend, double z0, double zend, int nr, int nz) : nr(nr), nz(nz) {
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