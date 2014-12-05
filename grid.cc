#include "grid.h"

Grid::Grid(double R0, double Rend, double z0, double zend, int nr, int nz)
   : nr_(nr),
     nz_(nz),
     z0_(z0),
     zend_(zend),
     R0_(R0),
     Rend_(Rend) {
    
    R_ = new double[nr_]();
    z_ = new double[nz_]();	
    dr_ = (Rend_ - R0_)/(nr_ - 1);
    dz_ = (zend_ - z0_)/(nz_ - 1);
    // FIXME make the grid
}

Grid::~Grid(){
    delete [] z_;
    delete [] R_;
}

// Given r returns i coordinate of containing grid cell
// returned as double so user can get position in cell
double Grid::celli(double r) {
    return (r - R0_)*(1.0/dr_);
}

// Given z returns j coordinate of containing grid cell
// returned as double so user can get position in cell
double Grid::cellj(double z) {
    return (z - z0_)*(1.0/dz_);
}
