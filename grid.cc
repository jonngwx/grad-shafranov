/*!
 * @file grid.cc
 * @author ???
 * @brief Implementation for class Grid.
 */

#include "include/grid.h"
#include "include/util.h"
#include <math.h>

Grid::Grid(double R0, double Rend, double z0, double zend, int nr, int nz)
   : nr_(nr),
     nz_(nz) {
    
    R_ = new double[nr_]();
    z_ = new double[nz_]();	
    dr_ = (Rend - R0)/(nr - 1.0);
    dz_ = (zend - z0)/(nz - 1.0);
    //Fill up R and Z
    linspace(R0,Rend,nr,R_);
    linspace(z0,zend,nz,z_);
}

Grid::~Grid(){
    delete [] z_;
    delete [] R_;
}

double Grid::celli(double r) {
    double i = (r - R_[0])*(1.0/dr_);
    //ensure that the returned position is between 0 and nr_-1.
    return fmin(nr_-1,fmax(i,0));
}

double Grid::cellj(double z) {
    double j = (z - z_[0])*(1.0/dz_);
    return fmin(nz_-1,fmax(j,0));
}
