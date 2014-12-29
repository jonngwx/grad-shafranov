#include "include/grid.h"
#include "include/util.h"

Grid::Grid(double R0, double Rend, double z0, double zend, int nr, int nz)
   : nr_(nr),
     nz_(nz) {
    
    R_ = new double[nr_]();
    z_ = new double[nz_]();	
    dr_ = (Rend - R0)/(nr - 1.0);
    dz_ = (zend - z0)/(nz - 1.0);
    linspace(R0,Rend,nr,R_);
    linspace(z0,zend,nz,z_);
}

Grid::~Grid(){
    delete [] z_;
    delete [] R_;
}

/*!
 * Given r returns i coordinate of containing grid cell
 * returned as double so user can get position in cell
 */
double Grid::celli(double r) {
    double i = (r - R_[0])*(1.0/dr_);
//    printf("i = %f\n\n", i);
    return i;
}

/*!
 * Given z returns j coordinate of containing grid cell
 * returned as double so user can get position in cell
 */
double Grid::cellj(double z) {
    double j = (z - z_[0])*(1.0/dz_);
//    printf("j = %f\n\n", j);
    return j;
}
