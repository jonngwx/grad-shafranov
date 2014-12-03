#include "field.h"

Field::Field(double R0, double Rend, double z0, double zend, int nr, int nz) : nr(nr), nz(nz) {
	R = new double[nr]();
	z = new double[nz]();
	f = new double*[nr];
	f_old = new double*[nr];
	dr = (Rend-R0)/(nr-1);
	dz = (zend-z0)/(nz-1);
	for (int i = 0; i < nr; ++i){
		f[i] = new double[nz]();
		f_old[i] = new double[nz]();
	}
}

Field::~Field(){
	delete [] R;
	delete [] z;
	for (int i = 0; i < nr; ++i){
		delete [] f[i];
		delete [] f_old[i];
	}
	delete [] f;
	delete [] f_old;
}


/**
 * @brief Evaluates the value of f at the boundary
 * @param i the index of the boundary value to return
 * @return value of f at the boundary
 */
double Field::f_boundary(int i)
{
	return 0.;
}
