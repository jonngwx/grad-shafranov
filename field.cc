#include "field.h"

Field::Field(int nr, int nz) : nr(nr), nz(nz) {
	f = new double*[nr];
	for (int i = 0; i < nr; ++i){
		f[i] = new double[nz]();
	}
}

Field::~Field(){
	for (int i = 0; i < nr; ++i){
		delete [] f[i];
	}
	delete [] f;
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
