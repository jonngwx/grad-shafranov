#include "field.h"

Field::Field(int nr, int nz) : nr_(nr), nz_(nz) {
	f_ = new double*[nr_];
	for (int i = 0; i < nr_; ++i){
		f_[i] = new double[nz_]();
	}
}

Field::~Field(){
	for (int i = 0; i < nr_; ++i){
		delete [] f_[i];
	}
	delete [] f_;
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
