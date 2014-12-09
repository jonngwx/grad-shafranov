#include "include/field.h"

Field::Field(const Grid &grid) : grid_(&grid) {
	f_ = new double*[grid.nr_];
	for (int i = 0; i < grid.nr_; ++i){
		f_[i] = new double[grid.nz_]();
	}
}

Field::~Field(){
	for (int i = 0; i < grid_->nr_; ++i){
		delete [] f_[i];
	}
	delete [] f_;
}

