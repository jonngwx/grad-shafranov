#ifndef FIELD_H
#define FIELD_H

#include <stdlib.h>

class Field(){
	public:
	Field(double R0, double Rend, double z0, double zend, int nr, int nz);
	~Field();
	
	const int nr;
	const int nz;
	double *R;
	double *z;
	double **f;
	double **f_old;
	double dr;
	double dz;
	
	private:
	double f_boundary(int i);
};


#endif