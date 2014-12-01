#ifndef FIELD_H
#define FIELD_H

class Field(){
	public:
	Field(int nr, int nz);
	
	
	const int nr;
	const int nz;
	double *R;
	double *z;
	double **f;
	double **f_old;
	
}


#endif