#ifndef FIELD_H
#define FIELD_H

#include <stdlib.h>

/**
 * \file Header file for Field class
 * */

class Field(){
	public:
	Field(double R0, double Rend, double z0, double zend, int nr, int nz);
	~Field();
	
	const int nr; /** < number of points in R direction */
	const int nz; /** < number of points in z direction */
	double *R; /** < pointer to array of radial grid points */
	double *z; /** < pointer to array of vertical grid points */
	double **f; /** < pointer to field data */
	double **f_old; /** < pointer to old field data */
	double dr; /** < grid spacing in radial direction */
	double dz; /** < grid spacing in vertical direction */
	
private:
	/**
	 * Returns a the value of the ith boundary cell
	 * TODO: do we really need this?
	 * 
	 * @param i index on the boundary?
	 * */
	double f_boundary(int i);
};


#endif