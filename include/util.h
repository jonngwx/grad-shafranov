/*! \file */ 

#include<stdio.h>

/**
Helper function to make an array from min to max with n points
replicates matlab's linspace
@param min first point
@param max last point
@param n number of points
@param *array pointer to array to store result
*/
void inline linspace(double min, double max, int n, double *array){
  array[0] = min;
  const double dn = (max - min)/(n - 1.0);
  for (int i = 1; i < n; i++){
    array[i] = array[i-1] + (dn);
  }
}


/**
 * @brief prints a 1d array
 * @param file file pointer to open output file
 * @param data 1d array to print
 * @param n number of array elements
 */
inline void print1d(FILE *file, double *data, int n){
    for (int i = 0; i < n; ++i){
        fprintf(file, "%15.8f ", data[i]);
    }
}

/**
 * @brief prints a 2d array to a single line
 * @param file file pointer to open output file
 * @param data 2d array to print data[x][y]
 * @param nx number of elements in x direction
 * @param ny number of elements in y direction
 */
inline void print2d(FILE *file, double **data, int nx, int ny){
    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            fprintf(file, "%15.8f ", data[i][j]);
        }
    }
}
