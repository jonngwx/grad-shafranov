/*!
 *  @file util.h
 *  @brief A few helper functions.
 */ 

#include<stdio.h>

/*!
 * Helper function to make an array from min to max with n points 
 * replicates matlab's linspace
 * @param[in] min first point
 * @param[in] max last point
 * @param[in] n number of points
 * @param[out] *array pointer to array to store result
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
 * @param[in] file file pointer to open output file
 * @param[in] data 1d array to print
 * @param[in] n number of array elements
 *
 * Should file be constant (for const correctness?)
*/
inline void print1d(FILE *file, const double *data, int n){
    for (int i = 0; i < n; ++i){
        fprintf(file, "%15.8f ", data[i]);
    }
}

/*!
 * @brief prints a 2d array to a single line
 * @param[in] file file pointer to open output file
 * @param[in] data 2d array to print data[x][y]
 * @param[in] nx number of elements in x direction
 * @param[in] ny number of elements in y direction
 * 
 * Should file be constant (for const correctness?)
 */
inline void print2d(FILE *file, const double * const * data, int nx, int ny){
    for (int i = 0; i < nx; ++i){
        for (int j = 0; j < ny; ++j){
            fprintf(file, "%15.8f ", data[i][j]);
        }
    }
}
