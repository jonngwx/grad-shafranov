/*! \file */ 
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
