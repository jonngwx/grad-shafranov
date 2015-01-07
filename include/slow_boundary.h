/*!
 * @file slow_boundary.h
 * @author ???
 * @brief Header declarations the SlowBoundary class.
 */

#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "field.h"
#include "grid.h"
#include "tsv_reader.h"
#include "boundary.h"

/*!
 * @brief This boundary-value calculating class iterates over the field of j_phi
 * and (should do ... not implemented yet) over the set of external coils) using
 * Green's functions.
 */
class SlowBoundary : public Boundary {
 public:
  /*!
   * Constructor of SlowBoundary.
   *
   * Initializes the Green's function array: a cache of values representing the
   * effect of psi at the boundary from a unit current at any point on the grid.
   *
   * @param grid The grid for which this class will calculate boundary
   *  conditions.
   * @param cond_data The external coil currents which
   *  influence the boundary conditions for psi.
   *  Currently the affect of these coils has not been implemented.
   */
  SlowBoundary(Grid* grid);
  SlowBoundary(Grid* grid, CoilData* cond_data);
  ~SlowBoundary();
  /*!
   * @brief Calculate the updated boundary values for psi based on the plasma
   * current.
   *
   * Should we add another parameter of the external coil currents???
   *
   * @param[out] psi The field to which these new boundary values will be
   *written.
   * @param[in] jphi The plasma current.
   * @return The number 0 to indicate success.
   */
  int CalcB(Field* psi, Field* jphi);

 private:
  double* R_;
  double* z_;
  double dr_;            // Size of a grid cell in the r direction
  double dz_;            // Size of a grid cell in the z direction
  CoilData* cond_data_;  // Contains data on the external coil currents
  int perim_;            // The number of perimeter cells of the grid
  double*** g_plasma_;   /*  A three dimensional array.
   * The value g_plasma[i][j][l]
   * represents how much psi is at a given boundary point l
   * due to a unit current at point [i][j] in the plasma.
   * Note that for cells on the boundary,
   * g_plasma[ce,ll][corresponding cell] = 0. */
  double** g_coils_;  
};

#endif  // SLOWBOUNDARY_H_
