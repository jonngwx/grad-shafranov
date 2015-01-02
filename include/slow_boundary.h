/*!
 * @file slow_boundary.h
 * @author ???
 * @brief Header declarations for SlowBoundary
 */

#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "field.h"
#include "grid.h"
#include "tsv_reader.h"
#include "boundary.h"

/*!
 * @brief something something
 */
class SlowBoundary : public Boundary {
  public:
    /*!
     * Constructor of SlowBoundary.
     * @param grid The grid for which this class will calculate boundary
     *  conditions.
     * @param cond_data The external coil currents which 
     *  influence the boundary conditions for psi.
     */
    SlowBoundary(Grid* grid, CoilData* cond_data); 
    ~SlowBoundary();
    /*! Which of these parameters are 'in' or 'out'?*
     * @return What integer does it return?*/
    int CalcB(Field* psi, Field* jphi);
  private:
    double *R_;
    double *z_;
    double dr_;          // Size of a grid cell in the r direction
    double dz_;          // Size of a grid cell in the z direction
    CoilData* cond_data; // Contains data on the external coil currents
    int perim_;          // The number of perimeter cells of the grid
    double ***g_plasma_; /*  A three dimensional array. 
     * The value g_plasma[i][j][l] 
     * represents how much psi is at a given boundary point l
     * due to a unit current at point [i][j] in the plasma.
     * Note that for cells on the boundary, 
     * g_plasma[ce,ll][corresponding cell] = 0. */
};

#endif // SLOWBOUNDARY_H_
