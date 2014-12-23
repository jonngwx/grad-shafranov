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

class SlowBoundary : public Boundary {
  public:
    SlowBoundary(Grid* grid, CoilData* cond_data); 
    ~SlowBoundary();
    /*! Which of these parameters are 'in' or 'out'?*
     * @return What integer does it return?*/
    int CalcB(Field* psi, Field* jphi);
  private:
    double *R_;
    double *z_;
    double dr_;
    double dz_;
    CoilData* cond_data_;
    int perim_;
    double ***g_plasma_;
};

#endif // SLOWBOUNDARY_H_

